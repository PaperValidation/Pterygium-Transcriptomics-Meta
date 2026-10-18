"""Verify every reference in the Paper A manuscript against Crossref (and PubMed).

For each reference we query Crossref (by DOI when cited, else by title), take the
best-matching published record, and compare journal, year, volume and first page with
what the manuscript cites.

Two subtleties that otherwise produce false alarms:

  * Journal names. The manuscript uses ISO abbreviations ("Sci Rep"); Crossref returns
    full titles ("Scientific Reports"). Tokens are compared after normalising the
    digraph "ph"->"p", which also covers Crossref's own misspelling of IOVS as
    "Investigative Opthalmology & Visual Science".
  * Year. Crossref's `issued` for a Frontiers article is the online date (2022-01-31)
    while PubMed's yearOfPublication is the volume year (2021). When Crossref and the
    citation disagree on the year, we ask Europe PMC; if PubMed agrees with the
    manuscript the citation stands and the discrepancy is recorded as a note.

Output: results/paper_a/reports/reference_verification.md
"""
import html
import json
import os
import re
import urllib.parse
import urllib.request

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(BASE, "submission", "Submission_Package_PaperA_FINAL")
OUT = os.path.join(BASE, "reports", "reference_verification.md")
UA = {"User-Agent": "PaperA-ref-check/1.0 (mailto:gaoweiz@163.com)",
      "Accept": "application/json"}


def fetch_json(url, timeout=30):
    req = urllib.request.Request(url, headers=UA)
    return json.loads(urllib.request.urlopen(req, timeout=timeout).read())


def norm_tok(t):
    return t.lower().replace("ph", "p")


def tok_close(a, b):
    a, b = norm_tok(a), norm_tok(b)
    if a == b:
        return True
    if len(a) < 2 or len(b) < 2:
        return a == b
    n = min(3, len(a), len(b))
    return b.startswith(a[:n]) or a.startswith(b[:n])


def journal_match(cited, crossref):
    if not cited or not crossref:
        return True
    # drop parenthetical place qualifiers used in some ISO abbreviations
    cited = re.sub(r"\s*\([^)]*\)", "", cited).strip()
    full = [t for t in re.split(r"[^a-z0-9]+", crossref.lower()) if t]
    abbr = [t for t in re.split(r"[^a-z0-9]+", cited.lower()) if t and len(t) > 1]
    if not abbr:
        return True
    skip = {"the", "of", "and", "for"}
    for a in abbr:
        if a in skip:
            continue
        if not any(tok_close(a, f) for f in full):
            return False
    return True


def sim(x, y):
    A = set(re.sub(r"[^a-z0-9 ]", " ", x.lower()).split())
    B = set(re.sub(r"[^a-z0-9 ]", " ", y.lower()).split())
    return len(A & B) / len(A | B) if A and B else 0.0


def parse_ref(txt):
    doi = ""
    dm = re.search(r"doi:\s*(10\.\S+)", txt)
    if dm:
        # greedy match, then drop the sentence-final period
        doi = dm.group(1).rstrip(".,;")

    m = re.search(r"\.\s+([A-Za-z][^.]{14,}?)\.\s+([A-Z][A-Za-z&:()\- ]{2,}?)\.\s+(\d{4})[;.]", txt)
    if not m:
        m = re.search(r"\.\s+([A-Za-z][^.]{14,}?)\.\s+([A-Z][A-Za-z&:()\- ]{2,}?)\.\s+(\d{4})", txt)
    title = m.group(1) if m else ""
    journal = m.group(2).strip() if m else ""
    year = m.group(3) if m else ""

    vm = re.search(r"(\d{4});\s*(\d+)(?:\((\d+)\))?\s*:\s*([0-9eE\-–]+)", txt)
    vol = vm.group(2) if vm else ""
    pages = vm.group(4) if vm else ""
    return dict(title=title, journal=journal, year=year, vol=vol, pages=pages, doi=doi)


def pubmed_year(doi):
    try:
        d = fetch_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
                       f"?query=DOI:%22{urllib.parse.quote(doi)}%22&format=json&resultType=core")
        res = d.get("resultList", {}).get("result", [])
        if res:
            return str(res[0].get("journalInfo", {}).get("yearOfPublication") or "")
    except Exception:
        pass
    return None


def main():
    md = open(os.path.join(PKG, "Manuscript_IOVS_v2.2.md"), encoding="utf-8").read()
    block = re.search(r"^## REFERENCES(.*?)^## FIGURE LEGENDS", md, re.S | re.M).group(1)
    refs = re.findall(r"^(\d+)\.\s+(.+)$", block, re.M)

    rows = []
    for num, txt in refs:
        p = parse_ref(txt)
        rec = dict(num=num, cited=p, status="", issues=[], notes=[], matched="")

        items = []
        if p["doi"]:
            try:
                items = [fetch_json("https://api.crossref.org/works/"
                                    + urllib.parse.quote(p["doi"]))["message"]]
            except Exception as e:
                rec["notes"].append(f"DOI lookup failed ({e})")
        if not items and p["title"]:
            try:
                q = urllib.parse.quote(p["title"])
                items = fetch_json("https://api.crossref.org/works?query.bibliographic="
                                   f"{q}&rows=5&select=DOI,title,container-title,issued,"
                                   "volume,issue,page")["message"]["items"]
            except Exception as e:
                rec["notes"].append(f"title query failed ({e})")

        if not items:
            rec["status"] = "UNVERIFIED"
            rows.append(rec)
            continue

        best, best_s = None, -9.0
        for it in items:
            t = html.unescape((it.get("title") or [""])[0])
            s = sim(p["title"], t)
            container = html.unescape((it.get("container-title") or [""])[0]).lower()
            if re.search(r"biorxiv|medrxiv|research square|ssrn|preprint", container):
                s -= 0.30
            if not container:
                s -= 0.20
            if not str(it.get("volume", "") or "").strip():
                s -= 0.10
            if not str(it.get("page", "") or "").strip():
                s -= 0.05
            if s > best_s:
                best, best_s = it, s

        if best is None or best_s < 0.45:
            rec["status"] = "NO CONFIDENT MATCH"
            rec["notes"].append(f"best title similarity {best_s:.2f}")
            rows.append(rec)
            continue

        cr_j = html.unescape((best.get("container-title") or [""])[0])
        cr_y = str((best.get("issued", {}).get("date-parts") or [[None]])[0][0] or "")
        cr_v = str(best.get("volume", "") or "")
        cr_p = str(best.get("page", "") or "")
        cr_doi = best.get("DOI", "")
        rec["matched"] = f"{cr_j} {cr_y};{cr_v}:{cr_p} (sim {best_s:.2f}) {cr_doi}"

        if p["vol"] and cr_v and p["vol"] != cr_v:
            rec["issues"].append(f"volume {p['vol']} vs {cr_v}")
        if p["pages"] and cr_p:
            a = re.split(r"[-–]", p["pages"])[0]
            b = re.split(r"[-–]", cr_p)[0]
            if a != b:
                rec["issues"].append(f"page {p['pages']} vs {cr_p}")
        if not journal_match(p["journal"], cr_j):
            rec["issues"].append(f"journal '{p['journal']}' vs '{cr_j}'")

        if p["year"] and cr_y and p["year"] != cr_y:
            pyr = pubmed_year(cr_doi) if cr_doi else None
            if pyr and pyr == p["year"]:
                rec["notes"].append(
                    f"Crossref issued {cr_y} (online date) but PubMed year {pyr} -> citation year correct")
            else:
                rec["issues"].append(f"year {p['year']} vs Crossref {cr_y}, PubMed {pyr}")

        rec["status"] = "OK" if not rec["issues"] else "MISMATCH"
        rows.append(rec)

    ok = [r for r in rows if r["status"] == "OK"]
    bad = [r for r in rows if r["status"] == "MISMATCH"]
    unv = [r for r in rows if r["status"] not in ("OK", "MISMATCH")]

    lines = ["# Reference verification — Paper A", "",
             f"Checked **{len(rows)}** references against Crossref; Europe PMC was used to "
             "adjudicate year conflicts between Crossref's online date and PubMed's volume year.",
             "", "| Result | Count |", "|---|---|",
             f"| Verified consistent | **{len(ok)}** |",
             f"| Requires correction | **{len(bad)}** |",
             f"| Could not verify | **{len(unv)}** |", "",
             "| # | Status | As cited | Crossref record | Notes |", "|---|---|---|---|---|"]
    for r in rows:
        c = r["cited"]
        cited = (f"{c['journal']} {c['year']};{c['vol']}:{c['pages']}"
                 if c["vol"] else f"{c['journal']} {c['year']} (no vol/pages)")
        note = "; ".join(r["issues"] + r["notes"])
        lines.append(f"| {r['num']} | {r['status']} | {cited} | {r['matched'] or '—'} | {note or '✓'} |")
    open(OUT, "w", encoding="utf-8").write("\n".join(lines) + "\n")

    for r in rows:
        flag = "OK  " if r["status"] == "OK" else ("BAD " if r["status"] == "MISMATCH" else "??  ")
        print(f"  {flag}[{r['num']:>2s}] {'; '.join(r['issues'] + r['notes'])[:118]}")
    print(f"\nOK {len(ok)}/{len(rows)}  mismatch {len(bad)}  unverified {len(unv)}")
    print(f"report -> {OUT}")


if __name__ == "__main__":
    main()
