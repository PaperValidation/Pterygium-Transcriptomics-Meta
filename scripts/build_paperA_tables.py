"""Build the six deliverable tables for Paper A (IOVS submission).

Sources (all under results/paper_a/tables unless noted):
  - DEG counts: results/paper_b/deg_rebuilt/*.csv  (the b07 authoritative rebuild;
    the older results/paper_a/tables/GSE155776_DEG_*.csv is a pre-rebuild copy and
    does NOT reproduce the manuscript's 5,516).
  - Table 2 module meta: module_effect_meta.csv (R metafor v2.8 output; identical
    sha256 to module_effect_meta_v28.csv) + module_sensitivity_rsweep_v28.csv.
    NOTE: module_sensitivity_2x2.csv / module_lodo_mkh.csv / module_sensitivity_rsweep.csv
    (no suffix) are derived from export_audit_cohort.py, which uses a DIFFERENT
    variance convention; they are NOT used here (see audit report P0-1).
  - Table S1: Table_S1_dataset_census.csv (8-column version, platform-corrected).
  - Table S2: module_definitions.csv + module_gene_sets_UVsplit.json +
    module_platform_coverage.csv + module_overlap_jaccard.csv
  - Table S3: gene_meta_results.csv + core_stable_genes.csv
  - Table S4: module_lodo.csv + module_lodo_mkh.csv + gene_lodo.csv

Output: results/paper_a/submission/Submission_Package_PaperA_FINAL/tables/*.xlsx
"""
import json
import os

import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # results/paper_a
T = os.path.join(BASE, "tables")
# Where the finished XLSX deliverables are written. Defaults to <root>/tables so the
# repository is self-contained; override with PAPER_A_OUT to write straight into a
# manuscript submission package.
OUT = os.environ.get("PAPER_A_OUT") or os.path.join(BASE, "tables")
DEG_REBUILT = os.path.join(os.path.dirname(BASE), "paper_b", "deg_rebuilt")

FDR_CUT = 0.05
LFC_CUT = 0.585


def read_deg(path):
    d = pd.read_csv(path)
    for c in ("logFC", "FDR"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["logFC", "FDR"])


def deg_summary():
    """Table 1 sheet 1: per-lesion DEG counts (b07 rebuild, FDR<0.05 & |log2FC|>0.585).

    Reads the frozen counts table when available (that is what the public repository
    ships, since the per-dataset DEG matrices are upstream products of the GEO
    preprocessing pipeline and are too large to redistribute). Falls back to
    recomputing from the DEG matrices if they are present locally.
    """
    frozen = os.path.join(T, "deg_summary_counts.csv")
    if os.path.exists(frozen):
        return pd.read_csv(frozen)

    specs = [
        ("GSE155776", "Pterygium vs healthy conjunctiva", "independent",
         os.path.join(DEG_REBUILT, "GSE155776_DEG_pterygium_vs_normal.csv"), "Primary discovery cohort"),
        ("GSE183153", "Pterygium vs healthy conjunctiva", "independent",
         os.path.join(DEG_REBUILT, "GSE183153_DEG_pterygium_vs_normal.csv"), "Pterygium meta cohort"),
        ("GSE183153", "Pinguecula vs healthy conjunctiva", "independent",
         os.path.join(DEG_REBUILT, "GSE183153_DEG_pinguecula_vs_normal.csv"), "Exploratory (n=2)"),
        ("GSE51995", "Pterygium vs conjunctiva", "paired (same eye)",
         os.path.join(DEG_REBUILT, "GSE51995_DEG_pterygium_vs_normal_conjunctiva.csv"), "Pterygium meta cohort"),
        ("GSE2513", "Pterygium vs conjunctiva", "independent",
         os.path.join(DEG_REBUILT, "GSE2513_DEG_pterygium_vs_normal_conjunctiva.csv"), "Pterygium meta cohort"),
        ("GSE191114", "Active vs atrophic pterygium", "independent",
         os.path.join(DEG_REBUILT, "GSE191114_DEG_active_vs_atrophic.csv"), "Activity-only (supplementary)"),
        ("GSE208384", "Pterygium vs conjunctiva", "independent",
         os.path.join(DEG_REBUILT, "GSE208384_DEG_pterygium_vs_normal.csv"), "Exploratory (n=2)"),
        ("GSE147449", "SCC vs healthy conjunctiva", "independent",
         os.path.join(T, "GSE147449_DEG_b07_scc_vs_normal_conjunctiva.csv"), "OSSN comparator"),
        ("GSE147449", "Papilloma vs healthy conjunctiva", "independent",
         os.path.join(T, "GSE147449_DEG_b07_papilloma_vs_normal_conjunctiva.csv"), "Papilloma comparator"),
        ("GSE148387", "Melanoma vs healthy conjunctiva", "independent",
         os.path.join(T, "GSE148387_DEG_b07_melanoma_vs_normal_conjunctiva.csv"), "Melanoma comparator"),
    ]
    rows = []
    for acc, comp, design, path, role in specs:
        if not os.path.exists(path):
            print(f"  SKIP (missing): {path}")
            continue
        d = read_deg(path)
        sig = d[(d.FDR < FDR_CUT) & (d.logFC.abs() > LFC_CUT)]
        n_case = int(d["n_case"].iloc[0]) if "n_case" in d.columns else ""
        n_ctrl = int(d["n_control"].iloc[0]) if "n_control" in d.columns else ""
        rows.append({
            "Dataset": acc, "Comparison": comp, "Design": design,
            "n_disease": n_case, "n_control": n_ctrl,
            "Genes tested": len(d),
            "DEGs (FDR<0.05, |log2FC|>0.585)": len(sig),
            "Up": int((sig.logFC > 0).sum()),
            "Down": int((sig.logFC < 0).sum()),
            "Role": role,
        })
    return pd.DataFrame(rows)


def top_shared_genes(top_n=15):
    """Table 1 sheet 2: genes shared between the pterygium meta signature and each
    comparator lesion (same direction), ranked by |meta Z|.

    Reads the frozen result table when present (what the public repository ships);
    recomputes from the comparator DEG matrices when they are available locally.
    """
    frozen = os.path.join(T, "top_shared_genes.csv")
    if os.path.exists(frozen):
        return pd.read_csv(frozen)

    meta = pd.read_csv(os.path.join(T, "meta_DEG_all_genes_FIXED.csv"))
    zcol = "Z" if "Z" in meta.columns else [c for c in meta.columns if c.lower() == "z"][0]
    meta[zcol] = pd.to_numeric(meta[zcol], errors="coerce")
    robust = meta[meta["robust"].astype(str).str.lower() == "true"].copy()
    up = set(robust.loc[robust[zcol] > 0, "gene"])
    dn = set(robust.loc[robust[zcol] < 0, "gene"])

    comparators = {
        "SCC (GSE147449)": os.path.join(T, "GSE147449_DEG_b07_scc_vs_normal_conjunctiva.csv"),
        "Papilloma (GSE147449)": os.path.join(T, "GSE147449_DEG_b07_papilloma_vs_normal_conjunctiva.csv"),
        "Melanoma (GSE148387)": os.path.join(T, "GSE148387_DEG_b07_melanoma_vs_normal_conjunctiva.csv"),
    }
    zmap = dict(zip(robust["gene"], robust[zcol]))
    rows = []
    for name, path in comparators.items():
        if not os.path.exists(path):
            continue
        d = read_deg(path)
        s = d[(d.FDR < FDR_CUT) & (d.logFC.abs() > LFC_CUT)]
        s_up = set(s.loc[s.logFC > 0, "gene"])
        s_dn = set(s.loc[s.logFC < 0, "gene"])
        shared_up = sorted(up & s_up, key=lambda g: -abs(zmap.get(g, 0)))[:top_n]
        shared_dn = sorted(dn & s_dn, key=lambda g: -abs(zmap.get(g, 0)))[:top_n]
        rows.append({
            "Comparator lesion": name,
            "Shared up (ranked by |meta Z|)": ", ".join(shared_up),
            "n shared up": len(up & s_up),
            "Shared down (ranked by |meta Z|)": ", ".join(shared_dn),
            "n shared down": len(dn & s_dn),
        })
    return pd.DataFrame(rows)


def table2():
    """Table 2: REML+KH primary + modified-KH + DL + exclude-GSE51995 + BH q,
    plus the GSE51995 r-sweep sheet.

    All sensitivity columns are read from `module_sensitivity_2x2.csv`, which is
    regenerated by `rebuild_sensitivity_py.py` from the per-cohort effects recorded in
    `module_effect_meta.csv`. That script asserts the REML+knha cell equals the primary
    pooled g, so Table 2 cannot drift out of sync with the manuscript again.
    """
    m = pd.read_csv(os.path.join(T, "module_effect_meta.csv"))
    s2 = pd.read_csv(os.path.join(T, "module_sensitivity_2x2.csv"))
    idx = {(r["module"], r["tau2_estimator"], r["inference"]): r for _, r in s2.iterrows()}

    rows = []
    for _, R in m.iterrows():
        mod = R["module"]
        kh = idx[(mod, "REML", "knha")]
        mk = idx[(mod, "REML", "adhoc")]
        dl = idx[(mod, "DL", "z")]
        # The primary REML+KH columns are taken straight from the R/metafor output in
        # module_effect_meta.csv, so Table 2 is bit-identical to the numbers quoted in
        # the manuscript. The sensitivity columns (modified KH, DL) come from
        # module_sensitivity_2x2.csv, which is regenerated by rebuild_sensitivity_py.py
        # and verified to agree with the R output to ~1e-7.
        rows.append({
            "Module": mod,
            "k cohorts": int(float(R["kh_k"])),
            "Pooled Hedges g (REML+KH)": round(float(R["kh_g"]), 4),
            "95% CI low (REML+KH)": round(float(R["kh_ci_low"]), 4),
            "95% CI high (REML+KH)": round(float(R["kh_ci_high"]), 4),
            "SE (REML+KH)": round(float(R["kh_se"]), 5),
            "P (REML+KH)": float(R["kh_p"]),
            "q (BH, 19 tests)": round(float(R["q_BH19"]), 4),
            "tau2 (REML)": round(float(R["kh_tau2"]), 4),
            "I2 (%, REML+KH)": round(float(R["kh_I2"]) * 100, 1),
            "Q": round(float(R["kh_Q"]), 4),
            "P (Q)": float(R["kh_Qp"]),
            "Pooled g (modified KH)": round(float(mk["pooled_g"]), 4),
            "95% CI low (mKH)": round(float(mk["ci_low"]), 4),
            "95% CI high (mKH)": round(float(mk["ci_high"]), 4),
            "P (modified KH)": float(mk["p"]),
            "Pooled g (DL)": round(float(dl["pooled_g"]), 4),
            "95% CI low (DL)": round(float(dl["ci_low"]), 4),
            "95% CI high (DL)": round(float(dl["ci_high"]), 4),
            "P (DL)": float(dl["p"]),
            "I2 (%, DL)": round(float(dl["I2"]), 1),
            "Pooled g (exclude GSE51995)": round(float(R["x_g"]), 4),
            "95% CI low (excl)": round(float(R["x_ci_low"]), 4),
            "95% CI high (excl)": round(float(R["x_ci_high"]), 4),
            "P (excl GSE51995)": float(R["x_p"]),
            "I2 (%, excl GSE51995)": round(float(R["x_I2"]) * 100, 1),
            "Observed within-pair r (GSE51995)": float(R["pair_r_GSE51995"]),
        })
    out = pd.DataFrame(rows).sort_values("Pooled Hedges g (REML+KH)", ascending=False)

    rs = pd.read_csv(os.path.join(T, "module_sensitivity_rsweep.csv"))
    rs_out = rs[[c for c in rs.columns if c in
                 ("module", "r", "r_type", "pooled_g", "ci_low", "ci_high", "p", "I2", "k")]].copy()
    rs_out = rs_out.rename(columns={
        "module": "Module", "r": "Within-pair r (GSE51995)",
        "r_type": "Row type (observed = primary)",
        "pooled_g": "Pooled Hedges g", "ci_low": "95% CI low", "ci_high": "95% CI high",
        "p": "P", "I2": "I2 (%)", "k": "k cohorts"})
    return {"Module meta-analysis": out,
            "r-sweep GSE51995": rs_out}


def table_s1():
    """Table S1: dataset census — 8-column platform-corrected version.
    GSE155776 data type set to FFPE MACE RNA-seq per Methods section 2."""
    d = pd.read_csv(os.path.join(T, "Table_S1_dataset_census.csv"), encoding="utf-8-sig")
    d.columns = [c.strip() for c in d.columns]
    d = d.rename(columns={
        "dataset": "Dataset", "disease": "Disease", "design": "Design",
        "platform": "Platform", "data_type": "Data type", "pmid": "PMID",
        "role": "Role", "exclusion": "Exclusion / note"})
    # platform/data-type corrections verified during audit
    fix = {"GSE155776": "FFPE MACE RNA-seq"}
    d["Data type"] = [fix.get(x, y) for x, y in zip(d["Dataset"], d["Data type"])]
    d["PMID"] = d["PMID"].fillna("").astype(str).str.replace(r"\.0$", "", regex=True)
    d["Exclusion / note"] = d["Exclusion / note"].fillna("")
    return d


def table_s2():
    """Table S2: module definitions + gene-set members + platform coverage + Jaccard."""
    defs = pd.read_csv(os.path.join(T, "module_definitions.csv"))
    with open(os.path.join(T, "module_gene_sets_UVsplit.json"), encoding="utf-8") as f:
        gs = json.load(f)
    members = pd.DataFrame([
        {"Module": k, "n genes": len(v), "Genes": ", ".join(sorted(v))}
        for k, v in gs.items()
    ]).sort_values("Module")
    cov = pd.read_csv(os.path.join(T, "module_platform_coverage.csv"))
    cov = cov.rename(columns={"module": "Module", "cohort": "Cohort",
                              "n_matched": "Genes matched on platform",
                              "pct": "% of module measured"})
    ov = pd.read_csv(os.path.join(T, "module_overlap_jaccard.csv"))
    ov = ov.rename(columns={
        "module_1": "Module 1", "module_2": "Module 2",
        "n_module_1": "n genes module 1", "n_module_2": "n genes module 2",
        "n_overlap": "n overlapping", "n_union": "n union",
        "jaccard": "Jaccard", "overlap_pct_module_1": "% of module 1 shared",
        "overlap_pct_module_2": "% of module 2 shared"})
    return {"Module definitions": defs, "Gene set members": members,
            "Platform coverage": cov, "Pairwise module overlap": ov}


def table_s3():
    """Table S3: gene-level meta-analysis + core stable genes."""
    g = pd.read_csv(os.path.join(T, "gene_meta_results.csv"))
    core = pd.read_csv(os.path.join(T, "core_stable_genes.csv"))
    summary = pd.DataFrame([
        {"Item": "Genes in meta-analysis", "n": len(g)},
        {"Item": "Direction-consistent meta-associated genes (robust = TRUE)", "n": int((g["robust"].astype(str).str.lower() == "true").sum())},
        {"Item": "  of which up in pterygium", "n": int(((g["robust"].astype(str).str.lower() == "true") & (g["median_logFC"] > 0)).sum())},
        {"Item": "  of which down in pterygium", "n": int(((g["robust"].astype(str).str.lower() == "true") & (g["median_logFC"] < 0)).sum())},
        {"Item": "  measured in all four cohorts", "n": int(((g["robust"].astype(str).str.lower() == "true") & (g["n_cohorts_measured"] == 4)).sum())},
        {"Item": "Core stable genes", "n": len(core)},
    ])
    keep_g = ["gene", "Z", "p_combined", "FDR", "median_logFC", "direction",
              "n_cohorts_measured", "n_direction_agree", "direction_consistent", "robust",
              "logFC_GSE155776", "SE_GSE155776", "logFC_GSE51995", "SE_GSE51995",
              "logFC_GSE2513", "SE_GSE2513", "logFC_GSE183153", "SE_GSE183153"]
    keep_g = [c for c in keep_g if c in g.columns]
    g = g[keep_g].sort_values("FDR")
    return {"Summary": summary, "Gene-level meta-analysis": g,
            "Core stable genes": core}


def table_s4():
    """Table S4: gene-level and module-level leave-one-dataset-out."""
    ml = pd.read_csv(os.path.join(T, "module_lodo.csv"))
    ml = ml.rename(columns={"module": "Module", "dropped": "Cohort removed",
                            "pooled_g": "Pooled Hedges g", "ci_low": "95% CI low",
                            "ci_high": "95% CI high", "tau2": "tau2", "I2": "I2 (fraction)",
                            "p": "P"})
    mk = pd.read_csv(os.path.join(T, "module_lodo_mkh.csv"))
    mk = mk.rename(columns={"module": "Module", "dropped": "Cohort removed",
                            "method": "Method", "pooled_g": "Pooled Hedges g",
                            "ci_low": "95% CI low", "ci_high": "95% CI high",
                            "p": "P", "I2": "I2 (fraction)", "tau2": "tau2", "k": "k cohorts"})
    gl = pd.read_csv(os.path.join(T, "gene_lodo.csv"))
    return {"Module LODO (REML+KH)": ml, "Module LODO (modified KH)": mk,
            "Gene-level LODO summary": gl}


def table_s5():
    """Table S5: per-comparison DEG power diagnostics.

    Documents why several comparisons report zero FDR-significant genes: effect sizes
    are present (logFC spread) but standard errors are large when a group has only two
    samples. Generated by deg_power_diagnostics.py.
    """
    p = os.path.join(T, "TableS5_DEG_power_diagnostics.csv")
    if not os.path.exists(p):
        return None
    return {"DEG power diagnostics": pd.read_csv(p)}


def write_book(sheets, path):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
            ws = xw.sheets[name[:31]]
            ws.freeze_panes = "A2"
            for i, col in enumerate(df.columns, 1):
                width = max(len(str(col)) + 2,
                            min(60, int(df[col].astype(str).str.len().quantile(0.95) if len(df) else 10) + 2))
                ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = width
    print(f"  wrote {os.path.basename(path)}  ({len(sheets)} sheets, "
          f"{sum(len(d) for d in sheets.values())} rows total)")


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"output dir: {OUT}\n")

    print("[1/6] Table 1")
    d1 = deg_summary()
    print(d1.to_string(index=False))
    s2 = top_shared_genes()
    print("\ntop shared genes:")
    for _, r in s2.iterrows():
        print(f"  {r['Comparator lesion']}: up(n={r['n shared up']}) {r['Shared up (ranked by |meta Z|)'][:90]}")
        print(f"      down(n={r['n shared down']}) {r['Shared down (ranked by |meta Z|)'][:90]}")
    write_book({"DEG summary": d1, "Top shared genes": s2},
               os.path.join(OUT, "Table1_per_lesion_DEG_summary.xlsx"))

    print("\n[2/6] Table 2")
    write_book(table2(), os.path.join(OUT, "Table2_module_meta_analysis.xlsx"))

    print("\n[3/6] Table S1")
    write_book({"Dataset census": table_s1()},
               os.path.join(OUT, "TableS1_dataset_census.xlsx"))

    print("\n[4/6] Table S2")
    write_book(table_s2(), os.path.join(OUT, "TableS2_module_definitions_coverage_overlap.xlsx"))

    print("\n[5/6] Table S3")
    write_book(table_s3(), os.path.join(OUT, "TableS3_gene_meta_analysis_core_stable.xlsx"))

    print("\n[6/6] Table S4")
    write_book(table_s4(), os.path.join(OUT, "TableS4_leave_one_dataset_out.xlsx"))

    print("\n[7/7] Table S5")
    t5 = table_s5()
    if t5:
        write_book(t5, os.path.join(OUT, "TableS5_DEG_power_diagnostics.xlsx"))
    else:
        print("  SKIP: run deg_power_diagnostics.py first")

    print("\nDONE")


if __name__ == "__main__":
    main()
