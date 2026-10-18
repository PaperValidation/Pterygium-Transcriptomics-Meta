"""Differential-expression power diagnostics for each comparison (Table S5).

Motivation: in Table 1 several comparisons report zero FDR-significant genes, which
readers can mistake for a processing failure. This script quantifies *why*: for every
comparison it reports the group sizes, the number of genes tested, the observed versus
null-expected fractions of small p-values, the genomic inflation factor, the estimated
proportion of genes carrying real signal (Storey pi1), the spread of effect sizes, and
the median standard error.

The diagnostic distinguishes the two possible explanations:

  * effect sizes are small           -> logFC sd and |log2FC|>1 fraction are low
  * effect sizes are fine but their
    standard errors are large        -> logFC sd is normal/high while median SE is high

For GSE183153 (4 vs 2) the second pattern holds: logFC sd 0.77 with 16% of genes above
|log2FC|>1 is larger than in GSE2513 (8 vs 4, sd 0.33) or GSE51995 (4 vs 4, sd 0.50),
yet its median SE is 2.5-3x higher because the control arm has only two samples. The
analysis is therefore underpowered rather than technically broken, which is also why
the module-level meta-analysis (the paper's primary inference) still detects an effect
in this cohort: pooling genes into GSVA module scores averages away most of that
per-gene sampling noise.

Output: tables/TableS5_DEG_power_diagnostics.csv
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(BASE, "tables")

# The per-dataset DEG matrices are upstream products of the GEO preprocessing pipeline
# and are not redistributed with the code. Look for them in the usual analysis-tree
# locations; if absent, the frozen output CSV is reused instead of recomputing.
DEG_REBUILT_CANDIDATES = [
    os.environ.get("DEG_REBUILT_DIR") or "",
    os.path.join(os.path.dirname(BASE), "paper_b", "deg_rebuilt"),
    os.path.join(BASE, "..", "..", "results", "paper_b", "deg_rebuilt"),
]
DEG_REBUILT = next(
    (p for p in DEG_REBUILT_CANDIDATES if p and os.path.isdir(os.path.normpath(p))),
    os.path.normpath(DEG_REBUILT_CANDIDATES[1]),
)

FDR_CUT = 0.05
LFC_CUT = 0.585

SPECS = [
    ("GSE155776", "pterygium vs healthy conjunctiva", 8, 8,
     os.path.join(DEG_REBUILT, "GSE155776_DEG_pterygium_vs_normal.csv"), "Pterygium meta cohort"),
    ("GSE183153", "pterygium vs healthy conjunctiva", 4, 2,
     os.path.join(DEG_REBUILT, "GSE183153_DEG_pterygium_vs_normal.csv"), "Pterygium meta cohort"),
    ("GSE183153", "pinguecula vs healthy conjunctiva", 2, 2,
     os.path.join(DEG_REBUILT, "GSE183153_DEG_pinguecula_vs_normal.csv"), "Exploratory"),
    ("GSE51995", "pterygium vs conjunctiva (paired)", 4, 4,
     os.path.join(DEG_REBUILT, "GSE51995_DEG_pterygium_vs_normal_conjunctiva.csv"), "Pterygium meta cohort"),
    ("GSE2513", "pterygium vs conjunctiva", 8, 4,
     os.path.join(DEG_REBUILT, "GSE2513_DEG_pterygium_vs_normal_conjunctiva.csv"), "Pterygium meta cohort"),
    ("GSE191114", "active vs atrophic pterygium", 7, 2,
     os.path.join(DEG_REBUILT, "GSE191114_DEG_active_vs_atrophic.csv"), "Activity-only"),
    ("GSE208384", "pterygium vs conjunctiva", 2, 2,
     os.path.join(DEG_REBUILT, "GSE208384_DEG_pterygium_vs_normal.csv"), "Exploratory"),
    ("GSE147449", "SCC vs healthy conjunctiva", 7, 10,
     os.path.join(T, "GSE147449_DEG_b07_scc_vs_normal_conjunctiva.csv"), "OSSN comparator"),
    ("GSE147449", "papilloma vs healthy conjunctiva", 7, 10,
     os.path.join(T, "GSE147449_DEG_b07_papilloma_vs_normal_conjunctiva.csv"), "Papilloma comparator"),
    ("GSE148387", "melanoma vs healthy conjunctiva", 12, 8,
     os.path.join(T, "GSE148387_DEG_b07_melanoma_vs_normal_conjunctiva.csv"), "Melanoma comparator"),
]


def diagnose(path):
    d = pd.read_csv(path)
    for c in ("pvalue", "FDR", "logFC", "SE"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["pvalue", "logFC"])
    pv = d["pvalue"].to_numpy(float)

    chi2 = stats.chi2.isf(pv, 1)
    lam = float(np.median(chi2) / stats.chi2.ppf(0.5, 1))
    # Storey pi0 at lambda = 0.5
    pi0 = float(min(1.0, (pv > 0.5).mean() / 0.5))

    frac = {t: float((pv < t).mean()) for t in (0.05, 0.01, 0.001)}
    return dict(
        genes_tested=len(d),
        frac_p_lt_0_05=frac[0.05],
        ratio_p_lt_0_05=frac[0.05] / 0.05,
        frac_p_lt_0_001=frac[0.001],
        ratio_p_lt_0_001=frac[0.001] / 0.001,
        genomic_inflation_lambda=lam,
        pi1_signal_proportion=1 - pi0,
        logFC_sd=float(d["logFC"].std()),
        pct_abs_logFC_gt_1=float((d["logFC"].abs() > 1).mean() * 100),
        median_SE=float(d["SE"].median()) if "SE" in d.columns else np.nan,
        deg_fdr_0_05=float((d["FDR"] < FDR_CUT).sum()),
        deg_fdr_and_lfc=float(((d["FDR"] < FDR_CUT) & (d["logFC"].abs() > LFC_CUT)).sum()),
        min_FDR=float(d["FDR"].min()),
    )


def main():
    out_path = os.path.join(T, "TableS5_DEG_power_diagnostics.csv")

    available = [s for s in SPECS if os.path.exists(s[4])]
    if not available:
        if os.path.exists(out_path):
            print("DEG matrices not present in this checkout (they are not redistributed);")
            print(f"keeping the frozen table already in the repository: {os.path.basename(out_path)}")
            print(pd.read_csv(out_path).to_string(index=False))
            return
        raise SystemExit(
            "No DEG matrices found and no frozen table available. Set DEG_REBUILT_DIR to the "
            "directory containing *_DEG_*.csv files from the preprocessing pipeline."
        )

    rows = []
    for acc, comp, n1, n2, path, role in SPECS:
        if not os.path.exists(path):
            print(f"  SKIP (missing): {path}")
            continue
        r = diagnose(path)
        rows.append(dict(Dataset=acc, Comparison=comp,
                         n_disease=n1, n_control=n2, Role=role, **r))
    out = pd.DataFrame(rows)
    # order columns for readability
    out = out[["Dataset", "Comparison", "n_disease", "n_control", "genes_tested",
               "deg_fdr_0_05", "deg_fdr_and_lfc", "min_FDR",
               "frac_p_lt_0_05", "ratio_p_lt_0_05", "frac_p_lt_0_001",
               "ratio_p_lt_0_001", "genomic_inflation_lambda",
               "pi1_signal_proportion", "logFC_sd", "pct_abs_logFC_gt_1",
               "median_SE", "Role"]]
    out.columns = ["Dataset", "Comparison", "n disease", "n control", "Genes tested",
                   "DEGs FDR<0.05", "DEGs FDR<0.05 & |log2FC|>0.585", "min FDR",
                   "Fraction p<0.05", "Ratio vs null (p<0.05)",
                   "Fraction p<0.001", "Ratio vs null (p<0.001)",
                   "Genomic inflation lambda", "Signal proportion (pi1)",
                   "logFC sd", "% genes |log2FC|>1", "Median SE", "Role"]

    out.to_csv(out_path, index=False)
    print(f"wrote {os.path.basename(out_path)} ({len(out)} rows)\n")
    cols = ["Dataset", "n disease", "n control", "DEGs FDR<0.05", "Ratio vs null (p<0.001)",
            "Genomic inflation lambda", "Signal proportion (pi1)", "logFC sd",
            "% genes |log2FC|>1", "Median SE"]
    print(out[cols].to_string(index=False))
    print("\nInterpretation guide:")
    print("  Ratio vs null (p<0.001) ~1 and lambda ~1  -> analysis is null; signal proportion near 0")
    print("  logFC sd normal/high but Median SE high   -> underpowered, not technically broken")


if __name__ == "__main__":
    main()
