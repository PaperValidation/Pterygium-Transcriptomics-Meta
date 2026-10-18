# Cross-cohort transcriptomic meta-analysis of pterygium

Analysis code and derived data for:

> **Cross-cohort transcriptomic meta-analysis supports a reproducible ECM-remodeling
> signal in pterygium** (manuscript under review)

This repository reproduces every statistic, table and figure reported in the paper.
No primary patient data were generated: all expression data are public GEO series
listed in `tables/Table_S1_dataset_census.csv`.

---

## What is here

```
scripts/
  metafor_py.py                  metafor-equivalent random-effects meta-analysis (validated)
  rebuild_sensitivity_py.py      regenerates the 2x2 / LODO / r-sweep sensitivity analyses
  build_paperA_tables.py         assembles Table 1, Table 2, Table S1-S5 (XLSX)
  deg_power_diagnostics.py       Table S5: why some comparisons yield no FDR-significant DEGs
  P1_figures_effectsize_v28.py   Figure 4 (forest plots) and Figure 5 (effect sizes)
  P1_figS2_composition.py        Figure S2 (GSE217707 single-cell composition)
  effect_size_meta_v28.R         original R/metafor implementation (reference)
tables/                          module scores, effect sizes, and all derived tables
figures/, figures_supplement/    output directories
```

## Requirements

Python 3.11+ with `numpy`, `scipy`, `pandas`, `openpyxl`, `matplotlib`, `Pillow`.
The R script additionally needs `metafor` and `jsonlite`, but **is not required** —
`metafor_py.py` reproduces its output (see validation below).

```bash
pip install -r requirements.txt
```

## Reproduce

```bash
# 1. Rebuild the sensitivity analyses on the primary analysis' basis (~2 s)
python scripts/rebuild_sensitivity_py.py

# 2. Build all deliverable tables (Table 1, 2, S1-S5)
python scripts/build_paperA_tables.py

# 3. Rebuild Figures 4 and 5
python scripts/P1_figures_effectsize_v28.py

# 4. Rebuild Figure S2
python scripts/P1_figS2_composition.py

# 5. (optional) Recompute the DEG power diagnostics behind Table S5,
#    needs the per-dataset DEG matrices from the preprocessing pipeline
python scripts/deg_power_diagnostics.py
```

All scripts locate their inputs relative to the repository root; set `PAPER_A_DIR`
to run them against a copy of `tables/` stored elsewhere.

---

## Statistical model

Primary inference is a random-effects meta-analysis of per-module Hedges' *g* across
four pterygium cohorts (GSE155776, GSE51995, GSE2513, GSE183153; 24 pterygium vs
18 control samples), with:

* **REML** estimation of the between-study variance τ²,
* **Knapp–Hartung** adjusted inference as the primary analysis,
* **modified Knapp–Hartung** (`adhoc`) and **DerSimonian–Laird** as sensitivity analyses,
* **Benjamini–Hochberg** correction across the 19 module-level tests.

For GSE51995 (paired same-eye design) the effect size is a paired standardized mean
difference using the raw-score SD denominator with the reported within-pair correlation.

### Headline results

| Module | Pooled *g* [95% CI] | I² | *P* (KH) | *q* (BH, 19) |
|---|---|---|---|---|
| ECM organization | 0.82 [0.35, 1.29] | 0% | 0.012 | 0.224 |
| EMT | 1.22 [0.07, 2.36] | 18% | 0.043 | 0.271 |
| senescence / SASP | 0.51 [0.06, 0.97] | 0% | 0.037 | 0.271 |

Exactly three modules reach nominal *P* < 0.05; **none** survives Benjamini–Hochberg
correction across the 19 module tests (minimum *q* = 0.224). Immune-evasion
(*g* = −0.20, *P* = 0.46) and genomic-instability (*g* = +0.20, *P* = 0.72) modules are
not reproducibly elevated.

At the gene level, signed Stouffer *Z* meta-analysis identifies **1,691**
direction-consistent meta-associated genes (844 up / 847 down), of which **1,358** were
measured in all four cohorts and **1,107** meet the objective core-stability criteria.

---

## Validation

`metafor_py.py` reimplements the pieces of `metafor` that this analysis uses, so the
pipeline can run without an R installation. `rebuild_sensitivity_py.py` asserts on
every run that it reproduces the primary analysis before proceeding:

```
matches (tolerance 1e-06): 19/19   worst |delta| = 2.8e-07
```

Agreement with the R `metafor` 5.0-1 output shipped in `tables/module_effect_meta.csv`,
across all 19 modules:

| Quantity | Max difference |
|---|---|
| pooled *g* | 2.8 × 10⁻⁷ |
| SE (Knapp–Hartung) | 1.4 × 10⁻⁷ |
| *P* (Knapp–Hartung) | 1.0 × 10⁻⁷ |
| *P* (modified KH) | 2.5 × 10⁻⁷ |
| τ² | ~10⁻⁶ |
| I² | < 10⁻¹³ percentage point |

The residual differences come from the optimiser tolerance in the stored R values, not
from the model: solving the REML score equation to machine precision places the optimum
~2 × 10⁻⁶ away from the value `metafor` wrote out, with a log-likelihood 10⁻¹¹ higher.
`I²` follows metafor's definition I² = τ²/(τ² + ṽ) with ṽ the typical sampling variance,
not the Cochran-Q form.

Two implementation details matter and are documented in the code:

1. **Knapp–Hartung uses the random-effects residual Q<sub>E</sub>**, not the fixed-effect
   Cochran Q. Using the latter gives visibly wrong SEs for low-heterogeneity modules.
2. **`effect_size_meta_v28.R` had a `cbind()` column-collision** on the UV-response
   modules: `module_scores_GSE2513.csv` and `module_scores_GSE183153.csv` already carry
   `UV_response_UP/DN` (rebuilt 2026-08-12), with the same names but different values
   from the older `module_scores_UVsplit_*.csv` (2026-08-11). The script now takes the
   newer BASE columns explicitly and asserts there are no duplicate names. This affects
   only the two UV modules and not the three headline results.

---

## Notes on the data

* **GSE83627 is excluded** from every analysis: it is a byte-level re-deposit of
  GSE51995 (expression correlation r = 1.0 across all eight sample pairs).
* Per-dataset differential expression was computed with limma at FDR < 0.05 and
  |log₂FC| > 0.585 from the public series matrices; `tables/deg_summary_counts.csv`
  holds the resulting counts. The per-dataset DEG matrices themselves are upstream
  products of the GEO preprocessing pipeline and are not redistributed here.
* **Several comparisons report zero FDR-significant genes** (GSE183153 4 vs 2,
  GSE191114 7 vs 2). This is a power limitation, not a processing failure, and
  `tables/TableS5_DEG_power_diagnostics.csv` documents it: in GSE183153 the spread of
  effect sizes is actually larger than in GSE2513 (8 vs 4; log₂FC SD 0.77 vs 0.33, with
  16% vs 1.8% of genes above |log₂FC| > 1), but the median standard error is ~2.9×
  higher (0.52 vs 0.18) because the control arm has only two samples, so nothing
  survives FDR correction. The genomic inflation factor is 1.05 and the estimated
  signal proportion (π1) is 1.6%, i.e. the test statistics behave as expected. This is
  also why cross-cohort inference in the paper rests on the module-level meta-analysis:
  pooling genes into GSVA scores averages away most of the per-gene sampling noise, and
  the ECM module still shows g = 0.46 in this cohort.
* Cohorts with two samples per group (GSE151872, GSE208384) and the pinguecula
  comparison are exploratory only; comparisons with n < 3 per group yield no
  FDR-significant genes, which is a power limitation rather than an absence of effect.
* GSE217707 is conjunctival melanoma, not pterygium. It is included only as a
  supplementary single-cell context figure and is **not** used to infer pterygium
  biology.

## Data availability

All expression data are public: GSE155776, GSE183153, GSE51995, GSE2513, GSE191114,
GSE147449, GSE148387, GSE217707 (accessions, platforms and PMIDs in
`tables/Table_S1_dataset_census.csv`). Code is in this repository.

## License

Code released under the MIT License (see `LICENSE`). Derived tables under CC BY 4.0.
