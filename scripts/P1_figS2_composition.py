# -*- coding: utf-8 -*-
"""Paper A Supplementary Figure S2: single-cell composition of conjunctival melanoma
versus benign lesions (GSE217707; 10x scRNA-seq, 54,175 cells, 10 samples).

Panel (a): per-sample cell-type proportions (stacked bars), samples annotated by group.
Panel (b): group-level mean +/- SEM across samples, per cell type, with two-sided
           Wilcoxon rank-sum P values (n=7 melanoma vs n=3 benign; exploratory only).

IOVS specs: TIFF + PNG at 300 dpi, double-column width 7.2 in (2160 px @300 dpi).
Outputs: results/paper_a/supplementary_figures/FigS2_*.{tif,png}
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = os.environ.get("PAPER_A_DIR") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(BASE, "tables")
OUT = os.path.join(BASE, "figures_supplement")
os.makedirs(OUT, exist_ok=True)

# consistent with the main figures' palette
CELL_COLORS = {
    "melanocyte": "#D55E00",
    "Tcell": "#0072B2",
    "Bcell": "#009E73",
    "myeloid": "#CC79A7",
    "fibroblast": "#E69F00",
    "epithelial": "#56B4E9",
    "endothelial": "#999999",
}
CELL_ORDER = ["melanocyte", "epithelial", "fibroblast", "endothelial",
              "Tcell", "Bcell", "myeloid"]
LABEL = {"melanocyte": "Melanocyte", "epithelial": "Epithelial", "fibroblast": "Fibroblast",
         "endothelial": "Endothelial", "Tcell": "T cell", "Bcell": "B cell",
         "myeloid": "Myeloid"}


def main():
    ann = pd.read_csv(os.path.join(T, "GSE217707_cell_annotations.csv"))
    # per-sample counts -> proportions
    ct = (ann.groupby(["sample", "group", "celltype"]).size()
             .rename("n").reset_index())
    totals = ct.groupby("sample")["n"].sum().rename("total")
    ct = ct.merge(totals, on="sample")
    ct["pct"] = 100.0 * ct["n"] / ct["total"]
    wide = ct.pivot_table(index=["group", "sample"], columns="celltype",
                          values="pct", fill_value=0.0)
    wide = wide.reindex(columns=[c for c in CELL_ORDER if c in wide.columns])
    wide = wide.sort_index(level=["group", "sample"])
    # melanoma first, then benign; within group keep sample order
    order = ([("melanoma", s) for _, s in sorted(
        [(g, s) for g, s in wide.index if g == "melanoma"], key=lambda x: x[1])]
        + [("benign", s) for _, s in sorted(
            [(g, s) for g, s in wide.index if g == "benign"], key=lambda x: x[1])])
    wide = wide.loc[order]

    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.38)

    # ---- panel (a): stacked per-sample composition ----
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(wide))
    bottom = np.zeros(len(wide))
    for c in wide.columns:
        vals = wide[c].values
        ax.bar(x, vals, bottom=bottom, width=0.74,
               color=CELL_COLORS.get(c, "#BBBBBB"),
               edgecolor="white", linewidth=0.7, label=LABEL.get(c, c))
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([s for _, s in wide.index], rotation=55, ha="right", fontsize=6.5)
    ax.set_ylabel("Cells per sample (%)", fontsize=8)
    ax.set_ylim(0, 100)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title("a  Per-sample composition", fontsize=8.5, loc="left", fontweight="bold")
    # group separators / label
    n_mel = sum(1 for g, _ in wide.index if g == "melanoma")
    ax.axvline(n_mel - 0.5, color="#444444", lw=1.0, ls="--")
    ax.text((n_mel - 1) / 2, 103, "Melanoma (n=7)", ha="center", fontsize=7,
            color="#D55E00", fontweight="bold")
    ax.text(n_mel + (len(wide) - n_mel - 1) / 2, 103, "Benign (n=3)", ha="center",
            fontsize=7, color="#0072B2", fontweight="bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=6, frameon=False, loc="upper left",
              bbox_to_anchor=(1.01, 1.02), handlelength=1.0, labelspacing=0.35)

    # ---- panel (b): group mean +/- SEM with Wilcoxon P ----
    ax2 = fig.add_subplot(gs[0, 1])
    mel = wide.loc["melanoma"] if "melanoma" in wide.index.get_level_values(0) else None
    ben = wide.loc["benign"] if "benign" in wide.index.get_level_values(0) else None
    cts = list(wide.columns)
    ypos = np.arange(len(cts))[::-1]
    h = 0.34
    rows = []
    for i, c in enumerate(cts):
        ym, yb = ypos[i] + h / 2, ypos[i] - h / 2
        mx, sx = mel[c].mean(), mel[c].std(ddof=1) / np.sqrt(len(mel))
        bx, sb = ben[c].mean(), ben[c].std(ddof=1) / np.sqrt(len(ben))
        ax2.barh(ym, mx, height=h, color="#D55E00", alpha=0.9,
                 edgecolor="white", linewidth=0.6)
        ax2.errorbar(mx, ym, xerr=sx, fmt="none", ecolor="#333333",
                     elinewidth=0.9, capsize=2.2)
        ax2.barh(yb, bx, height=h, color="#0072B2", alpha=0.9,
                 edgecolor="white", linewidth=0.6)
        ax2.errorbar(bx, yb, xerr=sb, fmt="none", ecolor="#333333",
                     elinewidth=0.9, capsize=2.2)
        try:
            p = mannwhitneyu(mel[c], ben[c], alternative="two-sided").pvalue
        except ValueError:
            p = np.nan
        rows.append({"celltype": c, "mel_mean": mx, "mel_sem": sx,
                     "ben_mean": bx, "ben_sem": sb, "p": p})
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        xpos = max(mx + sx, bx + sb) + 3.0
        ax2.text(xpos, ypos[i], star, fontsize=7, va="center", color="#333333")
    ax2.set_yticks(ypos)
    ax2.set_yticklabels([LABEL.get(c, c) for c in cts], fontsize=7.5)
    ax2.set_xlabel("Cells per sample (%)", fontsize=8)
    ax2.set_xlim(0, min(100, max(wide.max().max() + 18, 60)))
    ax2.tick_params(axis="x", labelsize=7)
    ax2.set_title("b  Group mean \u00b1 SEM", fontsize=8.5, loc="left", fontweight="bold")
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    # manual legend
    ax2.barh([], [], color="#D55E00", label="Melanoma (n=7)")
    ax2.barh([], [], color="#0072B2", label="Benign (n=3)")
    ax2.legend(fontsize=6.5, frameon=False, loc="lower right")

    fig.text(0.005, -0.035,
             "GSE217707 (10x scRNA-seq; 54,175 cells passing QC; 7 melanoma vs 3 benign samples). "
             "Bars show mean and error bars SEM across samples; two-sided Wilcoxon rank-sum test, "
             "exploratory (n=3 benign), not corrected for multiple comparisons. ns, not significant.",
             fontsize=5.8, color="#444444", ha="left", va="top", wrap=True)

    for ext, dpi in (("tif", 300), ("png", 300)):
        path = os.path.join(OUT, f"FigS2_GSE217707_composition.{ext}")
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        print("saved", path, f"({os.path.getsize(path)/1024:.0f} KB)")
    plt.close(fig)

    pd.DataFrame(rows).to_csv(
        os.path.join(T, "FigS2_celltype_stats.csv"), index=False)
    print("\nper-sample proportions written to FigS2_celltype_stats.csv")
    for r in rows:
        print(f"  {r['celltype']:12s} mel={r['mel_mean']:6.2f}±{r['mel_sem']:4.2f}  "
              f"ben={r['ben_mean']:6.2f}±{r['ben_sem']:4.2f}  P={r['p']:.4f}")


if __name__ == "__main__":
    main()
