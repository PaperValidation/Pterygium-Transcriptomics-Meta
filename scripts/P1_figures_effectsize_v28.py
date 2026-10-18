#!/usr/bin/env python3
"""Regenerate Fig4 (forest panel) and Fig5 (effect-size summary) with v2.8 numbers
(REML+KH primary, SMCRP paired, tau2/Q added to Fig4 annotations)."""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

BASE = os.environ.get("PAPER_A_DIR") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(BASE, "tables")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

df = pd.read_csv(os.path.join(T, "module_effect_meta.csv"))
df["module"] = df["module"].str.replace("_", " ")
COHORTS = ["GSE155776", "GSE51995", "GSE2513", "GSE183153"]
t3 = 3.182

def color_for(r):
    if r["q_BH19"] < 0.05:
        return "#B22222", "FDR q<0.05"
    if r["kh_p"] < 0.05:
        return "#C77D5E", "nominal P<0.05"
    return "#777777", "P>=0.05"

# ---------------------------------------------------------------- Fig 4
fig4_mods = ["ECM_organization", "focal_adhesion", "angiogenesis",
             "EMT", "immune_evasion", "genomic_instability"]

def forest_axes(ax, mod):
    row = df[df["module"] == mod.replace("_", " ")].iloc[0]
    labels, est, lo, hi = [], [], [], []
    for c in COHORTS:
        labels.append(c)
        g = row[f"{c}_g"]; se = row[f"{c}_se"]
        est.append(g); lo.append(g - 1.96*se); hi.append(g + 1.96*se)
    labels.append("Pooled (REML+KH)")
    est.append(row["kh_g"]); lo.append(row["kh_ci_low"]); hi.append(row["kh_ci_high"])
    y = np.arange(len(labels))[::-1]
    sig = row["kh_p"] < 0.05
    for i, (e, l, h) in enumerate(zip(est, lo, hi)):
        is_pool = (i == len(labels) - 1)
        col = "#111111" if is_pool else ("#B22222" if sig else "#555555")
        ax.plot([l, h], [y[i], y[i]], color=col, lw=1.6 if is_pool else 1.2, zorder=2)
        ax.scatter([e], [y[i]], color=col, s=34 if is_pool else 22, zorder=3)
    py = y[-1]; pe = est[-1]; pl = lo[-1]; ph = hi[-1]
    ax.add_patch(Polygon([(pl, py), (pe, py+0.18), (ph, py), (pe, py-0.18)],
                         fill=True, color="#111111", zorder=4))
    ax.axvline(0, color="#999999", lw=0.8, ls="--", zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(-3.5, 9.0)
    ax.tick_params(axis="x", labelsize=8)
    ax.set_xlabel("Hedges' g", fontsize=8.5)
    sub = (f"REML+KH g={row['kh_g']:.2f} [{row['kh_ci_low']:.2f},{row['kh_ci_high']:.2f}]\n"
           f"I\u00b2={row['kh_I2']*100:.0f}%, tau\u00b2={row['kh_tau2']:.3f}\n"
           f"P={row['kh_p']:.2e}, q={row['q_BH19']:.2f}")
    ax.set_title(mod.replace("_", " ") + ("  *" if sig else ""), fontsize=9.5, fontweight="bold")
    ax.text(0.99, 0.02, sub, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.0, color=("#B22222" if sig else "#555555"))
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.4))
for ax, mod in zip(axes.flat, fig4_mods):
    forest_axes(ax, mod)
fig.suptitle("Fig 4. Meta-analysis of per-cohort Hedges' g (19 modules; REML + Knapp\u2013Hartung)\n* nominal P<0.05; none reached BH q<0.05 across 19 modules",
             fontsize=10.5, fontweight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(FIG, "Fig4_module_forest_plots.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "Fig4_module_forest_plots.tif"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote Fig4")

# ---------------------------------------------------------------- Fig 5
d2 = df.sort_values("kh_g", ascending=True).reset_index(drop=True)
y = np.arange(len(d2))
fig, ax = plt.subplots(figsize=(7.6, 8.6))
for i, r in d2.iterrows():
    col, _ = color_for(r)
    lw = 2.0 if r["q_BH19"] < 0.05 else (1.8 if r["kh_p"] < 0.05 else 1.2)
    ms = 30 if r["q_BH19"] < 0.05 else (24 if r["kh_p"] < 0.05 else 15)
    ax.plot([r["kh_ci_low"], r["kh_ci_high"]], [y[i], y[i]], color=col, lw=lw, zorder=2)
    ax.scatter([r["kh_g"]], [y[i]], color=col, s=ms, zorder=3)
    ax.text(r["kh_ci_high"] + 0.08, y[i], f"I\u00b2={r['kh_I2']*100:.0f}%",
            va="center", fontsize=7, color="#444444")
ax.axvline(0, color="#999999", lw=0.9, ls="--", zorder=1)
ax.set_yticks(y)
ax.set_yticklabels(d2["module"], fontsize=8.6)
ax.set_xlabel("Pooled Hedges' g  (95% CI, REML + Knapp\u2013Hartung)", fontsize=9.5)
ax.set_xlim(-1.8, 7.4)
ax.set_title("Fig 5. Cross-cohort module effect sizes (19 modules, four-cohort meta-analysis)\n"
             "orange = nominal P<0.05 (none reached BH q<0.05); grey = P>=0.05; I\u00b2 at right",
             fontsize=10, fontweight="bold")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Fig5_module_effect_sizes.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "Fig5_module_effect_sizes.tif"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote Fig5")
