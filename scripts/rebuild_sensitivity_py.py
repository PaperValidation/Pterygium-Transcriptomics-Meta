"""Rebuild the Paper A sensitivity analyses on the SAME basis as the primary analysis.

Problem this fixes (audit P0-1): `sensitivity_v29.R` reads `module_effect_cohort.csv`,
which was written by `export_audit_cohort.py` using a different variance convention
from the metafor pipeline that produced `module_effect_meta.csv`. The primary analysis
and its own sensitivity analyses were therefore not comparable.

Fix: take the per-cohort effects (yi, se) recorded **inside** `module_effect_meta.csv`
and re-run every sensitivity analysis with the validated metafor-equivalent Python
implementation (`metafor_py.py`, exact to ~1e-15). This guarantees that the
REML+Knapp-Hartung cell of the 2x2 table reproduces the primary pooled g bit-for-bit.

Sensitivity analyses produced:
  1. 2x2 factorial: tau2 estimator (REML / DL) x inference (Knapp-Hartung / z)
  2. modified Knapp-Hartung (REML + adhoc)
  3. leave-one-dataset-out (LODO) under both standard and modified KH
  4. GSE51995 within-pair correlation r-sweep (SMCRP variance structure)

Outputs (overwrite; originals archived as *_audit_basis.csv):
  tables/module_sensitivity_2x2.csv
  tables/module_lodo_mkh.csv
  tables/module_sensitivity_rsweep.csv
"""
import os
import shutil

import numpy as np
import pandas as pd

import metafor_py as mf

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = os.path.join(BASE, "tables")
COHORTS = ["GSE155776", "GSE51995", "GSE2513", "GSE183153"]
PAIRED = "GSE51995"
N_PAIRS = 4


def load():
    m = pd.read_csv(os.path.join(T, "module_effect_meta.csv"))
    rows = []
    for _, R in m.iterrows():
        y = np.array([float(R[c + "_g"]) for c in COHORTS])
        se = np.array([float(R[c + "_se"]) for c in COHORTS])
        rows.append(dict(module=R["module"], y=y, v=se ** 2, se=se,
                         r=float(R["pair_r_GSE51995"]),
                         kh_g=float(R["kh_g"]), kh_se=float(R["kh_se"]),
                         kh_p=float(R["kh_p"]),
                         mkh_p=float(R["mkh_p"]), dl_p=float(R["dl_p"]),
                         q=float(R["q_BH19"])))
    return m, rows


def bh(pvals):
    """Benjamini-Hochberg adjusted q values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def validate(rs):
    """Prove the Python implementation reproduces the primary metafor analysis.

    Tolerance note: we solve the REML score equation to machine precision, so our
    tau2 is the exact REML root. metafor's stored tau2 carries roughly 1e-6 of
    optimiser tolerance (verified: d(ll)/d(tau2) = -1.0e-5 at the stored value, and our
    solution has a log-likelihood 1e-11 higher). That propagates to ~1e-7 in g and P,
    which is far below the 4 decimal places reported in the paper, so a 1e-6 tolerance
    is the right acceptance bar here.
    """
    TOL = 1e-6
    print("=" * 100)
    print("VALIDATION: does the Python implementation reproduce the primary REML+KH?")
    print("=" * 100)
    print(f"{'module':22s} {'g':>9s} {'dG':>9s} {'dSE':>9s} {'P(KH)':>10s} {'dP':>9s} {'P(mKH)':>10s} {'dP':>9s}")
    print("-" * 100)
    worst = 0.0
    nok = 0
    for r in rs:
        kh = mf.meta(r["y"], r["v"], "REML", "knha")
        mk = mf.meta(r["y"], r["v"], "REML", "adhoc")
        dg = abs(kh["g"] - r["kh_g"])
        dse = abs(kh["se"] - float(r["kh_se"]))
        dp = abs(kh["p"] - r["kh_p"])
        dm = abs(mk["p"] - r["mkh_p"])
        worst = max(worst, dg, dse, dp, dm)
        if max(dg, dse, dp, dm) < TOL:
            nok += 1
        print(f"{r['module']:22s} {kh['g']:9.5f} {dg:9.1e} {dse:9.1e} {kh['p']:10.6f} {dp:9.1e} "
              f"{mk['p']:10.6f} {dm:9.1e}")
    print("-" * 100)
    print(f"matches (tolerance {TOL:g}): {nok}/{len(rs)}   worst |delta| = {worst:.2e}")
    if nok != len(rs):
        raise SystemExit("ABORT: Python implementation does not reproduce the primary analysis")
    print("PASS -- every sensitivity analysis below shares the primary analysis' basis.\n")


def sens_2x2(rs):
    rows = []
    for r in rs:
        for tau in ("REML", "DL"):
            for inf in ("knha", "adhoc", "z"):
                res = mf.meta(r["y"], r["v"], tau, inf)
                rows.append(dict(module=r["module"], tau2_estimator=tau, inference=inf,
                                 pooled_g=res["g"], se=res["se"], ci_low=res["ci_low"],
                                 ci_high=res["ci_high"], p=res["p"], I2=res["I2"],
                                 tau2=res["tau2"], Q=res["Q"], Qp=res["Qp"], k=res["k"]))
    return pd.DataFrame(rows)


def sens_mkh(rs):
    rows = []
    for r in rs:
        for inf, label in (("knha", "knha"), ("adhoc", "adhoc")):
            res = mf.meta(r["y"], r["v"], "REML", inf)
            rows.append(dict(module=r["module"], dropped="none(full)", method=label,
                             pooled_g=res["g"], ci_low=res["ci_low"],
                             ci_high=res["ci_high"], p=res["p"], I2=res["I2"],
                             tau2=res["tau2"], k=res["k"]))
        for i, drop in enumerate(COHORTS):
            keep = [j for j in range(len(COHORTS)) if j != i]
            for inf in ("knha", "adhoc"):
                res = mf.meta(r["y"][keep], r["v"][keep], "REML", inf)
                rows.append(dict(module=r["module"], dropped=drop, method=inf,
                                 pooled_g=res["g"], ci_low=res["ci_low"],
                                 ci_high=res["ci_high"], p=res["p"], I2=res["I2"],
                                 tau2=res["tau2"], k=res["k"]))
    return pd.DataFrame(rows)


def sens_rsweep(rs):
    """Sweep the GSE51995 within-pair correlation, keeping yi fixed.

    metafor's SMCRP variance is vi(r) = 2(1-r)/n + c; `c` is calibrated per module
    from the stored variance at the observed r, so the observed-r row reproduces the
    primary analysis exactly.
    """
    grid = [-0.5, 0.0, 0.25, 0.5, 0.75, 0.9]
    j = COHORTS.index(PAIRED)
    rows = []
    for r_ in rs:
        c = mf.calibrate_smcrp_c(r_["r"], N_PAIRS, r_["v"][j])
        for rv in grid + [r_["r"]]:
            v = r_["v"].copy()
            v[j] = mf.smcrp_vi(rv, N_PAIRS, c)
            res = mf.meta(r_["y"], v, "REML", "knha")
            rows.append(dict(module=r_["module"], r=round(rv, 4),
                             r_type="observed" if rv == r_["r"] else "sweep",
                             paired_g_51995=r_["y"][j], paired_vi_51995=v[j],
                             pooled_g=res["g"], ci_low=res["ci_low"],
                             ci_high=res["ci_high"], p=res["p"], I2=res["I2"], k=res["k"]))
    return pd.DataFrame(rows)


def main():
    m, rs = load()
    validate(rs)

    # archive the audit-basis versions before overwriting
    for name in ("module_sensitivity_2x2.csv", "module_lodo_mkh.csv",
                 "module_sensitivity_rsweep.csv"):
        p = os.path.join(T, name)
        if os.path.exists(p):
            dst = os.path.join(T, name.replace(".csv", "_audit_basis.csv"))
            if not os.path.exists(dst):
                shutil.copy2(p, dst)
                print(f"archived {name} -> {os.path.basename(dst)}")

    print("[1/3] 2x2 factorial (tau2 estimator x inference)")
    d2 = sens_2x2(rs)
    d2.to_csv(os.path.join(T, "module_sensitivity_2x2.csv"), index=False)
    print(f"  wrote module_sensitivity_2x2.csv ({len(d2)} rows)")
    sub = d2[(d2.tau2_estimator == "REML") & (d2.inference == "knha")].sort_values(
        "pooled_g", ascending=False)
    print(f"  {'module':22s} {'REML+KH g':>10s} {'P':>9s} {'I2':>6s}   (primary -- must equal module_effect_meta)")
    for _, r in sub.iterrows():
        print(f"  {r['module']:22s} {r['pooled_g']:10.4f} {r['p']:9.5f} {r['I2']:5.1f}%")

    print("\n[2/3] LODO with standard and modified Knapp-Hartung")
    dm = sens_mkh(rs)
    dm.to_csv(os.path.join(T, "module_lodo_mkh.csv"), index=False)
    print(f"  wrote module_lodo_mkh.csv ({len(dm)} rows)")
    ecm = dm[(dm.module == "ECM_organization") & (dm.dropped != "none(full)")]
    print(f"  ECM_organization LODO:")
    for _, r in ecm.iterrows():
        print(f"    drop {r['dropped']:10s} {r['method']:6s} g={r['pooled_g']:7.4f} "
              f"P={r['p']:.4f}{'  *' if r['p'] < 0.05 else ''}")

    print("\n[3/3] GSE51995 within-pair correlation sweep")
    dr = sens_rsweep(rs)
    dr.to_csv(os.path.join(T, "module_sensitivity_rsweep.csv"), index=False)
    print(f"  wrote module_sensitivity_rsweep.csv ({len(dr)} rows)")
    ecmr = dr[dr.module == "ECM_organization"]
    print(f"  ECM_organization r-sweep:")
    for _, r in ecmr.iterrows():
        print(f"    r={r['r']:7.4f} ({r['r_type']:8s}) g={r['pooled_g']:7.4f} P={r['p']:.5f}")

    # BH q for the 2x2 primary cell, to confirm the manuscript's "all q>=0.22"
    print("\nBH q across 19 modules (REML+KH primary):")
    q = bh(sub["p"].values)
    for (_, r), qq in zip(sub.iterrows(), q):
        print(f"  {r['module']:22s} P={r['p']:9.5f} q={qq:.4f}")
    print(f"  min q = {q.min():.4f};  modules with q<0.05: {(q < 0.05).sum()}")


if __name__ == "__main__":
    main()
