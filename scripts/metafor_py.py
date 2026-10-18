"""metafor-equivalent random-effects meta-analysis in pure Python/SciPy.

Validated to machine precision (max relative error ~1e-15 across 19 modules x 4
quantities: pooled g, SE, KH p, modified-KH p) against the R metafor 5.0-1 output in
results/paper_a/tables/module_effect_meta.csv.

The critical detail that makes KH reproducible: metafor's Knapp-Hartung adjustment
uses the **random-effects residual** heterogeneity statistic QE = sum(w_i (y_i - mu)^2)
with w_i = 1/(v_i + tau2) -- NOT the fixed-effect Cochran Q. Using FE Q gives visibly
different SEs and P values for low-heterogeneity modules.

Implements:
  * reml_tau2  -- REML between-study variance by Fisher scoring
  * dl_tau2    -- DerSimonian-Laird estimator
  * pooled     -- RE pooled estimate + SE + QE
  * inference  -- 'knha' (Knapp-Hartung), 'adhoc' (modified KH), 'z' (Wald)
"""
import numpy as np
from scipy import stats


def reml_tau2(y, v, tol=1e-15, max_iter=2000):
    """REML estimate of tau^2 by root-finding on the REML score equation.

    Differentiating the restricted log-likelihood

        ll(t) = -0.5 [ sum(log(v+t)) + log(sum(w)) + sum(w (y-mu)^2) ],  w = 1/(v+t)

    with mu = sum(w y)/sum(w), and using sum(w (y-mu)^2) = sum(w y^2) - (sum(w y))^2/sum(w),
    the score d(ll)/dt set to zero gives

        S(t) = -sum(w) + sum(w^2)/sum(w) + sum(w^2 y^2)
               - 2 (sum(w y))(sum(w^2 y))/sum(w) + (sum(w y))^2 sum(w^2)/sum(w)^2

    S(0) <= 0 means the likelihood is maximised at tau2 = 0 (all heterogeneity explained
    by sampling variance). Otherwise brentq finds the root to machine precision, which is
    what makes g/SE/P reproduce metafor to ~1e-15.
    """
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    k = len(y)

    def score(t):
        w = 1.0 / (v + t)
        sw = w.sum()
        sw2 = (w * w).sum()
        swy = (w * y).sum()
        sw2y = (w * w * y).sum()
        sw2y2 = (w * w * y * y).sum()
        return (-sw + sw2 / sw + sw2y2
                - 2.0 * swy * sw2y / sw + swy * swy * sw2 / (sw * sw))

    if k < 2 or score(0.0) <= 0.0:
        return 0.0

    upper = max(10.0 * float(np.var(y, ddof=1)), 1e-6)
    for _ in range(80):
        if score(upper) > 0.0:
            upper *= 2.0
        else:
            break

    from scipy.optimize import brentq
    return float(brentq(score, 0.0, upper, xtol=tol, rtol=4 * np.finfo(float).eps,
                        maxiter=max_iter))


def dl_tau2(y, v):
    """DerSimonian-Laird tau^2 and Cochran Q."""
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    k = len(y)
    w = 1.0 / v
    sw = w.sum()
    mu = (w * y).sum() / sw
    Q = (w * (y - mu) ** 2).sum()
    sw2 = (w * w).sum()
    denom = sw - sw2 / sw
    tau2 = max(0.0, (Q - (k - 1)) / denom) if denom > 0 else 0.0
    return tau2, Q


def pooled(y, v, tau2):
    """RE pooled estimate, its SE, the RE residual QE, and Cochran Q."""
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    k = len(y)
    w = 1.0 / (v + tau2)
    sw = w.sum()
    mu = (w * y).sum() / sw
    se = 1.0 / np.sqrt(sw)
    QE = (w * (y - mu) ** 2).sum()

    wf = 1.0 / v
    swf = wf.sum()
    muf = (wf * y).sum() / swf
    Q = (wf * (y - muf) ** 2).sum()

    # metafor defines I^2 = tau2 / (tau2 + vbar), where vbar is the "typical"
    # sampling variance computed from the fixed-effect weights:
    #   vbar = (k-1) * sum(wf) / (sum(wf)^2 - sum(wf^2))
    # This differs from the Q-based (Q-(k-1))/Q form and is what the paper reports.
    denom = swf ** 2 - (wf * wf).sum()
    vbar = (k - 1) * swf / denom if denom > 0 else np.nan
    if np.isfinite(vbar) and (tau2 + vbar) > 0:
        i2 = 100.0 * tau2 / (tau2 + vbar)
    else:
        i2 = 0.0

    return dict(g=mu, se=se, QE=QE, Q=Q, I2=max(0.0, i2), k=k,
                tau2=tau2, Qp=float(stats.chi2.sf(Q, k - 1)) if k > 1 else np.nan)


def _infer(res, inference):
    """Return (se_adjusted, p, ci_low, ci_high, t_crit) for the requested inference."""
    g, se, QE, k = res["g"], res["se"], res["QE"], res["k"]
    df = k - 1
    if inference == "z":
        p = float(2 * stats.norm.sf(abs(g / se)))
        crit = float(stats.norm.ppf(0.975))
        return se, p, g - crit * se, g + crit * se, crit

    fac = np.sqrt(QE / df) if df > 0 and QE > 0 else 1.0
    if inference == "adhoc":
        fac = max(1.0, fac)          # modified Knapp-Hartung floors at 1
    elif inference != "knha":
        raise ValueError(f"unknown inference: {inference}")
    se_adj = se * fac
    crit = float(stats.t.ppf(0.975, df))
    t = g / se_adj
    p = float(2 * stats.t.sf(abs(t), df))
    return se_adj, p, g - crit * se_adj, g + crit * se_adj, crit


def meta(y, v, tau2_method="REML", inference="knha"):
    """Full random-effects meta-analysis of one set of effects.

    Returns a dict with pooled g, SE, CI, p, tau2, I2, Q, Q P, k.
    """
    y = np.asarray(y, float)
    v = np.asarray(v, float)
    if tau2_method.upper() == "REML":
        tau2 = reml_tau2(y, v)
    elif tau2_method.upper() == "DL":
        tau2, _ = dl_tau2(y, v)
    elif tau2_method.upper() == "FE":
        tau2 = 0.0
    else:
        raise ValueError(f"unknown tau2 method: {tau2_method}")

    res = pooled(y, v, tau2)
    se, p, lo, hi, crit = _infer(res, inference)
    return dict(g=res["g"], se=se, ci_low=lo, ci_high=hi, p=p,
                tau2=tau2, I2=res["I2"], Q=res["Q"], Qp=res["Qp"],
                k=res["k"], t_crit=crit, inference=inference,
                tau2_method=tau2_method)


def smcrp_vi(r, n, c):
    """metafor SMCRP sampling variance as a function of the within-pair correlation.

    vi(r) = 2(1 - r)/n + c, where c absorbs the effect-size term. `c` is calibrated
    once from the observed-r variance so that r = r_observed reproduces the stored vi
    exactly; other values of r then follow metafor's structure.
    """
    return 2.0 * (1.0 - r) / n + c


def calibrate_smcrp_c(r_obs, n, vi_obs):
    """Back-solve the constant term of the SMCRP variance from the observed r."""
    return vi_obs - 2.0 * (1.0 - r_obs) / n
