"""
Autocorrelation-aware statistics for MD time series.

* :func:`statistical_inefficiency` -- g = 1 + 2 Σ (1 − k/N) C(k), summed until the
  autocorrelation first drops to zero (Chodera et al., J. Chem. Theory Comput.
  2007, 3, 26; the estimator used by pymbar's ``timeseries``).
* :func:`describe_correlated` -- mean, SD, naive SEM and the corrected SEM with the
  number of effectively independent samples N/g.
* :func:`detect_equilibration` -- the start t₀ that maximises N_eff = (N − t₀)/g
  (Chodera, J. Chem. Theory Comput. 2016, 12, 1799).
* :func:`drift` -- linear trend with p-value, and a half-vs-half Welch test.
"""

from __future__ import annotations

import numpy as np


def statistical_inefficiency(x) -> float:
    x = np.asarray(x, float)
    n = len(x)
    if n < 4:
        return 1.0
    d = x - x.mean()
    var = float(np.mean(d * d))
    if var == 0:
        return 1.0
    g = 1.0
    for k in range(1, n - 1):
        c = float(np.sum(d[:n - k] * d[k:]) / ((n - k) * var))
        if c <= 0:
            break
        g += 2.0 * c * (1.0 - k / n)
    return max(1.0, g)


def describe_correlated(x) -> dict:
    x = np.asarray(x, float)
    n = len(x)
    if n == 0:
        return {"n": 0}
    sd = float(x.std(ddof=1)) if n > 1 else 0.0
    g = statistical_inefficiency(x)
    return {"n": n, "mean": float(x.mean()), "sd": sd,
            "sem_naive": sd / np.sqrt(n) if n else float("nan"),
            "stat_ineff": g, "n_eff": n / g, "sem": sd / np.sqrt(n / g) if n else float("nan")}


def detect_equilibration(x, step: int = 1, max_fraction: float = 0.8) -> dict:
    """Return ``{"t0": index, "g": g, "n_eff": N_eff}`` maximising N_eff = (N − t0)/g."""
    x = np.asarray(x, float)
    best = {"t0": 0, "g": 1.0, "n_eff": 0.0}
    for t0 in range(0, max(1, int(len(x) * max_fraction)), max(1, step)):
        g = statistical_inefficiency(x[t0:])
        ne = (len(x) - t0) / g
        if ne > best["n_eff"]:
            best = {"t0": t0, "g": g, "n_eff": ne}
    return best


def drift(t, x) -> dict:
    """Linear slope (per unit of t) with p-value, and a half-vs-half Welch t-test."""
    from scipy import stats
    t, x = np.asarray(t, float), np.asarray(x, float)
    ok = np.isfinite(t) & np.isfinite(x)
    t, x = t[ok], x[ok]
    if len(x) < 6 or np.ptp(t) == 0:
        return {"slope": None, "p_slope": None, "half_p": None}
    r = stats.linregress(t, x)
    h = len(x) // 2
    w = stats.ttest_ind(x[:h], x[h:], equal_var=False)
    return {"slope": float(r.slope), "p_slope": float(r.pvalue), "half_p": float(w.pvalue),
            "first_half_mean": float(x[:h].mean()), "second_half_mean": float(x[h:].mean())}
