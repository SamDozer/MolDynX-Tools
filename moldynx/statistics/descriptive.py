"""Descriptive statistics and bootstrap confidence intervals."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats as sps


def describe(x: Iterable[float], name: str = "value") -> dict:
    """Mean/median/std/SEM/95% t-CI/min/max for a 1-D series (NaNs dropped)."""
    x = np.asarray(list(x), dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    mean = float(np.mean(x)) if n else np.nan
    std = float(np.std(x, ddof=1)) if n > 1 else 0.0
    sem = std / np.sqrt(n) if n > 1 else 0.0
    if n > 1:
        t = sps.t.ppf(0.975, df=n - 1)
        lo, hi = mean - t * sem, mean + t * sem
    else:
        lo = hi = mean
    return {
        "metric": name, "n": n, "mean": mean,
        "median": float(np.median(x)) if n else np.nan,
        "std": std, "sem": sem,
        "ci95_low": float(lo), "ci95_high": float(hi),
        "min": float(np.min(x)) if n else np.nan,
        "max": float(np.max(x)) if n else np.nan,
    }


def bootstrap_ci(x: Iterable[float], statistic=np.mean, n_boot: int = 10000,
                 ci: float = 95.0, seed: int | None = 0) -> dict:
    """
    Bootstrap confidence interval for a statistic of a (correlated) series.

    Uses the percentile method. For MD time-series the samples are correlated, so
    this complements the t-based CI in :func:`describe`.
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(list(x), dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 2:
        v = float(statistic(x)) if n else np.nan
        return {"estimate": v, "ci_low": v, "ci_high": v, "n_boot": 0}
    idx = rng.integers(0, n, size=(n_boot, n))
    boot = np.array([statistic(x[i]) for i in idx])
    a = (100 - ci) / 2
    return {
        "estimate": float(statistic(x)),
        "ci_low": float(np.percentile(boot, a)),
        "ci_high": float(np.percentile(boot, 100 - a)),
        "n_boot": n_boot,
    }


def summary_frame(series_map: dict[str, Iterable[float]]) -> pd.DataFrame:
    """Tidy summary-statistics DataFrame from {name: series}."""
    return pd.DataFrame([describe(v, name=k) for k, v in series_map.items()]
                        ).set_index("metric")
