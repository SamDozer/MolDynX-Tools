"""Time-series helpers: moving/rolling/running averages, block averaging, plateau."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd


def moving_average(x: Iterable[float], window: int = 20) -> np.ndarray:
    s = pd.Series(np.asarray(x, dtype=float))
    return s.rolling(window=window, center=True, min_periods=1).mean().to_numpy()


def rolling_std(x: Iterable[float], window: int = 20) -> np.ndarray:
    s = pd.Series(np.asarray(x, dtype=float))
    return s.rolling(window=window, center=True, min_periods=1).std().to_numpy()


def running_mean(x: Iterable[float]) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.cumsum(x) / np.arange(1, x.size + 1)


def block_average(x: Iterable[float], n_blocks: int = 10) -> pd.DataFrame:
    """Contiguous-block means (correlation-aware uncertainty estimate)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    n_blocks = max(1, min(n_blocks, n))
    bs = n // n_blocks
    means = [x[b * bs:(b + 1) * bs].mean() for b in range(n_blocks) if x[b*bs:(b+1)*bs].size]
    return pd.DataFrame({"block": np.arange(1, len(means) + 1), "block_mean": means})


def block_average_sem(x: Iterable[float], n_blocks: int = 10) -> float:
    m = block_average(x, n_blocks)["block_mean"].to_numpy()
    return float(np.std(m, ddof=1) / np.sqrt(m.size)) if m.size >= 2 else 0.0


def plateau_detection(t, y, tail_fraction: float = 0.5, slope_tol: float | None = None) -> dict:
    """Assess convergence via tail-slope and first/second-half drift."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = y.size
    i0 = int(n * (1.0 - tail_fraction))
    tt, yy = t[i0:], y[i0:]
    slope, _ = np.polyfit(tt, yy, 1)
    span = (tt[-1] - tt[0]) or 1.0
    if slope_tol is None:
        slope_tol = np.std(yy, ddof=1) / span if span else np.inf
    half = n // 2
    mean_first, mean_second = np.mean(y[:half]), np.mean(y[half:])
    return {
        "tail_slope": float(slope),
        "slope_tolerance": float(slope_tol),
        "mean_first_half": float(mean_first),
        "mean_second_half": float(mean_second),
        "drift_second_minus_first": float(mean_second - mean_first),
        "tail_mean": float(np.mean(yy)),
        "tail_std": float(np.std(yy, ddof=1)),
        # floor guards the degenerate perfectly-flat-tail case (slope_tol == 0)
        "converged": bool(abs(slope) <= max(slope_tol, 1e-9)),
    }
