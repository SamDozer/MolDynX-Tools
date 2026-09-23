"""Correlation matrices (Pearson + Spearman)."""

from __future__ import annotations

import pandas as pd


def correlation_matrices(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (Pearson, Spearman) correlation matrices for numeric columns."""
    return df.corr(method="pearson"), df.corr(method="spearman")
