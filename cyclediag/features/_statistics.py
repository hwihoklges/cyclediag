"""Numerically safe descriptive statistics; undefined is not evidence."""

import numpy as np
import pandas as pd


def finite_pearson(x: pd.Series, y: pd.Series) -> float:
    """Pairwise-finite Pearson r, undefined for either constant paired series."""
    x, y = x.align(y, join="inner")
    a = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(y, errors="coerce").to_numpy(dtype=float)
    paired = np.isfinite(a) & np.isfinite(b)
    a, b = a[paired], b[paired]
    if len(a) < 2 or np.all(a == a[0]) or np.all(b == b[0]):
        return float("nan")
    # Scale before centering to avoid overflow for large finite magnitudes.
    a = a / np.max(np.abs(a))
    b = b / np.max(np.abs(b))
    a, b = a - a.mean(), b - b.mean()
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return float("nan")
    return float(np.clip(np.dot(a / na, b / nb), -1.0, 1.0))