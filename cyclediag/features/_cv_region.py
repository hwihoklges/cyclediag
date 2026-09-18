"""Detect CC / CV regions in a charge or discharge leg (current-based)."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class CvRegionInfo:
    """Per-leg CC/CV segmentation result."""

    has_cv: bool
    cv_start_index: int | None
    cc_end_index: int | None
    n_points: int
    n_cc: int
    n_cv: int
    cc_q_frac: float | None
    cv_time_s: float | None
    v_cv_start: float | None
    q_cv_start: float | None
    q_cc_end: float | None
    v_cc_end: float | None
    current_col: str | None
    method: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def resolve_current_column(df: pd.DataFrame) -> str | None:
    """Pick a current column for CC/CV detection (AvgCurrent preferred)."""
    priority: list[str] = []
    fallback: list[str] = []
    for col in df.columns:
        clean = str(col).lower().replace(" ", "").replace("_", "")
        if clean in ("avgcurrent", "avgcurrentma", "avgcurrenta"):
            priority.append(col)
        elif clean in ("current", "curr", "i", "currenta", "currentma"):
            fallback.append(col)
    # Normalization owns unit conversion; don't bypass canonical A with a raw
    # auxiliary AvgCurrent(mA) column that was intentionally left untouched.
    ordered = (["current"] if "current" in df.columns else []) + priority + fallback
    for col in dict.fromkeys(ordered):
        series = pd.to_numeric(df[col], errors="coerce").abs()
        valid = series[np.isfinite(series) & (series > 0)]
        if len(valid) >= 4:
            return col
    return None


def detect_cv_region(
    df: pd.DataFrame,
    *,
    v_col: str,
    q_col: str | None = None,
    t_col: str | None = None,
    current_col: str | None = None,
    cc_current_frac: float = 0.92,
    sustained_points: int = 5,
    min_q_frac_for_cv: float = 0.85,
) -> CvRegionInfo:
    """Find the first sustained CV region (current drop) in *df*.

    Returns indices into *df* (after any caller-side filtering). If no CV is
    detected, ``has_cv`` is False and the full segment is treated as CC.
    """
    n = len(df)
    empty = CvRegionInfo(
        has_cv=False,
        cv_start_index=None,
        cc_end_index=None,
        n_points=n,
        n_cc=n,
        n_cv=0,
        cc_q_frac=None,
        cv_time_s=None,
        v_cv_start=None,
        q_cv_start=None,
        q_cc_end=None,
        v_cc_end=None,
        current_col=None,
        method="too_short" if n < 4 else "no_current_col",
    )
    if n < 4:
        return empty

    cur_col = current_col or resolve_current_column(df)
    if cur_col is None:
        return empty

    try:
        curr = pd.to_numeric(df[cur_col], errors="coerce").abs().to_numpy(dtype=float)
        cc_slice = curr[: max(n // 2, 4)]
        cc_valid = cc_slice[np.isfinite(cc_slice) & (cc_slice > 0)]
        if len(cc_valid) < 4:
            return CvRegionInfo(
                **{**empty.to_dict(), "current_col": cur_col, "method": "weak_current"}
            )

        cc_median = float(np.median(cc_valid))
        if cc_median <= 0:
            return CvRegionInfo(
                **{**empty.to_dict(), "current_col": cur_col, "method": "zero_cc_current"}
            )

        threshold = cc_median * cc_current_frac
        need = max(3, int(sustained_points))

        v_arr = (
            pd.to_numeric(df[v_col], errors="coerce").to_numpy(dtype=float)
            if v_col in df.columns
            else None
        )
        q_arr = None
        q_max = None
        if q_col and q_col in df.columns:
            q_arr = pd.to_numeric(df[q_col], errors="coerce").to_numpy(dtype=float)
            finite_q = q_arr[np.isfinite(q_arr)]
            if finite_q.size:
                q_max = float(finite_q.max())

        t_arr = None
        if t_col and t_col in df.columns:
            t_arr = pd.to_numeric(df[t_col], errors="coerce").to_numpy(dtype=float)

        search_start = max(n // 10, 3)
        if q_max and q_max > 0 and q_arr is not None:
            high_q_idx = np.where(q_arr >= min_q_frac_for_cv * q_max)[0]
            if high_q_idx.size:
                search_start = max(search_start, int(high_q_idx[0]))

        def _sustained_low(start_idx: int) -> bool:
            if start_idx + need > n:
                return False
            window = curr[start_idx : start_idx + need]
            return bool(
                np.all(np.isfinite(window) & (window < threshold) | ~np.isfinite(window))
            )

        cv_start = None
        for i in range(search_start, n - need + 1):
            if not _sustained_low(i):
                continue
            if q_max and q_arr is not None:
                q_at = q_arr[i]
                if not np.isfinite(q_at) or q_at < min_q_frac_for_cv * q_max:
                    continue
            cv_start = i
            break

        if cv_start is None:
            return CvRegionInfo(
                has_cv=False,
                cv_start_index=None,
                cc_end_index=n - 1 if n else None,
                n_points=n,
                n_cc=n,
                n_cv=0,
                cc_q_frac=1.0 if n else None,
                cv_time_s=None,
                v_cv_start=None,
                q_cv_start=None,
                q_cc_end=_safe_idx(q_arr, n - 1),
                v_cc_end=_safe_idx(v_arr, n - 1),
                current_col=cur_col,
                method="cc_only",
            )

        cc_end = max(0, cv_start - 1)
        n_cc = cv_start  # rows 0 .. cv_start-1 are CC
        n_cv = n - cv_start

        cc_q_frac = None
        if q_max and q_max > 0 and q_arr is not None and cc_end >= 0:
            q_end = q_arr[cc_end]
            if np.isfinite(q_end):
                cc_q_frac = float(q_end / q_max)

        cv_time = None
        if t_arr is not None and cv_start < n:
            t0 = t_arr[cv_start]
            t1 = t_arr[-1]
            if np.isfinite(t0) and np.isfinite(t1):
                cv_time = float(max(0.0, t1 - t0))

        return CvRegionInfo(
            has_cv=True,
            cv_start_index=cv_start,
            cc_end_index=cc_end,
            n_points=n,
            n_cc=n_cc,
            n_cv=n_cv,
            cc_q_frac=cc_q_frac,
            cv_time_s=cv_time,
            v_cv_start=_safe_idx(v_arr, cv_start),
            q_cv_start=_safe_idx(q_arr, cv_start),
            q_cc_end=_safe_idx(q_arr, cc_end),
            v_cc_end=_safe_idx(v_arr, cc_end),
            current_col=cur_col,
            method="current_drop",
        )
    except Exception:
        return CvRegionInfo(
            **{**empty.to_dict(), "current_col": cur_col, "method": "error"}
        )


def _safe_idx(arr: np.ndarray | None, idx: int) -> float | None:
    if arr is None or idx < 0 or idx >= len(arr):
        return None
    val = arr[idx]
    return float(val) if np.isfinite(val) else None


def trim_to_cc_end(df: pd.DataFrame, info: CvRegionInfo, *, include_boundary: bool = True) -> pd.DataFrame:
    """Return CC portion of *df* (optionally including first CV point — legacy trim)."""
    if not info.has_cv or info.cv_start_index is None:
        return df
    end = info.cv_start_index + (1 if include_boundary else 0)
    if end <= 0:
        return df.iloc[:0].copy()
    return df.iloc[:end].copy()
