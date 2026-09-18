"""dQ/dV and dV/dQ peak indicator screening and cross-cell comparison."""

from __future__ import annotations

import re
from typing import Literal

import numpy as np
import pandas as pd

from cyclediag.features._statistics import finite_pearson
from cyclediag.analysis.indicator_screen import (
    _health_column,
    _late_early_split,
    compare_cells,
    screen_indicators,
    screen_indicators_by_file,
)
from cyclediag.features.lges_catalog import dqdv_peak_column_names

_DQDV_RE = re.compile(
    r"^(chg|dchg)_(dQdV_peak(\d+)|dVdQ_peak(\d+))(_V|_Q)?$"
)

PeakKind = Literal["dQdV", "dVdQ"]
LegKind = Literal["chg", "dchg"]
MetricKind = Literal["V", "H", "Q"]


def is_dqdv_column(col: str) -> bool:
    if col in dqdv_peak_column_names():
        return True
    return bool(_DQDV_RE.match(col) or "dQdV_peak" in col or "dVdQ_peak" in col)


def dqdv_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if is_dqdv_column(c)]


def parse_dqdv_column(col: str) -> dict | None:
    """Parse e.g. chg_dQdV_peak3_V → leg, kind, peak, metric."""
    m = re.match(r"^(chg|dchg)_(dQdV|dVdQ)_peak(\d+)$", col)
    if m:
        return {
            "leg": m.group(1),
            "kind": m.group(2),
            "peak": int(m.group(3)),
            "metric": "H",
            "column": col,
        }
    m = re.match(r"^(chg|dchg)_(dQdV|dVdQ)_peak(\d+)_([VQ])$", col)
    if m:
        return {
            "leg": m.group(1),
            "kind": m.group(2),
            "peak": int(m.group(3)),
            "metric": m.group(4),
            "column": col,
        }
    return None


def _dqdv_category(parsed: dict) -> str:
    kind = parsed["kind"]
    metric = parsed["metric"]
    if metric == "V":
        return f"{kind}_peak_voltage"
    if metric == "Q":
        return f"{kind}_peak_capacity"
    return f"{kind}_peak_height"


def _slope_per_100_cycles(cycles: np.ndarray, values: np.ndarray) -> float | None:
    mask = np.isfinite(cycles) & np.isfinite(values)
    if mask.sum() < 4:
        return None
    x = cycles[mask].astype(float)
    y = values[mask]
    try:
        coef = np.polyfit(x, y, 1)
        return float(coef[0] * 100.0)
    except (np.linalg.LinAlgError, ValueError):
        return None


def screen_dqdv_indicators(
    features: pd.DataFrame,
    *,
    reference_cycle: int = 1,
    health_col: str | None = None,
) -> pd.DataFrame:
    """Screen only dQ/dV peak columns; add peak metadata and drift metrics."""
    if features is None or features.empty:
        return pd.DataFrame()

    dq_cols = dqdv_columns(features)
    if not dq_cols:
        return pd.DataFrame()

    sub = features.copy()
    base_screen = screen_indicators(sub, reference_cycle=reference_cycle, health_col=health_col)
    if base_screen.empty:
        return pd.DataFrame()

    out = base_screen[base_screen["feature"].map(is_dqdv_column)].copy()
    if out.empty:
        return out

    df = features.sort_values("cycle")
    hcol = health_col or _health_column(df)
    early_m, late_m = _late_early_split(df["cycle"])
    base = df[df["cycle"] == reference_cycle]
    if base.empty:
        base = df.head(1)
    base_row = base.iloc[0]
    cycles = pd.to_numeric(df["cycle"], errors="coerce").to_numpy(dtype=float)

    meta_rows = []
    for _, row in out.iterrows():
        col = row["feature"]
        parsed = parse_dqdv_column(col)
        if not parsed:
            meta_rows.append({})
            continue

        s = pd.to_numeric(df[col], errors="coerce")
        b = base_row.get(col)
        slope = _slope_per_100_cycles(cycles, s.to_numpy(dtype=float))

        delta_mV = None
        if parsed["metric"] == "V" and b is not None and np.isfinite(b):
            last = s.iloc[-1]
            if np.isfinite(last):
                delta_mV = (float(last) - float(b)) * 1000.0

        height_fade_pct = None
        if parsed["metric"] == "H" and b is not None and np.isfinite(b) and abs(b) > 1e-12:
            last = s.iloc[-1]
            if np.isfinite(last):
                height_fade_pct = (float(last) - float(b)) / abs(float(b)) * 100.0

        corr_h = np.nan
        if hcol and hcol in df.columns:
            corr_h = finite_pearson(s, df[hcol])

        meta_rows.append({
            "leg": parsed["leg"],
            "peak_kind": parsed["kind"],
            "peak_num": parsed["peak"],
            "metric": parsed["metric"],
            "category": _dqdv_category(parsed),
            "slope_per_100cyc": slope,
            "delta_V_mV": delta_mV,
            "height_fade_pct": height_fade_pct,
            "corr_health": round(float(corr_h), 3) if corr_h is not None and pd.notna(corr_h) else row.get("corr_health"),
        })

    meta_df = pd.DataFrame(meta_rows)
    out = pd.concat([out.reset_index(drop=True), meta_df], axis=1)

    # Boost severity for large peak-V shift (material fingerprint drift)
    if "delta_V_mV" in out.columns:
        v_shift = out["delta_V_mV"].abs()
        if v_shift.notna().any():
            out.loc[v_shift.notna(), "severity"] = out.loc[v_shift.notna(), "severity"].combine(
                (v_shift / 50.0).clip(0, 1),
                max,
            )
            mask = v_shift >= 10
            out.loc[mask, "signal"] = out.loc[mask, "signal"].astype(str) + "; peak-V shift"

    return out.sort_values("severity", ascending=False).reset_index(drop=True)


def screen_dqdv_by_file(features: pd.DataFrame) -> pd.DataFrame:
    if features is None or features.empty:
        return pd.DataFrame()
    group_cols = [c for c in ("cell_id", "file") if c in features.columns]
    if not group_cols:
        return screen_dqdv_indicators(features)

    parts = []
    for _, grp in features.groupby(group_cols, sort=False):
        screened = screen_dqdv_indicators(grp)
        if screened.empty:
            continue
        for c in group_cols:
            screened[c] = grp[c].iloc[0]
        parts.append(screened)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def compare_cells_dqdv(features: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Cross-cell divergence for dQ/dV peak indicators only."""
    if features is None or features.empty:
        return pd.DataFrame()
    dq_cols = dqdv_columns(features)
    if not dq_cols:
        return pd.DataFrame()
    keep = [c for c in features.columns if c in dq_cols or c in ("cell_id", "file", "cycle")]
    sub = features[keep].copy()
    cmp = compare_cells(sub, **kwargs)
    if cmp.empty:
        return cmp
    parsed = cmp["feature"].map(parse_dqdv_column)
    cmp["leg"] = [p["leg"] if p else None for p in parsed]
    cmp["peak_kind"] = [p["kind"] if p else None for p in parsed]
    cmp["peak_num"] = [p["peak"] if p else None for p in parsed]
    cmp["metric"] = [p["metric"] if p else None for p in parsed]
    return cmp


def dqdv_trajectory_long(features: pd.DataFrame) -> pd.DataFrame:
    """Long-format table: one row per (cell, cycle, peak indicator)."""
    if features is None or features.empty:
        return pd.DataFrame()

    id_cols = [c for c in ("cell_id", "file", "cycle") if c in features.columns]
    dq_cols = dqdv_columns(features)
    if not id_cols or not dq_cols:
        return pd.DataFrame()

    rows: list[dict] = []
    for _, row in features.iterrows():
        for col in dq_cols:
            val = row.get(col)
            if val is None or (isinstance(val, float) and not np.isfinite(val)):
                continue
            parsed = parse_dqdv_column(col)
            if not parsed:
                continue
            rec = {c: row[c] for c in id_cols}
            rec.update({
                "indicator": col,
                "value": float(val),
                "leg": parsed["leg"],
                "peak_kind": parsed["kind"],
                "peak_num": parsed["peak"],
                "metric": parsed["metric"],
            })
            rows.append(rec)
    return pd.DataFrame(rows)


def top_dqdv_problems(screened: pd.DataFrame, *, n: int = 12) -> pd.DataFrame:
    if screened is None or screened.empty:
        return pd.DataFrame()
    # Prefer peak-V and height metrics with coverage
    pref = screened.copy()
    if "metric" in pref.columns:
        pref["_rank"] = pref["metric"].map({"V": 0, "H": 1, "Q": 2}).fillna(3)
        pref = pref.sort_values(["severity", "_rank"], ascending=[False, True])
        return pref.drop(columns=["_rank"], errors="ignore").head(n)
    return screened.head(n)
