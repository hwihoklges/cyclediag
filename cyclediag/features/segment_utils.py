"""Step-type classification and leg/rest segmentation for cycle extraction."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .cc_cv import resolve_current_column


def parse_step_labels(text: str) -> list[str]:
    return [p.strip().lower() for p in str(text).split(",") if p.strip()]


def classify_step_kind(
    step: str,
    *,
    charge_text: str = "charge",
    discharge_text: str = "discharge",
    rest_text: str = "rest",
) -> str:
    """Classify one StepType value. Discharge is checked before charge."""
    s = str(step).strip().lower()
    for lb in parse_step_labels(discharge_text):
        if s == lb or lb in s:
            return "discharge"
    for lb in parse_step_labels(charge_text):
        if s == lb or lb in s:
            return "charge"
    for lb in parse_step_labels(rest_text):
        if s == lb or lb in s:
            return "rest"
    return "other"


def _classify_with_current(
    kind: str,
    current: float | None,
    *,
    rest_current_max: float | None,
) -> str:
    if kind != "other":
        return kind
    if (
        rest_current_max is not None
        and current is not None
        and np.isfinite(current)
        and abs(current) <= rest_current_max
    ):
        return "rest"
    return "other"


def segment_kinds(
    steps: pd.Series,
    currents: pd.Series | None,
    *,
    charge_text: str,
    discharge_text: str,
    rest_text: str,
    rest_current_max: float | None,
) -> list[str]:
    kinds: list[str] = []
    for i, step in enumerate(steps.astype(str)):
        cur = None
        if currents is not None and i < len(currents):
            cur = pd.to_numeric(currents.iloc[i], errors="coerce")
            cur = float(cur) if np.isfinite(cur) else None
        base = classify_step_kind(
            step,
            charge_text=charge_text,
            discharge_text=discharge_text,
            rest_text=rest_text,
        )
        kinds.append(_classify_with_current(base, cur, rest_current_max=rest_current_max))
    return kinds


def leg_segment(
    cycle_df: pd.DataFrame,
    leg: str,
    *,
    charge_text: str,
    discharge_text: str,
    rest_text: str = "rest",
    rest_current_max: float | None = None,
) -> pd.DataFrame:
    st_col = "step_type"
    if cycle_df is None or cycle_df.empty or st_col not in cycle_df.columns:
        return pd.DataFrame()
    i_col = resolve_current_column(cycle_df)
    currents = pd.to_numeric(cycle_df[i_col], errors="coerce") if i_col else None
    kinds = segment_kinds(
        cycle_df[st_col],
        currents,
        charge_text=charge_text,
        discharge_text=discharge_text,
        rest_text=rest_text,
        rest_current_max=rest_current_max,
    )
    mask = [k == leg for k in kinds]
    result = cycle_df.loc[mask].copy()
    positions = np.flatnonzero(mask)
    result.attrs["acquisition_disjoint"] = bool(np.any(np.diff(positions) != 1))
    return result


def iter_rest_periods(
    cycle_df: pd.DataFrame,
    *,
    charge_text: str,
    discharge_text: str,
    rest_text: str = "rest",
    rest_current_max: float | None = None,
) -> list[tuple[str, pd.DataFrame]]:
    """Yield (after_leg, rest_df) for each rest block following charge/discharge."""
    st_col = "step_type"
    if cycle_df is None or cycle_df.empty or st_col not in cycle_df.columns:
        return []

    df = cycle_df.reset_index(drop=True)
    i_col = resolve_current_column(df)
    currents = pd.to_numeric(df[i_col], errors="coerce") if i_col else None
    kinds = segment_kinds(
        df[st_col],
        currents,
        charge_text=charge_text,
        discharge_text=discharge_text,
        rest_text=rest_text,
        rest_current_max=rest_current_max,
    )
    if not kinds:
        return []

    segments: list[tuple[str, int, int]] = []
    start = 0
    for i in range(1, len(kinds)):
        if kinds[i] != kinds[i - 1]:
            segments.append((kinds[i - 1], start, i - 1))
            start = i
    segments.append((kinds[-1], start, len(kinds) - 1))

    out: list[tuple[str, pd.DataFrame]] = []
    for idx, (kind, s_idx, e_idx) in enumerate(segments):
        if kind != "rest":
            continue
        prev_kind = segments[idx - 1][0] if idx > 0 else None
        if prev_kind not in ("charge", "discharge"):
            continue
        out.append((prev_kind, df.iloc[s_idx : e_idx + 1].copy()))
    return out
