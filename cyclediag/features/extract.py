"""Tier 1 feature extraction — one row per cycle × leg."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .cc_cv import detect_cv_region
from .lges_catalog import FEATURE_SET_LGES
from .segment_utils import leg_segment

FEATURE_SET = "vp_v1_basic"
FEATURE_SETS = (FEATURE_SET, FEATURE_SET_LGES)


@dataclass
class FeatureConfig:
    """Logical input units: A, Ah, V, s; f_Q_spec is mAh/g of active_mass_g."""
    charge_step: str = "charge"
    discharge_step: str = "discharge"
    active_mass_g: float | None = None
    cell_id: str | None = None
    feature_set: str = FEATURE_SET
    rest_labels: str = "rest"
    rest_current_max: float | None = 0.01


def _capacity_weighted_v_avg(seg: pd.DataFrame) -> float | None:
    if seg.empty or "voltage" not in seg.columns or "capacity" not in seg.columns:
        return None
    v = pd.to_numeric(seg["voltage"], errors="coerce").to_numpy(dtype=float)
    q = pd.to_numeric(seg["capacity"], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(v) & np.isfinite(q)
    # Do not bridge missing samples or integrate across a counter reset.
    if len(q) < 2 or not valid.all():
        return None
    dq = np.diff(q)
    if np.any(dq < 0) or not np.any(dq > 0):
        return None
    v_mid = (v[:-1] + v[1:]) / 2.0
    return float(np.sum(v_mid * dq) / np.sum(dq))


def extract_leg_features(seg: pd.DataFrame, *, leg: str, config: FeatureConfig) -> dict:
    """Scalar Tier 1 features for one charge or discharge leg."""
    row: dict = {"leg": leg, "feature_set": FEATURE_SET}
    if seg is None or seg.empty:
        return row

    seg = seg.copy()
    q_col = f"{leg}_capacity"
    if q_col not in seg.columns:
        q_col = "capacity"
    # One leg-specific series is shared by total capacity, weighted V and CC/CV.
    seg["capacity"] = pd.to_numeric(
        seg.get(q_col, pd.Series(np.nan, index=seg.index)), errors="coerce",
    )
    v = pd.to_numeric(seg.get("voltage", pd.Series(np.nan, index=seg.index)), errors="coerce")
    q = seg["capacity"]
    q_valid = bool(np.isfinite(q).all() and (q >= 0).all() and not (q.diff() < 0).any())
    row["capacity_status"] = "ok" if q_valid else "invalid_or_reset_counter"
    row["capacity_column"] = q_col
    row["unit_warnings"] = list(seg.attrs.get("unit_warnings", []))
    row["feature_units"] = {"f_Q_max": "Ah", "f_Q_spec": "mAh/g", "f_V_avg": "V"}
    q_max = float(q.max()) if q_valid else None
    row["f_Q_max"] = q_max
    if config.active_mass_g and q_max is not None and config.active_mass_g > 0:
        row["f_Q_spec"] = 1000.0 * q_max / config.active_mass_g
    else:
        row["f_Q_spec"] = None

    if v.notna().any():
        row["f_V_start"] = float(v.iloc[0]) if np.isfinite(v.iloc[0]) else None
        row["f_V_end"] = float(v.iloc[-1]) if np.isfinite(v.iloc[-1]) else None
    else:
        row["f_V_start"] = row["f_V_end"] = None
    row["f_V_avg"] = _capacity_weighted_v_avg(seg)
    if not q_valid:
        row.update(f_cc_Q_frac=None, f_cv_time_s=None, f_v_cc_end=None,
                   f_q_cc_end=None, has_cv=None, cv_method="invalid_capacity")
        return row

    cv = detect_cv_region(
        seg,
        v_col="voltage",
        q_col="capacity" if "capacity" in seg.columns else None,
        t_col="time" if "time" in seg.columns else None,
        current_col="current" if "current" in seg.columns else None,
    )
    row["f_cc_Q_frac"] = cv.cc_q_frac
    row["f_cv_time_s"] = cv.cv_time_s
    row["f_v_cc_end"] = cv.v_cc_end
    row["f_q_cc_end"] = cv.q_cc_end
    row["has_cv"] = cv.has_cv
    row["cv_method"] = cv.method
    return row


def extract_cycle_features(
    df: pd.DataFrame,
    cycle: int,
    *,
    config: FeatureConfig | None = None,
    filepath: str = "",
) -> list[dict]:
    cfg = config or FeatureConfig()
    if "cycle" not in df.columns:
        return []
    cycle_df = df[df["cycle"] == cycle]
    if cycle_df.empty:
        return []

    cell_id = cfg.cell_id or Path(filepath).stem
    rows: list[dict] = []
    for leg in ("charge", "discharge"):
        seg = leg_segment(
            cycle_df,
            leg,
            charge_text=cfg.charge_step,
            discharge_text=cfg.discharge_step,
        )
        if seg.empty:
            continue
        feats = extract_leg_features(seg, leg=leg, config=cfg)
        feats.update({
            "cell_id": cell_id,
            "file": filepath,
            "cycle": int(cycle),
        })
        rows.append(feats)
    return rows


def extract_features_table(
    df: pd.DataFrame,
    *,
    cycles: Iterable[int] | None = None,
    filepath: str = "",
    config: FeatureConfig | None = None,
    raw_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """All feature rows for one loaded CSV."""
    if "cycle" not in df.columns:
        raise ValueError("DataFrame must have logical 'cycle' column")
    cfg = config or FeatureConfig()
    if cfg.feature_set == FEATURE_SET_LGES:
        from .lges_extract import LgesExtractConfig, extract_lges_features_table

        lges_cfg = LgesExtractConfig(
            charge_step=cfg.charge_step,
            discharge_step=cfg.discharge_step,
            active_mass_g=cfg.active_mass_g,
            cell_id=cfg.cell_id,
            rest_labels=cfg.rest_labels,
            rest_current_max=cfg.rest_current_max,
        )
        return extract_lges_features_table(
            df,
            cycles=cycles,
            filepath=filepath,
            config=lges_cfg,
            raw_df=raw_df,
        )
    if cfg.cell_id is None and filepath:
        cfg.cell_id = Path(filepath).stem

    cycle_list = (
        list(cycles)
        if cycles is not None
        else sorted(df["cycle"].dropna().unique().astype(int))
    )
    rows: list[dict] = []
    for cyc in cycle_list:
        rows.extend(extract_cycle_features(df, int(cyc), config=cfg, filepath=filepath))
    return pd.DataFrame(rows)
