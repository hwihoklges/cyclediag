"""LGES cycle-level feature extraction (one row per cycle)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .dqdv_segment import prepare_leg_segment_for_dqdv
from .dqdv_peaks import (
    DEFAULT_DQDV_PEAK_CONFIG,
    DqdvPeakConfig,
    find_dqdv_peaks,
    find_dvdq_peaks,
    peaks_to_columns,
)
from .lges_extra_indicators import (
    correct_r_to_25c,
    dtw_distance,
    extract_absolute_dvdq_indicators,
    extract_shape_indicators,
    fit_rest_tau,
    rolling_slope,
    safe_diff,
    safe_ratio,
    vq_norm_curve,
)
from .extract import FeatureConfig
from .lges_catalog import (
    DELTA_ABS_COLS,
    DELTA_PCT_COLS,
    FEATURE_SET_LGES,
    RESISTANCE_OFFSETS_S,
    resistance_offset_label,
)
from .band_capacity import BandCapacityConfig, discharge_band_capacity
from .cell_meta import CellProtocolMeta, DEFAULT_Q_RATED_AH
from .segment_utils import iter_rest_periods, leg_segment
from .cc_cv import resolve_current_column
from .signal_cv import detect_cv_signal, signal_cv_to_row
from .units import capacity_to_ah, canonical_unit_factor
from cyclediag.io.rest_voltage import extract_cycle_rest_voltages


@dataclass
class LgesExtractConfig(FeatureConfig):
    rest_labels: str = "rest"
    rest_current_max: float | None = 0.5  # ~0.007*Q for 72 Ah; was 0.01 (too small)
    temperature_col: str | None = None
    baseline_cycle: int = 1
    auto_baseline: bool = True
    dqdv_peak_config: DqdvPeakConfig | None = None
    with_diagnosis: bool = True
    diagnosis_config_path: str | None = None
    enrich_assb: bool = True
    expected_pulse_current: float | None = None
    q_rated_ah: float = DEFAULT_Q_RATED_AH
    routine_c_rate: float = 0.5
    rpt_c_rate: float = 1.0 / 3.0
    dcir_c_rate: float = 1.0
    band_capacity: BandCapacityConfig | None = None
    capacity_unit: str | None = None  # "ah" | "mah" | None (header-driven)

    def protocol_meta(self) -> CellProtocolMeta:
        return CellProtocolMeta(
            q_rated_ah=float(self.q_rated_ah),
            routine_c_rate=float(self.routine_c_rate),
            rpt_c_rate=float(self.rpt_c_rate),
            dcir_c_rate=float(self.dcir_c_rate),
        )

    def resolved_pulse_current(self) -> float:
        if self.expected_pulse_current is not None:
            return float(self.expected_pulse_current)
        return self.protocol_meta().dcir_pulse_current_a

    def resolved_rest_current_max(self) -> float:
        if self.rest_current_max is not None:
            return float(self.rest_current_max)
        return self.protocol_meta().rest_current_max_a

    def resolved_dqdv_config(self) -> DqdvPeakConfig:
        return self.dqdv_peak_config or DEFAULT_DQDV_PEAK_CONFIG


def _capacity_to_ah(q_max: float | None, *, unit: str | None = None) -> float | None:
    if q_max is None or not np.isfinite(q_max):
        return None
    return capacity_to_ah(q_max, unit=unit)


def _resolve_temperature_col(df: pd.DataFrame, preferred: str | None) -> str | None:
    if preferred and preferred in df.columns:
        return preferred
    if "temperature" in df.columns:
        return "temperature"
    for name in ("Temperature", "Temp", "CellTemp", "Aux_Temperature"):
        if name in df.columns:
            return name
    return None


def _leg_capacity_ah(seg: pd.DataFrame, leg: str, *, unit: str | None = None) -> float | None:
    if seg is None or seg.empty:
        return None
    if leg == "charge":
        cols = ("charge_capacity", "capacity")
    else:
        cols = ("discharge_capacity", "capacity")
    for col in cols:
        if col not in seg.columns:
            continue
        q = pd.to_numeric(seg[col], errors="coerce")
        if not np.isfinite(q).all() or (q < 0).any() or (q.diff() < 0).any():
            return None
        if q.notna().any():
            return _capacity_to_ah(float(q.max()), unit=unit)
    return None


def _capacity_series(seg: pd.DataFrame, leg: str) -> np.ndarray | None:
    if seg is None or seg.empty:
        return None
    col = "charge_capacity" if leg == "charge" else "discharge_capacity"
    if col not in seg.columns:
        col = "capacity"
    if col not in seg.columns:
        return None
    q = pd.to_numeric(seg[col], errors="coerce").to_numpy(dtype=float)
    return q if np.isfinite(q).all() and np.all(q >= 0) and not np.any(np.diff(q) < 0) else None


def _relative_time_s(seg: pd.DataFrame) -> np.ndarray:
    if "step_time" in seg.columns:
        st = pd.to_numeric(seg["step_time"], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(st).sum() >= 2:
            st0 = st[np.isfinite(st)][0]
            rel = st - st0
            if np.nanmax(rel) > 0:
                return rel
        elif np.isfinite(st).any():
            return st - st[np.isfinite(st)][0]
    if "time" in seg.columns:
        t = pd.to_numeric(seg["time"], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(t).sum() >= 2:
            t0 = t[np.isfinite(t)][0]
            return t - t0
    return np.arange(len(seg), dtype=float)


def _avg_temperature(seg: pd.DataFrame, temp_col: str | None) -> float | None:
    if not temp_col or temp_col not in seg.columns or seg.empty:
        return None
    t = pd.to_numeric(seg[temp_col], errors="coerce")
    if t.notna().sum() == 0:
        return None
    return float(t.mean())


def _cutoff_v_i(seg: pd.DataFrame) -> tuple[float | None, float | None]:
    if seg.empty:
        return None, None
    v = pd.to_numeric(seg.get("voltage"), errors="coerce")
    i_col = resolve_current_column(seg) or "current"
    i = pd.to_numeric(seg.get(i_col), errors="coerce")
    v_end = float(v.iloc[-1]) if v.notna().any() else None
    i_end = float(i.iloc[-1]) if i.notna().any() else None
    return v_end, i_end


def _sample_v_i_at_offsets(
    seg: pd.DataFrame,
    offsets_s: Iterable[float],
) -> dict[float, tuple[float | None, float | None]]:
    if seg.empty:
        return {o: (None, None) for o in offsets_s}
    t_rel = _relative_time_s(seg)
    v = pd.to_numeric(seg["voltage"], errors="coerce").to_numpy(dtype=float)
    i_col = resolve_current_column(seg) or "current"
    i = pd.to_numeric(seg[i_col], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(t_rel) & np.isfinite(v)
    if not valid.any():
        return {o: (None, None) for o in offsets_s}
    t_arr = t_rel[valid]
    v_arr = v[valid]
    i_arr = i[valid] if len(i) == len(t_rel) else np.full_like(v_arr, np.nan)
    order = np.argsort(t_arr)
    t_arr, v_arr, i_arr = t_arr[order], v_arr[order], i_arr[order]
    v0 = float(v_arr[0])
    out: dict[float, tuple[float | None, float | None]] = {}
    for offset in offsets_s:
        if offset > t_arr[-1] + 1e-9:
            out[offset] = (None, None)
        elif valid.sum() == 1:
            out[offset] = (v0, float(i_arr[0]) if np.isfinite(i_arr[0]) else None)
        else:
            vt = float(np.interp(offset, t_arr, v_arr))
            it = float(np.interp(offset, t_arr, i_arr)) if np.isfinite(i_arr).any() else None
            out[offset] = (vt, it)
    out[0.0] = (v0, float(i_arr[0]) if np.isfinite(i_arr[0]) else None)
    return out


def _resistance_mohm(v0: float | None, vt: float | None, it: float | None) -> float | None:
    if v0 is None or vt is None or it is None:
        return None
    if not np.isfinite(v0) or not np.isfinite(vt) or not np.isfinite(it):
        return None
    if abs(it) < 1e-9:
        return None
    return abs(v0 - vt) / abs(it) * 1000.0


def _rest_voltages_for_leg(
    cycle_df: pd.DataFrame,
    *,
    charge_step: str,
    discharge_step: str,
    rest_labels: str,
    rest_current_max: float | None,
) -> dict[str, dict[str, float | None]]:
    """Return {'charge': {init, 60s, 30m, end}, 'discharge': {...}}."""
    work = cycle_df.copy()
    if "TotalCycle" not in work.columns and "cycle" in work.columns:
        work["TotalCycle"] = work["cycle"]

    offsets = sorted({0.0, 60.0, 1800.0})
    step_t_col = "step_time" if "step_time" in work.columns else None
    t_col = "time" if "time" in work.columns else None
    result = extract_cycle_rest_voltages(
        work,
        v_col="voltage",
        st_col="step_type",
        charge_label=charge_step,
        discharge_label=discharge_step,
        rest_labels_text=rest_labels,
        offsets_s=offsets,
        t_col=t_col,
        step_t_col=step_t_col,
        rest_current_max=rest_current_max,
    )

    out: dict[str, dict[str, float | None]] = {
        "charge": {"init": None, "60s": None, "30m": None, "end": None},
        "discharge": {"init": None, "60s": None, "30m": None, "end": None},
    }
    offset_key = {0.0: "init", 60.0: "60s", 1800.0: "30m"}

    for period in result.periods:
        leg = period.after_leg
        if leg not in out:
            continue
        for sample in period.samples:
            key = offset_key.get(sample.offset_s)
            if key and sample.voltage is not None and np.isfinite(sample.voltage):
                if sample.status in ("ok", "single_point") or key == "init":
                    out[leg][key] = float(sample.voltage)

    for after_leg, rest_df in iter_rest_periods(
        cycle_df,
        charge_text=charge_step,
        discharge_text=discharge_step,
        rest_text=rest_labels,
        rest_current_max=rest_current_max,
    ):
        v = pd.to_numeric(rest_df.get("voltage"), errors="coerce")
        if v.notna().any():
            out[after_leg]["end"] = float(v.iloc[-1])
        if out[after_leg]["init"] is None and v.notna().any():
            out[after_leg]["init"] = float(v.iloc[0])

    return out


def _empty_lges_row() -> dict:
    from .lges_catalog import all_lges_feature_columns

    row = {c: None for c in all_lges_feature_columns()}
    return row


def extract_lges_cycle_row(
    df: pd.DataFrame,
    cycle: int,
    *,
    config: LgesExtractConfig | None = None,
    filepath: str = "",
    raw_df: pd.DataFrame | None = None,
) -> dict:
    cfg = config or LgesExtractConfig()
    if cfg.capacity_unit:
        factor = canonical_unit_factor(cfg.capacity_unit, "ah")
        if df.attrs.get("unit_schema") and factor != 1.0:
            raise ValueError("Normalized input is already Ah; specify raw capacity units in ColumnMap.units")
        if factor != 1.0:
            # Legacy logical-frame override applies to every capacity consumer,
            # not just maxima: energy, peaks, CC/CV and band capacities too.
            df = df.copy()
            for col in ("capacity", "charge_capacity", "discharge_capacity", "charge_cv_capacity", "ChargeCVCapacity"):
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce") * factor
    cycle_df = df[df["cycle"] == cycle].copy()
    if cycle_df.empty:
        return {}

    cell_id = cfg.cell_id or Path(filepath).stem
    row: dict = {
        "cell_id": cell_id,
        "file": filepath,
        "cycle": int(cycle),
        "feature_set": FEATURE_SET_LGES,
    }
    row.update(_empty_lges_row())
    row["unit_warnings"] = list(df.attrs.get("unit_warnings", []))
    row["feature_units"] = {"chgCapa": "Ah", "dchgCapa": "Ah", "chg_E": "Wh",
                            "dchg_E": "Wh", "CE": "percent", "VE": "fraction", "EE": "fraction"}

    temp_col = _resolve_temperature_col(raw_df if raw_df is not None else cycle_df, cfg.temperature_col)
    seg_kw = dict(
        charge_text=cfg.charge_step,
        discharge_text=cfg.discharge_step,
        rest_text=cfg.rest_labels,
        rest_current_max=cfg.rest_current_max,
    )
    chg = leg_segment(cycle_df, "charge", **seg_kw)
    dchg = leg_segment(cycle_df, "discharge", **seg_kw)
    # Peak preprocessing may trim CV: never use it for total capacity/energy.
    peak_chg = prepare_leg_segment_for_dqdv(chg, "charge")
    peak_dchg = prepare_leg_segment_for_dqdv(dchg, "discharge")

    chg_v, chg_i = _cutoff_v_i(chg)
    dchg_v, _ = _cutoff_v_i(dchg)
    row["chg_V_cutoff"] = chg_v
    row["dchg_V_cutoff"] = dchg_v
    row["chg_I_cutoff"] = chg_i
    row["chg_temp_avg"] = _avg_temperature(chg, temp_col)
    row["dchg_temp_avg"] = _avg_temperature(dchg, temp_col)

    rest_v = _rest_voltages_for_leg(
        cycle_df,
        charge_step=cfg.charge_step,
        discharge_step=cfg.discharge_step,
        rest_labels=cfg.rest_labels,
        rest_current_max=cfg.rest_current_max,
    )
    for suffix, src in (("EoC", "charge"), ("EoD", "discharge")):
        for key, col_suffix in (("init", "init"), ("60s", "60s"), ("30m", "30m"), ("end", "end")):
            row[f"{suffix}_restV_{col_suffix}"] = rest_v[src].get(key)

    row["EoC_restV_relax"] = safe_diff(row.get("EoC_restV_end"), row.get("EoC_restV_init"))
    row["EoD_restV_relax"] = safe_diff(row.get("EoD_restV_end"), row.get("EoD_restV_init"))
    row["EoC_restV_relax_60s"] = safe_diff(row.get("EoC_restV_60s"), row.get("EoC_restV_init"))
    row["EoD_restV_relax_60s"] = safe_diff(row.get("EoD_restV_60s"), row.get("EoD_restV_init"))

    for prefix, seg in (("EoC_dchgR", dchg), ("EoD_chgR", chg)):
        samples = _sample_v_i_at_offsets(seg, RESISTANCE_OFFSETS_S)
        v0, _ = samples.get(0.0, (None, None))
        r_by_label: dict[str, float | None] = {}
        for off in RESISTANCE_OFFSETS_S:
            label = resistance_offset_label(off)
            vt, it = samples.get(off, (None, None))
            rval = _resistance_mohm(v0, vt, it)
            row[f"{prefix}_{label}"] = rval
            r_by_label[label] = rval
        r0 = r_by_label.get("0p1s")
        r10 = r_by_label.get("10s")
        r30 = r_by_label.get("30s")
        if r10 is not None and r0 is not None and np.isfinite(r10) and np.isfinite(r0):
            row[f"{prefix}_R10_minus_R0p1"] = float(r10) - float(r0)
        else:
            row[f"{prefix}_R10_minus_R0p1"] = None
        if r30 is not None and r0 is not None and np.isfinite(r30) and np.isfinite(r0):
            row[f"{prefix}_R30_minus_R0p1"] = float(r30) - float(r0)
        else:
            row[f"{prefix}_R30_minus_R0p1"] = None

    row["EoC_dchgR_10_60_ratio"] = safe_ratio(row.get("EoC_dchgR_10s"), row.get("EoC_dchgR_60s"))
    row["EoD_chgR_10_60_ratio"] = safe_ratio(row.get("EoD_chgR_10s"), row.get("EoD_chgR_60s"))
    row["EoC_dchgR_10s_T25"] = correct_r_to_25c(row.get("EoC_dchgR_10s"), row.get("dchg_temp_avg"))
    row["EoD_chgR_10s_T25"] = correct_r_to_25c(row.get("EoD_chgR_10s"), row.get("chg_temp_avg"))

    cap_unit = "ah"  # All capacity columns converted together above or by IO.
    chg_ah = _leg_capacity_ah(chg, "charge", unit=cap_unit)
    dchg_ah = _leg_capacity_ah(dchg, "discharge", unit=cap_unit)
    row["chgCapa"] = chg_ah
    row["dchgCapa"] = dchg_ah
    row["capacity_status"] = {
        "charge": "ok" if chg_ah is not None else "missing_or_invalid_counter",
        "discharge": "ok" if dchg_ah is not None else "missing_or_invalid_counter",
    }
    if chg_ah and dchg_ah and chg_ah > 0:
        row["CE"] = dchg_ah / chg_ah * 100.0
    else:
        row["CE"] = None

    # §5.14 signal-based CV (fixes chgCVcapa=0 when raw has ChargeCVCapacity)
    column_cv = None
    for cname in ("ChargeCVCapacity", "ChargeCVCapacity (Ah)", "charge_cv_capacity"):
        src = cycle_df if cycle_df is not None else chg
        if cname in src.columns:
            column_cv = float(pd.to_numeric(src[cname], errors="coerce").max())
            break
    cv = detect_cv_signal(chg, column_cv_ah=column_cv)
    row.update(signal_cv_to_row(cv))

    chg_q = _capacity_series(chg, "charge")
    dqdv_cfg = cfg.resolved_dqdv_config()
    peak_chg_q = _capacity_series(peak_chg, "charge")
    if peak_chg_q is not None and not peak_chg.empty and "voltage" in peak_chg.columns:
        v = pd.to_numeric(peak_chg["voltage"], errors="coerce").to_numpy(dtype=float)
        row.update(peaks_to_columns("chg", find_dqdv_peaks(v, peak_chg_q, config=dqdv_cfg), "dqdv"))
        row.update(peaks_to_columns("chg", find_dvdq_peaks(peak_chg_q, v, config=dqdv_cfg), "dvdq"))

    dchg_q = _capacity_series(dchg, "discharge")
    dchg_v = None
    if dchg_q is not None and not dchg.empty and "voltage" in dchg.columns:
        dchg_v = pd.to_numeric(dchg["voltage"], errors="coerce").to_numpy(dtype=float)
    peak_dchg_q = _capacity_series(peak_dchg, "discharge")
    if peak_dchg_q is not None and "voltage" in peak_dchg.columns:
        peak_v = pd.to_numeric(peak_dchg["voltage"], errors="coerce").to_numpy(dtype=float)
        row.update(peaks_to_columns("dchg", find_dqdv_peaks(peak_v, peak_dchg_q, config=dqdv_cfg), "dqdv"))
        row.update(peaks_to_columns("dchg", find_dvdq_peaks(peak_dchg_q, peak_v, config=dqdv_cfg), "dvdq"))

    # SOC bands, V_avg, energy, hysteresis, plateau, IC area, IR proxy, cliff/margin
    row.update(
        extract_shape_indicators(
            chg, dchg, chg_q, dchg_q,
            config=dqdv_cfg,
            dchg_v_cutoff=row.get("dchg_V_cutoff"),
        )
    )
    if dchg_q is not None and dchg_v is not None:
        row.update(extract_absolute_dvdq_indicators(dchg_q, dchg_v, config=dqdv_cfg))
    row.update(
        discharge_band_capacity(
            dchg,
            config=cfg.band_capacity or BandCapacityConfig(),
        )
    )

    # §4.3 energy efficiency / coulombic inefficiency
    e_chg = row.get("chg_E")
    e_dchg = row.get("dchg_E")
    ce = row.get("CE")
    if e_chg and e_dchg and e_chg > 0:
        row["EE"] = float(e_dchg) / float(e_chg)
        row["dE"] = float(e_chg) - float(e_dchg)
        if ce and ce > 0:
            row["VE"] = row["EE"] / (float(ce) / 100.0)
        else:
            row["VE"] = None
    else:
        row["EE"] = None
        row["VE"] = None
        row["dE"] = None
    # cycle duration for CI_per_hour
    dur_h = None
    if "time" in cycle_df.columns:
        tt = pd.to_numeric(cycle_df["time"], errors="coerce")
        if tt.notna().any():
            dur_h = float(tt.max() - tt.min()) / 3600.0
    row["cycle_duration_h"] = dur_h
    # Coulombic inefficiency per hour. The bare inefficiency (100 - CE) is not
    # emitted: it is an exact affine transform of CE, so it carried no extra
    # information while making CE count twice wherever both were consumed.
    if ce is not None and np.isfinite(ce) and dur_h and dur_h > 0:
        row["CI_per_hour"] = (100.0 - float(ce)) / dur_h
    else:
        row["CI_per_hour"] = None

    # Rest relaxation time constants
    row["EoC_restV_tau"] = None
    row["EoD_restV_tau"] = None
    for after_leg, rest_df in iter_rest_periods(
        cycle_df,
        charge_text=cfg.charge_step,
        discharge_text=cfg.discharge_step,
        rest_text=cfg.rest_labels,
        rest_current_max=cfg.rest_current_max,
    ):
        tau = fit_rest_tau(rest_df)
        if after_leg == "charge":
            row["EoC_restV_tau"] = tau
        elif after_leg == "discharge":
            row["EoD_restV_tau"] = tau

    # Discharge V(Q_norm) fingerprint for later DTW vs baseline
    if dchg_q is not None and not dchg.empty and "voltage" in dchg.columns:
        dv = pd.to_numeric(dchg["voltage"], errors="coerce").to_numpy(dtype=float)
        row["_dchg_vq_norm"] = vq_norm_curve(dchg_q, dv)
    else:
        row["_dchg_vq_norm"] = None

    return row


def apply_lges_delta_features(table: pd.DataFrame, *, baseline_cycle: int = 1) -> pd.DataFrame:
    """Add delta_* and *_inc columns vs baseline cycle (per cell_id + file)."""
    if table.empty:
        return table
    out = table.copy()
    group_cols = [c for c in ("cell_id", "file") if c in out.columns]
    if not group_cols:
        group_cols = ["file"] if "file" in out.columns else []

    for col in DELTA_ABS_COLS:
        out[f"delta_{col}"] = None
    for col in DELTA_PCT_COLS:
        out[f"{col}_inc"] = None
    out["delta_chgCapa_CCratio"] = None
    out["SoHQ"] = None
    out["CE_rev"] = None
    for col in ("dSoHQ_dN", "d2SoHQ", "EoC_dchgR_10s_growth_50", "CE_local_20"):
        if col not in out.columns:
            out[col] = None

    if not group_cols:
        groups = [("__all__", out)]
    else:
        groups = [(name, grp) for name, grp in out.groupby(group_cols, sort=False)]

    for _, grp in groups:
        idx = grp.index
        base_rows = grp[grp["cycle"] == baseline_cycle]
        if base_rows.empty:
            base_rows = grp.sort_values("cycle").head(1)
        base = base_rows.iloc[0]

        dchg_base = base.get("dchgCapa")
        if dchg_base and np.isfinite(dchg_base) and dchg_base > 0:
            out.loc[idx, "SoHQ"] = grp["dchgCapa"] / dchg_base * 100.0

        for col in DELTA_ABS_COLS:
            b = base.get(col)
            if b is None or not np.isfinite(b):
                continue
            out.loc[idx, f"delta_{col}"] = grp[col] - b

        if base.get("chgCapa_CCratio") is not None and np.isfinite(base.get("chgCapa_CCratio")):
            out.loc[idx, "delta_chgCapa_CCratio"] = (
                grp["chgCapa_CCratio"] - base["chgCapa_CCratio"]
            )

        for col in DELTA_PCT_COLS:
            b = base.get(col)
            if b is None or not np.isfinite(b) or abs(b) < 1e-12:
                continue
            out.loc[idx, f"{col}_inc"] = (grp[col] - b) / abs(b) * 100.0

        sorted_grp = grp.sort_values("cycle")
        cycles = sorted_grp["cycle"].to_numpy()
        chg_caps = sorted_grp["chgCapa"].to_numpy()
        dchg_caps = sorted_grp["dchgCapa"].to_numpy()
        ce_rev = np.full(len(cycles), np.nan)
        for i in range(len(cycles) - 1):
            if (
                dchg_caps[i] is not None
                and chg_caps[i + 1] is not None
                and np.isfinite(dchg_caps[i])
                and np.isfinite(chg_caps[i + 1])
                and dchg_caps[i] > 0
            ):
                ce_rev[i] = chg_caps[i + 1] / dchg_caps[i] * 100.0
        out.loc[sorted_grp.index, "CE_rev"] = ce_rev

        # Rolling trajectory features (need SoHQ first)
        sohq = pd.to_numeric(out.loc[sorted_grp.index, "SoHQ"], errors="coerce").to_numpy(dtype=float)
        cyc_x = pd.to_numeric(sorted_grp["cycle"], errors="coerce").to_numpy(dtype=float)
        dsoh = rolling_slope(sohq, cyc_x, window=21)
        out.loc[sorted_grp.index, "dSoHQ_dN"] = dsoh
        # second derivative on the slope series vs cycle
        d2 = np.full_like(dsoh, np.nan)
        ok = np.isfinite(dsoh) & np.isfinite(cyc_x)
        if ok.sum() >= 5:
            d2_full = np.gradient(np.where(ok, dsoh, np.nan), cyc_x)
            d2 = d2_full
        out.loc[sorted_grp.index, "d2SoHQ"] = d2

        if "EoC_dchgR_10s" in sorted_grp.columns:
            r10 = pd.to_numeric(sorted_grp["EoC_dchgR_10s"], errors="coerce").to_numpy(dtype=float)
            out.loc[sorted_grp.index, "EoC_dchgR_10s_growth_50"] = rolling_slope(
                r10, cyc_x, window=50,
            )
        if "CE" in sorted_grp.columns:
            ce = pd.to_numeric(sorted_grp["CE"], errors="coerce").to_numpy(dtype=float)
            ce_local = pd.Series(ce).rolling(window=20, min_periods=8).mean().to_numpy(dtype=float)
            out.loc[sorted_grp.index, "CE_local_20"] = ce_local

    return out


def extract_lges_features_table(
    df: pd.DataFrame,
    *,
    cycles: Iterable[int] | None = None,
    filepath: str = "",
    config: LgesExtractConfig | None = None,
    raw_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if "cycle" not in df.columns:
        raise ValueError("DataFrame must have logical 'cycle' column")
    cfg = config or LgesExtractConfig()
    if cfg.cell_id is None and filepath:
        cfg.cell_id = Path(filepath).stem

    cycle_list = (
        list(cycles)
        if cycles is not None
        else sorted(df["cycle"].dropna().unique().astype(int))
    )
    rows: list[dict] = []
    baseline_curve = None
    for cyc in cycle_list:
        row = extract_lges_cycle_row(
            df, int(cyc), config=cfg, filepath=filepath, raw_df=raw_df,
        )
        if not row:
            continue
        curve = row.pop("_dchg_vq_norm", None)
        if curve is not None:
            if baseline_curve is None and int(cyc) == int(cfg.baseline_cycle):
                baseline_curve = curve
                row["dchg_shape_DTW"] = 0.0
            elif baseline_curve is None:
                baseline_curve = curve
                row["dchg_shape_DTW"] = 0.0
            else:
                row["dchg_shape_DTW"] = dtw_distance(baseline_curve, curve)
        else:
            row["dchg_shape_DTW"] = None
        rows.append(row)
    table = pd.DataFrame(rows)
    if table.empty:
        return table

    # §9.1 ASSB enrichment (DCIR / Q_relax / quality / RCF)
    if cfg.enrich_assb and raw_df is not None:
        from cyclediag.features.enrich_assb import enrich_feature_table

        table, enrich_meta = enrich_feature_table(
            table,
            raw_df,
            rest_current_max=cfg.resolved_rest_current_max(),
            expected_pulse_current=cfg.resolved_pulse_current(),
            protocol_meta=cfg.protocol_meta(),
        )
        if cfg.auto_baseline and enrich_meta.get("baseline_cycle_auto"):
            cfg.baseline_cycle = int(enrich_meta["baseline_cycle_auto"])

    table = apply_lges_delta_features(table, baseline_cycle=cfg.baseline_cycle)

    # Protocol flags for the indicator-scoring track (routine vs RPT/DC-IR).
    # Causal diagnosis remains a separate optional step below.
    if raw_df is not None and not table.empty:
        from cyclediag.models.indicator_scoring import attach_protocol_flags

        table = attach_protocol_flags(table, raw_df)

    if cfg.with_diagnosis:
        from cyclediag.diagnosis import diagnose_feature_table

        table = diagnose_feature_table(
            table,
            config_path=cfg.diagnosis_config_path,
            baseline_cycle=cfg.baseline_cycle,
        )
    return table
