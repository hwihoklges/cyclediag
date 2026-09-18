"""dQ/dV and dV/dQ peak detection — Q-axis 500-pt grid (pne_studio parity)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter

from .lges_catalog import MAX_DQDV_PEAKS

try:
    from pne_studio2.core.dqdv_interp import build_dqdv_from_segment
except ImportError:
    try:
        from pne_studio.core.dqdv_interp import build_dqdv_from_segment
    except ImportError:
        from cyclediag.features._dqdv_interp import build_dqdv_from_segment


@dataclass
class DqdvPeakConfig:
    """Match pne_studio default: Q interpolation, 500 points, smooth-then-diff."""

    n_interp: int = 500
    interp_axis: str = "Q"
    deriv_mode: str = "smooth_then_diff"
    sg_window: int = 21
    sg_poly: int = 3
    prominence_frac: float = 0.02
    min_distance_frac: float = 0.04
    min_width_points: int = 5
    mad_prominence_factor: float = 4.0
    spike_ratio_max: float = 2.5
    merge_v_sep_v: float = 0.012


DEFAULT_DQDV_PEAK_CONFIG = DqdvPeakConfig()


# Voltage bands for shoulder-aware peak assign (Ch022 SJ900 charge/discharge).
DEFAULT_CHARGE_VOLTAGE_BANDS: tuple[tuple[float, float, str], ...] = (
    (3.48, 3.62, "P1_low"),
    (3.60, 3.78, "P2_shoulder"),
    (3.78, 3.94, "P3_main"),
    (3.94, 4.12, "P4_high"),
)
DEFAULT_DISCHARGE_VOLTAGE_BANDS: tuple[tuple[float, float, str], ...] = (
    (3.05, 3.28, "P1_low"),
    (3.55, 3.78, "P2_mid"),
    (3.82, 4.00, "P3_high"),
)


def find_dqdv_peaks_banded_prepared(
    vx: np.ndarray,
    y_smooth: np.ndarray,
    bands: tuple[tuple[float, float, str], ...],
    *,
    min_band_height_frac: float = 0.12,
) -> list[dict]:
    """One peak per voltage band using pre-interpolated SG-smoothed dQ/dV."""
    if len(vx) < 5:
        return []

    y_abs = np.abs(y_smooth)
    global_max = float(np.nanmax(y_abs)) if y_abs.size else 0.0
    floor = global_max * min_band_height_frac if global_max > 0 else 0.0

    peaks: list[dict] = []
    for v_min, v_max, label in bands:
        mask = (vx >= v_min) & (vx <= v_max) & np.isfinite(y_smooth)
        if not mask.any():
            continue
        idx_local = int(np.argmax(y_abs[mask]))
        indices = np.flatnonzero(mask)
        i = int(indices[idx_local])
        h = float(y_smooth[i])
        if abs(h) < floor:
            continue
        band_height_frac = abs(h) / global_max if global_max > 0 else 0.0
        assign_confidence = (
            min(1.0, band_height_frac / min_band_height_frac)
            if min_band_height_frac > 0
            else 1.0
        )
        peaks.append({
            "V": float(vx[i]),
            "H": h,
            "band": label,
            "band_v_min": v_min,
            "band_v_max": v_max,
            "band_height_frac": band_height_frac,
            "assign_confidence": assign_confidence,
        })

    peaks.sort(key=lambda p: p["V"])
    return peaks


def find_dqdv_peaks_banded(
    v: np.ndarray,
    q: np.ndarray,
    bands: tuple[tuple[float, float, str], ...],
    *,
    config: DqdvPeakConfig | None = None,
    min_band_height_frac: float = 0.12,
) -> list[dict]:
    """One peak per voltage band on SG-smoothed dQ/dV (shoulder split).

    Skips a band when its maximum |dQ/dV| is below ``min_band_height_frac`` of the
    global curve maximum (weak bump / empty band).
    """
    config = config or DqdvPeakConfig(sg_window=31)
    vx, dqdv, _, _ = prepare_dqdv_arrays(v, q, config)
    if len(vx) < 5:
        return []
    y_smooth = _smooth(dqdv, window=config.sg_window, poly=config.sg_poly)
    return find_dqdv_peaks_banded_prepared(
        vx, y_smooth, bands, min_band_height_frac=min_band_height_frac,
    )


def charge_discharge_bands(leg: str) -> tuple[tuple[float, float, str], ...]:
    if leg == "charge":
        return DEFAULT_CHARGE_VOLTAGE_BANDS
    return DEFAULT_DISCHARGE_VOLTAGE_BANDS


def _smooth(y: np.ndarray, *, window: int = 21, poly: int = 3) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    if len(y) < 5 or not np.isfinite(y).any():
        return y
    y = np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)
    if np.nanstd(y) < 1e-15:
        return y
    w = min(window, len(y) if len(y) % 2 else len(y) - 1)
    if w < poly + 2:
        w = poly + 2 if (poly + 2) % 2 else poly + 3
    if w < 5 or w > len(y):
        return y
    try:
        return savgol_filter(y, w, poly)
    except ValueError:
        return y


def _noise_mad(y: np.ndarray, y_smooth: np.ndarray) -> float:
    resid = np.asarray(y, dtype=float) - np.asarray(y_smooth, dtype=float)
    resid = resid[np.isfinite(resid)]
    if len(resid) < 3:
        return 0.0
    med = float(np.median(resid))
    mad = float(np.median(np.abs(resid - med)))
    return mad if np.isfinite(mad) else 0.0


def prepare_dqdv_arrays(
    v: np.ndarray,
    q: np.ndarray,
    config: DqdvPeakConfig | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Resample leg on Q (500 pts), return V, dQ/dV, Q, dV/dQ."""
    config = config or DEFAULT_DQDV_PEAK_CONFIG
    v = np.asarray(v, dtype=float)
    q = np.asarray(q, dtype=float)
    valid = np.isfinite(v) & np.isfinite(q)
    v, q = v[valid], q[valid]
    if len(v) < 4:
        return np.array([]), np.array([]), np.array([]), np.array([])

    raw = pd.DataFrame({"voltage": v, "capacity": q})
    axis = str(config.interp_axis).strip().upper()
    if axis not in ("V", "Q"):
        axis = "Q"
    n_pts = max(int(config.n_interp), 2)
    proc = build_dqdv_from_segment(
        raw,
        "voltage",
        "capacity",
        axis=axis,
        num_points=n_pts,
        use_interp=True,
        deriv_mode=config.deriv_mode,
        sg_window=config.sg_window,
        sg_poly=config.sg_poly,
    )
    if proc.empty or len(proc) < 5:
        return np.array([]), np.array([]), np.array([]), np.array([])

    vx = pd.to_numeric(proc["voltage"], errors="coerce").to_numpy(dtype=float)
    qx = pd.to_numeric(proc["capacity"], errors="coerce").to_numpy(dtype=float)
    dqdv = pd.to_numeric(proc["dQ/dV"], errors="coerce").to_numpy(dtype=float)
    dvdq = pd.to_numeric(proc["dV/dQ"], errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(vx) & np.isfinite(dqdv)
    return vx[mask], dqdv[mask], qx[mask], dvdq[mask]


def compute_dqdv(
    v: np.ndarray,
    q: np.ndarray,
    config: DqdvPeakConfig | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    vx, dqdv, _, _ = prepare_dqdv_arrays(v, q, config)
    return vx, dqdv


def compute_dvdq(
    q: np.ndarray,
    v: np.ndarray,
    config: DqdvPeakConfig | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    # dQ/dV's minimum voltage increment is a denominator guard. Applying it
    # to dV/dQ discards valid shallow slopes and changes the capacity domain.
    from ._dqdv_interp import (
        interpolate_dqdv_segment,
        normalize_deriv_mode,
        resolve_dqdv_interp_axis,
        smooth_state_columns,
    )

    cfg = config or DEFAULT_DQDV_PEAK_CONFIG
    q, v = np.asarray(q, dtype=float), np.asarray(v, dtype=float)
    valid = np.isfinite(q) & np.isfinite(v)
    q, v = q[valid], v[valid]
    if len(q) < 4:
        return np.array([]), np.array([])
    axis = resolve_dqdv_interp_axis(v, q, cfg.interp_axis)
    grid = interpolate_dqdv_segment(
        pd.DataFrame({"voltage": v, "capacity": q}), "voltage", "capacity",
        axis=axis, num_points=max(int(cfg.n_interp), 2),
    )
    if normalize_deriv_mode(cfg.deriv_mode) == "smooth_then_diff":
        grid = smooth_state_columns(
            grid, "voltage", "capacity", window=cfg.sg_window, poly=cfg.sg_poly,
        )
    qx = grid["capacity"].to_numpy(dtype=float)
    dq = np.diff(qx)
    dvdq = np.divide(np.diff(grid["voltage"].to_numpy(dtype=float)), dq,
                     out=np.full(dq.shape, np.nan), where=dq != 0)
    qx = qx[1:]
    mask = np.isfinite(qx) & np.isfinite(dvdq)
    return qx[mask], dvdq[mask]


def _robust_peak_indices(
    y: np.ndarray,
    config: DqdvPeakConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Return peak indices and SG-smoothed curve."""
    y_smooth = _smooth(y, window=config.sg_window, poly=config.sg_poly)
    y_use = np.abs(y_smooth)
    ymax = float(np.nanmax(y_use)) if y_use.size and np.nanmax(y_use) > 0 else 1e-9
    mad = _noise_mad(y, y_smooth)
    prom_floor = config.mad_prominence_factor * mad * 1.4826
    prom = max(config.prominence_frac * ymax, prom_floor, 1e-9)
    distance = max(config.min_width_points, int(config.n_interp * config.min_distance_frac))
    width = config.min_width_points

    idx, _props = find_peaks(
        y_use,
        prominence=prom,
        distance=distance,
        width=width,
    )
    if len(idx) == 0:
        return idx, y_smooth

    filtered: list[int] = []
    for i in idx:
        raw_h = abs(float(y[i])) if np.isfinite(y[i]) else 0.0
        smooth_h = abs(float(y_smooth[i])) if np.isfinite(y_smooth[i]) else 1e-12
        if smooth_h > 1e-15 and raw_h > config.spike_ratio_max * smooth_h:
            continue
        filtered.append(int(i))
    return np.asarray(filtered, dtype=int), y_smooth


def _merge_close_peaks(peaks: list[dict], min_v_sep: float) -> list[dict]:
    if len(peaks) <= 1:
        return peaks
    peaks = sorted(peaks, key=lambda p: p.get("V", p.get("Q", 0.0)))
    merged: list[dict] = []
    for pk in peaks:
        if not merged:
            merged.append(pk)
            continue
        prev = merged[-1]
        key = "V" if "V" in pk else "Q"
        if abs(float(pk[key]) - float(prev[key])) < min_v_sep:
            if abs(float(pk.get("H", 0))) >= abs(float(prev.get("H", 0))):
                merged[-1] = pk
        else:
            merged.append(pk)
    return merged


def find_dqdv_peaks_prepared(
    vx: np.ndarray,
    dqdv: np.ndarray,
    y_smooth: np.ndarray,
    max_peaks: int = MAX_DQDV_PEAKS,
    config: DqdvPeakConfig | None = None,
) -> list[dict]:
    """Peak find on pre-interpolated dQ/dV (avoids duplicate Q-grid work)."""
    config = config or DEFAULT_DQDV_PEAK_CONFIG
    if len(vx) < 5:
        return []

    idx, y_smooth = _robust_peak_indices(dqdv, config)
    if len(idx) == 0:
        return []

    prominences = np.abs(y_smooth[idx])
    order = np.argsort(prominences)[::-1][: max_peaks * 2]
    peaks: list[dict] = []
    for j in order:
        i = int(idx[j])
        peaks.append({"V": float(vx[i]), "H": float(y_smooth[i])})

    peaks = _merge_close_peaks(peaks, config.merge_v_sep_v)
    peaks.sort(key=lambda p: abs(float(p.get("H", 0))), reverse=True)
    peaks = peaks[:max_peaks]
    peaks.sort(key=lambda p: p["V"])
    return peaks


def find_dqdv_peaks(
    v: np.ndarray,
    q: np.ndarray,
    max_peaks: int = MAX_DQDV_PEAKS,
    config: DqdvPeakConfig | None = None,
) -> list[dict]:
    config = config or DEFAULT_DQDV_PEAK_CONFIG
    vx, dqdv, _, _ = prepare_dqdv_arrays(v, q, config)
    if len(vx) < 5:
        return []

    y_smooth = _smooth(dqdv, window=config.sg_window, poly=config.sg_poly)
    return find_dqdv_peaks_prepared(vx, dqdv, y_smooth, max_peaks=max_peaks, config=config)


def find_dvdq_peaks(
    q: np.ndarray,
    v: np.ndarray,
    max_peaks: int = MAX_DQDV_PEAKS,
    config: DqdvPeakConfig | None = None,
) -> list[dict]:
    config = config or DEFAULT_DQDV_PEAK_CONFIG
    qx, dvdq = compute_dvdq(q, v, config)
    if len(qx) < 5:
        return []

    y_smooth = _smooth(dvdq, window=config.sg_window, poly=config.sg_poly)
    idx, y_smooth = _robust_peak_indices(dvdq, config)
    if len(idx) == 0:
        return []

    prominences = np.abs(y_smooth[idx])

    order = np.argsort(prominences)[::-1][: max_peaks * 2]
    peaks: list[dict] = []
    for j in order:
        i = int(idx[j])
        peaks.append({"Q": float(qx[i]), "H": float(y_smooth[i])})

    peaks.sort(key=lambda p: abs(float(p.get("H", 0))), reverse=True)
    peaks = peaks[:max_peaks]
    peaks.sort(key=lambda p: p["Q"])
    return peaks


def peaks_to_columns(prefix: str, peaks: list[dict], kind: str) -> dict:
    """kind: 'dqdv' or 'dvdq'"""
    out: dict = {}
    for i in range(1, MAX_DQDV_PEAKS + 1):
        if kind == "dqdv":
            out[f"{prefix}_dQdV_peak{i}_V"] = None
            out[f"{prefix}_dQdV_peak{i}"] = None
        else:
            out[f"{prefix}_dVdQ_peak{i}_Q"] = None
            out[f"{prefix}_dVdQ_peak{i}"] = None
    for i, pk in enumerate(peaks[:MAX_DQDV_PEAKS], start=1):
        if kind == "dqdv":
            out[f"{prefix}_dQdV_peak{i}_V"] = pk.get("V")
            out[f"{prefix}_dQdV_peak{i}"] = pk.get("H")
        else:
            out[f"{prefix}_dVdQ_peak{i}_Q"] = pk.get("Q")
            out[f"{prefix}_dVdQ_peak{i}"] = pk.get("H")
    return out


def dvdq_intensity_at_soc(
    q: np.ndarray,
    v: np.ndarray,
    *,
    soc_target: float = 0.0,
    soc_window: float = 0.02,
    discharge: bool = True,
    config: DqdvPeakConfig | None = None,
    use_abs: bool = True,
) -> dict:
    """Sample dV/dQ intensity near a target SOC on one leg.

    SOC is normalized capacity along the leg:
    - discharge: SOC = 1 − Q_norm (Q↑ as SOC↓) → SOC0 ≈ end of discharge
    - charge:    SOC = Q_norm → SOC0 ≈ start of charge

    Intensity = mean(|dV/dQ|) (or signed mean) over SOC ∈ [target, target+window]
    clipped to [0, 1]. Falls back to nearest-point sample if the window is empty.
    """
    empty = {"intensity": None, "Q": None, "SOC": None, "n": 0}
    qx, dvdq = compute_dvdq(q, v, config)
    if len(qx) < 5:
        return empty

    qmin, qmax = float(np.nanmin(qx)), float(np.nanmax(qx))
    if not np.isfinite(qmin) or not np.isfinite(qmax) or qmax <= qmin:
        return empty

    q_norm = (qx - qmin) / (qmax - qmin)
    soc = (1.0 - q_norm) if discharge else q_norm

    lo = float(soc_target)
    hi = float(soc_target) + max(float(soc_window), 1e-6)
    # clamp window into [0, 1]
    lo = max(0.0, min(1.0, lo))
    hi = max(0.0, min(1.0, hi))
    if hi < lo:
        lo, hi = hi, lo

    mask = (soc >= lo) & (soc <= hi) & np.isfinite(dvdq) & np.isfinite(soc)
    if not mask.any():
        # nearest sample to target
        i = int(np.nanargmin(np.abs(soc - soc_target)))
        val = float(dvdq[i])
        return {
            "intensity": abs(val) if use_abs else val,
            "Q": float(qx[i]),
            "SOC": float(soc[i]),
            "n": 1,
        }

    vals = dvdq[mask]
    intensity = float(np.nanmean(np.abs(vals) if use_abs else vals))
    # representative Q/SOC: closest to target within window
    i_loc = int(np.nanargmin(np.abs(soc[mask] - soc_target)))
    idx = np.flatnonzero(mask)[i_loc]
    return {
        "intensity": intensity,
        "Q": float(qx[idx]),
        "SOC": float(soc[idx]),
        "n": int(mask.sum()),
    }


def dvdq_soc_columns(prefix: str, sample: dict, *, soc_label: str = "SOC0") -> dict:
    """Map ``dvdq_intensity_at_soc`` result to feature columns."""
    return {
        f"{prefix}_dVdQ_{soc_label}": sample.get("intensity"),
        f"{prefix}_dVdQ_{soc_label}_Q": sample.get("Q"),
    }
