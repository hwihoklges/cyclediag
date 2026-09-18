"""Extra LGES cycle indicators (SOC bands, shape, hysteresis, plateau)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cyclediag.features.dqdv_peaks import (
    DqdvPeakConfig,
    compute_dqdv,
    compute_dvdq,
    dvdq_intensity_at_soc,
)

_RGAS = 8.314  # J/mol/K
_EA_DEFAULT = 20_000.0  # J/mol — mild kinetic correction
_TREF = 298.15  # 25 °C


def dvdq_intensity_soc_band(
    q: np.ndarray,
    v: np.ndarray,
    *,
    soc_lo: float,
    soc_hi: float,
    discharge: bool = True,
    config: DqdvPeakConfig | None = None,
    use_abs: bool = True,
) -> dict:
    """Mean |dV/dQ| over SOC ∈ [soc_lo, soc_hi]."""
    empty = {"intensity": None, "Q": None, "SOC": None, "n": 0}
    qx, dvdq = compute_dvdq(q, v, config)
    if len(qx) < 5:
        return empty
    qmin, qmax = float(np.nanmin(qx)), float(np.nanmax(qx))
    if not np.isfinite(qmin) or not np.isfinite(qmax) or qmax <= qmin:
        return empty
    q_norm = (qx - qmin) / (qmax - qmin)
    soc = (1.0 - q_norm) if discharge else q_norm
    lo, hi = sorted((float(soc_lo), float(soc_hi)))
    lo = max(0.0, min(1.0, lo))
    hi = max(0.0, min(1.0, hi))
    mask = (soc >= lo) & (soc <= hi) & np.isfinite(dvdq)
    if not mask.any():
        return empty
    vals = dvdq[mask]
    mid = 0.5 * (lo + hi)
    i_loc = int(np.nanargmin(np.abs(soc[mask] - mid)))
    idx = np.flatnonzero(mask)[i_loc]
    return {
        "intensity": float(np.nanmean(np.abs(vals) if use_abs else vals)),
        "Q": float(qx[idx]),
        "SOC": float(soc[idx]),
        "n": int(mask.sum()),
    }


def capacity_weighted_v_avg(v: np.ndarray, q: np.ndarray) -> float | None:
    v = np.asarray(v, dtype=float)
    q = np.asarray(q, dtype=float)
    if v.shape != q.shape or v.ndim != 1 or len(v) < 2:
        return None
    if not (np.isfinite(v).all() and np.isfinite(q).all()):
        return None
    dq = np.diff(q)
    if np.any(dq < 0) or not np.any(dq > 0):
        return None
    # trapezoid weights on midpoints
    v_mid = 0.5 * (v[1:] + v[:-1])
    w = np.abs(dq)
    s = float(np.sum(w))
    if s <= 0:
        return float(np.nanmean(v))
    return float(np.sum(v_mid * w) / s)


def energy_wh(v: np.ndarray, q: np.ndarray, *, q_is_mah: bool = False) -> float | None:
    """Return Wh from a finite monotonic leg counter (default Ah; explicit mAh)."""
    v = np.asarray(v, dtype=float)
    q = np.asarray(q, dtype=float)
    if capacity_weighted_v_avg(v, q) is None:
        return None
    e = float(np.trapezoid(v, q) if hasattr(np, "trapezoid") else np.trapz(v, q))
    return e / 1000.0 if q_is_mah else e


def ir_drop_proxy(seg: pd.DataFrame, *, n_points: int = 5) -> float | None:
    """Early-leg |ΔV| over first n samples (ohmic-ish proxy)."""
    if seg is None or seg.empty or "voltage" not in seg.columns:
        return None
    v = pd.to_numeric(seg["voltage"], errors="coerce").dropna()
    if len(v) < max(3, n_points):
        return None
    v0 = float(v.iloc[0])
    v1 = float(v.iloc[min(n_points, len(v) - 1)])
    if not np.isfinite(v0) or not np.isfinite(v1):
        return None
    return abs(v1 - v0)


def hysteresis_metrics(
    chg_q: np.ndarray,
    chg_v: np.ndarray,
    dchg_q: np.ndarray,
    dchg_v: np.ndarray,
    *,
    n_grid: int = 500,
) -> dict:
    """Charge vs discharge V(Q_norm) loop area, max ΔV, and SOC-band areas (§5.11)."""
    empty = {
        "hyst_area": None,
        "hyst_max_dV": None,
        "hyst_area_low": None,
        "hyst_area_mid": None,
        "hyst_area_high": None,
        "hyst_frac_low": None,
        "hyst_frac_high": None,
        "hyst_max_dV_low": None,
        "hyst_max_dV_mid": None,
        "hyst_max_dV_high": None,
    }

    def _norm_vq(q, v):
        q = np.asarray(q, dtype=float)
        v = np.asarray(v, dtype=float)
        m = np.isfinite(q) & np.isfinite(v)
        q, v = q[m], v[m]
        if len(q) < 5:
            return None, None
        order = np.argsort(q)
        q, v = q[order], v[order]
        qmin, qmax = float(q[0]), float(q[-1])
        if qmax <= qmin:
            return None, None
        qn = (q - qmin) / (qmax - qmin)
        _, uid = np.unique(qn, return_index=True)
        return qn[uid], v[uid]

    qc, vc = _norm_vq(chg_q, chg_v)
    qd, vd = _norm_vq(dchg_q, dchg_v)
    if qc is None or qd is None:
        return empty
    # trim ends 2% per roadmap
    grid = np.linspace(0.02, 0.98, n_grid)
    try:
        vc_i = np.interp(grid, qc, vc)
        vd_i = np.interp(grid, qd, vd)
    except Exception:
        return empty
    dv = vc_i - vd_i
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    area = float(trap(np.abs(dv), grid))
    out = {
        "hyst_area": area,
        "hyst_max_dV": float(np.nanmax(np.abs(dv))),
        "hyst_area_low": None,
        "hyst_area_mid": None,
        "hyst_area_high": None,
        "hyst_frac_low": None,
        "hyst_frac_high": None,
        "hyst_max_dV_low": None,
        "hyst_max_dV_mid": None,
        "hyst_max_dV_high": None,
    }
    bands = {"low": (0.00, 0.20), "mid": (0.20, 0.80), "high": (0.80, 1.00)}
    for name, (a, b) in bands.items():
        m = (grid >= a) & (grid <= b)
        if m.sum() < 3:
            continue
        g, d = grid[m], dv[m]
        out[f"hyst_area_{name}"] = float(trap(np.abs(d), g))
        out[f"hyst_max_dV_{name}"] = float(np.nanmax(np.abs(d)))
    if area and area > 1e-15:
        if out["hyst_area_low"] is not None:
            out["hyst_frac_low"] = out["hyst_area_low"] / area
        if out["hyst_area_high"] is not None:
            out["hyst_frac_high"] = out["hyst_area_high"] / area
    return out


def plateau_metrics(
    q: np.ndarray,
    v: np.ndarray,
    *,
    discharge: bool = True,
    config: DqdvPeakConfig | None = None,
    flat_percentile: float = 30.0,
) -> dict:
    """Longest low-|dV/dQ| region → plateau_V, plateau_width (Q units)."""
    empty = {"plateau_V": None, "plateau_width": None}
    qx, dvdq = compute_dvdq(q, v, config)
    if len(qx) < 10:
        return empty
    mag = np.abs(dvdq)
    finite = np.isfinite(mag) & np.isfinite(qx)
    if finite.sum() < 10:
        return empty
    thr = float(np.nanpercentile(mag[finite], flat_percentile))
    if not np.isfinite(thr) or thr <= 0:
        return empty
    flat = finite & (mag <= thr)
    # longest True run
    best_a = best_b = -1
    a = None
    for i, flag in enumerate(flat):
        if flag and a is None:
            a = i
        elif not flag and a is not None:
            if best_a < 0 or (i - a) > (best_b - best_a):
                best_a, best_b = a, i
            a = None
    if a is not None and (best_a < 0 or (len(flat) - a) > (best_b - best_a)):
        best_a, best_b = a, len(flat)
    if best_a < 0 or best_b - best_a < 3:
        return empty
    sl = slice(best_a, best_b)
    v = np.asarray(v, dtype=float)
    q = np.asarray(q, dtype=float)
    m = np.isfinite(v) & np.isfinite(q)
    if m.sum() < 5:
        return empty
    order = np.argsort(q[m])
    qq, vv = q[m][order], v[m][order]
    _, uid = np.unique(qq, return_index=True)
    v_plat = np.interp(qx[sl], qq[uid], vv[uid])
    return {
        "plateau_V": float(np.nanmean(v_plat)),
        "plateau_width": float(np.nanmax(qx[sl]) - np.nanmin(qx[sl])),
    }


def dqdv_area_sum(
    v: np.ndarray,
    q: np.ndarray,
    *,
    config: DqdvPeakConfig | None = None,
) -> float | None:
    """∫ |dQ/dV| dV over the leg (IC area proxy)."""
    vx, dqdv = compute_dqdv(v, q, config)
    if len(vx) < 5:
        return None
    order = np.argsort(vx)
    vx, dqdv = vx[order], dqdv[order]
    return float(
        np.trapezoid(np.abs(dqdv), vx) if hasattr(np, "trapezoid") else np.trapz(np.abs(dqdv), vx)
    )


def correct_r_to_25c(r_mohm: float | None, temp_c: float | None, *, ea: float = _EA_DEFAULT) -> float | None:
    """Arrhenius map R(T) → equivalent at 25 °C."""
    if r_mohm is None or temp_c is None:
        return None
    if not np.isfinite(r_mohm) or not np.isfinite(temp_c):
        return None
    t_k = float(temp_c) + 273.15
    if t_k <= 0:
        return None
    return float(r_mohm * np.exp(-ea / _RGAS * (1.0 / t_k - 1.0 / _TREF)))


def safe_ratio(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    if not np.isfinite(a) or not np.isfinite(b) or abs(b) < 1e-15:
        return None
    return float(a / b)


def safe_diff(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    if not np.isfinite(a) or not np.isfinite(b):
        return None
    return float(a - b)


def dvdq_at_q_from_end(
    q: np.ndarray,
    v: np.ndarray,
    *,
    ah_from_end: float,
    config: DqdvPeakConfig | None = None,
) -> dict[str, float | None]:
    """|dV/dQ| at absolute Q = Qmax − ah_from_end (no SOC normalization)."""
    qx, dvdq = compute_dvdq(q, v, config)
    empty = {"intensity": None, "Q": None, "n": 0}
    if len(qx) < 5:
        return empty
    qmax = float(np.nanmax(qx))
    q_target = qmax - float(ah_from_end)
    if q_target < float(np.nanmin(qx)):
        return empty
    mag = np.abs(dvdq)
    fin = np.isfinite(qx) & np.isfinite(mag)
    if fin.sum() < 3:
        return empty
    qx, mag = qx[fin], mag[fin]
    order = np.argsort(qx)
    qx, mag = qx[order], mag[order]
    intensity = float(np.interp(q_target, qx, mag))
    return {"intensity": intensity, "Q": q_target, "n": int(len(qx))}


def soc0_cliff_width_abs(
    q: np.ndarray,
    v: np.ndarray,
    *,
    ah_window_from_end: float = 10.0,
    thr_factor: float = 2.0,
    config: DqdvPeakConfig | None = None,
) -> float | None:
    """Cliff width in absolute Ah within last ``ah_window_from_end`` Ah of discharge."""
    qx, dvdq = compute_dvdq(q, v, config)
    if len(qx) < 10:
        return None
    qmax = float(np.nanmax(qx))
    qmin = float(np.nanmin(qx))
    mag = np.abs(dvdq)
    mid_lo = qmin + 0.4 * (qmax - qmin)
    mid_hi = qmin + 0.6 * (qmax - qmin)
    mid = (qx >= mid_lo) & (qx <= mid_hi) & np.isfinite(mag)
    if mid.sum() < 5:
        return None
    thr = float(np.nanmedian(mag[mid])) * thr_factor
    if not np.isfinite(thr) or thr <= 0:
        return None
    end = (qx >= qmax - ah_window_from_end) & np.isfinite(mag) & (mag >= thr)
    if not end.any():
        return 0.0
    return float(np.nanmax(qx[end]) - np.nanmin(qx[end]))


def extract_absolute_dvdq_indicators(
    dchg_q: np.ndarray,
    dchg_v: np.ndarray,
    *,
    config: DqdvPeakConfig | None = None,
) -> dict:
    """Absolute-Ah dV/dQ landmarks + cliff metrics (additive; does not replace SOC-norm cols)."""
    from cyclediag.features.cliff_metrics import compute_cliff_metrics

    out: dict = {}
    if dchg_q is None or dchg_v is None or len(dchg_q) < 10:
        return out
    for ah in (2.0, 5.0, 10.0):
        s = dvdq_at_q_from_end(dchg_q, dchg_v, ah_from_end=ah, config=config)
        key = f"dchg_dVdQ_at_Qabs_{int(ah)}" if ah == int(ah) else f"dchg_dVdQ_at_Qabs_{ah}"
        out[key] = s.get("intensity")
    out["dchg_dVdQ_SOC0_cliff_width_abs"] = soc0_cliff_width_abs(
        dchg_q, dchg_v, config=config,
    )
    cliff = compute_cliff_metrics(dchg_q, dchg_v, config=config)
    for k, v in cliff.items():
        out[f"dchg_{k}"] = v
    return out


def soc0_cliff_width(
    q: np.ndarray,
    v: np.ndarray,
    *,
    soc_search: float = 0.15,
    thr_factor: float = 2.0,
    config: DqdvPeakConfig | None = None,
) -> float | None:
    """Q-width near SOC0 where |dV/dQ| exceeds thr_factor × mid-SOC median."""
    qx, dvdq = compute_dvdq(q, v, config)
    if len(qx) < 10:
        return None
    qmin, qmax = float(np.nanmin(qx)), float(np.nanmax(qx))
    if qmax <= qmin:
        return None
    soc = 1.0 - (qx - qmin) / (qmax - qmin)
    mag = np.abs(dvdq)
    mid = (soc >= 0.4) & (soc <= 0.6) & np.isfinite(mag)
    if mid.sum() < 5:
        return None
    thr = float(np.nanmedian(mag[mid])) * thr_factor
    if not np.isfinite(thr) or thr <= 0:
        return None
    end = (soc <= soc_search) & np.isfinite(mag) & (mag >= thr)
    if not end.any():
        return 0.0
    # contiguous from SOC=0 upward: take max Q - min Q among end mask
    return float(np.nanmax(qx[end]) - np.nanmin(qx[end]))


def _rest_time_seconds(rest_df: pd.DataFrame, n: int) -> np.ndarray:
    """Relative rest time (s). Prefer step_time / StepTime_sec*, then total time."""
    candidates: list[str] = []
    for col in rest_df.columns:
        key = str(col).split("(")[0].replace(" ", "").replace("_", "").lower()
        if key in ("steptime", "steptimesec", "step_time") or "steptime" in key:
            candidates.insert(0, col)  # prefer step-local clocks
        elif key in ("time", "totaltime", "totaltimesec") or key.endswith("totaltimesec"):
            candidates.append(col)
    # de-dupe preserving order
    seen: set[str] = set()
    ordered: list[str] = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    for col in ordered:
        t = pd.to_numeric(rest_df[col], errors="coerce").to_numpy(dtype=float)
        if np.isfinite(t).sum() >= max(2, n // 2):
            t0 = float(t[np.isfinite(t)][0])
            return t - t0
    return np.arange(n, dtype=float)


def fit_rest_tau(rest_df: pd.DataFrame, *, min_points: int = 8) -> float | None:
    """Exponential rest relaxation time constant τ (seconds).

    Model: V(t) ≈ V∞ + A·exp(−t/τ). V∞ = mean of last 20% samples.
    """
    if rest_df is None or rest_df.empty or "voltage" not in rest_df.columns:
        return None
    v = pd.to_numeric(rest_df["voltage"], errors="coerce").to_numpy(dtype=float)
    t = _rest_time_seconds(rest_df, len(v))
    mask = np.isfinite(v) & np.isfinite(t)
    v, t = v[mask], t[mask]
    if len(v) < min_points:
        return None
    order = np.argsort(t)
    t, v = t[order], v[order]
    n_tail = max(3, len(v) // 5)
    v_inf = float(np.nanmean(v[-n_tail:]))
    resid = v - v_inf
    # use points clearly away from asymptote
    use = np.abs(resid) > 1e-4 * max(1.0, abs(v[0] - v_inf))
    use[: max(1, len(use) // 10)] = True  # keep early points
    use[-n_tail:] = False
    if use.sum() < 5:
        return None
    y = np.log(np.abs(resid[use]))
    x = t[use]
    if not np.all(np.isfinite(y)):
        return None
    # y = log|A| - t/τ  → slope = -1/τ
    coef = np.polyfit(x, y, 1)
    slope = float(coef[0])
    if slope >= -1e-9:
        return None
    tau = -1.0 / slope
    if not np.isfinite(tau) or tau <= 0 or tau > 1e6:
        return None
    return float(tau)


def vq_norm_curve(
    q: np.ndarray,
    v: np.ndarray,
    *,
    n_grid: int = 128,
) -> np.ndarray | None:
    """Resample V onto uniform Q_norm ∈ [0,1]."""
    q = np.asarray(q, dtype=float)
    v = np.asarray(v, dtype=float)
    m = np.isfinite(q) & np.isfinite(v)
    q, v = q[m], v[m]
    if len(q) < 5:
        return None
    order = np.argsort(q)
    q, v = q[order], v[order]
    qmin, qmax = float(q[0]), float(q[-1])
    if qmax <= qmin:
        return None
    qn = (q - qmin) / (qmax - qmin)
    _, uid = np.unique(qn, return_index=True)
    grid = np.linspace(0.0, 1.0, n_grid)
    return np.interp(grid, qn[uid], v[uid])


def dtw_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Classic DTW distance (O(n²)); a,b 1-D same-ish length."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n, m = len(a), len(b)
    inf = 1e30
    dp = np.full((n + 1, m + 1), inf)
    dp[0, 0] = 0.0
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            cost = abs(ai - b[j - 1])
            dp[i, j] = cost + min(dp[i - 1, j], dp[i, j - 1], dp[i - 1, j - 1])
    return float(dp[n, m] / max(n, m))


def rolling_slope(y: np.ndarray, x: np.ndarray, window: int) -> np.ndarray:
    """Centered-ish trailing slope over ``window`` points (%/unit-x uses raw units)."""
    n = len(y)
    out = np.full(n, np.nan)
    w = max(3, int(window))
    for i in range(w - 1, n):
        sl = slice(i - w + 1, i + 1)
        xx = x[sl]
        yy = y[sl]
        ok = np.isfinite(xx) & np.isfinite(yy)
        if ok.sum() < max(3, w // 2):
            continue
        out[i] = float(np.polyfit(xx[ok], yy[ok], 1)[0])
    return out


def extract_shape_indicators(
    chg: pd.DataFrame,
    dchg: pd.DataFrame,
    chg_q: np.ndarray | None,
    dchg_q: np.ndarray | None,
    *,
    config: DqdvPeakConfig | None = None,
    dchg_v_cutoff: float | None = None,
) -> dict:
    """Bundle of Tier1/2 shape indicators for one cycle."""
    out: dict = {}
    chg_v = (
        pd.to_numeric(chg["voltage"], errors="coerce").to_numpy(dtype=float)
        if chg is not None and not chg.empty and "voltage" in chg.columns
        else None
    )
    dchg_v = (
        pd.to_numeric(dchg["voltage"], errors="coerce").to_numpy(dtype=float)
        if dchg is not None and not dchg.empty and "voltage" in dchg.columns
        else None
    )

    if chg_q is not None and chg_v is not None:
        out["chg_V_avg"] = capacity_weighted_v_avg(chg_v, chg_q)
        out["chg_E"] = energy_wh(chg_v, chg_q)
        out["chg_ir_drop_proxy"] = ir_drop_proxy(chg)
        out["chg_dQdV_area_sum"] = dqdv_area_sum(chg_v, chg_q, config=config)
        plat = plateau_metrics(chg_q, chg_v, discharge=False, config=config)
        out["chg_plateau_V"] = plat["plateau_V"]
        out["chg_plateau_width"] = plat["plateau_width"]
        s100 = dvdq_intensity_at_soc(
            chg_q, chg_v, soc_target=0.98, soc_window=0.02,
            discharge=False, config=config,
        )
        out["chg_dVdQ_SOC100"] = s100["intensity"]

    if dchg_q is not None and dchg_v is not None:
        out["dchg_V_avg"] = capacity_weighted_v_avg(dchg_v, dchg_q)
        out["dchg_E"] = energy_wh(dchg_v, dchg_q)
        out["dchg_ir_drop_proxy"] = ir_drop_proxy(dchg)
        out["dchg_dQdV_area_sum"] = dqdv_area_sum(dchg_v, dchg_q, config=config)
        plat = plateau_metrics(dchg_q, dchg_v, discharge=True, config=config)
        out["dchg_plateau_V"] = plat["plateau_V"]
        out["dchg_plateau_width"] = plat["plateau_width"]
        for label, tgt, win in (
            ("SOC0", 0.0, 0.02),
            ("SOC5", 0.05, 0.02),
            ("SOC10", 0.10, 0.02),
        ):
            s = dvdq_intensity_at_soc(
                dchg_q, dchg_v, soc_target=tgt, soc_window=win,
                discharge=True, config=config,
            )
            out[f"dchg_dVdQ_{label}"] = s["intensity"]
            if label == "SOC0":
                out["dchg_dVdQ_SOC0_Q"] = s["Q"]
        mid = dvdq_intensity_soc_band(
            dchg_q, dchg_v, soc_lo=0.40, soc_hi=0.60,
            discharge=True, config=config,
        )
        out["dchg_dVdQ_SOCmid"] = mid["intensity"]
        out["dchg_dVdQ_SOC0_cliff_width"] = soc0_cliff_width(
            dchg_q, dchg_v, config=config,
        )
        out["dchg_dVdQ_SOC0_to_mid_ratio"] = safe_ratio(
            out.get("dchg_dVdQ_SOC0"), out.get("dchg_dVdQ_SOCmid"),
        )
        # cutoff margin: V@SOC10 − lower cutoff (room before hitting wall)
        if dchg_v_cutoff is not None and np.isfinite(dchg_v_cutoff):
            qmin, qmax = float(np.nanmin(dchg_q)), float(np.nanmax(dchg_q))
            if qmax > qmin:
                q_at = qmin + 0.90 * (qmax - qmin)  # SOC≈0.10 on discharge
                order = np.argsort(dchg_q)
                v_at = float(np.interp(q_at, dchg_q[order], dchg_v[order]))
                out["dchg_V_cutoff_margin"] = safe_diff(v_at, float(dchg_v_cutoff))
            else:
                out["dchg_V_cutoff_margin"] = None
        else:
            out["dchg_V_cutoff_margin"] = None

    if (
        chg_q is not None and dchg_q is not None
        and chg_v is not None and dchg_v is not None
    ):
        hyst = hysteresis_metrics(chg_q, chg_v, dchg_q, dchg_v)
        out.update(hyst)

    return out
