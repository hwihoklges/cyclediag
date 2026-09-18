"""DC-IR R(t) three-component decomposition — IMPROVEMENT_ROADMAP §5.3 / §5.5."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


@dataclass
class DcirFitResult:
    R_ohmic: float | None = None
    R_ct: float | None = None
    tau_ct: float | None = None
    A_diff: float | None = None
    R_30s_total: float | None = None
    R_ohmic_frac: float | None = None
    R_ct_frac: float | None = None
    R_diff_frac: float | None = None
    dcir_fit_rmse: float | None = None
    dcir_fit_r2: float | None = None
    dcir_fit_cond: float | None = None
    dcir_fit_valid: bool = False
    R_recovery_tau1: float | None = None
    R_recovery_tau2: float | None = None
    relax_amp_ratio: float | None = None
    V_inf_est: float | None = None
    recovery_fit_r2: float | None = None
    tau_consistency_flag: bool | None = None
    pulse_current_A: float | None = None
    n_points: int = 0
    n_t_le_1s: int = 0
    flag: str = ""
    scientific_validity: dict[str, Any] = field(default_factory=lambda: {
        "model": "R_ohmic_plus_RC_plus_sqrt_time",
        "causally_identified": False,
        "validated_probability": False,
        "warnings": ["components_are_equivalent_circuit_proxies_not_unique_mechanisms",
                     "requires_stable_current_local_linearity_fixed_SOC_temperature",
                     "sqrt_time_diffusion_is_finite_window_approximation",
                     "fit_quality_does_not_validate_parameter_identifiability"],
    })


def _linreg(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return slope, intercept."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 2:
        return float("nan"), float("nan")
    A = np.vstack([x, np.ones(len(x))]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(slope), float(intercept)


def _detect_current_settle(
    t: np.ndarray,
    i: np.ndarray | None,
    *,
    settle_frac: float = 0.95,
) -> tuple[float, bool]:
    """Return (t_settle, ramp_ok). ramp_ok=False when rise time > 0.2 s."""
    if i is None or len(t) < 3:
        return 0.0, True
    i_abs = np.abs(np.asarray(i, dtype=float))
    t = np.asarray(t, dtype=float)
    active = np.isfinite(i_abs) & (i_abs > 1e-6)
    if not active.any():
        return 0.0, True
    i_target = float(np.nanmedian(i_abs[active]))
    if i_target < 1e-9:
        return 0.0, True
    settled = active & (i_abs >= settle_frac * i_target)
    if not settled.any():
        return 0.0, False
    t_settle = float(t[settled][0])
    return t_settle, t_settle <= 0.2


def fit_r_t_components(
    t: np.ndarray,
    r_mohm: np.ndarray,
    *,
    i: np.ndarray | None = None,
    refine_global: bool = True,
) -> DcirFitResult:
    """Stepwise fit R(t) = RΩ + Rct(1-exp(-t/τ)) + A√t  [mΩ]."""
    t = np.asarray(t, dtype=float)
    r = np.asarray(r_mohm, dtype=float)
    if t.shape != r.shape or t.ndim != 1:
        raise ValueError("t and r_mohm must be aligned one-dimensional traces")
    m = np.isfinite(t) & np.isfinite(r) & (t >= 0)
    if i is not None:
        i = np.asarray(i, dtype=float)
        if i.shape != t.shape:
            raise ValueError("current must align with time")
        m &= np.isfinite(i)
        i = i[m]
    t, r = t[m], r[m]
    out = DcirFitResult(n_points=int(len(t)), n_t_le_1s=int(np.sum(t <= 1.0)))
    if np.any(np.diff(t) <= 0) or np.any(r < 0):
        out.flag = "invalid_time_or_negative_resistance"
        return out
    if len(t) < 50:
        out.flag = "too_few_points"
        return out
    if out.n_t_le_1s < 8:
        out.flag = "insufficient_early_samples"
        # still try but mark invalid

    t_settle, ramp_ok = _detect_current_settle(t, i)
    if not ramp_ok:
        out.flag = (out.flag + "|current_ramp").strip("|")

    # 4a R_ohmic via √t extrapolate on [t_settle, t_settle+0.3]
    early = (t >= t_settle) & (t <= t_settle + 0.3)
    if early.sum() >= 3:
        slope_e, r_ohm = _linreg(np.sqrt(t[early]), r[early])
        out.R_ohmic = r_ohm if np.isfinite(r_ohm) else float(r[np.argmin(t)])
    else:
        out.R_ohmic = float(r[np.argmin(t)])
        out.flag = (out.flag + "|early_fallback").strip("|")

    # 4b A_diff on [10, 30]
    late = (t >= 10.0) & (t <= 30.0)
    if late.sum() >= 5:
        a_diff, _ = _linreg(np.sqrt(t[late]), r[late])
        out.A_diff = a_diff if np.isfinite(a_diff) and a_diff >= 0 else None
        if out.A_diff is None:
            out.flag = (out.flag + "|negative_or_invalid_diffusion_slope").strip("|")
    else:
        out.A_diff = None
        out.flag = (out.flag + "|late_fail").strip("|")

    r_ohm = float(out.R_ohmic or 0.0)
    a_diff = float(out.A_diff or 0.0)

    # 4c R_ct, tau on residual mid window
    mid = (t >= 0.3) & (t <= 10.0)
    resid = r - r_ohm - a_diff * np.sqrt(np.maximum(t, 0.0))
    tau_ct = None
    r_ct = None
    if mid.sum() >= 8:
        tm, rm = t[mid], resid[mid]

        def _exp_sat(x, rct, tau):
            return rct * (1.0 - np.exp(-x / np.maximum(tau, 1e-3)))

        try:
            rct0 = max(float(rm[np.argmin(np.abs(tm - 10.0))]), 1e-6)
            popt, _ = curve_fit(
                _exp_sat, tm, rm,
                p0=(rct0, 2.0),
                bounds=([0.0, 0.05], [max(5.0 * abs(r_ohm) + 1.0, rct0 * 5), 20.0]),
                maxfev=8000,
            )
            r_ct, tau_ct = float(popt[0]), float(popt[1])
        except Exception:
            r_ct = max(float(np.nanmax(rm)), 0.0)
            tau_ct = 2.0
            out.flag = (out.flag + "|ct_fit_fallback").strip("|")

    out.R_ct = r_ct
    out.tau_ct = tau_ct

    # optional global refine
    cond = None
    if refine_global and r_ct is not None and tau_ct is not None and a_diff is not None:

        def _full(x, ro, rct, tau, a):
            return ro + rct * (1.0 - np.exp(-x / np.maximum(tau, 1e-3))) + a * np.sqrt(x)

        try:
            popt, pcov = curve_fit(
                _full, t, r,
                p0=(r_ohm, r_ct, tau_ct, max(a_diff, 0.0)),
                bounds=(
                    [0.0, 0.0, 0.05, 0.0],
                    [max(r_ohm * 5, 1.0), max(r_ct * 5, 1.0), 20.0, max(a_diff * 10, 1.0)],
                ),
                maxfev=8000,
            )
            # condition number of covariance
            if pcov is not None and np.all(np.isfinite(pcov)):
                cond = float(np.linalg.cond(pcov))
            if cond is not None and cond > 1e8:
                out.flag = (out.flag + "|degenerate_keep_step").strip("|")
            else:
                out.R_ohmic, out.R_ct, out.tau_ct, out.A_diff = (float(x) for x in popt)
        except Exception:
            out.flag = (out.flag + "|refine_fail").strip("|")

    out.dcir_fit_cond = cond
    r_ohm = float(out.R_ohmic or 0.0)
    r_ct = float(out.R_ct or 0.0)
    a_diff = float(out.A_diff or 0.0)
    tau_ct = float(out.tau_ct or 2.0)

    r_fit = r_ohm + r_ct * (1.0 - np.exp(-t / max(tau_ct, 1e-3))) + a_diff * np.sqrt(t)
    resid_f = r - r_fit
    ss_res = float(np.nansum(resid_f ** 2))
    ss_tot = float(np.nansum((r - np.nanmean(r)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0
    rmse = float(np.sqrt(np.nanmean(resid_f ** 2)) / max(np.nanmean(r), 1e-9))
    out.dcir_fit_r2 = r2
    out.dcir_fit_rmse = rmse

    # R at 30s
    if np.any(t >= 29.0):
        out.R_30s_total = float(r[np.argmin(np.abs(t - 30.0))])
    else:
        out.R_30s_total = float(r[-1])
    r30 = float(out.R_30s_total or 0.0)
    if r30 > 1e-12:
        out.R_ohmic_frac = r_ohm / r30
        out.R_ct_frac = r_ct / r30
        out.R_diff_frac = (a_diff * np.sqrt(30.0)) / r30

    out.dcir_fit_valid = bool(
        rmse < 0.03
        and r2 > 0.98
        and (cond is None or cond < 1e8)
        and out.n_t_le_1s >= 8
        and ramp_ok
        and "sampling_too_sparse" not in out.flag
        and not any(tag in out.flag for tag in ("fallback", "fail", "degenerate", "invalid_diffusion"))
    )
    return out


def fit_recovery_tau(
    t: np.ndarray,
    v: np.ndarray,
    *,
    tau_ct_ref: float | None = None,
) -> dict[str, Any]:
    """Two-exp recovery fit on post-pulse rest (§5.5)."""
    out: dict[str, Any] = {
        "R_recovery_tau1": None,
        "R_recovery_tau2": None,
        "relax_amp_ratio": None,
        "V_inf_est": None,
        "recovery_fit_r2": None,
        "tau_consistency_flag": None,
    }
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    m = np.isfinite(t) & np.isfinite(v) & (t >= 0)
    t, v = t[m], v[m]
    if len(t) < 20 or float(np.nanmax(t)) < 60:
        return out

    def _f(x, vinf, b1, tau1, b2, tau2):
        return vinf - b1 * np.exp(-x / np.maximum(tau1, 0.5)) - b2 * np.exp(-x / np.maximum(tau2, 50.0))

    try:
        vinf0 = float(v[-1])
        p0 = (vinf0, float(v[0] - vinf0) * 0.5, 5.0, float(v[0] - vinf0) * 0.5, 300.0)
        popt, _ = curve_fit(
            _f, t, v, p0=p0,
            bounds=(
                [vinf0 - 1.0, -2.0, 0.5, -2.0, 50.0],
                [vinf0 + 1.0, 2.0, 50.0, 2.0, 3000.0],
            ),
            maxfev=8000,
        )
        vinf, b1, tau1, b2, tau2 = (float(x) for x in popt)
        if tau1 > tau2:
            b1, b2, tau1, tau2 = b2, b1, tau2, tau1
        yhat = _f(t, *popt)
        ss_res = float(np.nansum((v - yhat) ** 2))
        ss_tot = float(np.nansum((v - np.nanmean(v)) ** 2))
        out["V_inf_est"] = vinf
        out["R_recovery_tau1"] = tau1
        out["R_recovery_tau2"] = tau2
        amp = abs(b1) + abs(b2)
        out["relax_amp_ratio"] = (abs(b2) / amp) if amp > 1e-15 else None
        out["recovery_fit_r2"] = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else None
        if tau_ct_ref is not None and tau_ct_ref > 0 and tau1 > 0:
            out["tau_consistency_flag"] = abs(np.log10(tau_ct_ref / tau1)) <= 0.7
    except Exception:
        return out
    return out


def extract_pulse_trace(
    cycle_df: pd.DataFrame,
    *,
    rest_current_max: float = 0.5,
    pulse_duration_s: float = 30.0,
    expected_pulse_current: float | None = 70.0,
    meta: dict[str, Any] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float] | None:
    """Return (t, V, I_series, I_med, V0) for first long high-|I| pulse in cycle."""
    need = {"voltage", "current"}
    if not need.issubset(cycle_df.columns):
        return None
    v = pd.to_numeric(cycle_df["voltage"], errors="coerce").to_numpy(dtype=float)
    i = pd.to_numeric(cycle_df["current"], errors="coerce").to_numpy(dtype=float)
    if "step_time" in cycle_df.columns:
        st = pd.to_numeric(cycle_df["step_time"], errors="coerce").to_numpy(dtype=float)
    elif "time" in cycle_df.columns:
        st = pd.to_numeric(cycle_df["time"], errors="coerce").to_numpy(dtype=float)
        st = st - st[0]
    else:
        return None

    i_abs = np.abs(i)
    thr = rest_current_max
    if expected_pulse_current is not None:
        thr = max(thr, 0.5 * expected_pulse_current)

    # find rising edge into high current
    active = np.isfinite(i_abs) & (i_abs > thr)
    if not active.any():
        return None
    start = int(np.argmax(active))
    # walk until current drops or duration exceeded
    t0 = float(st[start]) if np.isfinite(st[start]) else 0.0
    end = start
    for j in range(start, len(i)):
        if not np.isfinite(st[j]):
            continue
        if float(st[j]) - t0 > pulse_duration_s + 1.0:
            break
        if j > start + 5 and i_abs[j] < rest_current_max:
            break
        end = j
    seg_len = end - start + 1
    dt_pts = np.diff(st[start : end + 1])
    dt_fin = dt_pts[np.isfinite(dt_pts) & (dt_pts > 0)]
    dt_med = float(np.median(dt_fin)) if len(dt_fin) else float("inf")
    min_pts = max(15, int(pulse_duration_s / max(dt_med, 1e-6) * 0.85))
    if meta is not None:
        meta["dt_median"] = dt_med
        meta["n_pulse_points"] = seg_len
    if dt_med > 0.6:
        if meta is not None:
            meta["flag"] = "sampling_too_sparse"
        return None
    if seg_len < min_pts:
        if meta is not None:
            meta["flag"] = "insufficient_pulse_points"
        return None

    sl = slice(start, end + 1)
    t = st[sl] - t0
    # clamp to ~30s
    keep = t <= pulse_duration_s + 0.15
    t, vv, ii = t[keep], v[sl][keep], i[sl][keep]
    i_med = float(np.nanmedian(ii))
    if expected_pulse_current and abs(i_med) < 0.5 * expected_pulse_current:
        return None
    if np.nanstd(ii[1:]) / max(abs(i_med), 1e-9) > 0.02:
        # unstable — still return but caller can flag
        pass
    v0 = float(v[start - 1]) if start > 0 and np.isfinite(v[start - 1]) else float(vv[0])
    return t, vv, ii, i_med, v0


def decompose_pulse_cycle(
    cycle_df: pd.DataFrame,
    *,
    rest_current_max: float = 0.5,
    expected_pulse_current: float | None = 70.0,
    fit_recovery: bool = True,
) -> DcirFitResult:
    """Full §5.3 (+ optional §5.5) for one DC-IR cycle DataFrame."""
    meta: dict[str, Any] = {}
    extracted = extract_pulse_trace(
        cycle_df,
        rest_current_max=rest_current_max,
        expected_pulse_current=expected_pulse_current,
        meta=meta,
    )
    if extracted is None:
        return DcirFitResult(flag=str(meta.get("flag") or "no_pulse"))
    t, v, i, i_med, v0 = extracted
    r = np.abs(v0 - v) / max(abs(i_med), 1e-9) * 1000.0  # mΩ
    fit = fit_r_t_components(t, r, i=i)
    if meta.get("flag"):
        fit.flag = (fit.flag + "|" + str(meta["flag"])).strip("|")
    fit.pulse_current_A = abs(i_med)

    if fit_recovery and "step_time" in cycle_df.columns:
        # post-pulse rest: after pulse current returns to rest
        i_abs = pd.to_numeric(cycle_df["current"], errors="coerce").abs().to_numpy()
        st = pd.to_numeric(cycle_df["step_time"], errors="coerce").to_numpy()
        v_all = pd.to_numeric(cycle_df["voltage"], errors="coerce").to_numpy()
        # find pulse end index approximately
        pulse_idx = np.where(i_abs > rest_current_max)[0]
        if len(pulse_idx):
            pe = int(pulse_idx[-1])
            # rest after
            rest_mask = np.arange(len(i_abs)) > pe
            rest_mask &= i_abs <= rest_current_max
            if rest_mask.sum() > 20:
                idx = np.flatnonzero(rest_mask)
                t_r = st[idx] - float(st[idx[0]])
                keep = t_r <= 1800.0
                rec = fit_recovery_tau(t_r[keep], v_all[idx][keep], tau_ct_ref=fit.tau_ct)
                for k, val in rec.items():
                    setattr(fit, k, val)
    return fit


def soc_ratio_features(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """From per-SOC DCIR rows build R_ratio_20_50, R_SOC_slope (§5.9 / §9.1#7)."""
    by = {int(r.get("soc")): r for r in rows if r.get("soc") is not None}
    out: dict[str, Any] = {}
    r20 = by.get(20, {}).get("R_30s_total")
    r50 = by.get(50, {}).get("R_30s_total")
    r80 = by.get(80, {}).get("R_30s_total")
    if r20 and r50 and r50 > 0:
        out["R_ratio_20_50"] = float(r20) / float(r50)
    if r80 and r50 and r50 > 0:
        out["R_ratio_80_50"] = float(r80) / float(r50)
    xs, ys = [], []
    for soc, key in ((20, r20), (50, r50), (80, r80)):
        if key is not None and np.isfinite(key):
            xs.append(float(soc))
            ys.append(float(key))
    if len(xs) >= 2:
        slope, intercept = _linreg(np.asarray(xs), np.asarray(ys))
        # R_SOC_slope in mOhm per %SOC. R20 - R80 is not reported separately:
        # for the three evenly spaced SOC points it equals -60 * slope.
        out["R_SOC_slope"] = slope
        if len(xs) == 3:
            coef = np.polyfit(xs, ys, 2)
            out["R_SOC_curvature"] = float(coef[0])
    return out


def result_to_dict(fit: DcirFitResult, *, suffix: str = "") -> dict[str, Any]:
    d = asdict(fit)
    if not suffix:
        return d
    return {f"{k}{suffix}": v for k, v in d.items()}
