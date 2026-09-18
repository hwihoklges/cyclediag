"""Cycle-level data-quality metrics — IMPROVEMENT_ROADMAP §5.13."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def _v_noise_sigma(v: np.ndarray, window: int = 50) -> float | None:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < window * 2:
        return None
    resid = []
    for i in range(0, len(v) - window, window // 2):
        sl = v[i : i + window]
        x = np.arange(len(sl), dtype=float)
        coef = np.polyfit(x, sl, 1)
        resid.append(sl - np.polyval(coef, x))
    if not resid:
        return None
    return float(np.nanstd(np.concatenate(resid)))


def _quant_step(v: np.ndarray) -> float | None:
    v = np.asarray(v, dtype=float)
    d = np.diff(v[np.isfinite(v)])
    d = np.abs(d[np.abs(d) > 0])
    if len(d) == 0:
        return None
    return float(np.min(d))


def cycle_quality_metrics(
    cycle_df: pd.DataFrame,
    *,
    expected_v_window: tuple[float, float] = (2.5, 4.2),
    rest_current_max: float = 0.5,
    tau_relax_est: float | None = None,
) -> dict[str, Any]:
    """Compute per-cycle quality fields from raw points."""
    out: dict[str, Any] = {
        "samples_per_mV": None,
        "v_noise_sigma": None,
        "quant_step_est": None,
        "dqdv_snr": None,
        "rest_sufficiency": None,
        "pulse_sample_count_1s": None,
        "pulse_current_stability": None,
        "leg_completeness": None,
        "temperature_available": False,
        "quality_score": None,
        "quality_gate_failed_groups": "",
        "quality_evidence_coverage": 0.0,
        "quality_status": "unknown",
        "scientific_validity": {
            "score_kind": "heuristic_data_quality_not_probability",
            "warnings": ["voltage_window_and_quality_targets_require_protocol_validation",
                         "voltage_span_not_per_leg_completeness",
                         "dqdv_snr_is_voltage_span_noise_proxy"],
        },
    }
    if cycle_df is None or cycle_df.empty or "voltage" not in cycle_df.columns:
        return out

    v = pd.to_numeric(cycle_df["voltage"], errors="coerce").to_numpy(dtype=float)
    out["v_noise_sigma"] = _v_noise_sigma(v)
    out["quant_step_est"] = _quant_step(v)

    finite_v = v[np.isfinite(v)]
    if len(finite_v) >= 2:
        v_span_mv = (float(np.nanmax(finite_v)) - float(np.nanmin(finite_v))) * 1000.0
        if v_span_mv > 1e-6:
            out["samples_per_mV"] = len(finite_v) / v_span_mv

    if "temperature" in cycle_df.columns:
        t = pd.to_numeric(cycle_df["temperature"], errors="coerce")
        out["temperature_available"] = bool(np.isfinite(t).any())

    lo, hi = expected_v_window
    if len(finite_v):
        covered = float(np.nanmax(finite_v) - np.nanmin(finite_v))
        expect = max(hi - lo, 1e-9)
        out["leg_completeness"] = min(1.0, covered / expect)

    if "current" in cycle_df.columns:
        i = pd.to_numeric(cycle_df["current"], errors="coerce").abs().to_numpy()
        if "step_time" in cycle_df.columns:
            st = pd.to_numeric(cycle_df["step_time"], errors="coerce").to_numpy()
        else:
            st = None
        pulse = np.isfinite(i) & (i > max(rest_current_max * 10, 20.0))
        if pulse.any() and st is not None:
            p0 = int(np.argmax(pulse))
            t0 = float(st[p0]) if np.isfinite(st[p0]) else 0.0
            early = pulse & np.isfinite(st) & ((st - t0) <= 1.0) & ((st - t0) >= 0)
            # only within first pulse step
            out["pulse_sample_count_1s"] = int(early.sum())
            ii = i[pulse]
            if len(ii) > 5:
                out["pulse_current_stability"] = float(np.nanstd(ii) / max(np.nanmedian(ii), 1e-9))

        rest = np.isfinite(i) & (i <= rest_current_max)
        if rest.any() and st is not None:
            # longest rest duration
            best = 0.0
            a = None
            for idx, flag in enumerate(rest):
                if flag and a is None:
                    a = idx
                elif not flag and a is not None:
                    dur = float(st[idx - 1] - st[a]) if np.isfinite(st[idx - 1]) and np.isfinite(st[a]) else 0.0
                    best = max(best, dur)
                    a = None
            if a is not None:
                dur = float(st[len(rest) - 1] - st[a]) if np.isfinite(st[-1]) and np.isfinite(st[a]) else 0.0
                best = max(best, dur)
            if tau_relax_est is not None and np.isfinite(tau_relax_est) and tau_relax_est > 0:
                out["rest_sufficiency"] = best / tau_relax_est
            else:
                out["scientific_validity"]["warnings"].append("rest_sufficiency_requires_measured_tau")

    # simple SNR proxy: v span / noise
    if out["v_noise_sigma"] and out["v_noise_sigma"] > 0 and len(finite_v):
        out["dqdv_snr"] = (float(np.nanmax(finite_v)) - float(np.nanmin(finite_v))) / out["v_noise_sigma"]

    # weighted geometric mean of clipped ratios
    targets = {
        "samples_per_mV": (out["samples_per_mV"], 0.5, 0.25),
        "dqdv_snr": (out["dqdv_snr"], 10.0, 0.25),
        "leg_completeness": (out["leg_completeness"], 0.9, 0.2),
        "rest_sufficiency": (out["rest_sufficiency"], 3.0, 0.15),
        "pulse_stability": (
            None if out["pulse_current_stability"] is None else max(0.0, 1.0 - out["pulse_current_stability"] / 0.02),
            1.0,
            0.15,
        ),
    }
    failed = []
    score = 1.0
    wsum = 0.0
    for name, (val, tgt, w) in targets.items():
        if val is None or not np.isfinite(val):
            continue
        q = float(np.clip(val / tgt, 0.0, 1.0))
        if q < 1.0:
            failed.append(name)
        score *= q ** w
        wsum += w
    out["quality_score"] = float(score ** (1.0 / wsum)) if wsum > 0 else None
    out["quality_evidence_coverage"] = wsum / sum(w for _, _, w in targets.values())
    out["quality_status"] = "failed" if failed else "partial" if 0 < wsum < 1 else "assessed" if wsum else "unknown"
    out["quality_gate_failed_groups"] = ",".join(failed)
    return out
