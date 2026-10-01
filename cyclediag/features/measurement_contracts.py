"""Acquisition and empirical-model contracts; no inferred calibration or endpoints."""

from dataclasses import asdict, dataclass
import math
from typing import Mapping

import numpy as np
import pandas as pd

VE_SEMANTICS = "observed_window_v2"


def paired_counter(q, v):
    """Validate acquisition order before removing only identical adjacent pairs."""
    if q is None or v is None:
        return None, None, "missing_counter_or_voltage"
    q, v = np.asarray(q, dtype=float), np.asarray(v, dtype=float)
    if q.ndim != 1 or q.shape != v.shape or len(q) < 2:
        return None, None, "insufficient_pairs"
    if not (np.isfinite(q).all() and np.isfinite(v).all()):
        return None, None, "nonfinite_counter_or_voltage"
    dq = np.diff(q)
    if np.any(dq < 0):
        return None, None, "counter_reset"
    if np.any((dq == 0) & (np.diff(v) != 0)):
        return None, None, "conflicting_voltage_at_same_counter"
    if not np.any(dq > 0):
        return None, None, "zero_observed_capacity"
    keep = np.r_[True, dq > 0]
    return q[keep], v[keep], None


def continuity_reason(seg):
    """Conservative contract for a single observed leg, not a full-cycle claim."""
    if seg is None or len(seg) < 2:
        return "insufficient_samples"
    if seg.attrs.get("acquisition_disjoint", False):
        return "disjoint_acquisition_segments"
    # Numeric original indexes preserve missing/intervening source rows.
    if pd.api.types.is_numeric_dtype(seg.index.dtype):
        idx = np.asarray(seg.index, dtype=float)
        if not np.isfinite(idx).all() or np.any(np.diff(idx) != 1):
            return "noncontiguous_source_index"
    if "data_point" in seg:
        points = pd.to_numeric(seg["data_point"], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(points).all() or np.any(np.diff(points) != 1):
            return "noncontiguous_data_point"
    clocks = [c for c in ("step_time", "time") if c in seg]
    if not clocks:
        return "missing_acquisition_time"
    for col in clocks:
        t = pd.to_numeric(seg[col], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(t).all():
            return "missing_acquisition_time"
        if np.any(np.diff(t) < 0):
            return "acquisition_time_reset"
    return None


def observed_leg(seg, q):
    v = pd.to_numeric(seg["voltage"], errors="coerce").to_numpy(dtype=float) if "voltage" in seg else None
    qq, vv, pair_reason = paired_counter(q, v)
    reason = pair_reason or continuity_reason(seg)
    meta = {"status": "unknown" if reason else "ok", "reason": reason,
            "window": "observed_endpoints_only", "full_cycle_verified": False,
            "counter_start_Ah": None, "counter_end_Ah": None}
    if q is not None and len(q):
        meta["counter_start_Ah"] = float(q[0]) if np.isfinite(q[0]) else None
        meta["counter_end_Ah"] = float(q[-1]) if np.isfinite(q[-1]) else None
    if reason:
        return None, None, None, meta
    dq = float(qq[-1] - qq[0])
    energy = float(np.sum((vv[1:] + vv[:-1]) * 0.5 * np.diff(qq)))
    if not math.isfinite(energy):
        meta.update(status="unknown", reason="nonfinite_energy")
        return None, None, None, meta
    return dq, energy, energy / dq, meta


def observed_ve_compatible(row, baseline=None):
    """Metadata-free legacy VE is deliberately not comparable to corrected VE."""
    def valid(r):
        sem = r.get("feature_semantics")
        obs = r.get("observed_window")
        return (isinstance(sem, Mapping) and sem.get("VE") == VE_SEMANTICS
                and isinstance(obs, Mapping)
                and all(isinstance(obs.get(k), Mapping) and obs[k].get("status") == "ok"
                        for k in ("charge", "discharge")))
    return valid(row) and (baseline is None or valid(baseline))


def feature_semantics_match(row, baseline, feature):
    """Compatibility for corrected baseline/delta fields, never infer legacy units."""
    if feature in {"VE", "VE_observed"}:
        return observed_ve_compatible(row, baseline)
    a, b = row.get("feature_semantics"), baseline.get("feature_semantics")
    return (isinstance(a, Mapping) and isinstance(b, Mapping)
            and a.get(feature) is not None and a.get(feature) == b.get(feature))


@dataclass(frozen=True)
class TemperatureCalibration:
    ea_j_mol: float
    reference_c: float
    domain_c: tuple[float, float]
    calibration_id: str
    measurement_definition: str
    soc: str
    protocol: str
    time_scale: str
    temperature_source: str


@dataclass(frozen=True)
class MeasurementTemperature:
    temperature_c: float
    measurement_definition: str
    soc: str
    protocol: str
    time_scale: str
    temperature_source: str


def temp_correct_ir(r_mohm, temp_c, reference_c=25.0, ea=None):
    """Pure R(T)=A exp(Ea/RgT) transform; Ea must be explicitly supplied.

    No calibration claim: extraction additionally requires a matched contract.
    """
    try:
        r, t, ref, activation = map(float, (r_mohm, temp_c, reference_c, ea))
        if not all(math.isfinite(x) for x in (r, t, ref, activation)) or min(t, ref) <= -273.15:
            return None
        value = r * math.exp(activation / 8.314 * (1 / (ref + 273.15) - 1 / (t + 273.15)))
        return value if math.isfinite(value) and (r == 0 or value != 0) else None
    except (TypeError, ValueError, OverflowError, ZeroDivisionError):
        return None


def calibrated_resistance(raw, calibration, measurement):
    meta = {"status": "not_requested", "reason": None,
            "model": "R(T)=A exp(Ea/RgT)", "calibration_verified": False}
    if calibration is None:
        return None, meta
    meta.update(status="unknown", reason="invalid_calibration")
    if not isinstance(calibration, TemperatureCalibration):
        return None, meta
    # Metadata is serialized through provenance.json_safe at export.
    from cyclediag.provenance import json_safe
    meta["calibration"] = json_safe(asdict(calibration))
    fields = ("measurement_definition", "soc", "protocol", "time_scale", "temperature_source")
    try:
        lo, hi = calibration.domain_c
        if (not all(math.isfinite(float(x)) for x in (lo, hi, calibration.reference_c, calibration.ea_j_mol))
                or lo <= -273.15 or lo > hi
                or not lo <= calibration.reference_c <= hi
                or not isinstance(calibration.calibration_id, str) or not calibration.calibration_id.strip()
                or any(not isinstance(getattr(calibration, f), str) or not getattr(calibration, f).strip() for f in fields)):
            return None, meta
        if not isinstance(measurement, MeasurementTemperature):
            meta["reason"] = "missing_measurement_temperature"
            return None, meta
        meta["measurement"] = json_safe(asdict(measurement))
        if any(getattr(calibration, f) != getattr(measurement, f) for f in fields):
            meta["reason"] = "measurement_definition_mismatch"
            return None, meta
        if not math.isfinite(measurement.temperature_c) or not lo <= measurement.temperature_c <= hi:
            meta["reason"] = "measurement_temperature_outside_domain"
            return None, meta
        value = temp_correct_ir(raw, measurement.temperature_c, calibration.reference_c, calibration.ea_j_mol)
    except (TypeError, ValueError, OverflowError):
        return None, meta
    meta.update(status="ok" if value is not None else "unknown",
                reason=None if value is not None else "invalid_or_overflow_transform")
    return value, meta