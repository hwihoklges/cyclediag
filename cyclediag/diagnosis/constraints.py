"""Physics/data constraints that cap causal-diagnosis confidence.

These do not change indicator scores (Track A). They only mark when a
physicochemical mode interpretation (Track B) is under-determined.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from .schema import ELIGIBILITY_CONTRACT_VERSION


def _affirmed(value: Any) -> bool:
    return (type(value) is bool and value is True) or (type(value).__name__ == "bool_" and bool(value))


def _finite_number(value: Any) -> float | None:
    if isinstance(value, (bool, str)) or value is None:
        return None
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def mode_eligibility(
    row: Mapping[str, Any], config: Mapping[str, Any], mode: str,
    baseline_row: Mapping[str, Any] | None, available_features: set[str],
) -> dict[str, Any]:
    """Conservative scientific interpretation contract; scores never fill missing gates.

    A configuration without declared requirements remains exploratory, not validated.
    Group evidence consists of *usable* terms, not merely present columns.
    """
    spec = (config.get("modes") or {}).get(mode, {}).get("eligibility")
    reasons: list[str] = []
    groups: dict[str, bool] = {}
    domain = config.get("chemistry")
    if not isinstance(domain, str) or not domain.strip() or not isinstance(row.get("chemistry"), str) or not row["chemistry"].strip():
        reasons.append("chemistry_unverified")
    elif row["chemistry"] != domain:
        reasons.append("chemistry_mismatch")

    if not isinstance(spec, Mapping) or not spec.get("required_groups"):
        reasons.append("eligibility_contract_unconfigured")
        return {"contract_version": ELIGIBILITY_CONTRACT_VERSION, "eligible": False,
                "reasons": reasons, "required_groups": groups}

    semantics = row.get("feature_semantics")
    if not isinstance(semantics, Mapping) or semantics.get("schema") != "lges_measurements_v2":
        reasons.append("feature_semantics_unverified")

    if (not _affirmed(row.get("protocol_comparable")) or row.get("protocol_kind") != "routine"
            or row.get("protocol_excluded") is not False):
        reasons.append("protocol_unverified_or_incompatible")
    quality = _finite_number(row.get("quality_score"))
    coverage = _finite_number(row.get("quality_evidence_coverage"))
    if (row.get("quality_status") not in ("assessed", "partial") or quality is None
            or quality <= 0 or coverage is None or coverage <= 0):
        reasons.append("quality_unverified")
    failed = row.get("quality_gate_failed_groups")
    if not isinstance(failed, str) or failed.strip().lower() in ("nan", "none"):
        reasons.append("quality_groups_unverified")
        failed_groups: set[str] = set()
    else:
        failed_groups = {s.strip() for s in failed.split(",") if s.strip()}

    for group in spec["required_groups"]:
        name = group["name"]
        features = group["any_of"]
        present = any(f in available_features for f in features)
        groups[name] = present
        if not present:
            reasons.append(f"required_group_missing:{name}")
        for metric in group.get("quality_metrics", []):
            gate = {"pulse_current_stability": "pulse_stability"}.get(metric, metric)
            value = _finite_number(row.get(metric))
            if gate in failed_groups or value is None or (metric != "pulse_current_stability" and value <= 0):
                reasons.append(f"required_quality_unverified:{name}:{metric}")
        if group.get("fit") and not _affirmed(row.get(group["fit"])):
            reasons.append(f"fit_unverified:{name}")

    terms = (config.get("modes") or {})[mode].get("evidence") or []
    used = [t for t in terms if t["feature"] in available_features]
    for term in used:
        feature = term["feature"]
        if feature.startswith(("R_ct", "R_ohmic", "A_diff", "tau_ct")):
            suffix = feature[feature.rfind("_soc"):] if "_soc" in feature else ""
            if not _affirmed(row.get(f"dcir_fit_valid{suffix}")):
                reasons.append(f"fit_unverified:{feature}")
    if any(t.get("direction") == "decrease_vs_baseline" or t["feature"].startswith("delta_") for t in used):
        if (not isinstance(baseline_row, Mapping) or baseline_row.get("protocol_kind") != "routine"
            or baseline_row.get("protocol_excluded") is not False
            or not _affirmed(baseline_row.get("protocol_comparable"))):
            reasons.append("baseline_protocol_unverified")
        if (isinstance(baseline_row, Mapping) and
            (baseline_row.get("chemistry") != row.get("chemistry") or
             not isinstance(row.get("protocol_id"), str) or not row["protocol_id"].strip() or
             row["protocol_id"] != baseline_row.get("protocol_id"))):
            reasons.append("baseline_protocol_mismatch_or_unknown")
        meta = row.get("baseline_semantics")
        if (not isinstance(meta, Mapping) or not isinstance(baseline_row, Mapping)
                or _finite_number(meta.get("cycle")) != _finite_number(baseline_row.get("cycle"))
                or _finite_number(meta.get("requested_cycle")) != _finite_number(meta.get("cycle"))
                or not isinstance(meta.get("feature_semantics"), Mapping)):
            reasons.append("baseline_unverified")
        else:
            sem = row.get("feature_semantics")
            if not isinstance(sem, Mapping) or sem != meta["feature_semantics"]:
                reasons.append("feature_semantics_mismatch")
    if any(t["feature"] in ("VE", "VE_observed") for t in used):
        from cyclediag.features.measurement_contracts import observed_ve_compatible
        if not observed_ve_compatible(row, baseline_row):
            reasons.append("observed_window_unverified")

    return {"contract_version": ELIGIBILITY_CONTRACT_VERSION,
            "eligible": not reasons, "reasons": list(dict.fromkeys(reasons)),
            "required_groups": groups}


def scientific_validity(row: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, Any]:
    """Machine-readable limits; proxy agreement is not causal identification."""
    warnings = ["LLI_LAM_kinetics_proxies_confounded", "correlated_evidence_not_independent",
                "matched_protocol_temperature_SOC_required", "absolute_losses_not_identified"]
    domain = config.get("chemistry")
    observed = row.get("chemistry")
    domain_status = "unverified" if not domain or not observed else "matched" if domain == observed else "mismatch"
    if domain_status != "matched":
        warnings.append(f"chemistry_domain_{domain_status}")
    try:
        quality_known = math.isfinite(float(row.get("quality_score")))
    except (TypeError, ValueError):
        quality_known = False
    if not quality_known:
        warnings.append("measurement_quality_unverified")
    unit_warnings = row.get("unit_warnings")
    if isinstance(unit_warnings, list):
        warnings.extend(unit_warnings)
    for key, value in row.items():
        if str(key).startswith(("quality_scientific_validity", "ecm_scientific_validity")) and isinstance(value, dict):
            warnings.extend(value.get("warnings", []))
    try:
        ve = float(row.get("VE"))
        if math.isfinite(ve) and not 0 < ve <= 1:
            warnings.append("VE_outside_nominal_fraction_range_check_protocol_and_units")
    except (TypeError, ValueError):
        pass
    for name in ("CE", "CE_local_20"):
        try:
            ce = float(row.get(name))
            if math.isfinite(ce) and (ce > 100 or ce < 0):
                warnings.append(f"{name}_outside_nominal_percent_range_not_causal_evidence")
        except (TypeError, ValueError):
            pass
    if any(str(k).startswith(("R_ct", "R_ohmic", "A_diff", "tau_ct")) for k in row):
        warnings.extend(["equivalent_circuit_parameters_not_unique_mechanisms",
                         "ECM_requires_stable_pulse_local_linearity_and_time_resolution"])
    warnings.extend(constraint_flags(row, config))
    return {
        "score_kind": "heuristic_pattern_support",
        "confidence_kind": "uncalibrated_heuristic_not_probability",
        "causally_identified": False, "validated_probability": False,
        "chemistry_domain": domain, "domain_status": domain_status,
        "feature_conventions": {"VE": "fraction", "CE": "percent"},
        "warnings": list(dict.fromkeys(warnings)),
    }


def constraint_flags(
    row: Mapping[str, Any],
    config: Mapping[str, Any] | None = None,
) -> list[str]:
    """Return active constraint tags for one cycle row."""
    flags: list[str] = []
    cfg = config or {}
    cons = cfg.get("constraints") or {}

    temp = row.get("temperature_available")
    if temp is False or temp == 0 or temp == "False":
        flags.append("no_temperature_log")

    if cons.get("stack_pressure_MPa") is None and row.get("stack_pressure_MPa") is None:
        flags.append("stack_pressure_unknown")

    if cons.get("halfcell_calibrated") is not True:
        flags.append("halfcell_uncalibrated")

    # protocol contamination — diagnosis should not trust these rows
    if row.get("protocol_excluded") is not False and row.get("protocol_excluded") is not None:
        flags.append("protocol_excluded")
    kind = str(row.get("protocol_kind") or "")
    if kind and kind not in ("routine", "unknown", "nan", ""):
        flags.append(f"protocol_{kind}")

    failed_groups = row.get("quality_gate_failed_groups")
    if isinstance(failed_groups, str) and failed_groups.strip() not in ("", "nan"):
        flags.append("quality_groups_failed")
    try:
        if float(row.get("quality_score")) <= 0:
            flags.append("quality_gate_failed")
    except (TypeError, ValueError):
        pass

    return flags


def confidence_multiplier(
    mode: str,
    flags: list[str],
) -> float:
    """Down-weight mode confidence when constraints apply."""
    m = 1.0
    if "quality_gate_failed" in flags:
        return 0.0
    if "protocol_excluded" in flags or any(f.startswith("protocol_") for f in flags):
        # engine should skip these rows; if scored anyway, collapse confidence
        return 0.0
    if "no_temperature_log" in flags and mode in (
        "impedance", "interface_R", "solid_diffusion", "transport", "contact_loss", "contact",
    ):
        m *= 0.85
    if "stack_pressure_unknown" in flags and mode in ("contact_loss", "contact"):
        m *= 0.7
    if "halfcell_uncalibrated" in flags:
        # Level-3 absolute estimates unavailable — slight caution on all modes
        m *= 0.95
    return float(m)
