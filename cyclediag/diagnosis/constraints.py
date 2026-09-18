"""Physics/data constraints that cap causal-diagnosis confidence.

These do not change indicator scores (Track A). They only mark when a
physicochemical mode interpretation (Track B) is under-determined.
"""

from __future__ import annotations

from typing import Any, Mapping
import math


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
    if bool(row.get("protocol_excluded")):
        flags.append("protocol_excluded")
    kind = str(row.get("protocol_kind") or "")
    if kind and kind not in ("routine", "unknown", "nan", ""):
        flags.append(f"protocol_{kind}")

    failed_groups = row.get("quality_gate_failed_groups")
    if isinstance(failed_groups, str) and failed_groups.strip() not in ("", "nan"):
        flags.append("quality_gate_failed")
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
