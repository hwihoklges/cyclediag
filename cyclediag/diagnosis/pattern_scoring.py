"""Method A — config-driven degradation-mode pattern scoring."""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .schema import (
    DIAGNOSIS_MODEL_VERSION,
    DIAGNOSIS_VERSION_FULLCELL,
    DiagnosisResult,
    PATTERN_MODES,
)

_CONFIG_DIR = Path(__file__).resolve().parent / "config"
_DEFAULT_CONFIG = _CONFIG_DIR / "mode_weights_fullcell_v1.json"


def load_mode_weights(path: str | Path | None = None) -> dict[str, Any]:
    """Load mode weight config (JSON). Optional YAML if PyYAML is installed.

    Default: ASSB Si-rich weights when present, else fullcell_v1.
    """
    if path is None:
        assb = _CONFIG_DIR / "mode_weights_assb_si_v1.json"
        cfg_path = assb if assb.exists() else _DEFAULT_CONFIG
    else:
        cfg_path = Path(path)
    text = cfg_path.read_text(encoding="utf-8")
    if cfg_path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "PyYAML required for .yaml config; use mode_weights_fullcell_v1.json"
            ) from exc
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    if not isinstance(data, dict) or "modes" not in data:
        raise ValueError(f"Invalid diagnosis weight config: {cfg_path}")
    return data


@lru_cache(maxsize=4)
def _cached_default_config() -> dict[str, Any]:
    assb = _CONFIG_DIR / "mode_weights_assb_si_v1.json"
    return load_mode_weights(assb if assb.exists() else _DEFAULT_CONFIG)


def _finite(x: Any) -> float | None:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return v


def _signed_evidence(
    value: float,
    *,
    direction: str,
    scale: float,
    baseline: float | None = None,
) -> float:
    """Map raw feature → signed evidence in [-1, 1] (positive = supports mode)."""
    s = abs(float(scale)) if scale and abs(scale) > 1e-15 else 1.0
    d = (direction or "increase").lower()

    if d == "increase":
        raw = value / s
    elif d == "decrease":
        raw = -value / s
    elif d == "either":
        raw = abs(value) / s
    elif d == "decrease_from_100":
        raw = (100.0 - value) / s
    elif d == "decrease_vs_baseline":
        if baseline is None or not math.isfinite(baseline) or abs(baseline) < 1e-15:
            return 0.0
        raw = (baseline - value) / (abs(baseline) * 0.1 + s * 0.1)
    else:
        raw = value / s

    # tanh soft clip
    return float(math.tanh(raw))


def _evidence_for_term(
    row: Mapping[str, Any],
    term: Mapping[str, Any],
    baseline_row: Mapping[str, Any] | None,
) -> tuple[str, float | None]:
    feat = str(term["feature"])
    val = _finite(row.get(feat))
    if val is None:
        return feat, None
    # Preserve raw CE in the input; out-of-range cycle accounting is not
    # evidence for or against a particular degradation mechanism.
    if feat in {"CE", "CE_local_20"} and not 0 <= val <= 100:
        return feat, None
    if feat == "VE":
        units = row.get("feature_units")
        if isinstance(units, Mapping) and units.get("VE", "fraction") != "fraction":
            return feat, None
        if not 0 < val <= 1:
            return feat, None
    if feat == "CI_per_hour" and (_finite(row.get("CE")) or 0) > 100:
        return feat, None
    if feat.startswith(("R_ct", "R_ohmic", "A_diff", "tau_ct")):
        suffix = feat[feat.rfind("_soc"):] if "_soc" in feat else ""
        fit_valid = row.get(f"dcir_fit_valid{suffix}")
        if fit_valid is not None and not bool(fit_valid):
            return feat, None
    baseline = None
    if term.get("direction") == "decrease_vs_baseline":
        ref = term.get("baseline_ref") or feat
        baseline = _finite((baseline_row or {}).get(ref))
        if baseline is None or baseline == 0:
            return feat, None
        if feat == "VE" and not 0 < baseline <= 1:
            return feat, None
    signed = _signed_evidence(
        val,
        direction=str(term.get("direction", "increase")),
        scale=float(term.get("scale", 1.0)),
        baseline=baseline,
    )
    return feat, signed


def score_mode_for_row(
    row: Mapping[str, Any],
    mode: str,
    config: dict[str, Any],
    *,
    baseline_row: Mapping[str, Any] | None = None,
    other_mode_scores: Mapping[str, float] | None = None,
) -> DiagnosisResult:
    """Compute Level-1 pattern DiagnosisResult for one mode on one cycle row."""
    modes = config.get("modes") or {}
    mode_cfg = modes.get(mode) or {}
    evidence_terms = list(mode_cfg.get("evidence") or [])
    support_thr = float(config.get("support_threshold", 0.25))
    conflict_thr = float(config.get("conflict_threshold", 0.25))
    min_ev = int(config.get("min_evidence_for_valid", 2))
    collision_thr = float(config.get("mode_collision_score", 0.65))
    from .constraints import scientific_validity
    validity = scientific_validity(row, config)
    validity["warnings"].extend(
        f"missing_baseline:{term.get('baseline_ref') or term['feature']}"
        for term in evidence_terms
        if term.get("direction") == "decrease_vs_baseline"
        and _finite((baseline_row or {}).get(term.get("baseline_ref") or term["feature"])) in (None, 0)
    )

    weights: list[float] = []
    signed_vals: list[float] = []
    names: list[str] = []
    for term in evidence_terms:
        w = float(term.get("weight", 1.0))
        if not math.isfinite(w) or w <= 0:
            continue
        name, signed = _evidence_for_term(row, term, baseline_row)
        if signed is None:
            continue
        names.append(name)
        signed_vals.append(signed)
        weights.append(max(0.0, w))

    n = len(signed_vals)
    data_quality = float(n / max(1, len(evidence_terms))) if evidence_terms else 0.0

    if n == 0 or sum(weights) <= 0:
        return DiagnosisResult(
            degradation_mode=mode,
            level=1,
            estimate=None,
            unit="pattern_score_0_1",
            confidence=0.0,
            evidence_count=0,
            supporting_features=[],
            conflicting_features=[],
            data_quality_score=data_quality,
            diagnosis_valid=False,
            scientific_validity=validity,
            diagnosis_version=str(config.get("diagnosis_version", DIAGNOSIS_VERSION_FULLCELL)),
            diagnosis_method=str(config.get("diagnosis_method", "rule_pattern")),
            diagnosis_model_version=str(
                config.get("diagnosis_model_version", DIAGNOSIS_MODEL_VERSION)
            ),
        )

    w_arr = np.asarray(weights, dtype=float)
    s_arr = np.asarray(signed_vals, dtype=float)
    # Pattern score: weighted mean of positive-part evidence, clipped to [0,1]
    pos = np.clip(s_arr, 0.0, None)
    score = float(np.average(pos, weights=w_arr))
    score = float(np.clip(score, 0.0, 1.0))

    supporting = [n for n, s in zip(names, signed_vals) if s >= support_thr]
    conflicting = [n for n, s in zip(names, signed_vals) if s <= -conflict_thr]

    agree = float(np.average((s_arr > 0).astype(float), weights=w_arr))
    # Prefer extract quality_score when present over evidence fill-rate alone.
    q_row = row.get("quality_score")
    try:
        q_row_f = float(q_row) if q_row is not None else float("nan")
    except (TypeError, ValueError):
        q_row_f = float("nan")
    dq = (
        0.6 * data_quality + 0.4 * max(0.0, min(1.0, q_row_f))
        if math.isfinite(q_row_f)
        else data_quality
    )
    conf = 0.45 * dq + 0.40 * agree + 0.15 * min(1.0, n / 5.0)

    # Mode collision: another mode also high → lower confidence
    if other_mode_scores:
        rivals = [
            v for k, v in other_mode_scores.items()
            if k != mode and v is not None and v >= collision_thr and score >= collision_thr
        ]
        if rivals:
            conf *= 0.7

    from .constraints import confidence_multiplier, constraint_flags

    cflags = constraint_flags(row, config)
    conf *= confidence_multiplier(mode, cflags)

    conf = float(np.clip(conf, 0.0, 1.0))
    valid = n >= min_ev and dq >= 0.3 and conf >= 0.25 and conf > 0.0
    valid = valid and validity["domain_status"] != "mismatch"
    if validity["domain_status"] == "mismatch":
        conf = 0.0

    return DiagnosisResult(
        degradation_mode=mode,
        level=1,
        estimate=score,
        unit="pattern_score_0_1",
        confidence=conf,
        evidence_count=n,
        supporting_features=supporting,
        conflicting_features=conflicting,
        data_quality_score=dq,
        diagnosis_valid=valid,
        status="heuristic_support" if valid else "insufficient_or_invalid_evidence",
        evidence_coverage=data_quality,
        scientific_validity=validity,
        diagnosis_version=str(config.get("diagnosis_version", DIAGNOSIS_VERSION_FULLCELL)),
        diagnosis_method=str(config.get("diagnosis_method", "rule_pattern")),
        diagnosis_model_version=str(
            config.get("diagnosis_model_version", DIAGNOSIS_MODEL_VERSION)
        ),
    )


def score_all_modes_for_row(
    row: Mapping[str, Any],
    config: dict[str, Any] | None = None,
    *,
    baseline_row: Mapping[str, Any] | None = None,
    modes: tuple[str, ...] | None = None,
) -> dict[str, DiagnosisResult]:
    """Two-pass scoring so mode-collision can use peer scores."""
    cfg = config or _cached_default_config()
    if modes is None:
        cfg_modes = tuple(cfg.get("modes", {}).keys())
        mode_list = cfg_modes or PATTERN_MODES
    else:
        mode_list = modes
    prelim: dict[str, float] = {}
    # first pass without collision
    first: dict[str, DiagnosisResult] = {}
    for mode in mode_list:
        r = score_mode_for_row(row, mode, cfg, baseline_row=baseline_row)
        first[mode] = r
        prelim[mode] = float(r.estimate or 0.0)
    # second pass with collision awareness
    out: dict[str, DiagnosisResult] = {}
    for mode in mode_list:
        out[mode] = score_mode_for_row(
            row, mode, cfg,
            baseline_row=baseline_row,
            other_mode_scores=prelim,
        )
    return out
