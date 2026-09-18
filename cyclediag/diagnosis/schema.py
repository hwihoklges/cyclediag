"""Diagnosis result schema and serialization helpers."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

DIAGNOSIS_VERSION_FULLCELL = "fullcell_v1"
DIAGNOSIS_VERSION_HC_CALIBRATED = "hc_calibrated_v1"
DIAGNOSIS_MODEL_VERSION = "pattern_rule_v1"


@dataclass
class DiagnosisResult:
    """One degradation-mode diagnosis for a cycle (or cell summary)."""

    degradation_mode: str
    level: int = 1
    estimate: float | None = None
    unit: str = "pattern_score_0_1"
    confidence: float = 0.0
    evidence_count: int = 0
    supporting_features: list[str] = field(default_factory=list)
    conflicting_features: list[str] = field(default_factory=list)
    data_quality_score: float = 0.0
    diagnosis_valid: bool = False
    diagnosis_version: str = DIAGNOSIS_VERSION_FULLCELL
    diagnosis_method: str = "rule_pattern"
    diagnosis_model_version: str = DIAGNOSIS_MODEL_VERSION
    status: str = "unknown"
    evidence_coverage: float = 0.0
    scientific_validity: dict[str, Any] = field(default_factory=lambda: {
        "score_kind": "heuristic_pattern_support",
        "confidence_kind": "uncalibrated_heuristic_not_probability",
        "causally_identified": False,
        "validated_probability": False,
    })

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# Modes emitted as wide columns (Level 1) — liquid-cell default
PATTERN_MODES: tuple[str, ...] = (
    "LLI",
    "LAM_PE",
    "LAM_NE",
    "impedance",
    "transport",
    "plating",
    "contact",
)

# ASSB Si-rich default (IMPROVEMENT_ROADMAP §3)
ASSB_PATTERN_MODES: tuple[str, ...] = (
    "contact_loss",
    "interface_R",
    "SE_decomposition",
    "microshort",
    "LAM_PE",
    "LLI",
    "solid_diffusion",
)

# Map internal mode → column prefix
MODE_COLUMN_PREFIX: dict[str, str] = {
    "LLI": "LLI",
    "LAM_PE": "LAM_PE",
    "LAM_NE": "LAM_NE",
    "impedance": "impedance",
    "transport": "transport_limitation",
    "plating": "plating_risk",
    "contact": "contact_loss",
    "contact_loss": "contact_loss",
    "interface_R": "interface_R",
    "SE_decomposition": "SE_decomposition",
    "microshort": "microshort",
    "solid_diffusion": "solid_diffusion",
}


def score_column_name(mode: str) -> str:
    prefix = MODE_COLUMN_PREFIX.get(mode, mode)
    if mode in ("LLI", "LAM_PE", "LAM_NE", "impedance"):
        return f"{prefix}_pattern_score"
    if mode in ("contact_loss", "interface_R", "SE_decomposition", "microshort", "solid_diffusion"):
        return f"{prefix}_score"
    return f"{prefix}_score"


def confidence_column_name(mode: str) -> str:
    prefix = MODE_COLUMN_PREFIX.get(mode, mode)
    return f"{prefix}_confidence"
