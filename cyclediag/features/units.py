"""Capacity / current unit helpers — prefer header units over value heuristics."""

from __future__ import annotations

import re
import math
from typing import Any


_UNIT_RE = re.compile(
    r"\((?P<u>m?ah|a|ma|v|mv|ohm|mohm|s|sec|celsius|degc|°c)\)",
    re.IGNORECASE,
)


def parse_unit_from_header(header: str | None) -> str | None:
    if not header:
        return None
    h = str(header)
    bracket = re.search(r"[\[(]\s*([^\])]+)\s*[\])]", h)
    if bracket:
        unit = bracket.group(1).strip().lower()
        if unit in {"a", "ma", "ah", "mah", "v", "mv", "s", "sec", "min", "h", "ms"}:
            return unit
    if h.lower().endswith("_sec"):
        return "s"
    m = _UNIT_RE.search(h)
    if m:
        return m.group("u").lower().replace("°", "deg")
    hl = h.lower().replace(" ", "")
    if hl.endswith("mah") or "capacity(mah)" in hl:
        return "mah"
    if hl.endswith("(ah)") or hl.endswith("_ah"):
        return "ah"
    if "current(ma)" in hl or hl.endswith("ma"):
        return "ma"
    if "current(a)" in hl:
        return "a"
    return None


def canonical_unit_factor(unit: str, canonical: str) -> float:
    """Supported dimensional conversions only; no magnitude-based inference."""
    factors = {
        "a": {"a": 1.0, "ma": 0.001},
        "ah": {"ah": 1.0, "mah": 0.001, "a*h": 1.0, "a·h": 1.0},
        "v": {"v": 1.0, "mv": 0.001},
        "s": {"s": 1.0, "sec": 1.0, "ms": 0.001, "min": 60.0, "h": 3600.0},
    }
    try:
        return factors[canonical.lower()][unit.strip().lower()]
    except KeyError as exc:
        raise ValueError(f"Unsupported {canonical} unit: {unit}") from exc


def capacity_to_ah(
    q: Any,
    *,
    unit: str | None = None,
    header: str | None = None,
) -> float | None:
    """Convert capacity to Ah using explicit unit / header; never guess from magnitude."""
    try:
        v = float(q)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    u = (unit or parse_unit_from_header(header) or "").lower()
    if u:
        return v * canonical_unit_factor(u, "ah")
    # Logical columns after normalize are already Ah for Studio Ah exports.
    # Do NOT divide large values — that heuristic breaks ~72 Ah cells.
    return v


def current_to_a(
    i: Any,
    *,
    unit: str | None = None,
    header: str | None = None,
) -> float | None:
    try:
        v = float(i)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    u = (unit or parse_unit_from_header(header) or "").lower()
    return v * canonical_unit_factor(u or "a", "a")
