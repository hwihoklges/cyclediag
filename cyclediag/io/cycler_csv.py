"""Load PNE-style cycler CSV into a normalized DataFrame."""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Any, Mapping

import pandas as pd

# Logical name -> default column name (pne_studio presets aligned)
PNE_DEFAULT_COLUMNS: dict[str, str] = {
    "cycle": "CycleIndex",
    "voltage": "Voltage(V)",
    "capacity": "ChargeCapacity",
    "step_type": "StepType",
    "current": "Current(A)",
    "time": "TotalTime_sec",
    "data_point": "Data_Point",
}


@dataclass
class ColumnMap:
    """Maps logical fields to CSV columns; ``units`` maps logical names to units.

    Explicit units must agree with labelled headers (conflicts raise ValueError).
    Bare columns retain the historical canonical-unit assumption, recorded in
    attrs['unit_metadata'] and attrs['unit_warnings']; never infer magnitude.
    Normalized frames use A, Ah, V, s and may be normalized again unchanged.
    """

    cycle: str = "CycleIndex"
    voltage: str = "Voltage(V)"
    capacity: str = "ChargeCapacity"
    step_type: str = "StepType"
    current: str = "Current(A)"
    time: str = "TotalTime_sec"
    step_time: str = "StepTime_sec"
    data_point: str = "Data_Point"
    discharge_capacity: str = "DischargeCapacity"
    temperature: str = "Temperature"
    units: dict[str, str] = field(default_factory=dict)

    @classmethod
    def pne_default(cls) -> ColumnMap:
        return cls()

    @classmethod
    def studio_default(cls) -> ColumnMap:
        """PNE Studio / pne_studio2 UI column names."""
        return cls(
            cycle="TotalCycle",
            voltage="Voltage",
            capacity="Capacity",
            step_type="StepType",
            current="Current",
            time="TotalTime_sec",
            step_time="StepTime_sec",
            data_point="Data_Point",
            discharge_capacity="DischargeCapacity",
            temperature="Temperature",
        )


def _resolve_column(
    df: pd.DataFrame,
    name: str,
    aliases: tuple[str, ...] = (),
    *,
    exclude: set[str] | None = None,
) -> str | None:
    """Resolve a logical column, preferring alias order over CSV column order.

    Example: prefer ``StepTime_sec (sec)`` over empty ``StepTime`` when both exist.
    """
    exclude = exclude or set()
    if name in df.columns and name not in exclude:
        return name

    def _norm(label: str) -> str:
        return re.split(r"[\[(]", str(label))[0].replace(" ", "").replace("_", "").lower()

    # Alias priority first (more specific names listed before bare names).
    for alias in (name, *aliases):
        a = _norm(alias)
        for col in df.columns:
            if col in exclude:
                continue
            if _norm(col) == a:
                return col
    return None


_CYCLE_ALIASES = ("TotalCycle", "CycleIndex", "CycleNum", "Cycle")
_VOLTAGE_ALIASES = ("Voltage", "Voltage(V)")
_CHARGE_CAP_ALIASES = ("Capacity", "ChargeCapacity")
_DISCHARGE_CAP_ALIASES = ("DischargeCapacity",)
_STEP_ALIASES = ("StepType", "Step")
_CURRENT_ALIASES = ("Current", "Current(A)", "AvgCurrent", "AvgCurrent(A)", "Current(mA)")
_TIME_ALIASES = ("TotalTime_sec", "TotalTime")
_STEP_TIME_ALIASES = ("StepTime_sec", "StepTime")
_TEMP_ALIASES = ("Temperature", "Temp", "CellTemp", "Aux_Temperature", "AuxTemp")


def normalize_cycler_dataframe(
    df: pd.DataFrame,
    column_map: ColumnMap | Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    """Rename and convert labelled quantities to A, Ah, V, s; preserve provenance."""
    if column_map is None:
        cmap = ColumnMap.pne_default()
    elif isinstance(column_map, ColumnMap):
        cmap = column_map
    else:
        cmap = ColumnMap(**dict(column_map))

    rename: dict[str, str] = {}
    used: set[str] = set()

    def _map(logical: str, preferred: str, aliases: tuple[str, ...]) -> None:
        resolved = logical if logical in df.columns else _resolve_column(df, preferred, aliases, exclude=used)
        if resolved:
            rename[resolved] = logical
            used.add(resolved)

    _map("cycle", cmap.cycle, _CYCLE_ALIASES)
    _map("voltage", cmap.voltage, _VOLTAGE_ALIASES)
    _map("step_type", cmap.step_type, _STEP_ALIASES)
    _map("current", cmap.current, _CURRENT_ALIASES)
    _map("time", cmap.time, _TIME_ALIASES)
    _map("step_time", cmap.step_time, _STEP_TIME_ALIASES)
    _map("data_point", cmap.data_point, ("Data_Point", "DataPoint"))
    _map("temperature", cmap.temperature, _TEMP_ALIASES)
    _map("charge_cv_capacity", "ChargeCVCapacity", ("ChargeCVCapacity",))

    charge_cap = "charge_capacity" if "charge_capacity" in df.columns else _resolve_column(df, cmap.capacity, _CHARGE_CAP_ALIASES, exclude=used)
    if charge_cap:
        rename[charge_cap] = "charge_capacity"
        used.add(charge_cap)

    discharge_cap = "discharge_capacity" if "discharge_capacity" in df.columns else _resolve_column(
        df, cmap.discharge_capacity, _DISCHARGE_CAP_ALIASES, exclude=used,
    )
    if discharge_cap:
        rename[discharge_cap] = "discharge_capacity"
        used.add(discharge_cap)

    out = df.rename(columns=rename).copy()
    from cyclediag.features.units import canonical_unit_factor, parse_unit_from_header

    dimensions = {
        "current": "a", "voltage": "v", "time": "s", "step_time": "s",
        "capacity": "ah", "charge_capacity": "ah", "discharge_capacity": "ah",
        "charge_cv_capacity": "ah",
    }
    unknown = set(cmap.units) - set(dimensions)
    if unknown:
        raise ValueError(f"Unknown logical unit fields: {sorted(unknown)}")
    metadata = dict(df.attrs.get("unit_metadata", {}))
    warnings = list(df.attrs.get("unit_warnings", []))
    sources = {logical: source for source, logical in rename.items()}
    for logical, canonical in dimensions.items():
        if logical not in out.columns:
            continue
        source = sources.get(logical, logical)
        base = re.split(r"[\[(]", str(source))[0].replace(" ", "").replace("_", "").lower()
        candidates = [str(c) for c in df.columns if re.split(r"[\[(]", str(c))[0].replace(" ", "").replace("_", "").lower() == base]
        if len(candidates) > 1:
            warnings.append(f"ambiguous_columns:{logical}:selected={source}:candidates={','.join(candidates)}")
        explicit = cmap.units.get(logical)
        if explicit is None and canonical == "ah":
            explicit = cmap.units.get("capacity")
        if logical in metadata:
            # Overrides describe input values, not the original pre-conversion CSV.
            if explicit and canonical_unit_factor(explicit, canonical) != 1.0:
                raise ValueError(f"{logical} is already canonical; set raw units before normalization")
            continue
        header_unit = parse_unit_from_header(source)
        if header_unit is None and re.search(r"[\[(].+[\])]", source):
            raise ValueError(f"Unsupported unit-bearing header: {source}")
        if explicit and header_unit:
            if canonical_unit_factor(explicit, canonical) != canonical_unit_factor(header_unit, canonical):
                raise ValueError(f"Conflicting units for {source}: {explicit} vs {header_unit}")
        unit = explicit or header_unit or canonical
        factor = canonical_unit_factor(unit, canonical)
        out[logical] = pd.to_numeric(out[logical], errors="coerce") * factor
        status = "explicit" if explicit else "header" if header_unit else "assumed_canonical"
        metadata[logical] = {"source_column": source, "source_unit": unit,
                             "unit": canonical, "factor": factor, "status": status}
        if status == "assumed_canonical":
            warnings.append(f"unit_assumed:{logical}:{canonical}")

    if "charge_capacity" in out.columns:
        out["capacity"] = out["charge_capacity"]
        metadata["capacity"] = dict(metadata["charge_capacity"])
    elif "discharge_capacity" in out.columns:
        out["capacity"] = out["discharge_capacity"]
        metadata["capacity"] = dict(metadata["discharge_capacity"])

    numeric_cols = (
        "cycle", "voltage", "capacity", "charge_capacity", "discharge_capacity",
        "current", "time", "step_time", "data_point", "temperature",
    )
    for col in numeric_cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    if "cycle" in out.columns:
        out = out.dropna(subset=["cycle"])
        out["cycle"] = out["cycle"].astype(int)
    out.attrs["unit_metadata"] = metadata
    out.attrs["unit_warnings"] = list(dict.fromkeys(warnings))
    out.attrs["unit_schema"] = "canonical_A_Ah_V_s_v1"
    out.attrs["unit_missing_fields"] = [c for c in ("current", "voltage", "capacity", "time") if c not in out.columns]
    return out


def load_cycler_csv(
    path: str,
    column_map: ColumnMap | Mapping[str, Any] | None = None,
) -> pd.DataFrame:
    """Load CSV and rename to logical column names where possible."""
    df = pd.read_csv(path, on_bad_lines="skip")
    return normalize_cycler_dataframe(df, column_map)
