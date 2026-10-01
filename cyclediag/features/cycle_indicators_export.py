"""Export cycle-level degradation indicators for offline inspection (no GUI)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from cyclediag.features.cycle_indicators_plots import save_cycle_indicator_pngs
from cyclediag.features.lges_extract import LgesExtractConfig, extract_lges_features_table
from cyclediag.io.classification_pairs import (
    TaggedCycle,
    baseline_raw_cycle_for_tagged,
    resolve_tagged_cycles_for_raw,
)
from cyclediag.io.cycler_csv import ColumnMap, load_cycler_csv

# Compact review columns (rest V + R + capacity health)
INSPECT_COLS: tuple[str, ...] = (
    "cell_id",
    "tagged_cycle",
    "pair_label",
    "cycle",
    "dchgCapa",
    "chgCapa",
    "SoHQ",
    "CE",
    "chgCapa_CCratio",
    "chgCVtime",
    "EoC_restV_init",
    "EoC_restV_60s",
    "EoC_restV_30m",
    "EoC_restV_end",
    "EoD_restV_init",
    "EoD_restV_60s",
    "EoD_restV_30m",
    "EoD_restV_end",
    "EoC_dchgR_0p1s",
    "EoC_dchgR_10s",
    "EoC_dchgR_30s",
    "EoC_dchgR_60s",
    "EoC_dchgR_R10_minus_R0p1",
    "EoC_dchgR_R30_minus_R0p1",
    "EoD_chgR_0p1s",
    "EoD_chgR_10s",
    "EoD_chgR_30s",
    "EoD_chgR_60s",
    "EoD_chgR_R10_minus_R0p1",
    "EoD_chgR_R30_minus_R0p1",
    "dchg_Q_low_frac",
    "dchg_Q_high_frac",
    "dchg_f_graphite_proxy",
    "delta_EoC_restV_end",
    "delta_EoD_restV_end",
    "EoC_dchgR_10s_inc",
    "EoD_chgR_10s_inc",
    "dchg_dVdQ_SOC0",
    "dchg_dVdQ_SOC5",
    "dchg_dVdQ_SOC10",
    "dchg_dVdQ_SOCmid",
    "dchg_dVdQ_SOC0_cliff_width",
    "dchg_dVdQ_SOC0_to_mid_ratio",
    "chg_dVdQ_SOC100",
    "EoC_restV_relax",
    "EoD_restV_relax",
    "EoC_restV_tau",
    "EoD_restV_tau",
    "EoC_dchgR_10_60_ratio",
    "EoD_chgR_10_60_ratio",
    "chg_V_avg",
    "dchg_V_avg",
    "chg_E",
    "dchg_E",
    "hyst_area",
    "hyst_max_dV",
    "chg_ir_drop_proxy",
    "dchg_ir_drop_proxy",
    "dchg_plateau_V",
    "dchg_plateau_width",
    "dchg_dQdV_area_sum",
    "dchg_V_cutoff_margin",
    "dchg_shape_DTW",
    "dSoHQ_dN",
    "d2SoHQ",
    "EoC_dchgR_10s_growth_50",
    "CE_local_20",
    "delta_dchg_dQdV_peak1_V",
    "EoC_dchgR_10s_T25",
    "LLI_pattern_score",
    "LAM_PE_pattern_score",
    "LAM_NE_pattern_score",
    "impedance_pattern_score",
    "transport_limitation_score",
    "plating_risk_score",
    "contact_loss_score",
    "LLI_confidence",
    "LAM_PE_confidence",
    "LAM_NE_confidence",
    "diagnosis_quality_score",
    "diagnosis_valid",
    "diagnosis_heuristic_support_valid",
    "diagnosis_status",
    "scientific_validity",
    "diagnosis_version",
)


@dataclass
class CycleIndicatorExportResult:
    features: pd.DataFrame
    inspect: pd.DataFrame
    out_xlsx: Path | None = None
    out_csv: Path | None = None
    out_pngs: list[Path] = field(default_factory=list)


def _select_cols(df: pd.DataFrame, cols: Sequence[str]) -> pd.DataFrame:
    present = [c for c in cols if c in df.columns]
    return df[present].copy()


def discover_raw_csvs(path: Path) -> list[Path]:
    """Return one or more raw CSV paths from a file or folder."""
    path = Path(path)
    if path.is_file():
        return [path]
    if not path.is_dir():
        raise FileNotFoundError(path)
    found = sorted(path.glob("*_raw.csv"))
    if not found:
        found = sorted(path.glob("*.csv"))
    return found


def _attach_tagged_metadata(
    features: pd.DataFrame,
    tagged: list[TaggedCycle],
) -> pd.DataFrame:
    if features.empty or not tagged:
        return features.iloc[0:0].copy()
    by_raw = {int(t.raw_cycle): t for t in tagged}
    out = features.copy()
    out["cycle"] = pd.to_numeric(out["cycle"], errors="coerce")
    out["tagged_cycle"] = out["cycle"].map(
        lambda c: by_raw[int(c)].tagged_cycle if pd.notna(c) and int(c) in by_raw else None
    )
    out["pair_label"] = out["cycle"].map(
        lambda c: by_raw[int(c)].pair_label if pd.notna(c) and int(c) in by_raw else None
    )
    out["tagged_source"] = out["cycle"].map(
        lambda c: by_raw[int(c)].source if pd.notna(c) and int(c) in by_raw else None
    )
    out = out[out["tagged_cycle"].notna()].copy()
    return out.sort_values(["cell_id", "tagged_cycle"], na_position="last")


def extract_cycle_indicators(
    csv_path: Path | str,
    *,
    cell_id: str | None = None,
    cycles: Iterable[int] | None = None,
    tagged_only: bool = False,
    column_map: ColumnMap | None = None,
    config: LgesExtractConfig | None = None,
) -> pd.DataFrame:
    """Load one cycler CSV and return one row per cycle (LGES indicators)."""
    csv_path = Path(csv_path)
    cmap = column_map or ColumnMap.studio_default()
    cfg = config or LgesExtractConfig()
    if cell_id:
        cfg.cell_id = cell_id
    elif cfg.cell_id is None:
        cfg.cell_id = csv_path.stem.replace("_raw", "")

    tagged_rows: list[TaggedCycle] = []
    cycle_list = list(cycles) if cycles is not None else None
    if tagged_only:
        tagged_rows = resolve_tagged_cycles_for_raw(csv_path)
        if not tagged_rows:
            return pd.DataFrame()
        cycle_list = [t.raw_cycle for t in tagged_rows]
        bl = baseline_raw_cycle_for_tagged(csv_path)
        if bl is not None:
            cfg.baseline_cycle = bl

    df = load_cycler_csv(csv_path, column_map=cmap)
    table = extract_lges_features_table(
        df,
        cycles=cycle_list,
        filepath=str(csv_path),
        config=cfg,
        raw_df=df,
    )
    if tagged_only and not table.empty:
        table = _attach_tagged_metadata(table, tagged_rows)
    return table


def extract_cycle_indicators_many(
    paths: Sequence[Path | str],
    *,
    cycles: Iterable[int] | None = None,
    tagged_only: bool = False,
    column_map: ColumnMap | None = None,
    config: LgesExtractConfig | None = None,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for p in paths:
        table = extract_cycle_indicators(
            p,
            cycles=cycles,
            tagged_only=tagged_only,
            column_map=column_map,
            config=config,
        )
        if not table.empty:
            frames.append(table)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def write_cycle_indicator_workbook(
    features: pd.DataFrame,
    out_path: Path | str,
    *,
    inspect_cols: Sequence[str] = INSPECT_COLS,
) -> Path:
    """Write Excel with Inspect / Full / Meta sheets."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    inspect = _select_cols(features, inspect_cols)

    meta = pd.DataFrame(
        [
            {"item": "n_rows", "value": len(features)},
            {"item": "n_cells", "value": int(features["cell_id"].nunique()) if "cell_id" in features else 0},
            {"item": "cycle_min", "value": features["cycle"].min() if "cycle" in features and not features.empty else None},
            {"item": "cycle_max", "value": features["cycle"].max() if "cycle" in features and not features.empty else None},
            {"item": "inspect_cols", "value": ", ".join(inspect.columns)},
        ]
    )

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        meta.to_excel(writer, sheet_name="Meta", index=False)
        inspect.to_excel(writer, sheet_name="Inspect", index=False)
        features.to_excel(writer, sheet_name="Full", index=False)

        # Per-cell quick sheets (cap to avoid huge workbooks)
        if "cell_id" in inspect.columns:
            for i, (cid, grp) in enumerate(inspect.groupby("cell_id", sort=False)):
                if i >= 20:
                    break
                safe = str(cid)[:28].replace("/", "_").replace("\\", "_")
                grp.to_excel(writer, sheet_name=safe or f"cell{i}", index=False)

    return out_path


def summarize_cycle_indicators(features: pd.DataFrame) -> str:
    """Short text summary for console inspection."""
    if features is None or features.empty:
        return "No features extracted."

    lines = [
        f"rows={len(features)}",
        f"cells={features['cell_id'].nunique() if 'cell_id' in features else '?'}",
    ]
    if "tagged_cycle" in features.columns and features["tagged_cycle"].notna().any():
        lines.append(
            f"tagged_cycles={int(features['tagged_cycle'].min())}-"
            f"{int(features['tagged_cycle'].max())}"
        )
        if "tagged_source" in features.columns:
            src = features["tagged_source"].dropna().iloc[0] if features["tagged_source"].notna().any() else "?"
            lines.append(f"tagged_source={src}")
    elif "cycle" in features.columns:
        lines.append(
            f"raw_cycles={int(features['cycle'].min())}-{int(features['cycle'].max())}"
        )
    for col in ("SoHQ", "EoC_restV_end", "EoD_restV_end", "EoC_dchgR_10s", "EoD_chgR_10s", "CE"):
        if col not in features.columns:
            continue
        s = pd.to_numeric(features[col], errors="coerce")
        if s.notna().any():
            lines.append(
                f"{col}: first={s.dropna().iloc[0]:.4g}, "
                f"last={s.dropna().iloc[-1]:.4g}, "
                f"coverage={100 * s.notna().mean():.0f}%"
            )
    for col in (
        "LLI_pattern_score", "LAM_PE_pattern_score", "LAM_NE_pattern_score",
        "impedance_pattern_score",
    ):
        if col not in features.columns:
            continue
        s = pd.to_numeric(features[col], errors="coerce")
        if s.notna().any():
            lines.append(
                f"{col}: mean={s.mean():.3f}, last={s.dropna().iloc[-1]:.3f}"
            )
    if "diagnosis_version" in features.columns and features["diagnosis_version"].notna().any():
        lines.append(f"diagnosis_version={features['diagnosis_version'].dropna().iloc[0]}")
    return "\n".join(x for x in lines if x)


def export_cycle_indicators(
    input_path: Path | str,
    out_dir: Path | str | None = None,
    *,
    stem: str | None = None,
    cycles: Iterable[int] | None = None,
    tagged_only: bool = True,
    write_csv: bool = True,
    write_xlsx: bool = True,
    write_png: bool = True,
    per_cell_png: bool = True,
    column_map: ColumnMap | None = None,
    config: LgesExtractConfig | None = None,
) -> CycleIndicatorExportResult:
    """Extract indicators from file/folder and write Inspect Excel / CSV / PNG."""
    input_path = Path(input_path)
    paths = discover_raw_csvs(input_path)
    if not paths:
        raise FileNotFoundError(f"No CSV found under {input_path}")

    features = extract_cycle_indicators_many(
        paths,
        cycles=cycles,
        tagged_only=tagged_only,
        column_map=column_map,
        config=config,
    )
    inspect = _select_cols(features, INSPECT_COLS)

    out_dir = Path(out_dir) if out_dir else (input_path.parent if input_path.is_file() else input_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = stem or (
        input_path.stem if input_path.is_file() else f"{input_path.name}_cycle_indicators"
    )
    base = base.replace("_raw", "") + "_cycle_indicators"
    if tagged_only:
        base += "_tagged"

    out_xlsx = out_csv = None
    out_pngs: list[Path] = []
    if write_xlsx and not features.empty:
        out_xlsx = write_cycle_indicator_workbook(features, out_dir / f"{base}.xlsx")
    if write_csv and not features.empty:
        out_csv = out_dir / f"{base}_inspect.csv"
        inspect.to_csv(out_csv, index=False)
    if write_png and not features.empty:
        out_pngs = save_cycle_indicator_pngs(
            features, out_dir, stem=base, per_cell=per_cell_png,
        )

    return CycleIndicatorExportResult(
        features=features,
        inspect=inspect,
        out_xlsx=out_xlsx,
        out_csv=out_csv,
        out_pngs=out_pngs,
    )
