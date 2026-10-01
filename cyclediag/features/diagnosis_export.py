"""Unified offline diagnosis export (tables + PNG report bundle)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from cyclediag.analysis.indicator_screen import screen_indicators_by_file
from cyclediag.analysis.indicator_screen_plots import plot_sohq_correlation_report
from cyclediag.analysis.sohq_inflection import plot_sohq_inflection_report
from cyclediag.diagnosis import diagnose_feature_table
from cyclediag.diagnosis.schema import PATTERN_MODES, confidence_column_name, score_column_name
from cyclediag.features.cycle_indicators_export import (
    INSPECT_COLS,
    _select_cols,
    discover_raw_csvs,
    extract_cycle_indicators_many,
    write_cycle_indicator_workbook,
)
from cyclediag.features.cycle_indicators_plots import (
    plot_sohq_rest_v_linear_proxy,
    save_cycle_indicator_pngs,
)
from cyclediag.features.lges_extract import LgesExtractConfig
from cyclediag.io.cycler_csv import ColumnMap


@dataclass
class DiagnosisExportResult:
    """Full diagnosis bundle: tables, screening, and PNG reports."""

    features: pd.DataFrame
    inspect: pd.DataFrame
    screened: pd.DataFrame
    out_xlsx: Path | None = None
    out_csv: Path | None = None
    out_screen_csv: Path | None = None
    out_diagnosis_csv: Path | None = None
    out_diagnosis_json: Path | None = None
    out_pngs: list[Path] = field(default_factory=list)


def _write_degradation_diagnosis_outputs(
    features: pd.DataFrame,
    out_dir: Path,
    stem: str,
) -> tuple[pd.DataFrame, Path | None, Path | None]:
    """Ensure Level-1 diagnosis columns, write score CSV + JSON sidecar."""
    if features is None or features.empty:
        return features, None, None
    json_path = out_dir / f"{stem}_degradation_diagnosis.json"
    work = diagnose_feature_table(features, write_json_sidecar=json_path)
    cols = [
        c for c in (
            ["cell_id", "tagged_cycle", "cycle", "SoHQ"]
            + [score_column_name(m) for m in PATTERN_MODES]
            + [confidence_column_name(m) for m in PATTERN_MODES]
            + [c for c in work.columns if c.endswith(("_scientific_eligible", "_eligibility_reasons"))]
            + [
                "diagnosis_quality_score", "diagnosis_valid",
                "diagnosis_heuristic_support_valid", "diagnosis_status", "scientific_validity",
                "diagnosis_method", "diagnosis_model_version", "diagnosis_version",
            ]
        )
        if c in work.columns
    ]
    csv_path = out_dir / f"{stem}_degradation_diagnosis.csv"
    work[cols].to_csv(csv_path, index=False)
    return work, csv_path, json_path


def save_diagnosis_pngs(
    features: pd.DataFrame,
    out_dir: Path | str,
    *,
    stem: str,
    screened: pd.DataFrame | None = None,
    per_cell: bool = True,
    include_overview: bool = True,
    include_sohq_proxy: bool = True,
    include_sohq_corr: bool = True,
    include_sohq_inflection: bool = True,
) -> list[Path]:
    """Write overview / SoHQ proxy / corr / inflection PNGs for one export stem."""
    if features is None or features.empty:
        return []

    out_dir = Path(out_dir)
    saved: list[Path] = []

    if include_overview:
        saved.extend(
            save_cycle_indicator_pngs(
                features, out_dir, stem=stem, per_cell=per_cell,
            )
        )

    if include_sohq_proxy and {"EoC_restV_end", "EoD_restV_end", "SoHQ"}.issubset(features.columns):
        try:
            proxy_path = out_dir / f"{stem}_sohq_rest_v_linear_proxy.png"
            plot_sohq_rest_v_linear_proxy(features, proxy_path)
            saved.append(proxy_path)
        except ValueError:
            pass

    if include_sohq_corr and "SoHQ" in features.columns:
        try:
            corr_path = out_dir / f"{stem}_sohq_correlation_report.png"
            plot_sohq_correlation_report(
                features,
                corr_path,
                screened=screened,
            )
            if corr_path.is_file():
                saved.append(corr_path)
        except ValueError:
            pass

    if include_sohq_inflection and "SoHQ" in features.columns:
        try:
            infl_path = out_dir / f"{stem}_sohq_inflection_regimes.png"
            out, _ = plot_sohq_inflection_report(features, infl_path)
            if out is not None and Path(out).is_file():
                saved.append(Path(out))
        except ValueError:
            pass

    return saved


def export_diagnosis_bundle(
    input_path: Path | str,
    out_dir: Path | str | None = None,
    *,
    stem: str | None = None,
    cycles: Iterable[int] | None = None,
    tagged_only: bool = True,
    write_csv: bool = True,
    write_xlsx: bool = True,
    write_screen_csv: bool = True,
    write_png: bool = True,
    per_cell_png: bool = True,
    include_sohq_proxy: bool = True,
    include_sohq_corr: bool = True,
    column_map: ColumnMap | None = None,
    config: LgesExtractConfig | None = None,
) -> DiagnosisExportResult:
    """Extract indicators, screen vs SoHQ, and write tables + PNG report bundle."""
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
    out_dir = Path(out_dir) if out_dir else (input_path.parent if input_path.is_file() else input_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    base = stem or (
        input_path.stem if input_path.is_file() else f"{input_path.name}_diagnosis"
    )
    base = base.replace("_raw", "")
    if tagged_only:
        base += "_tagged"

    features, out_diagnosis_csv, out_diagnosis_json = _write_degradation_diagnosis_outputs(
        features, out_dir, base,
    )
    inspect = _select_cols(features, INSPECT_COLS)
    screened = screen_indicators_by_file(features) if not features.empty else pd.DataFrame()

    out_xlsx = out_csv = out_screen_csv = None
    out_pngs: list[Path] = []

    if write_xlsx and not features.empty:
        out_xlsx = write_cycle_indicator_workbook(features, out_dir / f"{base}_indicators.xlsx")
    if write_csv and not features.empty:
        out_csv = out_dir / f"{base}_indicators_inspect.csv"
        inspect.to_csv(out_csv, index=False)
    if write_screen_csv and not screened.empty:
        out_screen_csv = out_dir / f"{base}_indicator_screen.csv"
        screened.to_csv(out_screen_csv, index=False)

    if write_png and not features.empty:
        out_pngs = save_diagnosis_pngs(
            features,
            out_dir,
            stem=f"{base}_indicators",
            screened=screened,
            per_cell=per_cell_png,
            include_overview=True,
            include_sohq_proxy=include_sohq_proxy,
            include_sohq_corr=include_sohq_corr,
            include_sohq_inflection=True,
        )

    return DiagnosisExportResult(
        features=features,
        inspect=inspect,
        screened=screened,
        out_xlsx=out_xlsx,
        out_csv=out_csv,
        out_screen_csv=out_screen_csv,
        out_diagnosis_csv=out_diagnosis_csv,
        out_diagnosis_json=out_diagnosis_json,
        out_pngs=out_pngs,
    )


def export_diagnosis_from_features(
    features: pd.DataFrame,
    out_dir: Path | str,
    *,
    stem: str = "diagnosis",
    write_csv: bool = True,
    write_xlsx: bool = True,
    write_screen_csv: bool = True,
    write_png: bool = True,
    per_cell_png: bool = True,
    include_sohq_proxy: bool = True,
    include_sohq_corr: bool = True,
) -> DiagnosisExportResult:
    """Write diagnosis bundle from an already-extracted feature table (GUI path)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    features, out_diagnosis_csv, out_diagnosis_json = _write_degradation_diagnosis_outputs(
        features, out_dir, stem,
    )
    inspect = _select_cols(features, INSPECT_COLS)
    screened = screen_indicators_by_file(features) if not features.empty else pd.DataFrame()

    out_xlsx = out_csv = out_screen_csv = None
    out_pngs: list[Path] = []

    if write_xlsx and not features.empty:
        out_xlsx = write_cycle_indicator_workbook(features, out_dir / f"{stem}_indicators.xlsx")
    if write_csv and not features.empty:
        out_csv = out_dir / f"{stem}_indicators_inspect.csv"
        inspect.to_csv(out_csv, index=False)
    if write_screen_csv and not screened.empty:
        out_screen_csv = out_dir / f"{stem}_indicator_screen.csv"
        screened.to_csv(out_screen_csv, index=False)
    if write_png and not features.empty:
        out_pngs = save_diagnosis_pngs(
            features,
            out_dir,
            stem=f"{stem}_indicators",
            screened=screened,
            per_cell=per_cell_png,
            include_sohq_proxy=include_sohq_proxy,
            include_sohq_corr=include_sohq_corr,
            include_sohq_inflection=True,
        )

    return DiagnosisExportResult(
        features=features,
        inspect=inspect,
        screened=screened,
        out_xlsx=out_xlsx,
        out_csv=out_csv,
        out_screen_csv=out_screen_csv,
        out_diagnosis_csv=out_diagnosis_csv,
        out_diagnosis_json=out_diagnosis_json,
        out_pngs=out_pngs,
    )
