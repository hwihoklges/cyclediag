"""High-level diagnosis API (GUI-free; replaces Studio Diagnosis tab workflow)."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

from cyclediag.analysis.dqdv_screen import screen_dqdv_by_file, top_dqdv_problems
from cyclediag.analysis.indicator_screen import (
    compare_cells,
    screen_indicators_by_file,
    top_problem_indicators,
)
from cyclediag.features.lges_extract import LgesExtractConfig, extract_lges_features_table
from cyclediag.io.cycler_csv import ColumnMap, load_cycler_csv, normalize_cycler_dataframe
from cyclediag.models.indicator_scoring import score_indicators, top_scored_indicators
from cyclediag.models.predict import predict_features
from cyclediag.provenance import attach_provenance, build_provenance, config_hash


def extract_features(
    source: str | Path | pd.DataFrame,
    *,
    filepath: str = "",
    cycles: Iterable[int] | None = None,
    column_map: ColumnMap | None = None,
    config: LgesExtractConfig | None = None,
) -> pd.DataFrame:
    """Load (if needed) and extract LGES cycle indicators."""
    cmap = column_map or ColumnMap.studio_default()
    cfg = deepcopy(config) if config is not None else LgesExtractConfig()
    cycles = list(cycles) if cycles is not None else None
    if isinstance(source, pd.DataFrame):
        df = normalize_cycler_dataframe(source.copy(), column_map=cmap)
        path = filepath or str(cfg.cell_id or "dataframe")
    else:
        path = str(source)
        df = load_cycler_csv(path, column_map=cmap)
        if cfg.cell_id is None:
            cfg.cell_id = Path(path).stem
    provenance = build_provenance(source, df, config=cfg, column_map=cmap, cycles=cycles,
                                  selectors={"entry_point": "extract_features"})
    table = extract_lges_features_table(
        df,
        cycles=cycles,
        filepath=path,
        config=cfg,
        raw_df=df,
    )
    return attach_provenance(table, df, provenance)


def score_dataframe(
    features: pd.DataFrame,
    *,
    reference: pd.DataFrame | None = None,
    raw_df: pd.DataFrame | None = None,
    routine_only: bool = True,
    top_n: int = 15,
) -> dict[str, Any]:
    """Indicator scoring track — *how much* each indicator moved, not *why*.

    Separate from physicochemical causal diagnosis (``cyclediag.diagnosis`` /
    ``diagnose_feature_table``). Default ``routine_only=True`` excludes RPT /
    post-RPT / DC-IR spikes from the score.
    """
    from cyclediag.analysis.indicator_layers import split_by_layer

    result = score_indicators(
        features,
        reference=reference,
        raw_df=raw_df,
        routine_only=routine_only,
        grain="both",
    )
    layers = split_by_layer(result.indicator_summary)
    return {
        "score_layer": "indicator",
        "causal_track": "separate",
        "cycle_scores": result.cycle_scores,
        "cycle_contributions": result.cycle_contributions,
        "indicator_summary": result.indicator_summary,
        "top_indicators": top_scored_indicators(result.indicator_summary, n=top_n),
        "by_layer": layers,
        "meta": result.meta,
    }


def diagnose_dataframe(
    features: pd.DataFrame,
    *,
    reference: pd.DataFrame | None = None,
    with_screen: bool = True,
    routine_only: bool = True,
) -> dict[str, Any]:
    """Indicator scoring + optional descriptive screens.

    Name kept for compatibility. This is **not** the causal diagnosis track
    (LLI/LAM mode scores) — that lives in ``cyclediag.diagnosis`` and is
    invoked from extract via ``with_diagnosis=True``.
    """
    scored = predict_features(
        features, reference=reference, routine_only=routine_only,
    )
    indicator = score_dataframe(
        features, reference=reference, routine_only=routine_only,
    )
    out: dict[str, Any] = {
        "features": features,
        "scored": scored,
        "score_layer": "indicator",
        "indicator_summary": indicator["indicator_summary"],
        "indicator_top": indicator["top_indicators"],
        "indicator_screen": pd.DataFrame(),
        "top_indicators": indicator["top_indicators"],
        "dqdv_screen": pd.DataFrame(),
        "top_dqdv": pd.DataFrame(),
        "compare_cells": pd.DataFrame(),
        "meta": indicator["meta"],
    }
    if with_screen and features is not None and not features.empty:
        screened = screen_indicators_by_file(features)
        out["indicator_screen"] = screened
        # Prefer scored summary when available; fall back to descriptive screen.
        if out["top_indicators"] is None or out["top_indicators"].empty:
            out["top_indicators"] = top_problem_indicators(screened)
        dq = screen_dqdv_by_file(features)
        out["dqdv_screen"] = dq
        out["top_dqdv"] = top_dqdv_problems(dq)
        if "cell_id" in features.columns and features["cell_id"].nunique() >= 2:
            out["compare_cells"] = compare_cells(features)
    return out


def diagnose_csv(
    path: str | Path,
    *,
    cycles: Iterable[int] | None = None,
    column_map: ColumnMap | None = None,
    config: LgesExtractConfig | None = None,
    with_screen: bool = True,
) -> dict[str, Any]:
    """End-to-end: CSV → features → anomaly / screens."""
    feats = extract_features(
        path,
        cycles=cycles,
        column_map=column_map,
        config=config,
    )
    result = diagnose_dataframe(feats, with_screen=with_screen)
    provenance = deepcopy(feats.attrs["provenance"])
    provenance["extraction_config_sha256"] = provenance["config_sha256"]
    provenance["config_sha256"] = config_hash({
        "extraction": provenance["extraction_config_sha256"],
        "with_screen": with_screen, "routine_only": True, "entry_point": "diagnose_csv",
    })
    result["provenance"] = provenance
    return result


def diagnose_folder(
    input_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    encoding: str = "cp949",
    top_n: int = 12,
    write_pngs: bool = True,
) -> Mapping[str, Any]:
    """Batch report over a folder (StepEnd or cycler CSVs)."""
    from cyclediag.analysis.batch_report import run_batch_report

    return run_batch_report(
        input_dir,
        output_dir=output_dir,
        encoding=encoding,
        top_n=top_n,
        write_pngs=write_pngs,
    )


def screen_problems(features: pd.DataFrame, *, n: int = 15) -> pd.DataFrame:
    """Convenience: ranked problem indicators for one feature table."""
    screened = screen_indicators_by_file(features)
    return top_problem_indicators(screened, n=n)


def compare_doe(
    *,
    doe: str = "DOE2",
    fixtures_root: str | Path | None = None,
    out_dir: str | Path | None = None,
    early_cycles: int = 30,
    run_diagnosis: bool = True,
    write_plots: bool = True,
) -> Mapping[str, Any]:
    """Compare DOE arms (DOE2 = SJ900 vs SJ1300, same cathode / different anode)."""
    from cyclediag.analysis.doe_compare import DoeCompareConfig, run_doe_compare

    return run_doe_compare(
        DoeCompareConfig(
            doe_id=doe,
            fixtures_root=Path(fixtures_root) if fixtures_root else None,
            out_dir=Path(out_dir) if out_dir else None,
            early_cycles=early_cycles,
            run_diagnosis=run_diagnosis,
            write_plots=write_plots,
        )
    )
