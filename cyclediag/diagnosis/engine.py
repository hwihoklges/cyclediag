"""Apply degradation-mode diagnosis to a cycle feature table.

Track B — physicochemical causal interpretation. Separate from the indicator
scoring track in ``cyclediag.models.indicator_scoring``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from cyclediag.models.indicator_scoring import filter_scoring_rows

from .constraints import constraint_flags, scientific_validity
from .pattern_scoring import load_mode_weights, score_all_modes_for_row
from .schema import (
    DIAGNOSIS_MODEL_VERSION,
    DIAGNOSIS_VERSION_FULLCELL,
    PATTERN_MODES,
    confidence_column_name,
    score_column_name,
)

META_COLS = (
    "diagnosis_quality_score",
    "diagnosis_valid",
    "diagnosis_method",
    "diagnosis_model_version",
    "diagnosis_version",
    "diagnosis_constraints",
    "diagnosis_scored_row",
    "diagnosis_status",
    "scientific_validity",
)


def diagnosis_wide_columns(modes: tuple[str, ...] = PATTERN_MODES) -> list[str]:
    cols: list[str] = []
    for mode in modes:
        cols.append(score_column_name(mode))
        cols.append(confidence_column_name(mode))
        cols.append(f"{mode}_supporting_features")
        cols.append(f"{mode}_conflicting_features")
        cols.append(f"{mode}_evidence_count")
        cols.append(f"{mode}_evidence_coverage")
        cols.append(f"{mode}_status")
    cols.extend(META_COLS)
    # Level-2 placeholders (null until validated models exist)
    cols.extend(["LLI_est", "LAM_PE_est", "LAM_NE_est", "electrode_slippage_est"])
    # Level-3 placeholders
    cols.extend([
        "LLI_est_hc_calibrated",
        "LAM_PE_est_hc_calibrated",
        "LAM_NE_est_hc_calibrated",
    ])
    return cols


def _baseline_row(grp: pd.DataFrame, baseline_cycle: int | None) -> dict[str, Any]:
    if grp.empty:
        return {}
    if baseline_cycle is not None and "cycle" in grp.columns:
        hit = grp[grp["cycle"] == baseline_cycle]
        if not hit.empty:
            return hit.iloc[0].to_dict()
        return {}  # An explicitly requested reference must not be fabricated.
    return grp.sort_values("cycle").iloc[0].to_dict()


def diagnose_feature_table(
    table: pd.DataFrame,
    *,
    config_path: str | Path | None = None,
    baseline_cycle: int | None = 1,
    write_json_sidecar: Path | str | None = None,
    routine_only: bool = True,
) -> pd.DataFrame:
    """Append Level-1 pattern diagnosis columns to a cycle feature DataFrame.

    Does **not** require half-cell data. Level-2/3 estimate columns are left null.

    By default only routine (non-excluded) cycles are scored so RPT / DC-IR
    spikes do not drive heuristic mode support. Non-scored rows keep NaN scores and
    ``diagnosis_scored_row=False``.
    """
    if table is None or table.empty:
        out = table.copy() if table is not None else pd.DataFrame()
        cfg0 = load_mode_weights(config_path) if config_path else load_mode_weights()
        modes0 = tuple(cfg0.get("modes", {}).keys()) or PATTERN_MODES
        for c in diagnosis_wide_columns(modes0):
            if c not in out.columns:
                out[c] = None
        return out

    cfg = load_mode_weights(config_path) if config_path else load_mode_weights()
    mode_list = tuple(cfg.get("modes", {}).keys()) or PATTERN_MODES
    out = table.copy()
    n = len(out)
    for c in diagnosis_wide_columns(mode_list):
        if c.endswith(("_features", "_status")) or c in (
            "diagnosis_method", "diagnosis_model_version", "diagnosis_version",
            "diagnosis_constraints", "scientific_validity",
        ):
            out[c] = pd.Series([None] * n, index=out.index, dtype=object)
        elif c in ("diagnosis_valid", "diagnosis_scored_row"):
            out[c] = False
        else:
            out[c] = np.nan

    scored_pool = filter_scoring_rows(out, routine_only=routine_only)
    scored_idx = set(scored_pool.index)

    group_cols = [c for c in ("cell_id", "file", "source_file") if c in out.columns]
    if group_cols:
        groups = list(out.groupby(group_cols, sort=False, dropna=False))
    else:
        groups = [(("__all__",), out)]

    sidecar_rows: list[dict[str, Any]] = []

    for _, grp in groups:
        # baseline from routine/scored rows when possible
        base_src = grp.loc[grp.index.isin(scored_idx)] if scored_idx else grp
        base = _baseline_row(base_src if not base_src.empty else grp, baseline_cycle)
        for idx, row in grp.iterrows():
            row_dict = row.to_dict()
            cflags = constraint_flags(row_dict, cfg)
            out.at[idx, "scientific_validity"] = scientific_validity(row_dict, cfg)
            out.at[idx, "diagnosis_status"] = "unknown"
            out.at[idx, "diagnosis_constraints"] = ",".join(cflags) if cflags else ""

            if idx not in scored_idx:
                out.at[idx, "diagnosis_scored_row"] = False
                out.at[idx, "diagnosis_valid"] = False
                out.at[idx, "diagnosis_method"] = str(cfg.get("diagnosis_method", "rule_pattern"))
                out.at[idx, "diagnosis_model_version"] = str(
                    cfg.get("diagnosis_model_version", DIAGNOSIS_MODEL_VERSION)
                )
                out.at[idx, "diagnosis_version"] = str(
                    cfg.get("diagnosis_version", DIAGNOSIS_VERSION_FULLCELL)
                )
                for c in (
                    "LLI_est", "LAM_PE_est", "LAM_NE_est", "electrode_slippage_est",
                    "LLI_est_hc_calibrated", "LAM_PE_est_hc_calibrated", "LAM_NE_est_hc_calibrated",
                ):
                    out.at[idx, c] = None
                continue

            results = score_all_modes_for_row(
                row_dict, cfg, baseline_row=base, modes=mode_list,
            )
            qualities = []
            valids = []
            for mode, res in results.items():
                out.at[idx, score_column_name(mode)] = res.estimate
                out.at[idx, confidence_column_name(mode)] = res.confidence
                out.at[idx, f"{mode}_supporting_features"] = ",".join(res.supporting_features)
                out.at[idx, f"{mode}_conflicting_features"] = ",".join(res.conflicting_features)
                out.at[idx, f"{mode}_evidence_count"] = res.evidence_count
                out.at[idx, f"{mode}_evidence_coverage"] = res.evidence_coverage
                out.at[idx, f"{mode}_status"] = res.status
                qualities.append(res.data_quality_score)
                valids.append(res.diagnosis_valid)
                sidecar_rows.append({
                    "cycle": row_dict.get("cycle"),
                    "tagged_cycle": row_dict.get("tagged_cycle"),
                    "cell_id": row_dict.get("cell_id"),
                    "constraints": cflags,
                    **res.to_dict(),
                })

            out.at[idx, "diagnosis_quality_score"] = float(np.nanmean(qualities)) if qualities else 0.0
            validity = out.at[idx, "scientific_validity"]
            validity["warnings"] = list(dict.fromkeys(
                warning for result in results.values()
                for warning in result.scientific_validity["warnings"]
            ))
            out.at[idx, "diagnosis_valid"] = bool(any(valids))
            out.at[idx, "diagnosis_status"] = "heuristic_support" if any(valids) else "unknown"
            out.at[idx, "diagnosis_scored_row"] = True
            out.at[idx, "diagnosis_method"] = str(cfg.get("diagnosis_method", "rule_pattern"))
            out.at[idx, "diagnosis_model_version"] = str(
                cfg.get("diagnosis_model_version", DIAGNOSIS_MODEL_VERSION)
            )
            out.at[idx, "diagnosis_version"] = str(
                cfg.get("diagnosis_version", DIAGNOSIS_VERSION_FULLCELL)
            )
            # Level 2/3 intentionally null (not validated absolute estimates)
            for c in (
                "LLI_est", "LAM_PE_est", "LAM_NE_est", "electrode_slippage_est",
                "LLI_est_hc_calibrated", "LAM_PE_est_hc_calibrated", "LAM_NE_est_hc_calibrated",
            ):
                out.at[idx, c] = None

    if write_json_sidecar is not None:
        path = Path(write_json_sidecar)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(sidecar_rows, indent=2, ensure_ascii=False), encoding="utf-8")

    return out
