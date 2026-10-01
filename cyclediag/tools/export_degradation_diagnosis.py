"""CLI: full-cell LLI/LAM pattern diagnosis from cycle feature table or raw CSV."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cyclediag.diagnosis import diagnose_feature_table  # noqa: E402
from cyclediag.diagnosis.schema import (  # noqa: E402
    PATTERN_MODES,
    confidence_column_name,
    score_column_name,
)
from cyclediag.features.cycle_indicators_export import (  # noqa: E402
    export_cycle_indicators,
    summarize_cycle_indicators,
)
from cyclediag.features.lges_extract import LgesExtractConfig  # noqa: E402
from cyclediag.io.cycler_csv import ColumnMap  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Full-cell LLI/LAM/impedance pattern diagnosis (Level 1). "
            "Half-cell data is NOT required."
        ),
    )
    p.add_argument("--input", type=Path, required=True, help="Raw CSV, folder, or feature CSV")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--stem", type=str, default=None)
    p.add_argument(
        "--from-features",
        action="store_true",
        help="Treat --input as an existing cycle-feature CSV (skip re-extract)",
    )
    p.add_argument("--all-cycles", action="store_true", help="Use all raw cycles (default: tagged)")
    p.add_argument("--config", type=Path, default=None, help="mode_weights JSON path")
    p.add_argument("--no-json", action="store_true", help="Skip diagnosis JSON sidecar")
    args = p.parse_args()

    out_dir = args.out_dir or (
        args.input.parent if args.input.is_file() else args.input
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.stem or args.input.stem.replace("_raw", "").replace(
        "_cycle_indicators_tagged_full", ""
    ).replace("_cycle_indicators_tagged_inspect", "")

    json_path = None if args.no_json else out_dir / f"{stem}_diagnosis.json"

    if args.from_features:
        import pandas as pd

        features = pd.read_csv(args.input)
        features = diagnose_feature_table(
            features,
            config_path=args.config,
            write_json_sidecar=json_path,
        )
    else:
        cfg = LgesExtractConfig(
            with_diagnosis=True,
            diagnosis_config_path=str(args.config) if args.config else None,
        )
        result = export_cycle_indicators(
            args.input,
            out_dir,
            stem=stem,
            tagged_only=not args.all_cycles,
            write_csv=True,
            write_xlsx=False,
            write_png=False,
            column_map=ColumnMap.studio_default(),
            config=cfg,
        )
        features = result.features
        if json_path is not None and not features.empty:
            diagnose_feature_table(
                features,
                config_path=args.config,
                write_json_sidecar=json_path,
            )

    diag_csv = out_dir / f"{stem}_diagnosis_scores.csv"
    ordered = [
        c for c in (
            ["cell_id", "tagged_cycle", "cycle", "SoHQ"]
            + [score_column_name(m) for m in PATTERN_MODES]
            + [confidence_column_name(m) for m in PATTERN_MODES]
            + [c for c in features.columns if c.endswith(("_scientific_eligible", "_eligibility_reasons"))]
            + [
                "diagnosis_quality_score", "diagnosis_valid",
                "diagnosis_heuristic_support_valid", "diagnosis_status", "scientific_validity",
                "diagnosis_method", "diagnosis_model_version", "diagnosis_version",
            ]
        )
        if c in features.columns
    ]

    if not features.empty:
        features[ordered].to_csv(diag_csv, index=False)
        features.to_csv(out_dir / f"{stem}_cycle_indicators_with_diagnosis.csv", index=False)

    print(summarize_cycle_indicators(features))
    print(f"diagnosis_csv: {diag_csv}")
    if json_path is not None:
        print(f"diagnosis_json: {json_path}")

    if not features.empty and "LLI_pattern_score" in features.columns:
        tail = features.tail(max(1, len(features) // 10))
        summary = {
            m: float(tail[score_column_name(m)].mean())
            for m in PATTERN_MODES
            if score_column_name(m) in tail.columns
        }
        print("late_life_mean_scores:", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
