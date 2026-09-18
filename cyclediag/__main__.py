"""CLI: python -m cyclediag …"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from cyclediag import __version__


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="cyclediag",
        description="Standalone cycle / voltage-profile diagnosis (cyclediag replacement)",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command")

    sub.add_parser("planning", help="Print path to planning docs")

    extract = sub.add_parser("extract", help="Raw CSV → feature table")
    extract.add_argument("--input", required=True)
    extract.add_argument("--out", required=True)
    extract.add_argument("--charge-step", default="charge")
    extract.add_argument("--discharge-step", default="discharge")
    extract.add_argument("--column-map", choices=("studio", "pne"), default="studio")
    extract.add_argument("--cv-only", action="store_true", help="CV regions table only")
    extract.add_argument(
        "--feature-set",
        choices=("vp_v1_basic", "vp_lges_cycle_v1", "vp_lges_cycle_v2"),
        default="vp_lges_cycle_v2",
        help="Feature catalog (default: LGES cycle v2)",
    )

    diagnose = sub.add_parser("diagnose", help="CSV → features + anomaly + screens")
    diagnose.add_argument("--input", required=True)
    diagnose.add_argument("--out-dir", required=True)
    diagnose.add_argument("--charge-step", default="charge")
    diagnose.add_argument("--discharge-step", default="discharge")
    diagnose.add_argument("--column-map", choices=("studio", "pne"), default="studio")
    diagnose.add_argument("--no-screen", action="store_true")

    peaks = sub.add_parser("peaks", help="dQ/dV peak feature pipeline")
    peaks_sub = peaks.add_subparsers(dest="peaks_command", required=True)

    peaks_export = peaks_sub.add_parser("export", help="Raw CSV → peak feature tables")
    peaks_export.add_argument("--input", required=True)
    peaks_export.add_argument("--out-dir", required=True)
    peaks_export.add_argument("--cell-id", required=True)
    peaks_export.add_argument(
        "--good-cycles",
        default="10,50,80,163,210,255,283,327,426",
        help="Reference cycles for delta features",
    )
    peaks_export.add_argument("--relaxed", action="store_true")
    peaks_export.add_argument("--usable-mad-factor", type=float, default=2.0)
    peaks_export.add_argument("--max-noise", type=float, default=0.008)
    peaks_export.add_argument("--max-charge-hf", type=float, default=0.68)
    peaks_export.add_argument("--max-discharge-hf", type=float, default=0.58)
    peaks_export.add_argument("--max-band-gap", type=int, default=0)
    peaks_export.add_argument("--min-usable-score", type=float, default=0.0)

    peaks_evo = peaks_sub.add_parser("evolution", help="2D evolution map + Viterbi ridge tracking")
    peaks_evo.add_argument("--input", required=True)
    peaks_evo.add_argument("--stepend", default="", help="Step-end CSV for protocol/rate filtering")
    peaks_evo.add_argument("--out-dir", required=True)
    peaks_evo.add_argument("--domain", choices=("Q", "V"), default="Q")
    peaks_evo.add_argument("--leg", choices=("discharge", "charge"), default="discharge")
    peaks_evo.add_argument("--normalize", choices=("none", "area", "capacity", "local_contrast"), default="local_contrast")
    peaks_evo.add_argument("--sg-window", type=int, default=7)
    peaks_evo.add_argument("--n-interp", type=int, default=2500)
    peaks_evo.add_argument("--n-grid", type=int, default=1000)
    peaks_evo.add_argument("--lam-scale", type=float, default=0.1)
    peaks_evo.add_argument("--max-tracks", type=int, default=12)
    peaks_evo.add_argument("--no-roi", action="store_true", help="Disable ROI-independent extraction")
    peaks_evo.add_argument("--no-preflight", action="store_true")
    peaks_evo.add_argument("--abort-on-preflight-fail", action="store_true")

    train = sub.add_parser("train", help="Train peak ML model (Isolation Forest)")
    train.add_argument("--features", required=True)
    train.add_argument("--out", required=True)
    train.add_argument("--good-cycles", default="10,50,80,163,210,255,283,327")
    train.add_argument("--train-on", choices=("good_cycles", "usable", "all_complete"), default="usable")
    train.add_argument("--no-require-usable", action="store_true")

    predict = sub.add_parser("predict", help="Feature table → anomaly scores")
    predict.add_argument("--features", required=True)
    predict.add_argument("--out", required=True)
    predict.add_argument("--reference", default="")
    predict.add_argument("--model", default="", help="Trained peak ML model directory")

    report = sub.add_parser("report", help="Batch diagnosis report (folder)")
    report.add_argument("--input-dir", required=True)
    report.add_argument("--output-dir", default="")
    report.add_argument("--encoding", default="cp949")
    report.add_argument("--top-n", type=int, default=12)

    compare_doe = sub.add_parser(
        "compare-doe",
        help="Compare DOE arms (DOE2: SJ900 vs SJ1300, same cathode / different anode)",
    )
    compare_doe.add_argument("--doe", default="DOE2", help="DOE id (default DOE2)")
    compare_doe.add_argument(
        "--fixtures-root",
        default="",
        help="Path to example/fixtures (default: auto-detect)",
    )
    compare_doe.add_argument(
        "--out",
        default="",
        help="Output directory (default: example/output/<DOE>_compare)",
    )
    compare_doe.add_argument("--early-cycles", type=int, default=30)
    compare_doe.add_argument("--no-diagnosis", action="store_true")
    compare_doe.add_argument("--no-plots", action="store_true")
    compare_doe.add_argument(
        "--cycles",
        default="",
        help="Optional cycle list e.g. 1,5,10,20,50 (smoke/faster). Default: all",
    )

    return p


def _load_table(path: str):
    import pandas as pd

    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    return pd.read_csv(p)


def _write_table(df, path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if p.suffix.lower() == ".parquet":
        df.to_parquet(p, index=False)
    else:
        df.to_csv(p, index=False)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command is None:
        build_parser().print_help()
        return 0

    if args.command == "planning":
        root = Path(__file__).resolve().parent / "planning"
        print(root)
        print("  ROADMAP.md  - phases")
        print("  FEATURES.md - feature catalog")
        print("  NOTES.md    - notes")
        return 0

    if args.command == "extract":
        from cyclediag.features.cv_extract import extract_cv_regions_table
        from cyclediag.features.extract import FeatureConfig, extract_features_table
        from cyclediag.features.lges_catalog import FEATURE_SET_LGES
        from cyclediag.features.lges_extract import LgesExtractConfig, extract_lges_features_table
        from cyclediag.io.cycler_csv import ColumnMap, load_cycler_csv
        from cyclediag.provenance import attach_provenance, build_provenance, save_features_csv

        cmap = (
            ColumnMap.studio_default()
            if args.column_map == "studio"
            else ColumnMap.pne_default()
        )
        df = load_cycler_csv(args.input, column_map=cmap)
        cfg = {"charge_step": args.charge_step, "discharge_step": args.discharge_step}
        if args.cv_only:
            table = extract_cv_regions_table(
                df,
                charge_step=args.charge_step,
                discharge_step=args.discharge_step,
            )
        elif args.feature_set in (FEATURE_SET_LGES, "vp_lges_cycle_v1", "vp_lges_cycle_v2"):
            cfg = LgesExtractConfig(
                charge_step=args.charge_step,
                discharge_step=args.discharge_step,
                cell_id=Path(args.input).stem,
            )
            table = extract_lges_features_table(df, filepath=args.input, config=cfg, raw_df=df)
        else:
            cfg = FeatureConfig(
                charge_step=args.charge_step,
                discharge_step=args.discharge_step,
                feature_set=args.feature_set,
            )
            table = extract_features_table(df, filepath=args.input, config=cfg)
        provenance = build_provenance(args.input, df, config=cfg, column_map=cmap,
                                      selectors={"entry_point": "cli_extract",
                                                 "feature_set": args.feature_set, "cv_only": args.cv_only})
        attach_provenance(table, df, provenance)
        if Path(args.out).suffix.lower() == ".parquet":
            if Path(args.out).resolve() == Path(args.input).resolve():
                raise ValueError("Output must differ from input")
            _write_table(table, args.out)
        else:
            save_features_csv(table, args.out, source_path=args.input)
        print(f"Wrote {len(table)} row(s) → {args.out}")
        return 0

    if args.command == "diagnose":
        from cyclediag.api import diagnose_csv
        from cyclediag.features.lges_extract import LgesExtractConfig
        from cyclediag.io.cycler_csv import ColumnMap

        cmap = (
            ColumnMap.studio_default()
            if args.column_map == "studio"
            else ColumnMap.pne_default()
        )
        cfg = LgesExtractConfig(
            charge_step=args.charge_step,
            discharge_step=args.discharge_step,
            cell_id=Path(args.input).stem,
        )
        result = diagnose_csv(
            args.input,
            column_map=cmap,
            config=cfg,
            with_screen=not args.no_screen,
        )
        out = Path(args.out_dir)
        out.mkdir(parents=True, exist_ok=True)
        from cyclediag.provenance import save_features_csv

        export_features = result["features"].copy()
        export_features.attrs["provenance"] = result["provenance"]
        save_features_csv(export_features, out / "features.csv", source_path=args.input)
        _write_table(result["scored"], str(out / "diagnosis_scores.csv"))
        if not result["indicator_screen"].empty:
            _write_table(result["indicator_screen"], str(out / "indicator_screen.csv"))
        if not result["top_indicators"].empty:
            _write_table(result["top_indicators"], str(out / "top_indicators.csv"))
        if not result["dqdv_screen"].empty:
            _write_table(result["dqdv_screen"], str(out / "dqdv_screen.csv"))
        print(f"Diagnosis written → {out}")
        print(f"  cycles: {len(result['scored'])}")
        return 0

    if args.command == "peaks" and args.peaks_command == "export":
        from cyclediag.features.peak_export import (
            export_peak_feature_table,
            peak_trajectory_config_from_args,
        )

        good = [int(x.strip()) for x in args.good_cycles.split(",") if x.strip()]
        cfg = peak_trajectory_config_from_args(
            relaxed=args.relaxed,
            usable_mad_factor=args.usable_mad_factor,
            max_noise=args.max_noise,
            max_charge_hf=args.max_charge_hf,
            max_discharge_hf=args.max_discharge_hf,
            max_band_gap=args.max_band_gap,
            min_usable_score=args.min_usable_score,
        )
        out_dir = Path(args.out_dir) / args.cell_id
        meta = export_peak_feature_table(
            Path(args.input),
            out_dir,
            cell_id=args.cell_id,
            good_cycles=good,
            config=cfg,
        )
        print(f"Cell: {meta['cell_id']}")
        print(f"Total cycles: {meta['n_cycles_total']}")
        print(f"Usable: {meta['n_cycles_usable']} → {meta['outputs']['usable']}")
        return 0

    if args.command == "peaks" and args.peaks_command == "evolution":
        from cyclediag.features.dqdv_peaks import DqdvPeakConfig
        from cyclediag.features.peak_evolution import PeakEvolutionConfig, track_peaks_pipeline
        from cyclediag.io.cycler_csv import ColumnMap, load_cycler_csv
        from cyclediag.io.stepemd_csv import load_stepemd_csv

        cmap = ColumnMap.studio_default()
        raw = load_cycler_csv(args.input, column_map=cmap)
        step = load_stepemd_csv(args.stepend) if args.stepend else None
        cfg = PeakEvolutionConfig(
            domain=args.domain,
            leg=args.leg,
            normalize=args.normalize,
            n_grid=args.n_grid,
            dqdv_config=DqdvPeakConfig(
                n_interp=args.n_interp,
                sg_window=args.sg_window,
                sg_poly=3,
                merge_v_sep_v=0.003,
                min_distance_frac=0.015,
            ),
            lam_scale=args.lam_scale,
            max_tracks=args.max_tracks,
            use_roi_extract=not args.no_roi,
        )
        result = track_peaks_pipeline(
            raw,
            step,
            config=cfg,
            run_preflight=not args.no_preflight,
            abort_on_preflight_fail=args.abort_on_preflight_fail,
        )
        out = Path(args.out_dir)
        result.to_csv(out)
        plots = result.plot(out)
        print(f"Tracks: {len(result.tracks)}")
        print(f"Events: {len(result.events)}")
        print(f"Validation: {json.dumps(result.validation, default=str)[:2000]}")
        if result.preflight is not None:
            print(result.preflight.to_string(index=False))
        print(f"Wrote -> {out}")
        for p in plots:
            print(f"  plot: {p}")
        return 0

    if args.command == "predict":
        feats = _load_table(args.features)
        if args.model:
            from cyclediag.models.peak_ml import PeakMlBundle, predict_peak_model

            bundle = PeakMlBundle.load(args.model)
            out = predict_peak_model(feats, bundle)
        else:
            from cyclediag.models.predict import predict_features

            ref = _load_table(args.reference) if args.reference else None
            out = predict_features(feats, reference=ref)
        _write_table(out, args.out)
        print(f"Wrote {len(out)} row(s) → {args.out}")
        return 0

    if args.command == "train":
        from cyclediag.models.peak_ml import PeakMlConfig, train_peak_model

        feats = _load_table(args.features)
        good = [int(x.strip()) for x in args.good_cycles.split(",") if x.strip()]
        cfg = PeakMlConfig(
            train_on=args.train_on,
            require_usable=not args.no_require_usable,
        )
        bundle = train_peak_model(feats, good_cycles=good, config=cfg)
        out_dir = bundle.save(args.out)
        print(f"Peak ML model saved → {out_dir}")
        print(f"  train rows: {bundle.train_rows}")
        print(f"  train cycles: {bundle.train_cycles}")
        print(f"  features: {len(bundle.feature_columns)}")
        return 0

    if args.command == "report":
        from cyclediag.analysis.batch_report import run_batch_report

        out_dir = args.output_dir or None
        summary = run_batch_report(
            args.input_dir,
            output_dir=out_dir,
            encoding=args.encoding,
            top_n=args.top_n,
        )
        odir = summary["output_dir"]
        print(f"Report written → {odir}")
        print(f"  HTML: {odir}/diagnosis_report.html")
        print(f"  Cells: {summary['n_cells_ok']}/{summary['n_files']}")
        return 0

    if args.command == "compare-doe":
        from cyclediag.analysis.doe_compare import DoeCompareConfig, run_doe_compare

        cycles = None
        if getattr(args, "cycles", ""):
            cycles = [int(x.strip()) for x in args.cycles.split(",") if x.strip()]
        cfg = DoeCompareConfig(
            doe_id=args.doe,
            fixtures_root=Path(args.fixtures_root) if args.fixtures_root else None,
            out_dir=Path(args.out) if args.out else None,
            early_cycles=args.early_cycles,
            run_diagnosis=not args.no_diagnosis,
            write_plots=not args.no_plots,
            cycles=cycles,
        )
        summary = run_doe_compare(cfg)
        print(f"DOE compare written → {summary['out_dir']}")
        print(f"  arms: {[a['arm_id'] for a in summary['arms']]}")
        print(f"  feature rows: {summary['n_feature_rows']}")
        print(f"  narrative: {summary['out_dir']}/narrative.txt")
        for line in summary.get("narrative", [])[:12]:
            print(f"  · {line}")
        return 0

    print("Unknown command — run with --help")
    return 1


if __name__ == "__main__":
    sys.exit(main())
