"""Peak-feature ML — one-class anomaly model from good/reference cycles."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

META_COLS = frozenset({
    "cell_id", "source_file", "cycle", "exclude_flags",
    "quality_median", "quality_threshold", "good_cycle_ref",
    "usable", "usable_auto", "usable_charge", "usable_discharge",
})

DEFAULT_PEAK_FEATURE_PREFIXES = ("cha_", "dis_", "d_cha_", "d_dis_")
DEFAULT_PEAK_FEATURE_SUFFIXES = ("_V", "_H")


@dataclass
class PeakMlConfig:
    contamination: float = 0.05
    n_estimators: int = 200
    random_state: int = 42
    train_on: str = "usable"  # usable | good_cycles | all_complete
    require_usable: bool = True
    alert_quantile: float = 0.90
    watch_quantile: float = 0.75


@dataclass
class PeakMlBundle:
    pipeline: Pipeline
    feature_columns: list[str]
    config: PeakMlConfig
    train_cycles: list[int] = field(default_factory=list)
    train_rows: int = 0
    score_thresholds: dict[str, float] = field(default_factory=dict)
    version: str = "peak_ml_v1"

    def save(self, out_dir: str | Path) -> Path:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        joblib_path = out_dir / "model.joblib"
        meta_path = out_dir / "model_meta.json"
        joblib.dump(
            {
                "pipeline": self.pipeline,
                "feature_columns": self.feature_columns,
                "config": asdict(self.config),
                "train_cycles": self.train_cycles,
                "train_rows": self.train_rows,
                "score_thresholds": self.score_thresholds,
                "version": self.version,
            },
            joblib_path,
        )
        meta = {
            "version": self.version,
            "train_rows": self.train_rows,
            "train_cycles": self.train_cycles,
            "n_features": len(self.feature_columns),
            "feature_columns": self.feature_columns,
            "config": asdict(self.config),
            "score_thresholds": self.score_thresholds,
            "model_path": str(joblib_path),
        }
        meta_path.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
        return out_dir

    @classmethod
    def load(cls, model_dir: str | Path) -> PeakMlBundle:
        model_dir = Path(model_dir)
        raw: dict[str, Any] = joblib.load(model_dir / "model.joblib")
        cfg_dict = dict(raw.get("config", {}))
        config = PeakMlConfig(**{k: v for k, v in cfg_dict.items() if k in PeakMlConfig.__dataclass_fields__})
        return cls(
            pipeline=raw["pipeline"],
            feature_columns=list(raw["feature_columns"]),
            config=config,
            train_cycles=list(raw.get("train_cycles", [])),
            train_rows=int(raw.get("train_rows", 0)),
            score_thresholds=dict(raw.get("score_thresholds", {})),
            version=str(raw.get("version", "peak_ml_v1")),
        )


def peak_feature_columns(df: pd.DataFrame) -> list[str]:
    cols: list[str] = []
    for c in df.columns:
        if c in META_COLS:
            continue
        if not pd.api.types.is_numeric_dtype(df[c]):
            continue
        if c.startswith(("noise_", "hf_", "quality_", "band_", "n_pts", "usable_score")):
            continue
        if any(c.startswith(p) for p in DEFAULT_PEAK_FEATURE_PREFIXES) and any(
            c.endswith(s) for s in DEFAULT_PEAK_FEATURE_SUFFIXES
        ):
            cols.append(c)
            continue
        if c.startswith("d_cha_") or c.startswith("d_dis_"):
            cols.append(c)
    return list(dict.fromkeys(cols))


def _select_train_rows(
    df: pd.DataFrame,
    *,
    good_cycles: list[int] | None,
    config: PeakMlConfig,
) -> pd.DataFrame:
    work = df.copy()
    if "cycle" not in work.columns:
        raise ValueError("feature table must include 'cycle' column")

    if config.train_on == "good_cycles":
        if not good_cycles:
            raise ValueError("good_cycles required when train_on='good_cycles'")
        work = work[work["cycle"].isin(good_cycles)]
    elif config.train_on == "usable":
        if "usable" not in work.columns:
            raise ValueError("usable column required when train_on='usable'")
        work = work[work["usable"]]
    elif config.train_on == "all_complete":
        if "band_gap_total" in work.columns:
            work = work[work["band_gap_total"] == 0]
    else:
        raise ValueError(f"Unknown train_on: {config.train_on}")

    if config.require_usable and "usable" in work.columns:
        work = work[work["usable"]]

    feat_cols = peak_feature_columns(work)
    if not feat_cols:
        raise ValueError("no peak feature columns found for training")
    complete = work[feat_cols].notna().sum(axis=1) >= max(4, int(0.6 * len(feat_cols)))
    work = work[complete]
    if work.empty:
        raise ValueError("no training rows after filtering")
    return work


def train_peak_model(
    features: pd.DataFrame,
    *,
    good_cycles: list[int] | None = None,
    config: PeakMlConfig | None = None,
) -> PeakMlBundle:
    config = config or PeakMlConfig()
    train_df = _select_train_rows(features, good_cycles=good_cycles, config=config)
    feat_cols = peak_feature_columns(train_df)
    x = train_df[feat_cols].to_numpy(dtype=float)

    pipeline = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        (
            "model",
            IsolationForest(
                n_estimators=config.n_estimators,
                contamination=config.contamination,
                random_state=config.random_state,
            ),
        ),
    ])
    pipeline.fit(x)

    scores = _raw_scores(pipeline, x)
    thresholds = {
        "watch": float(np.quantile(scores, min(0.99, config.watch_quantile))),
        "alert": float(np.quantile(scores, min(0.995, config.alert_quantile))),
    }

    return PeakMlBundle(
        pipeline=pipeline,
        feature_columns=feat_cols,
        config=config,
        train_cycles=sorted(int(c) for c in train_df["cycle"].unique()),
        train_rows=len(train_df),
        score_thresholds=thresholds,
    )


def _raw_scores(pipeline: Pipeline, x: np.ndarray) -> np.ndarray:
    model: IsolationForest = pipeline.named_steps["model"]
    # Higher = more anomalous
    return -model.decision_function(pipeline[:-1].transform(x))


def _normalize_scores(raw: np.ndarray) -> np.ndarray:
    if len(raw) == 0:
        return raw
    lo, hi = float(np.min(raw)), float(np.max(raw))
    if hi - lo < 1e-12:
        return np.zeros_like(raw)
    return (raw - lo) / (hi - lo)


def predict_peak_model(
    features: pd.DataFrame,
    bundle: PeakMlBundle,
) -> pd.DataFrame:
    out = features.copy()
    feat_cols = bundle.feature_columns
    missing = [c for c in feat_cols if c not in out.columns]
    for c in missing:
        out[c] = np.nan

    x = out[feat_cols].to_numpy(dtype=float, copy=True)
    x[~np.isfinite(x)] = np.nan
    valid = np.isfinite(x).any(axis=1)
    raw = np.full(len(out), np.nan, dtype=float)
    if valid.any():
        raw[valid] = _raw_scores(bundle.pipeline, x[valid])

    out["ml_raw_score"] = raw
    normalized = np.full(len(out), np.nan)
    if valid.any():
        normalized[valid] = _normalize_scores(raw[valid])
    out["ml_anomaly_score"] = normalized
    out["ml_evidence_coverage"] = np.isfinite(x).mean(axis=1)

    watch = bundle.score_thresholds.get("watch", np.nan)
    alert = bundle.score_thresholds.get("alert", np.nan)
    flags: list[str] = []
    for s in raw:
        if not np.isfinite(s):
            flags.append("unknown")
        elif np.isfinite(alert) and s >= alert:
            flags.append("alert")
        elif np.isfinite(watch) and s >= watch:
            flags.append("watch")
        else:
            flags.append("ok")
    out["ml_flag"] = flags
    out["ml_is_outlier"] = out["ml_flag"].isin(["watch", "alert"])
    # v1 bundles store cycle numbers, not source/cell identity. A number match is
    # not reference membership, and must never suppress an alert or unknown.
    out["ml_reference_membership"] = "unknown_legacy_cycle_only"
    out["ml_train_cycle_number_match"] = out["cycle"].isin(bundle.train_cycles) if "cycle" in out else False
    out["ml_score_kind"] = "heuristic_anomaly_not_probability"
    return out


def top_deviating_features(
    row: pd.Series,
    bundle: PeakMlBundle,
    *,
    train_median: pd.Series,
    top_n: int = 3,
) -> str:
    """Human-readable top |z|-like deviations vs training median."""
    parts: list[tuple[float, str]] = []
    for col in bundle.feature_columns:
        val = row.get(col)
        ref = train_median.get(col)
        if pd.isna(val) or pd.isna(ref):
            continue
        parts.append((abs(float(val) - float(ref)), f"{col}={float(val):.4g}"))
    parts.sort(reverse=True)
    return ", ".join(p for _, p in parts[:top_n])
