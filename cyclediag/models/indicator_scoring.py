"""Indicator scoring track — how much each indicator moved, not *why*.

This module is intentionally separate from ``cyclediag.diagnosis``.

* **Indicator scoring** (here): protocol-aware, one column per indicator family,
  z-scored against a reference window. Outputs per-cycle rollups and per-indicator
  summaries. No LLI / LAM / contact_loss labels.
* **Causal interpretation** (``cyclediag.diagnosis``): mode-weight pattern scores
  that map evidence onto physicochemical degradation modes. Run later, on purpose,
  and never mixed into these tables.

Default row selection is ``routine_only=True`` so RPT / post-RPT / DC-IR spikes
do not inflate the scores. Anchor-only quantities (Q_relax, OCV drift, DC-IR
decomposition) simply have low coverage on the routine subset; they are scored
when present rather than forced in.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Literal

import numpy as np
import pandas as pd

from cyclediag.features._statistics import finite_pearson
from cyclediag.features.indicator_registry import (
    ROLE_INDICATOR,
    family_of,
    primary_indicator_columns,
    role_of,
)

FLAG_WATCH = 0.55
FLAG_ALERT = 0.75

SCORE_META_COLS = (
    "indicator_score",
    "indicator_flag",
    "indicator_top",
    "indicator_n_scored",
    "scoring_row",
    "score_layer",
)


@dataclass
class IndicatorScoreResult:
    """Outputs of the indicator scoring track."""

    cycle_scores: pd.DataFrame
    """Wide feature rows that were scored, plus rollup columns."""

    cycle_contributions: pd.DataFrame
    """Long table: one row per (cycle, family) with |z| contribution."""

    indicator_summary: pd.DataFrame
    """Per-indicator (family) summary for the scored subset — no mode labels."""

    meta: dict[str, Any] = field(default_factory=dict)

    @property
    def empty(self) -> bool:
        return self.cycle_scores is None or self.cycle_scores.empty


def attach_protocol_flags(
    features: pd.DataFrame,
    raw_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """Attach ``protocol_kind`` / ``protocol_excluded`` from the raw cycler table.

    No-op when flags are already present or ``raw_df`` is missing. Detection uses
    step-end grouping — the same rules as peak export / presentation tools.
    """
    if features is None or features.empty:
        return features
    out = features.copy()
    if "protocol_kind" in out.columns and "protocol_excluded" in out.columns:
        return out
    if (
        raw_df is None
        or raw_df.empty
        or "cycle" not in raw_df.columns
        or "cycle" not in out.columns
    ):
        if "protocol_kind" not in out.columns:
            out["protocol_kind"] = "unknown"
        if "protocol_excluded" not in out.columns:
            out["protocol_excluded"] = False
        return out

    from cyclediag.io.cycle_protocol import build_protocol_exclusion, detect_protocol_flags

    step_keys = [c for c in ("cycle", "StepNo", "step_no", "step") if c in raw_df.columns]
    if "cycle" in step_keys and len(step_keys) >= 2:
        se = raw_df.groupby(step_keys[:2], as_index=False).tail(1)
    else:
        se = raw_df
    try:
        prot = build_protocol_exclusion(se)
        flags = prot.flags if prot.flags is not None and not prot.flags.empty else detect_protocol_flags(se)
    except Exception:
        out["protocol_kind"] = "unknown"
        out["protocol_excluded"] = False
        return out

    if flags is None or flags.empty or "cycle" not in flags.columns:
        out["protocol_kind"] = "unknown"
        out["protocol_excluded"] = False
        return out

    kind_map = dict(zip(flags["cycle"].astype(int), flags["protocol_kind"].astype(str)))
    excl = set(int(c) for c in (prot.excluded or set()))
    cyc = pd.to_numeric(out["cycle"], errors="coerce")
    out["protocol_kind"] = cyc.map(
        lambda c: kind_map.get(int(c), "unknown") if pd.notna(c) else "unknown"
    )
    out["protocol_excluded"] = cyc.map(
        lambda c: int(c) in excl if pd.notna(c) else False
    )
    return out


def filter_scoring_rows(
    features: pd.DataFrame,
    *,
    routine_only: bool = True,
    include_kinds: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Select rows that enter the indicator score.

    ``routine_only`` (default) drops RPT, post-RPT, capa-check and DC-IR-like
    cycles. Pass ``routine_only=False`` only for debug / full-table views.
    """
    if features is None or features.empty:
        return pd.DataFrame() if features is None else features.iloc[0:0].copy()
    df = features
    if not routine_only and not include_kinds:
        return df.copy()
    if "protocol_excluded" in df.columns:
        excl = df["protocol_excluded"].fillna(False).astype(bool)
        df = df.loc[~excl]
    if include_kinds and "protocol_kind" in df.columns:
        return df.loc[df["protocol_kind"].astype(str).isin(set(include_kinds))].copy()
    if routine_only and "protocol_kind" in df.columns:
        kind = df["protocol_kind"].astype(str)
        # When flags were never detected every row is "unknown" — keep those
        # rather than scoring an empty table. Once real kinds exist, routine only.
        if (kind == "routine").any():
            df = df.loc[kind == "routine"]
        elif not kind.isin(["unknown", "nan"]).all():
            df = df.loc[kind == "routine"]
    return df.copy()


def _reference_stats(
    ref: pd.DataFrame,
    cols: list[str],
    *,
    early_frac: float = 0.2,
) -> tuple[pd.Series, pd.Series]:
    """Median / std from the early portion of the reference rows."""
    if ref.empty:
        med = pd.Series({c: np.nan for c in cols})
        std = med.copy()
        return med, std
    if "cycle" in ref.columns and len(ref) >= 5:
        cyc = pd.to_numeric(ref["cycle"], errors="coerce")
        thr = float(cyc.quantile(early_frac))
        window = ref.loc[cyc <= thr]
        if len(window) < 3:
            window = ref
    else:
        window = ref
    values = window.reindex(columns=cols).apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    med = values.median()
    # Translation preserves variance and makes identical non-binary decimals
    # exactly zero, rather than treating mean-roundoff as calibrated variance.
    std = (values - med).std()
    return med, std


def _health_col(df: pd.DataFrame) -> str | None:
    for c in ("SoHQ", "dchgCapa", "capacity", "f_Q_max"):
        if c in df.columns and pd.to_numeric(df[c], errors="coerce").notna().sum() >= 3:
            return c
    return None


def score_indicators(
    features: pd.DataFrame,
    *,
    reference: pd.DataFrame | None = None,
    raw_df: pd.DataFrame | None = None,
    routine_only: bool = True,
    early_frac: float = 0.2,
    grain: Literal["cycle", "summary", "both"] = "both",
    scale_floors: dict[str, float] | None = None,
) -> IndicatorScoreResult:
    """Score indicators without assigning a degradation mode.

    Parameters
    ----------
    features
        Cycle feature table (output of ``extract_lges_features_table``).
    reference
        Optional reference population for median/std. Defaults to the early
        routine window of ``features`` itself.
    raw_df
        Optional raw cycler table used to attach protocol flags when missing.
    routine_only
        If True (default), score only non-excluded routine cycles.
    early_frac
        Fraction of earliest cycles used as the within-cell baseline.
    grain
        ``cycle`` / ``summary`` / ``both``.
    """
    floors = scale_floors or {}
    if any(not np.isfinite(v) or v <= 0 for v in floors.values()):
        raise ValueError("scale_floors must be positive finite calibrated feature scales")
    empty_meta = {
        "score_kind": "heuristic_deviation_not_probability",
        "zero_variance_policy": "unknown_without_calibrated_floor",
        "scale_floors": dict(floors),
        "score_layer": "indicator",
        "routine_only": routine_only,
        "n_input_rows": 0 if features is None else int(len(features)),
        "n_scored_rows": 0,
        "n_families": 0,
        "causal_track": "separate — see cyclediag.diagnosis",
    }
    if features is None or features.empty:
        return IndicatorScoreResult(
            cycle_scores=pd.DataFrame(),
            cycle_contributions=pd.DataFrame(),
            indicator_summary=pd.DataFrame(),
            meta=empty_meta,
        )

    tagged = attach_protocol_flags(features, raw_df)
    scored_rows = filter_scoring_rows(tagged, routine_only=routine_only)
    cols = primary_indicator_columns(scored_rows if not scored_rows.empty else tagged)
    empty_meta["n_families"] = len(cols)
    empty_meta["n_scored_rows"] = int(len(scored_rows))

    if scored_rows.empty or not cols:
        out = tagged.copy()
        out["scoring_row"] = False
        out["indicator_score"] = np.nan
        out["indicator_flag"] = "unknown"
        out["indicator_top"] = ""
        out["indicator_n_scored"] = 0
        out["score_layer"] = "indicator"
        return IndicatorScoreResult(
            cycle_scores=out if grain in ("cycle", "both") else pd.DataFrame(),
            cycle_contributions=pd.DataFrame(),
            indicator_summary=pd.DataFrame(),
            meta=empty_meta,
        )

    ref_src = reference if reference is not None else scored_rows
    if reference is None:
        ref_src = filter_scoring_rows(
            attach_protocol_flags(ref_src, raw_df), routine_only=routine_only,
        )
    # Implicit baselines are per cell/source, never pooled across cells.
    values = scored_rows[cols].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    z_mat = pd.DataFrame(np.nan, index=scored_rows.index, columns=cols)
    constant_changed = pd.Series(0, index=scored_rows.index)
    constant_unchanged = pd.Series(0, index=scored_rows.index)
    group_cols = [c for c in ("cell_id", "file", "source_file") if c in scored_rows.columns]
    groups = scored_rows.groupby(group_cols, sort=False, dropna=False) if group_cols and reference is None else [(None, scored_rows)]
    for _, group in groups:
        ref_group = group if reference is None else ref_src
        med, std = _reference_stats(ref_group, cols, early_frac=early_frac)
        delta = values.loc[group.index] - med
        zero = std.eq(0)
        constant_changed.loc[group.index] = (delta.ne(0) & delta.notna() & zero).sum(axis=1)
        constant_unchanged.loc[group.index] = (delta.eq(0) & zero).sum(axis=1)
        for col, floor in floors.items():
            if col in std and pd.notna(std[col]):
                std[col] = max(std[col], floor)
        z_mat.loc[group.index] = (delta / std.replace(0, np.nan)).abs()
    z_mat = z_mat.replace([np.inf, -np.inf], np.nan)
    row_score = (z_mat.mean(axis=1, skipna=True) / 3.0).clip(0.0, 1.0)
    top_labels: list[str] = []
    n_scored: list[int] = []
    for idx in z_mat.index:
        row = z_mat.loc[idx].dropna()
        n_scored.append(int(row.shape[0]))
        if row.empty:
            top_labels.append("")
            continue
        top = row.sort_values(ascending=False).head(3)
        top_labels.append(", ".join(f"{k}={v:.2g}" for k, v in top.items()))

    cycle_out = tagged.copy()
    cycle_out["scoring_row"] = False
    cycle_out["indicator_score"] = np.nan
    cycle_out["indicator_flag"] = "unknown"
    cycle_out["indicator_constant_changed"] = constant_changed.reindex(cycle_out.index).fillna(0).astype(int)
    cycle_out["indicator_constant_unchanged"] = constant_unchanged.reindex(cycle_out.index).fillna(0).astype(int)
    cycle_out["indicator_top"] = ""
    cycle_out["indicator_n_scored"] = 0
    cycle_out["score_layer"] = "indicator"
    cycle_out.loc[scored_rows.index, "scoring_row"] = True
    cycle_out.loc[scored_rows.index, "indicator_score"] = row_score.to_numpy()
    cycle_out.loc[scored_rows.index, "indicator_flag"] = [
        "unknown" if not np.isfinite(s) else "alert" if s >= FLAG_ALERT else ("watch" if s >= FLAG_WATCH else "ok")
        for s in row_score.to_numpy()
    ]
    cycle_out.loc[scored_rows.index, "indicator_top"] = top_labels
    cycle_out.loc[scored_rows.index, "indicator_n_scored"] = n_scored
    cycle_out["indicator_evidence_coverage"] = cycle_out["indicator_n_scored"] / len(cols)
    incomplete_ok = (cycle_out["indicator_flag"] == "ok") & (cycle_out["indicator_evidence_coverage"] < 1)
    cycle_out.loc[incomplete_ok, "indicator_flag"] = "unknown"

    # --- long contributions ---
    contrib_rows: list[dict[str, Any]] = []
    if grain in ("cycle", "both"):
        cyc = (
            pd.to_numeric(scored_rows["cycle"], errors="coerce")
            if "cycle" in scored_rows.columns
            else pd.Series(np.nan, index=scored_rows.index)
        )
        cell = scored_rows["cell_id"] if "cell_id" in scored_rows.columns else None
        for col in cols:
            series = z_mat[col]
            for idx, z in series.items():
                if not np.isfinite(z):
                    continue
                cyc_val = cyc.loc[idx]
                contrib_rows.append({
                    "cell_id": None if cell is None else cell.loc[idx],
                    "cycle": float(cyc_val) if np.isfinite(cyc_val) else np.nan,
                    "feature": col,
                    "family": family_of(col),
                    "role": role_of(col),
                    "abs_z": float(z),
                    "contribution": float(min(z / 3.0, 1.0)),
                    "score_layer": "indicator",
                })
    contributions = pd.DataFrame(contrib_rows)

    # --- per-indicator summary (no mechanism labels) ---
    summary = pd.DataFrame()
    if grain in ("summary", "both"):
        summary = _summarize_indicators(scored_rows, cols, z_mat)

    meta = {
        **empty_meta,
        "n_scored_rows": int(len(scored_rows)),
        "n_families": len(cols),
        "early_frac": early_frac,
        "reference": "provided" if reference is not None else "early_routine_window",
        "protocol_kinds_scored": (
            sorted(scored_rows["protocol_kind"].astype(str).unique())
            if "protocol_kind" in scored_rows.columns else []
        ),
    }
    return IndicatorScoreResult(
        cycle_scores=cycle_out if grain in ("cycle", "both") else pd.DataFrame(),
        cycle_contributions=contributions,
        indicator_summary=summary,
        meta=meta,
    )


def _summarize_indicators(
    df: pd.DataFrame,
    cols: list[str],
    z_mat: pd.DataFrame,
) -> pd.DataFrame:
    hcol = _health_col(df)
    cyc = pd.to_numeric(df["cycle"], errors="coerce") if "cycle" in df.columns else None
    rows: list[dict[str, Any]] = []
    if cyc is not None and cyc.notna().any():
        late_thr = float(cyc.quantile(0.8))
        early_thr = float(cyc.quantile(0.2))
        late_m = cyc >= late_thr
        early_m = cyc <= early_thr
    else:
        late_m = pd.Series(False, index=df.index)
        early_m = late_m

    for col in cols:
        s = pd.to_numeric(df[col], errors="coerce")
        z = z_mat[col]
        cov = float(s.notna().mean())
        if cov < 0.05:
            continue
        corr_cycle = finite_pearson(s, cyc) if cyc is not None else np.nan
        corr_health = finite_pearson(s, df[hcol]) if hcol else np.nan
        early = float(s[early_m].mean()) if early_m.any() else np.nan
        late = float(s[late_m].mean()) if late_m.any() else np.nan
        delta = late - early if np.isfinite(early) and np.isfinite(late) else np.nan
        med_abs_z = float(z.median(skipna=True)) if z.notna().any() else np.nan
        late_z = z[late_m].dropna()
        late_abs_z = float(late_z.median()) if not late_z.empty else np.nan

        # Indicator score in [0, 1]: how strongly this column moved.
        # Deliberately descriptive — not a mode probability.
        parts = []
        if np.isfinite(corr_health):
            parts.append(abs(corr_health))
        elif np.isfinite(corr_cycle):
            parts.append(abs(corr_cycle) * 0.85)
        if np.isfinite(late_abs_z):
            parts.append(float(np.clip(late_abs_z / 3.0, 0.0, 1.0)))
        elif np.isfinite(med_abs_z):
            parts.append(float(np.clip(med_abs_z / 3.0, 0.0, 1.0)))
        score = float(np.clip(max(parts), 0.0, 1.0)) if parts and z.notna().any() else np.nan

        rows.append({
            "feature": col,
            "family": family_of(col),
            "role": ROLE_INDICATOR,
            "coverage": cov,
            "scored_evidence_coverage": float(z.notna().mean()),
            "status": "scored" if z.notna().any() else "unknown",
            "corr_cycle": corr_cycle,
            "corr_health": corr_health,
            "early_mean": early,
            "late_mean": late,
            "delta_late_early": delta,
            "median_abs_z": med_abs_z,
            "late_median_abs_z": late_abs_z,
            "indicator_score": score,
            "score_layer": "indicator",
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    from cyclediag.analysis.indicator_layers import annotate_layers

    out = annotate_layers(out)
    return out.sort_values("indicator_score", ascending=False).reset_index(drop=True)


def top_scored_indicators(
    summary: pd.DataFrame,
    *,
    n: int = 15,
    min_score: float = 0.35,
) -> pd.DataFrame:
    """Highest-scoring indicators, already one per family."""
    if summary is None or summary.empty:
        return pd.DataFrame()
    filt = summary[summary["indicator_score"] >= min_score]
    if filt.empty:
        filt = summary
    if "family" in filt.columns:
        filt = filt.drop_duplicates(subset="family", keep="first")
    return filt.head(n)
