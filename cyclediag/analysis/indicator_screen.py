"""Per-file indicator screening and cross-cell comparison."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cyclediag.features._statistics import finite_pearson
from cyclediag.features.indicator_registry import (
    ROLE_COVARIATE,
    ROLE_INDICATOR,
    ROLE_TARGET,
    annotate_indicators,
    indicator_columns,
    primary_indicator_columns,
)

# The descriptive screen reports measurable quantities, targets and protocol
# covariates included, and tags each row with its registry role and family.
# Anything that ranks or aggregates uses one representative per family.
_SCREEN_ROLES = (ROLE_INDICATOR, ROLE_TARGET, ROLE_COVARIATE)
_HEALTH_COLS = ("SoHQ", "dchgCapa", "capacity", "f_Q_max")


def _health_column(df: pd.DataFrame) -> str | None:
    for c in _HEALTH_COLS:
        if c in df.columns and pd.to_numeric(df[c], errors="coerce").notna().sum() >= 3:
            return c
    return None


def _late_early_split(cycles: pd.Series, late_frac: float = 0.2) -> tuple[pd.Series, pd.Series]:
    cyc = pd.to_numeric(cycles, errors="coerce").dropna()
    if cyc.empty:
        mask = pd.Series(False, index=cycles.index)
        return mask, mask
    hi = float(cyc.quantile(1.0 - late_frac))
    lo = float(cyc.quantile(late_frac))
    late = pd.to_numeric(cycles, errors="coerce") >= hi
    early = pd.to_numeric(cycles, errors="coerce") <= lo
    return early, late


def screen_indicators(
    features: pd.DataFrame,
    *,
    reference_cycle: int = 1,
    health_col: str | None = None,
) -> pd.DataFrame:
    """Rank indicators for one cell: drift, health correlation, late vs early change."""
    if features is None or features.empty:
        return pd.DataFrame()

    df = features.sort_values("cycle").copy()
    hcol = health_col or _health_column(df)
    cols = indicator_columns(df, roles=_SCREEN_ROLES)
    if not cols:
        return pd.DataFrame()

    base = df[df["cycle"] == reference_cycle]
    if base.empty:
        base = df.head(1)
    base_row = base.iloc[0]

    early_m, late_m = _late_early_split(df["cycle"])
    rows: list[dict] = []

    for col in cols:
        s = pd.to_numeric(df[col], errors="coerce")
        cov = float(s.notna().mean() * 100.0)
        if cov < 5:
            continue
        std = float(s.std(skipna=True))
        if not np.isfinite(std) or std < 1e-12:
            continue

        corr_cycle = finite_pearson(s, df["cycle"])
        corr_health = np.nan
        if hcol:
            corr_health = finite_pearson(s, df[hcol])

        b = base_row.get(col)
        delta_vs_ref = None
        if b is not None and np.isfinite(b):
            delta_vs_ref = float(s.iloc[-1] - b) if np.isfinite(s.iloc[-1]) else None

        early_mean = float(s[early_m].mean()) if early_m.any() else np.nan
        late_mean = float(s[late_m].mean()) if late_m.any() else np.nan
        delta_late_early = (
            late_mean - early_mean
            if np.isfinite(early_mean) and np.isfinite(late_mean)
            else np.nan
        )

        # Severity: prefer |health corr|, else |cycle corr|, else |late-early| normalized
        sev = 0.0
        if corr_health is not None and np.isfinite(corr_health):
            sev = abs(float(corr_health))
        elif corr_cycle is not None and np.isfinite(corr_cycle):
            sev = abs(float(corr_cycle)) * 0.85
        if np.isfinite(delta_late_early):
            sev = max(sev, min(1.0, abs(delta_late_early) / (std + 1e-9)) * 0.7)

        reason = []
        if corr_health is not None and abs(corr_health) >= 0.7:
            reason.append(f"|r| vs {hcol}={corr_health:.2f}")
        if corr_cycle is not None and abs(corr_cycle) >= 0.5:
            reason.append(f"|r| vs cycle={corr_cycle:.2f}")
        if np.isfinite(delta_late_early) and abs(delta_late_early) > 0.5 * std:
            reason.append("late≠early")

        rows.append({
            "feature": col,
            "coverage_pct": round(cov, 1),
            "corr_cycle": round(float(corr_cycle), 3),
            "corr_health": round(float(corr_health), 3),
            "delta_vs_ref": delta_vs_ref,
            "delta_late_early": delta_late_early,
            "severity": round(sev, 3),
            "signal": "; ".join(reason) if reason else "weak",
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = annotate_indicators(out)
    return out.sort_values("severity", ascending=False).reset_index(drop=True)


def screen_indicators_by_file(
    features: pd.DataFrame,
    *,
    file_col: str = "file",
    cell_col: str = "cell_id",
) -> pd.DataFrame:
    """Run screen_indicators per file/cell; add grouping columns."""
    if features is None or features.empty:
        return pd.DataFrame()
    group_cols = [c for c in (cell_col, file_col) if c in features.columns]
    if not group_cols:
        group_cols = [file_col] if file_col in features.columns else []

    parts = []
    if group_cols:
        for _, grp in features.groupby(group_cols, sort=False):
            screened = screen_indicators(grp)
            if screened.empty:
                continue
            for c in group_cols:
                screened[c] = grp[c].iloc[0]
            parts.append(screened)
    else:
        parts.append(screen_indicators(features))

    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def compare_cells(
    features: pd.DataFrame,
    *,
    cycle_col: str = "cycle",
    cell_col: str = "cell_id",
    late_frac: float = 0.2,
) -> pd.DataFrame:
    """Find indicators that diverge most across cells (late-cycle window)."""
    if features is None or features.empty:
        return pd.DataFrame()
    if cell_col not in features.columns:
        return pd.DataFrame()

    cells = features[cell_col].dropna().unique()
    if len(cells) < 2:
        return pd.DataFrame()

    cols = primary_indicator_columns(features)
    cyc = pd.to_numeric(features[cycle_col], errors="coerce")
    late_thr = float(cyc.quantile(1.0 - late_frac))
    late = features[cyc >= late_thr].copy()
    if late.empty:
        late = features.copy()

    rows: list[dict] = []
    for col in cols:
        pivot = late.pivot_table(index=cycle_col, columns=cell_col, values=col, aggfunc="mean")
        if pivot.shape[1] < 2:
            continue
        # per-cycle spread across cells
        spread = pivot.std(axis=1, skipna=True)
        mean_spread = float(spread.mean()) if spread.notna().any() else np.nan
        if not np.isfinite(mean_spread):
            continue
        # trend of spread (growing divergence?)
        spread_vals = spread.dropna()
        spread_corr = np.nan
        if len(spread_vals) >= 3:
            spread_corr = finite_pearson(spread_vals,
                pd.Series(spread_vals.index.astype(float), index=spread_vals.index)
            )
        cell_means = pivot.mean(skipna=True)
        pair_diff = float(cell_means.max() - cell_means.min())
        rows.append({
            "feature": col,
            "mean_spread_late": round(mean_spread, 6),
            "max_cell_gap_late": round(pair_diff, 6),
            "spread_trend_r": round(float(spread_corr), 3),
            "n_cells": int(pivot.shape[1]),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["divergence_score"] = (
        out["mean_spread_late"].rank(pct=True) * 0.5
        + out["max_cell_gap_late"].rank(pct=True) * 0.5
    )
    return out.sort_values("divergence_score", ascending=False).reset_index(drop=True)


def top_problem_indicators(
    screened: pd.DataFrame,
    *,
    n: int = 15,
    min_severity: float = 0.35,
    one_per_family: bool = True,
) -> pd.DataFrame:
    """Most severe indicators, by default at most one per indicator family.

    Without the family filter a single physical signal fills the list with its
    own aliases — end-of-charge rest voltage alone had six columns.
    """
    if screened is None or screened.empty:
        return pd.DataFrame()
    filt = screened[screened["severity"] >= min_severity]
    if filt.empty:
        filt = screened
    if one_per_family and "family" in filt.columns:
        # already sorted by severity, so the first row of each family is its best
        filt = filt.drop_duplicates(subset="family", keep="first")
    return filt.head(n)
