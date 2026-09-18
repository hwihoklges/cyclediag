"""PNG overview plots for cycle-level indicators."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ._statistics import finite_pearson

PANEL_SPECS: tuple[tuple[str, tuple[str, ...], str], ...] = (
    ("Capacity / SoH", ("SoHQ", "CE"), "%"),
    ("EoC rest V (after charge)", ("EoC_restV_init", "EoC_restV_60s", "EoC_restV_end"), "V"),
    ("EoD rest V (after discharge)", ("EoD_restV_init", "EoD_restV_60s", "EoD_restV_end"), "V"),
    ("EoC dchg R (after charge rest)", ("EoC_dchgR_10s", "EoC_dchgR_30s", "EoC_dchgR_60s"), "mΩ"),
    ("EoD chg R (after discharge rest)", ("EoD_chgR_10s", "EoD_chgR_30s", "EoD_chgR_60s"), "mΩ"),
    ("CC ratio / CV time", ("chgCapa_CCratio", "chgCVtime"), "% / s"),
)

_METRIC_COLORS = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
)
_CELL_LINESTYLES = ("-", "--", "-.", ":")
_CELL_MARKERS = ("o", "s", "^", "D")


def _safe_stem(name: str) -> str:
    return (
        str(name)
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
        [:48]
        or "cell"
    )


def _metric_color(col_idx: int) -> str:
    return _METRIC_COLORS[col_idx % len(_METRIC_COLORS)]


def _x_axis_column(work: pd.DataFrame) -> tuple[str, str]:
    if "tagged_cycle" in work.columns and work["tagged_cycle"].notna().any():
        return "tagged_cycle", "Tagged cycle #"
    return "cycle", "Cycle"


def plot_cycle_indicator_overview(
    features: pd.DataFrame,
    out_path: Path | str,
    *,
    title: str | None = None,
    dpi: int = 140,
    max_cells: int = 8,
) -> Path | None:
    """Save a 2×3 overview PNG (rest V, R, SoHQ…) vs cycle.

    Each metric line uses a distinct color within a panel; multiple cells
    differ by linestyle/marker.
    """
    if features is None or features.empty or "cycle" not in features.columns:
        return None

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    work = features.copy()
    work["cycle"] = pd.to_numeric(work["cycle"], errors="coerce")
    work = work.dropna(subset=["cycle"]).sort_values("cycle")
    if work.empty:
        return None

    if "cell_id" not in work.columns:
        work["cell_id"] = "cell"

    cells = list(work["cell_id"].astype(str).unique())[:max_cells]
    multi_cell = len(cells) > 1
    x_col, x_label = _x_axis_column(work)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), facecolor="white")
    axes = axes.ravel()

    for ax, (panel_title, cols, ylabel) in zip(axes, PANEL_SPECS):
        plotted = False
        for ci, cell in enumerate(cells):
            sub = work[work["cell_id"].astype(str) == cell].sort_values(x_col)
            x = pd.to_numeric(sub[x_col], errors="coerce").to_numpy(dtype=float)
            ls = _CELL_LINESTYLES[ci % len(_CELL_LINESTYLES)]
            mk = _CELL_MARKERS[ci % len(_CELL_MARKERS)]
            for mi, col in enumerate(cols):
                if col not in sub.columns:
                    continue
                y = pd.to_numeric(sub[col], errors="coerce").to_numpy(dtype=float)
                if not np.isfinite(y).any():
                    continue
                label = f"{cell} · {col}" if multi_cell else col
                ax.plot(
                    x, y,
                    color=_metric_color(mi),
                    linestyle=ls,
                    linewidth=1.5,
                    marker=mk,
                    markersize=3 if not multi_cell else 2.5,
                    markevery=max(1, len(x) // 25),
                    label=label,
                )
                plotted = True
        ax.set_title(panel_title, fontsize=10)
        ax.set_xlabel(x_label, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        if plotted:
            ax.legend(fontsize=6, loc="best", framealpha=0.9, ncol=1)
        else:
            ax.text(
                0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="gray",
            )

    n_cells = work["cell_id"].nunique()
    if x_col == "tagged_cycle":
        cyc0, cyc1 = int(work["tagged_cycle"].min()), int(work["tagged_cycle"].max())
        n_cyc = int(work["tagged_cycle"].nunique())
        src = ""
        if "tagged_source" in work.columns and work["tagged_source"].notna().any():
            src = f" ({work['tagged_source'].iloc[0]})"
        title_default = (
            f"Cycle indicators  ·  {n_cyc} tagged cycles ({cyc0}–{cyc1}){src}  ·  "
            f"{n_cells} cell(s)"
        )
    else:
        cyc0, cyc1 = int(work["cycle"].min()), int(work["cycle"].max())
        n_cyc = work["cycle"].nunique()
        title_default = (
            f"Cycle indicators  ·  {n_cyc} cycles ({cyc0}–{cyc1})  ·  {n_cells} cell(s)"
        )
    fig.suptitle(
        title or title_default,
        fontsize=12,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def save_cycle_indicator_pngs(
    features: pd.DataFrame,
    out_dir: Path | str,
    *,
    stem: str = "cycle_indicators",
    per_cell: bool = True,
    max_cells: int = 12,
    dpi: int = 140,
) -> list[Path]:
    """Write overview PNG (+ optional one PNG per cell)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    overview = plot_cycle_indicator_overview(
        features,
        out_dir / f"{stem}_overview.png",
        title=None,
        dpi=dpi,
        max_cells=max_cells,
    )
    if overview:
        saved.append(overview)

    if not per_cell or features is None or features.empty:
        return saved
    if "cell_id" not in features.columns:
        return saved

    cells = list(features["cell_id"].astype(str).unique())
    if len(cells) <= 1:
        return saved

    for cell in cells[:max_cells]:
        sub = features[features["cell_id"].astype(str) == cell]
        path = plot_cycle_indicator_overview(
            sub,
            out_dir / f"{stem}_{_safe_stem(cell)}.png",
            title=f"Cycle indicators — {cell}",
            dpi=dpi,
            max_cells=1,
        )
        if path:
            saved.append(path)
    return saved


@dataclass(frozen=True)
class SohqRestVLinearFit:
    """Linear proxy: SoHQ ≈ a·EoC_restV_end + b·EoD_restV_end + c."""

    a_eoc: float
    b_eod: float
    intercept: float
    r2: float
    rmse_pct: float
    mae_pct: float
    pearson_r: float
    n_points: int

    def predict(self, eoc_end: np.ndarray, eod_end: np.ndarray) -> np.ndarray:
        return self.a_eoc * eoc_end + self.b_eod * eod_end + self.intercept

    def formula_text(self) -> str:
        return (
            f"SoHQ_hat = {self.a_eoc:.4g}·EoC_restV_end "
            f"+ {self.b_eod:.4g}·EoD_restV_end "
            f"+ {self.intercept:.4g}"
        )


def fit_sohq_from_rest_v_end(
    features: pd.DataFrame,
    *,
    eoc_col: str = "EoC_restV_end",
    eod_col: str = "EoD_restV_end",
    sohq_col: str = "SoHQ",
) -> tuple[SohqRestVLinearFit, pd.DataFrame]:
    """Fit SoHQ ~ linear combo of charge/discharge rest V_end."""
    need = [eoc_col, eod_col, sohq_col]
    for col in need:
        if col not in features.columns:
            raise ValueError(f"Missing column: {col}")

    work = features.copy()
    work[sohq_col] = pd.to_numeric(work[sohq_col], errors="coerce")
    work[eoc_col] = pd.to_numeric(work[eoc_col], errors="coerce")
    work[eod_col] = pd.to_numeric(work[eod_col], errors="coerce")
    work = work.dropna(subset=need)
    if len(work) < 3:
        raise ValueError("Need at least 3 valid rows for linear fit")

    y = work[sohq_col].to_numpy(dtype=float)
    x_eoc = work[eoc_col].to_numpy(dtype=float)
    x_eod = work[eod_col].to_numpy(dtype=float)
    design = np.column_stack([x_eoc, x_eod, np.ones(len(work))])
    coef, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    a, b, c = coef
    pred = design @ coef

    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    mae = float(np.mean(np.abs(y - pred)))
    pearson = finite_pearson(pd.Series(y), pd.Series(pred))

    fit = SohqRestVLinearFit(
        a_eoc=float(a),
        b_eod=float(b),
        intercept=float(c),
        r2=r2,
        rmse_pct=rmse,
        mae_pct=mae,
        pearson_r=pearson,
        n_points=len(work),
    )

    out = work.copy()
    out["SoHQ_hat"] = pred
    out["SoHQ_residual"] = y - pred
    return fit, out


def plot_sohq_rest_v_linear_proxy(
    features: pd.DataFrame,
    out_path: Path | str,
    *,
    title: str | None = None,
    dpi: int = 140,
) -> tuple[Path, SohqRestVLinearFit]:
    """PNG: how well EoC/EoD rest V_end linear combo mimics SoHQ."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fit, work = fit_sohq_from_rest_v_end(features)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    x_col = "tagged_cycle" if "tagged_cycle" in work.columns else "cycle"
    x_label = "Tagged cycle #" if x_col == "tagged_cycle" else "Cycle"
    x = pd.to_numeric(work[x_col], errors="coerce").to_numpy(dtype=float)
    y = work["SoHQ"].to_numpy(dtype=float)
    pred = work["SoHQ_hat"].to_numpy(dtype=float)
    resid = work["SoHQ_residual"].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), facecolor="white")

    ax = axes[0]
    ax.plot(x, y, color="#1f77b4", marker="o", markersize=3, linewidth=1.4, label="SoHQ (actual)")
    ax.plot(
        x, pred, color="#ff7f0e", marker="s", markersize=2.5, linewidth=1.4,
        label="SoHQ_hat (linear)",
    )
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel("SoHQ (%)", fontsize=9)
    ax.set_title("SoHQ vs linear proxy", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")

    ax = axes[1]
    lim_lo = min(y.min(), pred.min()) - 1.0
    lim_hi = max(y.max(), pred.max()) + 1.0
    ax.scatter(y, pred, s=18, alpha=0.65, color="#2ca02c", edgecolors="none")
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", linewidth=1, label="y = x")
    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_xlabel("SoHQ actual (%)", fontsize=9)
    ax.set_ylabel("SoHQ_hat (%)", fontsize=9)
    ax.set_title("Predicted vs actual", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")

    ax = axes[2]
    ax.axhline(0.0, color="k", linewidth=0.8, alpha=0.5)
    ax.plot(x, resid, color="#d62728", marker=".", markersize=4, linewidth=1.2)
    ax.set_xlabel(x_label, fontsize=9)
    ax.set_ylabel("Residual (%)", fontsize=9)
    ax.set_title("SoHQ - SoHQ_hat", fontsize=10)
    ax.grid(True, alpha=0.3)

    cell = ""
    if "cell_id" in work.columns and work["cell_id"].notna().any():
        cell = f" · {work['cell_id'].iloc[0]}"

    fig.suptitle(
        title or (
            f"SoHQ linear proxy from rest V_end{cell}\n"
            f"{fit.formula_text()}\n"
            f"R²={fit.r2:.4f}  RMSE={fit.rmse_pct:.3f}%  MAE={fit.mae_pct:.3f}%  "
            f"r={fit.pearson_r:.4f}  n={fit.n_points}"
        ),
        fontsize=10,
        fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path, fit
