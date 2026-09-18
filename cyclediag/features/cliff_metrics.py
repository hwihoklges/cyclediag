"""Absolute-Ah cliff (depletion / Gr–Si transition) metrics for Si/Gr mechanism discrimination."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from cyclediag.features.dqdv_peaks import DqdvPeakConfig, DEFAULT_DQDV_PEAK_CONFIG, compute_dvdq


def _effective_sg_width(q: np.ndarray, config: DqdvPeakConfig) -> float | None:
    """Approximate Savitzky–Golay smoothing width in Ah (for QC)."""
    if len(q) < 3:
        return None
    q_span = float(np.nanmax(q) - np.nanmin(q))
    if q_span <= 0:
        return None
    n = config.n_interp
    dq = q_span / max(n - 1, 1)
    return float(dq * config.sg_window)


def compute_cliff_metrics(
    q: np.ndarray,
    v: np.ndarray,
    *,
    thr_factor: float = 2.0,
    mid_q_range: tuple[float, float] = (0.4, 0.6),
    config: DqdvPeakConfig | None = None,
    min_continuous: int = 3,
) -> dict[str, Any]:
    """Extract cliff onset in **absolute Ah** (no per-cycle SOC normalization).

    Algorithm
    ---------
    1. ``compute_dvdq`` on absolute cumulative discharge Q.
    2. Threshold = ``thr_factor`` × median |dV/dQ| in mid-Q band (40–60 % of Qmax).
     3. Primary: strongest curvature peak in the interior 20–92 % capacity band.
     4. Secondary: first sustained |dV/dQ| ≥ thr in the rear 45–92 % band,
         for at least ``min_continuous`` consecutive points.
    5. Fallback: argmax |dV/dQ| in [0.15·Qmax, 0.92·Qmax] (excludes final spike only).

    Returns
    -------
    Q_cliff_abs, Q_cliff_norm, Q_tail_abs, V_at_cliff, dVdQ_post_cliff,
    cliff_thr_used, cliff_valid, cliff_sg_width_ah
    """
    cfg = config or DEFAULT_DQDV_PEAK_CONFIG
    empty: dict[str, Any] = {
        "Q_cliff_abs": None,
        "Q_cliff_norm": None,
        "Q_tail_abs": None,
        "V_at_cliff": None,
        "dVdQ_post_cliff": None,
        "cliff_thr_used": None,
        "cliff_valid": False,
        "cliff_sg_width_ah": _effective_sg_width(np.asarray(q, dtype=float), cfg),
    }
    qx, dvdq = compute_dvdq(np.asarray(q, dtype=float), np.asarray(v, dtype=float), cfg)
    if len(qx) < 10:
        return empty

    qmin = float(np.nanmin(qx))
    qmax = float(np.nanmax(qx))
    q_span = qmax - qmin
    if q_span <= 1e-9:
        return empty

    mag = np.abs(dvdq)
    fin = np.isfinite(mag) & np.isfinite(qx)
    if fin.sum() < 10:
        return empty

    qx, mag, dvdq = qx[fin], mag[fin], dvdq[fin]
    qmax = float(np.nanmax(qx))

    # mid-Q band threshold (absolute Q — not normalized SOC)
    lo_mid = qmin + mid_q_range[0] * (qmax - qmin)
    hi_mid = qmin + mid_q_range[1] * (qmax - qmin)
    mid = (qx >= lo_mid) & (qx <= hi_mid)
    if mid.sum() < 5:
        return empty

    med = float(np.nanmedian(mag[mid]))
    if not math.isfinite(med) or med <= 0:
        return empty
    thr = med * thr_factor
    empty["cliff_thr_used"] = thr

    # Primary: strongest interior curvature before final depletion. A fixed
    # 70 % preference can select roundoff when the transition moves later.
    dmag = np.abs(np.gradient(mag, qx))
    q_cut_hi = qmax - 0.08 * max(qmax - qmin, 1e-9)
    band = (
        (qx >= qmin + 0.20 * (qmax - qmin))
        & (qx <= q_cut_hi)
        & np.isfinite(dmag)
    )
    cliff_idx: int | None = None
    if band.sum() >= 7:
        idx_band = np.flatnonzero(band)
        thr_d = float(np.nanpercentile(dmag[band], 70))
        local_peaks: list[int] = []
        for pos, gi in enumerate(idx_band):
            if dmag[gi] < thr_d:
                continue
            left = dmag[idx_band[pos - 1]] if pos > 0 else 0.0
            right = dmag[idx_band[pos + 1]] if pos < len(idx_band) - 1 else 0.0
            if dmag[gi] >= left and dmag[gi] >= right:
                local_peaks.append(int(gi))
        candidates = local_peaks or [int(idx_band[int(np.nanargmax(dmag[band]))])]
        cliff_idx = max(candidates, key=lambda i: dmag[i])

    q_rear_lo = qmin + 0.45 * (qmax - qmin)
    q_rear_hi = qmax - 0.08 * max(qmax - qmin, 1e-9)
    rear = (qx >= q_rear_lo) & (qx <= q_rear_hi)

    # Secondary: forward scan in rear half for sustained |dV/dQ| exceedance
    if cliff_idx is None:
        rear_idx = np.flatnonzero(rear)
        if len(rear_idx) >= min_continuous:
            run = 0
            for k, ri in enumerate(rear_idx):
                if mag[ri] >= thr:
                    run += 1
                    if run >= min_continuous:
                        cliff_idx = int(rear_idx[k - min_continuous + 1])
                        break
                else:
                    run = 0

    # fallback: dominant steepness in upper-mid band (exclude top 8 %)
    if cliff_idx is None:
        search = (
            (qx >= qmin + 0.15 * (qmax - qmin))
            & (qx <= qmax - 0.08 * max(qmax - qmin, 1e-9))
        )
        if search.any():
            local = np.where(search)[0]
            cliff_idx = int(local[int(np.nanargmax(mag[search]))])

    if cliff_idx is None:
        return empty

    q_cliff = float(qx[cliff_idx])
    q_total = qmax - qmin
    q_tail = float(qmax - q_cliff) if qmax > q_cliff else 0.0

    # V at cliff via original arrays
    v_arr = np.asarray(v, dtype=float)
    q_arr = np.asarray(q, dtype=float)
    m2 = np.isfinite(v_arr) & np.isfinite(q_arr)
    v_at = None
    if m2.sum() >= 2:
        order = np.argsort(q_arr[m2])
        qa, va = q_arr[m2][order], v_arr[m2][order]
        v_at = float(np.interp(q_cliff, qa, va))

    post = mag[cliff_idx : min(cliff_idx + 20, len(mag))]
    d_post = float(np.nanmean(post)) if post.size else None

    return {
        "Q_cliff_abs": q_cliff,
        "Q_cliff_norm": q_cliff / q_total if q_total > 1e-9 else None,
        "Q_tail_abs": q_tail,
        "V_at_cliff": v_at,
        "dVdQ_post_cliff": d_post,
        "cliff_thr_used": thr,
        "cliff_valid": True,
        "cliff_sg_width_ah": _effective_sg_width(q_arr, cfg),
    }


def cliff_metrics_to_row(metrics: dict[str, Any], *, prefix: str = "") -> dict[str, Any]:
    """Prefix keys for wide-table merge (default: dchg_)."""
    p = prefix or ""
    return {f"{p}{k}" if not k.startswith(p) else k: v for k, v in metrics.items()}
