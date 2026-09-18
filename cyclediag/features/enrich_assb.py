"""Post-extract enrichment for ASSB / roadmap §9.1 metrics."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from cyclediag.features.cell_meta import CellProtocolMeta
from cyclediag.features.dcir_decompose import (
    decompose_pulse_cycle,
    result_to_dict,
    soc_ratio_features,
)
from cyclediag.features.ocv_drift import attach_ocv_drift_to_features
from cyclediag.features.quality import cycle_quality_metrics
from cyclediag.features.rpt_metrics import (
    Q_RELAX_NOISE_FLOOR_PCT,
    attach_per,
    attach_rcf,
    compute_q_relax_for_blocks,
)
from cyclediag.features.self_discharge import self_discharge_for_cycle

_SOC_ORDER = (80, 50, 20)


def _find_pulse_cycles(
    raw_df: pd.DataFrame,
    *,
    expected_pulse_current: float,
) -> list[int]:
    out: list[int] = []
    for cyc, g in raw_df.groupby("cycle"):
        i = pd.to_numeric(g.get("current"), errors="coerce")
        if i is None or i.isna().all():
            continue
        if float(i.abs().max()) >= 0.5 * expected_pulse_current:
            out.append(int(cyc))
    return sorted(set(out))


def _group_consecutive(cycles: list[int], *, min_len: int = 1) -> list[list[int]]:
    if not cycles:
        return []
    blocks: list[list[int]] = []
    cur = [cycles[0]]
    for c in cycles[1:]:
        if c == cur[-1] + 1:
            cur.append(c)
        else:
            if len(cur) >= min_len:
                blocks.append(cur)
            cur = [c]
    if len(cur) >= min_len:
        blocks.append(cur)
    return blocks


def _q(row: pd.Series) -> float | None:
    for col in ("dchgCapa", "dchg_capa"):
        if col in row.index and pd.notna(row[col]):
            try:
                v = float(row[col])
            except (TypeError, ValueError):
                continue
            if np.isfinite(v):
                return v
    return None


def _write_q_relax_pair(
    out: pd.DataFrame,
    c1: int,
    c2: int,
    *,
    prefix: str,
    noise_floor_pct: float,
) -> bool:
    """Write ``{prefix}`` / ``{prefix}_pct`` / ``{prefix}_significant`` for c1,c2."""
    by = {int(r.cycle): r for r in out.itertuples()}
    if c1 not in by or c2 not in by:
        return False
    q1, q2 = _q(pd.Series(by[c1]._asdict())), _q(pd.Series(by[c2]._asdict()))
    if q1 is None or q2 is None or q2 == 0 or min(q1, q2) < 10:
        return False
    q_relax = q2 - q1
    q_relax_pct = q_relax / q2 * 100.0
    for col, default in (
        (prefix, np.nan),
        (f"{prefix}_pct", np.nan),
        (f"{prefix}_significant", None),
    ):
        if col not in out.columns:
            out[col] = default
    for c in (c1, c2):
        mask = out["cycle"] == c
        out.loc[mask, prefix] = q_relax
        out.loc[mask, f"{prefix}_pct"] = q_relax_pct
        out.loc[mask, f"{prefix}_significant"] = abs(q_relax_pct) > noise_floor_pct
    return True


def attach_q_relax_from_dcir_blocks(
    features: pd.DataFrame,
    dcir_blocks: list[list[int]],
    *,
    noise_floor_pct: float = Q_RELAX_NOISE_FLOOR_PCT,
    rpt_blocks: list[list[int]] | None = None,
) -> pd.DataFrame:
    """Attach Q_relax with explicit source aliases.

    * ``Q_relax_dcir*`` — two high-Q cycles immediately before each DC-IR block
    * ``Q_relax_rpt*`` — first two cycles of each RPT block (when provided)
    * ``Q_relax*`` (primary) — prefers RPT definition, else DCIR-pre
    * ``Q_relax_source`` — ``rpt_block`` | ``dcir_pre``
    """
    out = features.copy()
    for col in (
        "Q_relax", "Q_relax_pct", "Q_relax_significant",
        "Q_relax_dcir", "Q_relax_dcir_pct", "Q_relax_dcir_significant",
        "Q_relax_rpt", "Q_relax_rpt_pct", "Q_relax_rpt_significant",
        "Q_relax_source",
    ):
        if col not in out.columns:
            out[col] = (
                None if col.endswith("significant") or col == "Q_relax_source" else np.nan
            )

    for block in dcir_blocks or []:
        if not block:
            continue
        start = int(block[0])
        _write_q_relax_pair(
            out, start - 2, start - 1,
            prefix="Q_relax_dcir", noise_floor_pct=noise_floor_pct,
        )

    for block in rpt_blocks or []:
        if len(block) < 2:
            continue
        _write_q_relax_pair(
            out, int(block[0]), int(block[1]),
            prefix="Q_relax_rpt", noise_floor_pct=noise_floor_pct,
        )

    # Primary: RPT when present on the row, else DCIR-pre.
    for idx in out.index:
        rpt = out.at[idx, "Q_relax_rpt_pct"]
        dcir = out.at[idx, "Q_relax_dcir_pct"]
        if pd.notna(rpt):
            out.at[idx, "Q_relax"] = out.at[idx, "Q_relax_rpt"]
            out.at[idx, "Q_relax_pct"] = rpt
            out.at[idx, "Q_relax_significant"] = out.at[idx, "Q_relax_rpt_significant"]
            out.at[idx, "Q_relax_source"] = "rpt_block"
        elif pd.notna(dcir):
            out.at[idx, "Q_relax"] = out.at[idx, "Q_relax_dcir"]
            out.at[idx, "Q_relax_pct"] = dcir
            out.at[idx, "Q_relax_significant"] = out.at[idx, "Q_relax_dcir_significant"]
            out.at[idx, "Q_relax_source"] = "dcir_pre"
    return out


def first_capa_baseline(features: pd.DataFrame, dcir_blocks: list[list[int]]) -> int | None:
    """First post-formation capa cycle2 (before first DC-IR)."""
    if not dcir_blocks:
        return None
    start = int(dcir_blocks[0][0])
    c2 = start - 1
    if (features["cycle"] == c2).any():
        return c2
    return None


def _set_cell(out: pd.DataFrame, mask, col: str, val) -> None:
    if col not in out.columns:
        out[col] = pd.Series(None, index=out.index, dtype=object)
    if isinstance(val, (dict, list)):
        from copy import deepcopy

        out[col] = out[col].astype(object)
        for idx in out.index[mask]:
            out.at[idx, col] = deepcopy(val)
        return
    if out[col].dtype != object and not (
        isinstance(val, (int, float, np.floating, np.integer)) and not isinstance(val, (bool, np.bool_))
    ):
        out[col] = out[col].astype(object)
    if isinstance(val, (bool, np.bool_)):
        out[col] = out[col].astype(object)
        out.loc[mask, col] = bool(val)
    elif val is None or (isinstance(val, float) and not np.isfinite(val)):
        out.loc[mask, col] = np.nan if out[col].dtype != object else None
    elif isinstance(val, str):
        out[col] = out[col].astype(object)
        out.loc[mask, col] = val
    else:
        try:
            out.loc[mask, col] = val
        except (TypeError, ValueError):
            out[col] = out[col].astype(object)
            out.loc[mask, col] = val


def enrich_feature_table(
    features: pd.DataFrame,
    raw_df: pd.DataFrame | None,
    *,
    rest_current_max: float | None = None,
    expected_pulse_current: float | None = None,
    protocol_meta: CellProtocolMeta | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Attach Q_relax, DCIR decompose, self-discharge, quality, RCF, PER."""
    pm = protocol_meta or CellProtocolMeta()
    if rest_current_max is None:
        rest_current_max = pm.rest_current_max_a
    if expected_pulse_current is None:
        expected_pulse_current = pm.dcir_pulse_current_a
    meta: dict[str, Any] = {
        "baseline_cycle_auto": None,
        "dcir_blocks": [],
        "capa_anchor_cycles": [],
        "ocv_drift_blocks": [],
    }
    if features is None or features.empty:
        return features, meta

    out = features.copy()
    if raw_df is None or "cycle" not in raw_df.columns:
        return out, meta

    pulse = _find_pulse_cycles(raw_df, expected_pulse_current=expected_pulse_current)
    blocks = _group_consecutive(pulse, min_len=1)
    # prefer length-3 SOC blocks
    dcir_blocks = [b[:3] for b in blocks if len(b) >= 3] or blocks
    meta["dcir_blocks"] = dcir_blocks
    meta["baseline_cycle_auto"] = first_capa_baseline(out, dcir_blocks)

    # RPT 2-cycle pairs sit immediately before each DC-IR block on SJ900.
    # Prefer protocol-derived blocks when available; else use the pre-DCIR pair.
    rpt_blocks: list[list[int]] = []
    try:
        from cyclediag.io.cycle_protocol import build_protocol_exclusion

        se = raw_df.groupby(["cycle", "StepNo"], as_index=False).tail(1) if "StepNo" in raw_df.columns else raw_df
        prot = build_protocol_exclusion(se)
        rpt_blocks = list(prot.rpt_blocks or [])
    except Exception:
        rpt_blocks = []
    if not rpt_blocks:
        rpt_blocks = [
            [int(b[0]) - 2, int(b[0]) - 1]
            for b in dcir_blocks if b
        ]
    out = attach_q_relax_from_dcir_blocks(out, dcir_blocks, rpt_blocks=rpt_blocks)
    meta["rpt_blocks"] = rpt_blocks

    capa_anchors: set[int] = set()
    for block in dcir_blocks:
        start = int(block[0])
        for c in (start - 2, start - 1):
            if (out["cycle"] == c).any():
                capa_anchors.add(c)
    meta["capa_anchor_cycles"] = sorted(capa_anchors)
    # routine ≈ not pulse and not tiny Q
    pulse_set = set(pulse)
    routine_mask = ~out["cycle"].isin(pulse_set)
    out = attach_rcf(out, routine_mask=routine_mask, rpt_cycles=capa_anchors)
    out = attach_per(out, dI_A=pm.per_delta_i_a)

    for soc in _SOC_ORDER:
        for k in (
            "R_ohmic", "R_ct", "tau_ct", "A_diff", "R_30s_total",
            "R_ohmic_frac", "dcir_fit_valid",
            "self_discharge_rate", "sd_fit_valid",
            "R_recovery_tau1", "R_recovery_tau2", "V_inf_est",
        ):
            col = f"{k}_soc{soc}"
            if col not in out.columns:
                out[col] = np.nan

    for block in dcir_blocks:
        soc_rows: list[dict[str, Any]] = []
        for i, cyc in enumerate(block[:3]):
            soc = _SOC_ORDER[i] if len(block) >= 3 and i < 3 else None
            g = raw_df[raw_df["cycle"] == int(cyc)]
            fit = decompose_pulse_cycle(
                g,
                rest_current_max=rest_current_max,
                expected_pulse_current=expected_pulse_current,
            )
            sd = self_discharge_for_cycle(
                g,
                rest_current_max=rest_current_max,
                expected_pulse_current=expected_pulse_current,
            )
            payload = result_to_dict(fit)
            if "scientific_validity" in payload:
                payload["ecm_scientific_validity"] = payload.pop("scientific_validity")
            mask = out["cycle"] == int(cyc)
            if soc is None:
                for k, val in payload.items():
                    _set_cell(out, mask, k, val)
                _set_cell(out, mask, "self_discharge_rate", sd.get("self_discharge_rate"))
                continue
            suf = f"_soc{soc}"
            for k, val in payload.items():
                _set_cell(out, mask, f"{k}{suf}", val)
            _set_cell(out, mask, f"self_discharge_rate_soc{soc}", sd.get("self_discharge_rate"))
            _set_cell(out, mask, f"sd_fit_valid_soc{soc}", sd.get("sd_fit_valid"))
            _set_cell(out, mask, f"V_inf_rest_soc{soc}", sd.get("V_inf_rest"))
            _set_cell(out, mask, f"relax_completeness_soc{soc}", sd.get("relax_completeness"))
            soc_rows.append({"soc": soc, "R_30s_total": fit.R_30s_total})
        ratios = soc_ratio_features(soc_rows)
        for cyc in block[:3]:
            mask = out["cycle"] == int(cyc)
            for k, val in ratios.items():
                _set_cell(out, mask, k, val)

    out, ocv_blocks = attach_ocv_drift_to_features(
        out,
        dcir_blocks,
        raw_df,
        rest_current_max=rest_current_max,
        expected_pulse_current=expected_pulse_current,
    )
    meta["ocv_drift_blocks"] = ocv_blocks.to_dict(orient="records") if not ocv_blocks.empty else []

    for col in (
        "samples_per_mV", "v_noise_sigma", "quant_step_est", "dqdv_snr",
        "rest_sufficiency", "pulse_sample_count_1s", "pulse_current_stability",
        "leg_completeness", "temperature_available", "quality_score",
        "quality_gate_failed_groups",
    ):
        if col not in out.columns:
            out[col] = pd.Series([None] * len(out), dtype=object)

    for cyc in out["cycle"].astype(int).unique():
        g = raw_df[raw_df["cycle"] == int(cyc)]
        qm = cycle_quality_metrics(g, rest_current_max=rest_current_max)
        mask = out["cycle"] == int(cyc)
        for k, val in qm.items():
            key = "quality_scientific_validity" if k == "scientific_validity" else k
            _set_cell(out, mask, key, val)

    return out, meta
