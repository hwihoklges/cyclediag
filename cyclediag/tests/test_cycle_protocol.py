"""Tests for RPT / capacheck protocol exclusion."""

from __future__ import annotations

import pandas as pd

from cyclediag.io.cycle_protocol import build_protocol_exclusion


def _step_row(cycle: int, dchg_mah: float, n_steps: int = 4, i_abs: float = 1000.0) -> pd.DataFrame:
    """Minimal step-end table: routine = 4 rows, discharge cap on discharge step."""
    rows = []
    for i in range(n_steps):
        st = ["charge", "rest", "discharge", "rest"][i % 4] if n_steps == 4 else "charge"
        dchg = dchg_mah if st == "discharge" else 0.0
        if st == "charge":
            cur = i_abs
        elif st == "discharge":
            cur = -i_abs
        else:
            cur = 0.0
        rows.append({
            "cycle": cycle,
            "step_type": st,
            "discharge_capacity": dchg / 1000.0,
            "current": cur,
        })
    return pd.DataFrame(rows)


def test_routine_cycle_not_excluded():
    flags = build_protocol_exclusion(_step_row(10, 72000.0))
    assert 10 not in flags.excluded
    assert flags.flags.loc[flags.flags["cycle"] == 10, "protocol_kind"].iloc[0] == "routine"


def test_capacheck_low_dchg_excluded():
    rows = [
        {"cycle": 4, "step_type": "charge", "discharge_capacity": 0, "current": 1000},
        {"cycle": 4, "step_type": "rest", "discharge_capacity": 0, "current": 0},
        {"cycle": 4, "step_type": "discharge", "discharge_capacity": 14440, "current": -1000},
        {"cycle": 4, "step_type": "rest", "discharge_capacity": 0, "current": 0},
        {"cycle": 4, "step_type": "discharge", "discharge_capacity": 14440, "current": -1000},
        {"cycle": 4, "step_type": "rest", "discharge_capacity": 0, "current": 0},
        {"cycle": 4, "step_type": "charge", "discharge_capacity": 0, "current": 1000},
        {"cycle": 4, "step_type": "rest", "discharge_capacity": 0, "current": 0},
    ]
    step_df = pd.DataFrame(rows)
    step_df["discharge_capacity"] /= 1000.0  # Explicit raw mAh boundary.
    excl = build_protocol_exclusion(step_df)
    assert 4 in excl.excluded
    assert 4 in excl.capacheck_cycles


def test_faded_routine_not_excluded():
    excl = build_protocol_exclusion(_step_row(481, 49900.0, n_steps=4))
    assert 481 not in excl.excluded


def test_post_rpt_buffer_after_block():
    parts = [
        _step_row(108, 72000.0),
        _step_row(109, 14000.0, n_steps=8),
        _step_row(110, 14000.0, n_steps=6),
    ]
    step_df = pd.concat(parts, ignore_index=True)
    excl = build_protocol_exclusion(step_df, post_rpt_exclude=2)
    assert 109 in excl.excluded
    assert 111 in excl.post_rpt_cycles
    assert 112 in excl.post_rpt_cycles
    assert 113 not in excl.excluded


def test_capa_full_0p33c_before_rpt_detected():
    """Two 0.33C full-capa cycles before DC-IR block (M01Ch022 107/108 pattern)."""
    parts = []
    for c in range(100, 107):
        parts.append(_step_row(c, 66000.0, i_abs=38.7))
    parts.append(_step_row(107, 68900.0, i_abs=25.8))
    parts.append(_step_row(108, 68900.0, i_abs=25.8))
    parts.append(_step_row(109, 14000.0, n_steps=8, i_abs=25.8))
    parts.append(_step_row(110, 21000.0, n_steps=6, i_abs=25.8))
    parts.append(_step_row(111, 21000.0, n_steps=6, i_abs=25.8))
    excl = build_protocol_exclusion(pd.concat(parts, ignore_index=True))
    assert 107 in excl.capa_full_cycles
    assert 108 in excl.capa_full_cycles
    assert 107 in excl.excluded
    assert 105 not in excl.capa_full_cycles
    assert excl.flags.loc[excl.flags.cycle == 107, "protocol_kind"].iloc[0] == "capa_full"
