"""Tests for LGES cycle feature extraction."""

from __future__ import annotations

import pandas as pd
import pytest

from cyclediag.features.extract import FeatureConfig, extract_features_table
from cyclediag.features.lges_catalog import FEATURE_SET_LGES, all_lges_feature_columns
from cyclediag.features.lges_extract import LgesExtractConfig, extract_lges_features_table
from cyclediag.io.cycler_csv import ColumnMap, normalize_cycler_dataframe


def _cycle_df(cycle: int = 1, q_scale: float = 1.0) -> pd.DataFrame:
    rows = []
    t = 0.0
    for st, n, v0, dv, cur in (
        ("charge", 30, 3.0, 0.04, 1.0),
        ("rest", 65, 4.2, 0.0, 0.0),
        ("discharge", 30, 4.1, -0.04, -1.0),
        ("rest", 65, 3.2, 0.0, 0.0),
    ):
        for i in range(n):
            rows.append({
                "cycle": cycle,
                "step_type": st,
                "voltage": v0 + i * dv if st != "rest" else v0,
                "charge_capacity": abs(i * 10.0 * q_scale) if st == "charge" else 100.0 * q_scale,
                "discharge_capacity": abs(i * 9.5 * q_scale) if st == "discharge" else 0.0,
                "capacity": abs(i * 10.0 * q_scale) if st == "charge" else abs(i * 9.5 * q_scale),
                "current": cur,
                "time": t,
                "step_time": float(i),
            })
            t += 1.0
    return pd.DataFrame(rows)


def test_lges_one_row_per_cycle():
    df = _cycle_df()
    out = extract_lges_features_table(df, filepath="cell_a.csv", config=LgesExtractConfig())
    assert len(out) == 1
    assert out["feature_set"].iloc[0] == FEATURE_SET_LGES
    assert out["cell_id"].iloc[0] == "cell_a"
    assert out["chgCapa"].notna().iloc[0]
    assert out["dchgCapa"].notna().iloc[0]
    assert out["EoC_restV_60s"].notna().iloc[0]
    assert out["EoD_restV_end"].notna().iloc[0]


def test_lges_via_feature_set_dispatch():
    df = _cycle_df()
    cfg = FeatureConfig(feature_set=FEATURE_SET_LGES)
    out = extract_features_table(df, filepath="x.csv", config=cfg)
    assert len(out) == 1
    assert "EoC_restV_init" in out.columns


def test_lges_delta_and_sohq_two_cycles():
    df = pd.concat([_cycle_df(1), _cycle_df(2, q_scale=0.95)], ignore_index=True)
    out = extract_lges_features_table(df, config=LgesExtractConfig())
    assert len(out) == 2
    assert out.loc[out["cycle"] == 2, "SoHQ"].iloc[0] == pytest.approx(95.0, rel=0.02)
    assert "delta_EoC_restV_60s" in out.columns


def test_lges_catalog_columns_present():
    df = _cycle_df()
    out = extract_lges_features_table(df, config=LgesExtractConfig())
    for col in all_lges_feature_columns():
        assert col in out.columns


def test_lges_discharge_capacity_from_studio_columns():
    raw = pd.DataFrame({
        "TotalCycle": [1] * 95,
        "StepType": ["charge"] * 30 + ["rest"] * 15 + ["discharge"] * 30 + ["rest"] * 20,
        "Voltage": [3.5] * 95,
        "Capacity": [3000.0] * 95,
        "DischargeCapacity": [0.0] * 45 + [2850.0] * 50,
        "Current": [1.0] * 30 + [0.0] * 15 + [-1.0] * 30 + [0.0] * 20,
        "TotalTime_sec": list(range(95)),
        "StepTime_sec": list(range(95)),
    })
    df = normalize_cycler_dataframe(raw, ColumnMap(units={"capacity": "mAh"}))
    out = extract_lges_features_table(df, config=LgesExtractConfig())
    assert out["dchgCapa"].iloc[0] == pytest.approx(2.85, rel=0.01)
    assert out["chgCapa"].iloc[0] == pytest.approx(3.0, rel=0.01)
