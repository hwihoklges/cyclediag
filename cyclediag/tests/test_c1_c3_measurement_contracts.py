"""Analytical contracts only; not empirical calibration or causal validation."""

from dataclasses import asdict, replace
import json

import numpy as np
import pandas as pd
import pytest

from cyclediag.api import extract_features
from cyclediag.diagnosis.pattern_scoring import _evidence_for_term
from cyclediag.features.lges_extract import (
    LgesExtractConfig, apply_lges_delta_features, extract_lges_cycle_row,
    extract_lges_features_table,
)
from cyclediag.features.lges_extra_indicators import hysteresis_metrics, correct_r_to_25c
from cyclediag.features.measurement_contracts import (
    MeasurementTemperature, TemperatureCalibration, calibrated_resistance,
    observed_leg, temp_correct_ir, VE_SEMANTICS,
)
from cyclediag.provenance import canonical_json, config_hash, read_feature_csv, save_features_csv


def cycle(offset=0.0):
    q = np.linspace(0, 1, 21)
    return pd.DataFrame({
        "cycle": 1, "step_type": ["charge"] * 21 + ["discharge"] * 21,
        "voltage": np.r_[np.full(21, 4.0), np.full(21, 3.5)],
        "charge_capacity": np.r_[offset + q, np.full(21, offset + 1)],
        "discharge_capacity": np.r_[np.zeros(21), offset + .8 * q],
        "current": np.r_[np.ones(21), -np.ones(21)],
        "step_time": np.r_[np.arange(21), np.arange(21)].astype(float),
        "time": np.arange(42, dtype=float), "temperature": 40.0,
    })


def cfg(**kwargs):
    return LgesExtractConfig(with_diagnosis=False, enrich_assb=False, **kwargs)


def calibration(**kwargs):
    base = TemperatureCalibration(20000., 25., (0., 60.), "synthetic-fit",
        "early-leg abs(deltaV)/abs(I) 10s proxy", "endpoint_proxy", "synthetic",
        "10s", "independent sensor at measurement")
    return replace(base, **kwargs)


def measurement(cal=None, temperature_c=40., **kwargs):
    cal = cal or calibration()
    fields = {k: getattr(cal, k) for k in (
        "measurement_definition", "soc", "protocol", "time_scale", "temperature_source")}
    return MeasurementTemperature(temperature_c, **(fields | kwargs))


@pytest.mark.parametrize("scale,offset", [(1, 0), (72, 7), (.02, -3)])
def test_hysteresis_reversible_offset_and_scale(scale, offset):
    q = np.linspace(0, 1, 31)
    out = hysteresis_metrics(scale*q + offset, 3 + q, 2*scale*q + offset, 4 - q)
    assert out["hyst_area"] == pytest.approx(0, abs=2e-12)
    out = hysteresis_metrics(scale*q + offset, 3.1 + q, 2*scale*q + offset, 4 - q)
    assert out["hyst_area"] == pytest.approx(.096, abs=2e-12)
    assert [out[f"hyst_area_{b}"] for b in ("low", "mid", "high")] == pytest.approx([.018, .06, .018])
    assert out["hysteresis_metadata"]["area_unit"] == "V"
    assert out["hysteresis_metadata"]["physical_SOC_verified"] is False
    assert out["hysteresis_metadata"]["equilibrium_verified"] is False


@pytest.mark.parametrize("center,band", [(.1, "low"), (.5, "mid"), (.9, "high")])
def test_localized_triangle_and_exact_boundaries(center, band):
    q = np.sort(np.r_[0., .02, .2, .8, .98, 1., center-.03, center, center+.03])
    delta = .1 * np.maximum(0, 1 - np.abs(q-center)/.03)
    out = hysteresis_metrics(q, 3+q+delta, 1-q[::-1], (3+q)[::-1], n_grid=3)
    assert out["hyst_area"] == pytest.approx(.003)
    assert out[f"hyst_area_{band}"] == pytest.approx(.003)
    assert sum(out[f"hyst_area_{b}"] for b in ("low", "mid", "high")) == pytest.approx(.003)


@pytest.mark.parametrize("q,v,reason", [
    ([0, .5, .2, 1], [3, 3.5, 3.2, 4], "counter_reset"),
    ([0, .5, .5, 1], [3, 3.5, 3.6, 4], "conflicting_voltage_at_same_counter"),
    ([0, np.nan, 1], [3, 3.5, 4], "nonfinite_counter_or_voltage"),
    ([0, .5, 1], [3, np.inf, 4], "nonfinite_counter_or_voltage"),
])
def test_hysteresis_rejects_acquisition_defects(q, v, reason):
    out = hysteresis_metrics(q, v, [0, 1], [4, 3])
    assert out["hyst_area"] is None
    assert out["hysteresis_metadata"]["reason"]["charge"] == reason


def test_identical_duplicate_pairs_are_deduplicated():
    out = hysteresis_metrics([0, .5, .5, 1], [3, 3.5, 3.5, 4], [0, 1], [4, 3])
    assert out["hyst_area"] == pytest.approx(0, abs=1e-12)


def test_temperature_default_is_unavailable_and_raw_unchanged():
    assert correct_r_to_25c(10, 40) is None
    row = extract_lges_cycle_row(cycle(), 1, config=cfg())
    for feature in ("EoC_dchgR_10s", "EoD_chgR_10s"):
        assert row[feature] == 0
        assert row[f"{feature}_T25"] is None
        assert row["temperature_correction"][feature]["status"] == "not_requested"


@pytest.mark.parametrize("ea", [0., 20000., -12000.])
def test_explicit_temperature_identity_roundtrip(ea):
    assert temp_correct_ir(10, 25, 25, ea) == pytest.approx(10)
    mapped = temp_correct_ir(10, 40, 25, ea)
    assert temp_correct_ir(mapped, 25, 40, ea) == pytest.approx(10)
    if ea == 0:
        assert mapped == pytest.approx(10)


@pytest.mark.parametrize("r,t,ref,ea", [
    (10, -273.15, 25, 20000), (10, 25, -274, 20000),
    (10, np.nan, 25, 20000), (np.inf, 25, 25, 0),
    (10, 40, 25, np.inf), (10, 40, 25, 1e308), (10, 40, 25, -1e308),
])
def test_invalid_temperature_and_overflow_are_unknown(r, t, ref, ea):
    assert temp_correct_ir(r, t, ref, ea) is None


def test_calibration_requires_definition_match_and_both_temperatures_in_domain():
    cal = calibration()
    value, meta = calibrated_resistance(10, cal, measurement())
    assert value > 10 and meta["status"] == "ok"
    assert meta["calibration_verified"] is False
    for bad in (replace(cal, reference_c=70), replace(cal, ea_j_mol=np.nan),
                replace(cal, calibration_id=""), replace(cal, domain_c=(-274, 60))):
        value, meta = calibrated_resistance(10, bad, measurement())
        assert value is None and meta["status"] == "unknown"
        json.loads(canonical_json(meta))
    for m in (None, measurement(temperature_c=70), measurement(soc="different")):
        value, meta = calibrated_resistance(10, cal, m)
        assert value is None and meta["status"] == "unknown"


def test_extraction_never_substitutes_average_leg_temperature():
    key = "EoC_dchgR_10s"
    cal = calibration()
    row = extract_lges_cycle_row(cycle(), 1, config=cfg(temperature_calibrations={key: cal}))
    assert row["dchg_temp_avg"] == 40
    assert row[f"{key}_T25"] is None
    assert row["temperature_correction"][key]["reason"] == "missing_measurement_temperature"


def test_observed_counter_offset_and_energy_identity():
    a, b = [extract_lges_cycle_row(cycle(offset), 1, config=cfg()) for offset in (0, 5)]
    assert a["chgCapa"] == 1 and b["chgCapa"] == 6
    assert a["CE"] == pytest.approx(80)
    assert b["CE"] == pytest.approx(5.8/6*100)
    for row in (a, b):
        assert row["chgCapa_observed"] == pytest.approx(1)
        assert row["dchgCapa_observed"] == pytest.approx(.8)
        assert row["chg_E_observed"] == pytest.approx(4)
        assert row["dchg_E_observed"] == pytest.approx(2.8)
        assert row["VE"] == row["VE_observed"] == pytest.approx(3.5/4)
        assert row["CE_observed"]/100 * row["VE_observed"] == pytest.approx(row["EE_observed"])
        assert row["observed_window"]["charge"]["full_cycle_verified"] is False
    assert b["observed_window"]["charge"]["counter_start_Ah"] == 5


@pytest.mark.parametrize("defect,reason", [
    ("intervening", "disjoint_acquisition_segments"),
    ("index_gap", "noncontiguous_source_index"),
    ("reset", "acquisition_time_reset"),
    ("missing_time", "missing_acquisition_time"),
    ("nan_counter", "nonfinite_counter_or_voltage"),
    ("counter_reset", "counter_reset"),
])
def test_extraction_does_not_integrate_disjoint_or_invalid_legs(defect, reason):
    frame = cycle()
    if defect == "intervening":
        frame.loc[10, "step_type"] = "rest"
    elif defect == "index_gap":
        frame = frame.drop(index=10)
    elif defect == "reset":
        frame.loc[10:, "step_time"] -= 10
    elif defect == "missing_time":
        frame.loc[10, "time"] = np.nan
    elif defect == "nan_counter":
        frame.loc[10, "charge_capacity"] = np.nan
    elif defect == "counter_reset":
        frame.loc[10, "charge_capacity"] = 0
    row = extract_lges_cycle_row(frame, 1, config=cfg())
    assert row["chg_E"] is None and row["chgCapa_observed"] is None
    assert row["VE"] is None and row["hyst_area"] is None
    assert row["observed_window"]["charge"]["reason"] == reason
    assert row["hysteresis_metadata"]["status"] == "unknown"


def test_observed_leg_rejects_missing_clock_and_data_point_gap():
    frame = cycle().iloc[:21].drop(columns=["time", "step_time"])
    assert observed_leg(frame, frame.charge_capacity.to_numpy())[3]["reason"] == "missing_acquisition_time"
    frame["time"] = np.arange(21)
    frame["data_point"] = np.r_[np.arange(10), np.arange(11, 22)]
    assert observed_leg(frame, frame.charge_capacity.to_numpy())[3]["reason"] == "noncontiguous_data_point"


def test_ve_diagnosis_requires_matching_current_and_baseline_semantics():
    row = extract_lges_cycle_row(cycle(), 1, config=cfg())
    term = {"feature": "VE", "direction": "decrease_vs_baseline", "scale": .05}
    assert _evidence_for_term(row, term, row)[1] == pytest.approx(0)
    assert _evidence_for_term(row, term, {"VE": .95})[1] is None
    assert _evidence_for_term({"VE": .8}, term, row)[1] is None
    wrong = dict(row, feature_semantics={"VE": "legacy"})
    assert _evidence_for_term(row, term, wrong)[1] is None
    assert row["feature_semantics"]["VE"] == VE_SEMANTICS


def test_corrected_deltas_do_not_mix_legacy_baselines():
    frame = pd.concat([cycle(), cycle().assign(cycle=2)], ignore_index=True)
    table = extract_lges_features_table(frame, config=cfg())
    assert table.loc[1, "delta_chg_V_avg"] == pytest.approx(0)
    table.at[0, "feature_semantics"] = {}
    out = apply_lges_delta_features(table)
    assert pd.isna(out.loc[1, "delta_chg_V_avg"])
    assert pd.isna(out.loc[1, "delta_hyst_area"])


def test_high_level_api_sidecar_and_calibration_provenance(tmp_path):
    key = "EoC_dchgR_10s"
    cal = calibration(reference_c=30.)
    config = cfg(temperature_calibrations={key: cal},
                 measurement_temperatures={"1": {key: measurement(cal)}})
    assert config_hash(config) == config_hash(asdict(config))
    assert config_hash(config) != config_hash(cfg())
    table = extract_features(cycle(5), config=config)
    default = extract_features(cycle(5), config=cfg())
    assert table.attrs["provenance"]["algorithm_schema"] == "cyclediag_scientific_v2"
    assert table.attrs["provenance"]["config_sha256"] != default.attrs["provenance"]["config_sha256"]
    assert table.loc[0, "temperature_correction"][key]["status"] == "ok"
    assert pd.isna(table.loc[0, f"{key}_T25"])
    assert table.loc[0, f"{key}_Tref"] == 0
    path = tmp_path / "features.csv"
    save_features_csv(table, path)
    restored = read_feature_csv(path)
    for col in ("feature_semantics", "observed_window", "hysteresis_metadata", "temperature_correction", "feature_units", "baseline_semantics"):
        assert restored.loc[0, col] == table.loc[0, col]
    assert restored.attrs["provenance"] == table.attrs["provenance"]
    assert restored.loc[0, "VE_observed"] == pytest.approx(.875)