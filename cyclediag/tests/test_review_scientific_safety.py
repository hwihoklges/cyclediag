"""Synthetic invariants, not validation of unique degradation mechanisms."""

from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cyclediag.io.cycler_csv import ColumnMap, load_cycler_csv, normalize_cycler_dataframe
from cyclediag.features.extract import FeatureConfig, extract_features_table, extract_leg_features
from cyclediag.features.lges_extract import LgesExtractConfig, extract_lges_cycle_row
from cyclediag.features.quality import cycle_quality_metrics
from cyclediag.features.dcir_decompose import fit_r_t_components
from cyclediag.features.units import capacity_to_ah, current_to_a
from cyclediag.features.lges_extra_indicators import capacity_weighted_v_avg, energy_wh
from cyclediag.models.indicator_scoring import score_indicators
from cyclediag.models.peak_ml import PeakMlBundle, PeakMlConfig, predict_peak_model
from cyclediag.diagnosis.engine import diagnose_feature_table
from cyclediag.diagnosis.pattern_scoring import load_mode_weights, score_mode_for_row


def raw_cycle(milli=False):
    n = 20
    scale = 1000 if milli else 1
    return pd.DataFrame({
        "CycleIndex": [1] * (2 * n),
        "StepType": ["charge"] * n + ["discharge"] * n,
        f"Voltage({'mV' if milli else 'V'})": np.r_[np.linspace(3, 4.2, n), np.linspace(4.1, 3, n)] * scale,
        f"Current({'mA' if milli else 'A'})": np.r_[np.ones(n), -np.ones(n)] * scale,
        f"ChargeCapacity({'mAh' if milli else 'Ah'})": np.r_[np.linspace(0, .1, n), np.full(n, .1)] * scale,
        f"DischargeCapacity({'mAh' if milli else 'Ah'})": np.r_[np.zeros(n), np.linspace(0, .095, n)] * scale,
        "TotalTime(sec)": np.arange(2 * n, dtype=float),
    })


def test_units_equivalent_and_normalization_idempotent(tmp_path):
    a = normalize_cycler_dataframe(raw_cycle())
    ma = normalize_cycler_dataframe(raw_cycle(True))
    cols = ["voltage", "current", "charge_capacity", "discharge_capacity", "time"]
    np.testing.assert_allclose(a[cols], ma[cols])
    pd.testing.assert_frame_equal(ma, normalize_cycler_dataframe(ma))
    assert ma.attrs["unit_metadata"]["current"]["source_unit"] == "ma"
    assert ma["current"].iloc[-1] == -1
    path = tmp_path / "raw.csv"
    raw_cycle(True).to_csv(path, index=False)
    pd.testing.assert_frame_equal(load_cycler_csv(str(path)), ma)
    cfg = FeatureConfig(active_mass_g=2)
    fa, fm = [extract_features_table(d, config=cfg).set_index("leg") for d in (a, ma)]
    numeric = ["f_Q_max", "f_Q_spec", "f_V_avg", "f_q_cc_end", "f_cc_Q_frac"]
    np.testing.assert_allclose(fa[numeric].astype(float), fm[numeric].astype(float))
    assert fm.loc["charge", "f_Q_max"] == pytest.approx(.1)
    assert fm.loc["discharge", "f_Q_max"] == pytest.approx(.095)
    assert fm.loc["discharge", "f_q_cc_end"] == pytest.approx(.095)
    assert fm.loc["discharge", "f_Q_spec"] == pytest.approx(47.5)
    assert fm.loc["discharge", "f_V_avg"] == pytest.approx(3.55)


def test_unit_overrides_and_unknowns():
    bare = pd.DataFrame({"Capacity": [72000], "Current": [-1000]})
    assumed = normalize_cycler_dataframe(bare)
    assert assumed.capacity.iloc[0] == 72000  # no magnitude inference
    assert "unit_assumed:charge_capacity:ah" in assumed.attrs["unit_warnings"]
    cmap = ColumnMap(units={"capacity": "mAh", "current": "mA"})
    out = normalize_cycler_dataframe(bare, asdict(cmap))
    assert out.capacity.iloc[0] == 72
    assert out.current.iloc[0] == -1
    with pytest.raises(ValueError, match="Conflicting"):
        normalize_cycler_dataframe(raw_cycle(), cmap)
    with pytest.raises(ValueError, match="Unsupported"):
        normalize_cycler_dataframe(pd.DataFrame({"Voltage(kV)": [3]}))
    with pytest.raises(ValueError):
        normalize_cycler_dataframe(bare, ColumnMap(units={"current": "Ah"}))
    assert capacity_to_ah(float("inf")) is None
    assert current_to_a(float("inf")) is None


@pytest.mark.parametrize("unit,scale", [("min", 60), ("h", 3600), ("ms", .001)])
def test_time_units(unit, scale):
    frame = normalize_cycler_dataframe(pd.DataFrame({f"TotalTime({unit})": [0, 2]}))
    assert frame.time.iloc[-1] == pytest.approx(2 * scale)


def test_lges_units_and_explicit_override_no_double_conversion():
    cfg = LgesExtractConfig(with_diagnosis=False, enrich_assb=False)
    a, ma = [normalize_cycler_dataframe(raw_cycle(m)) for m in (False, True)]
    ra, rm = [extract_lges_cycle_row(d, 1, config=cfg) for d in (a, ma)]
    for key in ("chgCapa", "dchgCapa", "CE", "VE", "chg_E", "dchg_E"):
        assert ra[key] == pytest.approx(rm[key])
    assert rm["CE"] == pytest.approx(95)
    assert 0 < rm["VE"] < 1  # fraction, not percent
    with pytest.raises(ValueError, match="already Ah"):
        extract_lges_cycle_row(ma, 1, config=LgesExtractConfig(capacity_unit="mah"))
    legacy = a.copy()
    legacy.attrs = {}
    for col in ("capacity", "charge_capacity", "discharge_capacity"):
        legacy[col] *= 1000
    legacy_row = extract_lges_cycle_row(legacy, 1, config=LgesExtractConfig(capacity_unit="mah"))
    assert legacy_row["chg_E"] == pytest.approx(ra["chg_E"])


@pytest.mark.parametrize("capacity", [[0, 1, .2, .8], [0, np.nan, 1, 2], [0, 1, np.inf, 2]])
def test_invalid_counters_do_not_generate_capacity_or_weighted_results(capacity):
    seg = pd.DataFrame({"capacity": capacity, "voltage": [3, 3.5, 4, 4.2], "current": [1] * 4})
    row = extract_leg_features(seg, leg="charge", config=FeatureConfig())
    assert row["f_Q_max"] is None
    assert row["f_V_avg"] is None
    assert row["f_cc_Q_frac"] is None
    assert row["capacity_status"] == "invalid_or_reset_counter"


def test_flat_counter_and_nonfinite_voltage_have_no_weighted_mean():
    for q, v in [([1, 1], [3, 4]), ([0, 1], [3, np.inf])]:
        row = extract_leg_features(pd.DataFrame({"capacity": q, "voltage": v}), leg="charge", config=FeatureConfig())
        assert row["f_V_avg"] is None


def test_lges_weighted_energy_reject_resets_and_convert_explicit_mah():
    assert capacity_weighted_v_avg([3, 4, 3.5], [0, 1, .2]) is None
    assert energy_wh([3, 4, 3.5], [0, 1, .2]) is None
    assert energy_wh([3, np.nan, 4], [0, .5, 1]) is None
    assert energy_wh([3, 4], [0, 1]) == pytest.approx(3.5)
    assert energy_wh([3, 4], [0, 1000], q_is_mah=True) == pytest.approx(3.5)
    frame = normalize_cycler_dataframe(raw_cycle())
    frame.loc[35, "discharge_capacity"] = 0
    row = extract_lges_cycle_row(frame, 1)
    assert row["dchgCapa"] is None
    assert row["CE"] is None
    assert row["capacity_status"]["discharge"] == "missing_or_invalid_counter"


def test_lges_cv_does_not_reduce_total_charge_capacity():
    n = 100
    frame = pd.DataFrame({"cycle": [1] * n, "step_type": ["charge"] * n,
                          "voltage": np.r_[np.linspace(3, 4.2, 80), np.full(20, 4.2)],
                          "charge_capacity": np.linspace(0, 1, n),
                          "current": np.r_[np.ones(80), np.linspace(.5, .05, 20)],
                          "time": np.arange(n, dtype=float)})
    row = extract_lges_cycle_row(frame, 1, config=LgesExtractConfig())
    assert row["chgCapa"] == pytest.approx(1)
    assert row["chg_E"] == pytest.approx(energy_wh(frame.voltage, frame.charge_capacity))


def test_constant_reference_is_unknown_without_calibrated_scale():
    ref = pd.DataFrame({"cycle": [1, 2, 3], "VE": [.9] * 3})
    data = pd.DataFrame({"cycle": [4, 5, 6], "VE": [.9, .7, np.nan]})
    res = score_indicators(data, reference=ref)
    out = res.cycle_scores
    assert out.indicator_flag.tolist() == ["unknown"] * 3
    assert out.indicator_n_scored.tolist() == [0] * 3
    assert out.indicator_score.isna().all()
    assert out.indicator_constant_unchanged.tolist() == [1, 0, 0]
    assert out.indicator_constant_changed.tolist() == [0, 1, 0]
    assert res.indicator_summary.indicator_score.isna().all()
    calibrated = score_indicators(data, reference=ref, scale_floors={"VE": .01}).cycle_scores
    assert calibrated.indicator_flag.tolist() == ["ok", "alert", "unknown"]


def test_all_missing_indicators_are_unknown():
    out = score_indicators(pd.DataFrame({"cycle": [1], "VE": [np.nan]})).cycle_scores
    assert out.indicator_flag.iloc[0] == "unknown"
    assert out.indicator_n_scored.iloc[0] == 0
    empty_reference = score_indicators(pd.DataFrame({"cycle": [1, 2, 3], "VE": [.9, .8, .7]}),
                                      reference=pd.DataFrame()).cycle_scores
    assert empty_reference.indicator_n_scored.eq(0).all()
    assert empty_reference.indicator_flag.eq("unknown").all()


def test_implicit_indicator_baselines_separate_cells():
    frame = pd.DataFrame({"cell_id": ["a"] * 3 + ["b"] * 3,
                          "cycle": [1, 2, 3] * 2, "VE": [.9] * 3 + [.7] * 3})
    out = score_indicators(frame).cycle_scores
    assert out.indicator_n_scored.eq(0).all()
    assert out.indicator_constant_unchanged.eq(1).all()


def test_peak_cycle_match_never_suppresses_alert_or_missing(monkeypatch, tmp_path):
    import cyclediag.models.peak_ml as ml
    # Prediction regression independent of stochastic fitting behavior.
    bundle = PeakMlBundle(pipeline=None, feature_columns=["cha_P1_V"], config=PeakMlConfig(),
                          train_cycles=[1], score_thresholds={"watch": .1, "alert": .2})
    bundle.save(tmp_path)
    bundle = PeakMlBundle.load(tmp_path)  # v1 schema stays readable
    monkeypatch.setattr(ml, "_raw_scores", lambda pipeline, x: np.full(len(x), .5))
    data = pd.DataFrame({"cell_id": ["different_cell", "different_cell"],
                         "cycle": [1, 1], "cha_P1_V": [4.2, np.nan]})
    out = predict_peak_model(data, bundle)
    assert out.ml_flag.tolist() == ["alert", "unknown"]
    assert out.ml_is_outlier.tolist() == [True, False]
    assert pd.isna(out.ml_anomaly_score.iloc[1])
    assert out.ml_reference_membership.eq("unknown_legacy_cycle_only").all()
    missing = predict_peak_model(data.drop(columns="cha_P1_V"), bundle)
    assert missing.ml_flag.eq("unknown").all()


CONFIG_DIR = Path(__file__).resolve().parents[1] / "diagnosis" / "config"


def ve_metadata():
    return {"feature_semantics": {"VE": "observed_window_v2"},
            "observed_window": {leg: {"status": "ok"} for leg in ("charge", "discharge")}}


@pytest.mark.parametrize("filename", ["mode_weights_fullcell_v1.json", "mode_weights_assb_si_v1.json"])
def test_ve_loss_monotonic_and_missing_reference_unscorable(filename):
    cfg = load_mode_weights(CONFIG_DIR / filename)
    for mode, mode_cfg in cfg["modes"].items():
        for term in mode_cfg["evidence"]:
            if term["feature"] != "VE":
                continue
            assert term["direction"] == "decrease_vs_baseline"
            assert term["baseline_ref"] == "VE"
            only_ve = {**cfg, "modes": {mode: {"evidence": [term]}}}
            scores = [score_mode_for_row({"VE": v, **ve_metadata()}, mode, only_ve,
                                        baseline_row={"VE": .95, **ve_metadata()}).estimate
                      for v in (.95, .9, .8)]
            assert scores[0] == 0 < scores[1] < scores[2]
            missing = score_mode_for_row({"VE": .8}, mode, only_ve)
            assert missing.estimate is None
            assert missing.evidence_count == 0
            assert missing.status == "unknown"
            assert not missing.diagnosis_valid
            percent_input = score_mode_for_row({"VE": 80}, mode, only_ve, baseline_row={"VE": 95})
            assert percent_input.evidence_count == 0
            assert percent_input.estimate is None


def test_diagnosis_unknown_quality_ce_and_scientific_semantics():
    cfg = {"min_evidence_for_valid": 1, "modes": {"LLI": {"evidence": [
        {"feature": "CE", "direction": "decrease_from_100", "scale": 2}]}}}
    for row in ({}, {"CE": 105}, {"CE": np.inf}):
        result = score_mode_for_row(row, "LLI", cfg)
        assert result.estimate is None
        assert result.evidence_coverage == 0
        assert not result.scientific_validity["validated_probability"]
        assert not result.scientific_validity["causally_identified"]
    invalid = score_mode_for_row({"CE": 95, "quality_score": 0}, "LLI", cfg)
    assert not invalid.diagnosis_valid
    assert invalid.confidence == 0
    row = {"CE": 105}
    result = score_mode_for_row(row, "LLI", cfg)
    assert row["CE"] == 105  # retain measured accounting, don't clip
    assert any("outside_nominal" in w for w in result.scientific_validity["warnings"])


def test_diagnosis_missing_baseline_and_cell_separation():
    table = pd.DataFrame({"cell_id": ["a", "a", "b", "b"], "cycle": [1, 2, 1, 2],
                          "VE": [.95, .8, .8, .8]}, index=[10, 11, 12, 13])
    path = CONFIG_DIR / "mode_weights_fullcell_v1.json"
    for key in ve_metadata():
        table[key] = [ve_metadata()[key] for _ in range(len(table))]
    out = diagnose_feature_table(table, config_path=path)
    assert out.loc[11, "LLI_pattern_score"] > 0
    assert out.loc[13, "LLI_pattern_score"] == 0
    unknown = diagnose_feature_table(table, config_path=path, baseline_cycle=99)
    assert unknown.LLI_evidence_count.eq(0).all()
    assert unknown.diagnosis_status.eq("unknown").all()


def test_quality_missing_tau_and_ecm_metadata():
    out = cycle_quality_metrics(pd.DataFrame({"voltage": [3, 3, 3], "current": [0, 0, 0],
                                             "step_time": [0, 60, 120], "temperature": [0, 0, 0]}))
    assert out["rest_sufficiency"] is None
    assert out["temperature_available"]
    assert out["quality_evidence_coverage"] < 1
    assert "requires_measured_tau" in " ".join(out["scientific_validity"]["warnings"])
    empty = cycle_quality_metrics(pd.DataFrame())
    assert empty["quality_status"] == "unknown"
    assert empty["quality_score"] is None
    fit = fit_r_t_components(np.array([0., 1., 0.]), np.ones(3))
    assert not fit.dcir_fit_valid
    assert not fit.scientific_validity["causally_identified"]


def test_ecm_negative_diffusion_not_reflected_positive():
    t = np.linspace(0, 30, 301)
    fit = fit_r_t_components(t, 2 - .1 * np.sqrt(t), refine_global=False)
    assert fit.A_diff is None
    assert not fit.dcir_fit_valid