"""C4: evidence support is not scientific eligibility or causal validation."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cyclediag.diagnosis.engine import diagnose_feature_table
from cyclediag.diagnosis.pattern_scoring import load_mode_weights, score_mode_for_row
from cyclediag.features.diagnosis_export import _write_degradation_diagnosis_outputs

CONFIG = Path(__file__).resolve().parents[1] / "diagnosis" / "config"


def row_for(chemistry="Liion_NMC_Gr"):
        return {"chemistry": chemistry, "cycle": 2, "protocol_kind": "routine", "protocol_excluded": False,
            "protocol_comparable": True, "protocol_id": "p1", "quality_score": 0.9,
            "quality_evidence_coverage": 0.7, "quality_status": "partial",
            "quality_gate_failed_groups": "rest_sufficiency", "leg_completeness": 1.,
            "feature_semantics": {"schema": "lges_measurements_v2"},
            "samples_per_mV": 1., "dqdv_snr": 20., "pulse_sample_count_1s": 8,
            "pulse_current_stability": .005, "rest_sufficiency": None}


def baseline(chemistry="Liion_NMC_Gr"):
        return {"chemistry": chemistry, "cycle": 1, "protocol_kind": "routine", "protocol_excluded": False,
            "protocol_comparable": True, "protocol_id": "p1", "feature_semantics": {"schema": "lges_measurements_v2"}}


def with_baseline(row):
    row.update(feature_semantics={"schema": "lges_measurements_v2"},
               baseline_semantics={"cycle": 1, "requested_cycle": 1,
                                   "feature_semantics": {"schema": "lges_measurements_v2"}})
    return row


def test_high_support_cannot_fill_required_group_and_irrelevant_failure_does_not_block():
    cfg = load_mode_weights(CONFIG / "mode_weights_fullcell_v1.json")
    row = with_baseline(row_for())
    row.update(delta_dchg_V_cutoff_margin=-.5, delta_EoD_restV_end=.5,
               ocv_parallel_shift=.5, delta_EoC_restV_end=-.5, CE=95)
    missing = score_mode_for_row(row, "LLI", cfg, baseline_row=baseline())
    assert missing.estimate > .5 and missing.heuristic_support_valid
    assert not missing.diagnosis_valid
    assert "required_group_missing:independent_shape_or_efficiency" in missing.eligibility_reasons
    row["delta_dchg_dVdQ_SOC0"] = 8.
    good = score_mode_for_row(row, "LLI", cfg, baseline_row=baseline())
    assert good.scientific_eligible and good.diagnosis_valid
    assert not good.scientific_validity["validated_probability"]
    assert not good.scientific_validity["causally_identified"]
    row["quality_gate_failed_groups"] = "samples_per_mV"
    assert "required_quality_unverified:independent_shape_or_efficiency:samples_per_mV" in score_mode_for_row(
        row, "LLI", cfg, baseline_row=baseline()).eligibility_reasons


@pytest.mark.parametrize("change,reason", [
    ({"chemistry": "ASSB_SJ900_Si_rich"}, "chemistry_mismatch"),
    ({"chemistry": None}, "chemistry_unverified"),
    ({"protocol_comparable": "False"}, "protocol_unverified_or_incompatible"),
    ({"protocol_comparable": np.nan}, "protocol_unverified_or_incompatible"),
    ({"protocol_excluded": "False"}, "protocol_unverified_or_incompatible"),
    ({"quality_status": "unknown"}, "quality_unverified"),
    ({"quality_status": "failed"}, "quality_unverified"),
    ({"quality_score": np.nan}, "quality_unverified"),
    ({"quality_evidence_coverage": np.nan}, "quality_unverified"),
    ({"feature_semantics": {}}, "feature_semantics_unverified"),
    ({"leg_completeness": None}, "required_quality_unverified:plateau:leg_completeness"),
    ({"baseline_semantics": {"cycle": 1, "requested_cycle": 99}}, "baseline_unverified"),
    ({"protocol_id": "other"}, "baseline_protocol_mismatch_or_unknown"),
])
def test_bad_metadata_cannot_pass(change, reason):
    cfg = load_mode_weights(CONFIG / "mode_weights_fullcell_v1.json")
    row = with_baseline(row_for())
    row.update(delta_dchg_plateau_V=.1, dchg_plateau_width=2.,
               delta_chg_dQdV_peak1_V=.1)
    ref = baseline()
    ref["dchg_plateau_width"] = 4.
    assert score_mode_for_row(row, "LAM_PE", cfg, baseline_row=ref).scientific_eligible
    row.update(change)
    result = score_mode_for_row(row, "LAM_PE", cfg, baseline_row=ref)
    assert not result.diagnosis_valid and reason in result.eligibility_reasons


@pytest.mark.parametrize("bad", ["False", "True", 1, np.nan, None])
def test_assb_fit_boolean_must_be_explicit_true(bad):
    cfg = load_mode_weights(CONFIG / "mode_weights_assb_si_v1.json")
    row = row_for(cfg["chemistry"])
    row.update(quality_gate_failed_groups="", R_ct_soc50=2., tau_ct_soc50=4.,
               EoC_dchgR_60s_inc=80., dcir_fit_valid_soc50=bad)
    result = score_mode_for_row(row, "interface_R", cfg)
    assert not result.scientific_eligible
    assert "fit_unverified:interface_fit" in result.eligibility_reasons
    row["dcir_fit_valid_soc50"] = True
    assert score_mode_for_row(row, "interface_R", cfg).scientific_eligible


def test_legacy_exploratory_sidecar_and_csv_retain_reasons(tmp_path):
    row = {"cycle": 1, "chemistry": "Liion_NMC_Gr", "delta_dchg_V_cutoff_margin": -.5,
           "delta_EoD_restV_end": .5, "delta_dchg_dVdQ_SOC0": 8., "CE": 95}
    cfg = load_mode_weights(CONFIG / "mode_weights_fullcell_v1.json")
    result = score_mode_for_row(row, "LLI", cfg)
    assert result.heuristic_support_valid and not result.diagnosis_valid
    assert "protocol_unverified_or_incompatible" in result.eligibility_reasons
    sidecar = tmp_path / "modes.json"
    out = diagnose_feature_table(pd.DataFrame([row]), config_path=CONFIG / "mode_weights_fullcell_v1.json",
                                 write_json_sidecar=sidecar)
    assert out.loc[0, "diagnosis_heuristic_support_valid"]
    assert not out.loc[0, "diagnosis_valid"]
    assert not out.loc[0, "LLI_scientific_eligible"]
    assert out.loc[0, "diagnosis_status"] == "heuristic_unverified"
    records = json.loads(sidecar.read_text(encoding="utf8"))
    assert records[0]["scientific_validity"]["eligibility"]["reasons"]
    _, csv, _ = _write_degradation_diagnosis_outputs(out, tmp_path, "sample")
    assert "LLI_eligibility_reasons" in pd.read_csv(csv).columns


def test_configured_modes_have_nonempty_declarations():
    for name in ("mode_weights_fullcell_v1.json", "mode_weights_assb_si_v1.json"):
        cfg = load_mode_weights(CONFIG / name)
        for mode, spec in cfg["modes"].items():
            terms = {t["feature"] for t in spec["evidence"]}
            groups = spec["eligibility"]["required_groups"]
            assert len(groups) >= 2, mode
            assert all(set(g["any_of"]) <= terms for g in groups)