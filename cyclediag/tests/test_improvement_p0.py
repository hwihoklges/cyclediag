"""Tests for roadmap §9.1 new modules."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cyclediag.features.dcir_decompose import fit_r_t_components, _detect_current_settle
from cyclediag.features.self_discharge import fit_self_discharge_rest
from cyclediag.features.signal_cv import detect_cv_signal
from cyclediag.features.units import capacity_to_ah
from cyclediag.features.lges_extra_indicators import hysteresis_metrics
from cyclediag.features.rpt_metrics import compute_q_relax_for_blocks
from cyclediag.features.quality import cycle_quality_metrics
from cyclediag.diagnosis.pattern_scoring import load_mode_weights


def test_capacity_to_ah_no_mah_heuristic():
    # 72 Ah must stay Ah (old heuristic wrongly /1000)
    assert capacity_to_ah(72.0) == pytest.approx(72.0)
    assert capacity_to_ah(72000.0, unit="mah") == pytest.approx(72.0)
    assert capacity_to_ah(72.0, header="ChargeCapacity (Ah)") == pytest.approx(72.0)


def test_fit_r_t_components_recovers_synthetic():
    rng = np.random.default_rng(0)
    t = np.arange(0.0, 30.0 + 1e-9, 0.1)
    r_ohm, r_ct, tau, a = 1.2, 0.8, 2.0, 0.15
    r = r_ohm + r_ct * (1 - np.exp(-t / tau)) + a * np.sqrt(t)
    r = r + rng.normal(0, 0.005, size=t.shape)
    fit = fit_r_t_components(t, r)
    assert fit.R_ohmic == pytest.approx(r_ohm, rel=0.15)
    assert fit.R_ct == pytest.approx(r_ct, rel=0.25)
    assert fit.A_diff == pytest.approx(a, rel=0.35)
    assert fit.dcir_fit_r2 is not None and fit.dcir_fit_r2 > 0.95


def test_self_discharge_late_linear():
    t = np.linspace(0, 3600, 500)
    # mild leak 2 mV/h = 2e-3/3600 V/s
    k = 2e-3 / 3600.0
    v = 3.9 - 0.02 * np.exp(-t / 200.0) - k * t
    out = fit_self_discharge_rest(t, v)
    assert out["self_discharge_rate"] == pytest.approx(2.0, rel=0.5)


def test_signal_cv_detects_taper():
    n_cc, n_cv = 60, 40
    v = np.concatenate([np.linspace(3.5, 4.2, n_cc), np.full(n_cv, 4.2)])
    i = np.concatenate([np.full(n_cc, 38.0), np.linspace(30.0, 2.0, n_cv)])
    q = np.concatenate([np.linspace(0, 60, n_cc), np.linspace(60, 70, n_cv)])
    t = np.arange(n_cc + n_cv) * 30.0
    df = pd.DataFrame({"voltage": v, "current": i, "charge_capacity": q, "step_time": t})
    res = detect_cv_signal(df, column_cv_ah=10.0)
    assert res.has_cv
    assert res.chgCVcapa is not None and res.chgCVcapa > 0


def test_hysteresis_soc_bands():
    q = np.linspace(0, 1, 200)
    vc = 3.5 + 0.5 * q
    vd = 3.9 - 0.5 * q
    h = hysteresis_metrics(q, vc, q, vd)
    assert h["hyst_area"] is not None and h["hyst_area"] > 0
    assert h["hyst_area_low"] is not None
    assert h["hyst_frac_low"] is not None


def test_q_relax_block():
    feats = pd.DataFrame({
        "cycle": [107, 108, 109],
        "dchgCapa": [68.95, 68.90, 14.0],
    })
    out = compute_q_relax_for_blocks(feats, [[107, 108, 109]])
    assert len(out) == 1
    assert out.iloc[0]["Q_relax"] == pytest.approx(-0.05, abs=0.01)


def test_quality_metrics_basic():
    v = np.linspace(4.2, 2.5, 200) + np.random.default_rng(0).normal(0, 0.0005, 200)
    df = pd.DataFrame({
        "voltage": v,
        "current": np.full(200, -38.0),
        "step_time": np.arange(200) * 30.0,
        "temperature": np.zeros(200),
    })
    q = cycle_quality_metrics(df)
    assert q["temperature_available"] is True  # zero degrees C is a valid measurement
    assert q["samples_per_mV"] is not None
    assert q["quality_score"] is not None


def test_assb_mode_weights_load():
    cfg = load_mode_weights()
    assert "contact_loss" in cfg["modes"]
    assert "LAM_NE" not in cfg["modes"]


def test_fit_r_t_components_rejects_current_ramp():
    t = np.arange(0.0, 30.0 + 1e-9, 0.1)
    r_ohm, r_ct, tau, a = 1.2, 0.8, 2.0, 0.15
    r_true = r_ohm + r_ct * (1 - np.exp(-t / tau)) + a * np.sqrt(t)
    i_target = 38.0
    i = np.minimum(i_target, i_target * t / 1.0)
    r_meas = r_true * (i_target / np.maximum(i, 1e-6))
    fit = fit_r_t_components(t, r_meas, i=i)
    assert "current_ramp" in fit.flag
    assert fit.dcir_fit_valid is False


def test_detect_current_settle_flags_slow_ramp():
    t = np.linspace(0, 2, 200)
    i = np.minimum(38.0, 38.0 * t / 1.0)
    t_settle, ok = _detect_current_settle(t, i)
    assert ok is False
    assert t_settle == pytest.approx(1.0, abs=0.05)


def test_ocv_drift_classify():
    from cyclediag.features.ocv_drift import _classify_drift

    mode, par = _classify_drift(
        d80=-0.05, d50=-0.05, d20=-0.05, d_spread_20_80=0.0,
    )
    assert mode == "parallel_shift"
    assert par == pytest.approx(-0.05, abs=0.001)

    mode0, par0 = _classify_drift(d80=0.0, d50=0.0, d20=0.0, d_spread_20_80=0.0)
    assert mode0 == "stable"
    assert par0 == pytest.approx(0.0, abs=1e-9)

    mode_sp, par_sp = _classify_drift(
        d80=0.0, d50=0.0, d20=-0.03, d_spread_20_80=-0.03,
    )
    assert mode_sp in ("spread_change", "local_soc20", "spread_and_shift")
    assert abs(par_sp) < 0.02


def test_band_capacity_high_low_split():
    from cyclediag.features.band_capacity import BandCapacityConfig, discharge_band_capacity

    v = np.linspace(4.0, 2.5, 400)
    q = np.linspace(0, 70, 400)
    seg = pd.DataFrame({"voltage": v, "discharge_capacity": q})
    out = discharge_band_capacity(seg, config=BandCapacityConfig(v_high=3.5, v_low=3.0))
    assert out["dchg_Q_high_frac"] is not None
    assert out["dchg_Q_low_frac"] is not None
    assert out["dchg_Q_high_frac"] + out["dchg_Q_low_frac"] <= 1.0 + 1e-6


def test_cell_meta_per_delta_i():
    from cyclediag.features.cell_meta import CellProtocolMeta

    pm = CellProtocolMeta(q_rated_ah=72.0, routine_c_rate=0.5, rpt_c_rate=1.0 / 3.0)
    assert pm.per_delta_i_a == pytest.approx(12.0, abs=0.5)
    assert pm.dcir_pulse_current_a == pytest.approx(72.0, abs=0.1)
