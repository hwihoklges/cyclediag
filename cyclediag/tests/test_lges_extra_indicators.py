"""Tests for extra LGES shape / SOC-band indicators."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cyclediag.features.lges_extra_indicators import (
    capacity_weighted_v_avg,
    correct_r_to_25c,
    dtw_distance,
    extract_shape_indicators,
    fit_rest_tau,
    hysteresis_metrics,
    plateau_metrics,
    rolling_slope,
    safe_ratio,
    soc0_cliff_width,
    vq_norm_curve,
)


def test_capacity_weighted_v_avg():
    q = np.linspace(0, 10, 50)
    v = np.linspace(3.0, 4.0, 50)
    avg = capacity_weighted_v_avg(v, q)
    assert avg is not None
    assert 3.4 < avg < 3.6


def test_hysteresis_positive_area():
    q = np.linspace(0, 1, 100)
    chg_v = 3.0 + q
    dchg_v = 3.9 - q
    out = hysteresis_metrics(q, chg_v, q, dchg_v)
    assert out["hyst_area"] is not None and out["hyst_area"] > 0
    assert out["hyst_max_dV"] == pytest.approx(0.1, abs=0.02)


def test_correct_r_to_25c_hotter_is_lower_raw_maps_up():
    r45 = 10.0
    r25 = correct_r_to_25c(r45, 45.0, ea=20_000.0)
    assert r25 is not None and r25 > r45


def test_extract_shape_has_soc_bands():
    n = 200
    q = np.linspace(0, 100, n)
    chg = pd.DataFrame({
        "step_time": np.arange(n, dtype=float),
        "voltage": 3.0 + 1.1 * (q / 100),
        "capacity": q,
        "charge_capacity": q,
        "current": np.ones(n),
    })
    dchg = pd.DataFrame({
        "step_time": np.arange(n, dtype=float),
        "voltage": 4.1 - 0.9 * (q / 100) - 0.5 * np.exp(-((q - 97) ** 2) / 4),
        "capacity": q,
        "discharge_capacity": q,
        "current": -np.ones(n),
    })
    out = extract_shape_indicators(chg, dchg, q, q)
    for key in (
        "dchg_dVdQ_SOC0", "dchg_dVdQ_SOC5", "dchg_dVdQ_SOC10", "dchg_dVdQ_SOCmid",
        "chg_dVdQ_SOC100", "chg_V_avg", "dchg_V_avg", "hyst_area", "dchg_E",
    ):
        assert key in out
        assert out[key] is not None


def test_plateau_on_flat_midsection():
    n = 300
    q = np.linspace(0, 100, n)
    # steep, flat middle, steep
    v = np.where(q < 30, 4.1 - 0.02 * q, np.where(q < 70, 3.5, 3.5 - 0.02 * (q - 70)))
    out = plateau_metrics(q, v, discharge=True)
    assert out["plateau_V"] is not None
    assert out["plateau_width"] is not None and out["plateau_width"] > 10


def test_safe_ratio():
    assert safe_ratio(10.0, 5.0) == 2.0
    assert safe_ratio(1.0, 0.0) is None


def test_soc0_cliff_width_positive_near_end():
    n = 250
    q = np.linspace(0, 100, n)
    # steep |dV/dQ| only near discharge end (high Q)
    v = 4.1 - 0.008 * q - 0.8 * np.exp(-((q - 98) ** 2) / 2)
    w = soc0_cliff_width(q, v)
    assert w is not None and w >= 0


def test_fit_rest_tau_recovers_order():
    t = np.linspace(0, 600, 120)
    tau_true = 80.0
    v = 4.0 + 0.05 * np.exp(-t / tau_true)
    rest = pd.DataFrame({"voltage": v, "step_time": t})
    tau = fit_rest_tau(rest)
    assert tau is not None
    assert 40 < tau < 160


def test_dtw_identical_is_zero():
    a = np.linspace(3.0, 4.0, 64)
    assert dtw_distance(a, a) == pytest.approx(0.0, abs=1e-9)


def test_vq_norm_and_rolling_slope():
    q = np.linspace(0, 10, 50)
    v = 3.0 + 0.1 * q
    curve = vq_norm_curve(q, v)
    assert curve is not None and len(curve) == 128
    x = np.arange(40, dtype=float)
    y = 2.0 * x + 1.0
    slopes = rolling_slope(y, x, window=10)
    assert np.isfinite(slopes[-1])
    assert slopes[-1] == pytest.approx(2.0, abs=0.05)


def test_extract_shape_new_cols():
    n = 200
    q = np.linspace(0, 100, n)
    chg = pd.DataFrame({
        "voltage": 3.0 + 1.1 * (q / 100),
        "capacity": q,
        "charge_capacity": q,
        "current": np.ones(n),
    })
    dchg = pd.DataFrame({
        "voltage": 4.1 - 0.9 * (q / 100) - 0.5 * np.exp(-((q - 97) ** 2) / 4),
        "capacity": q,
        "discharge_capacity": q,
        "current": -np.ones(n),
    })
    out = extract_shape_indicators(chg, dchg, q, q, dchg_v_cutoff=2.5)
    for key in (
        "dchg_dVdQ_SOC0_cliff_width",
        "dchg_dVdQ_SOC0_to_mid_ratio",
        "dchg_V_cutoff_margin",
    ):
        assert key in out
        assert out[key] is not None
