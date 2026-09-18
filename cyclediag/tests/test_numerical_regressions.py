"""Numerical regressions without relaxing scientific acceptance criteria."""

import numpy as np
import pandas as pd
import pytest

from cyclediag.analysis.indicator_screen import compare_cells, screen_indicators
from cyclediag.features._statistics import finite_pearson
from cyclediag.features.cliff_metrics import compute_cliff_metrics
from cyclediag.features.dqdv_peaks import DqdvPeakConfig, compute_dvdq
from cyclediag.models.indicator_scoring import score_indicators
from cyclediag.tests.test_cliff_metrics import _synth_curve


@pytest.mark.parametrize("constant", [0.0, 0.7, 0.9])
def test_constant_correlation_is_undefined_on_either_side(constant):
    flat = pd.Series([constant] * 3)
    varying = pd.Series([1.0, 2.0, 3.0])
    assert np.isnan(finite_pearson(flat, varying))
    assert np.isnan(finite_pearson(varying, flat))


def test_correlation_checks_only_aligned_finite_pairs():
    x = pd.Series([0.7, 0.7, 0.8, np.inf], index=[1, 2, 3, 4])
    y = pd.Series([2.0, 1.0, np.nan, 4.0], index=[2, 1, 3, 4])
    assert np.isnan(finite_pearson(x, y))
    assert finite_pearson(pd.Series([1., 2., 3.]), pd.Series([6., 4., 2.])) == pytest.approx(-1)


def test_small_real_reference_variation_is_not_erased():
    ref = pd.DataFrame({"cycle": [1, 2, 3], "VE": [.7 - 1e-10, .7, .7 + 1e-10]})
    out = score_indicators(ref, reference=ref).cycle_scores
    assert out.indicator_n_scored.eq(1).all()
    assert out.indicator_constant_unchanged.eq(0).all()


@pytest.mark.parametrize("n_interp", [250, 500, 1000])
def test_dvdq_retains_shallow_slope_independent_of_voltage_step(n_interp):
    q = np.linspace(0, 56, 1200)
    v = 4.1 - .00875 * q
    qx, dvdq = compute_dvdq(q, v, DqdvPeakConfig(n_interp=n_interp))
    assert len(qx) == n_interp - 1
    assert qx.min() < .01 * q.max()
    assert qx.max() == pytest.approx(q.max())
    np.testing.assert_allclose(dvdq, -.00875, rtol=1e-9, atol=1e-12)


def test_cliff_tracks_transition_not_numerical_noise():
    # Stronger localization assertion complements, not replaces, the CV test.
    for c_si in (32, 24, 16, 8):
        result = compute_cliff_metrics(*_synth_curve(40, c_si))
        assert result["cliff_valid"]
        assert result["Q_cliff_abs"] == pytest.approx(40, abs=.2)


def test_constant_health_and_spread_are_not_correlation_evidence():
    frame = pd.DataFrame({"cycle": range(1, 11), "VE": np.linspace(.9, .8, 10),
                          "SoHQ": [.7] * 10})
    screen = screen_indicators(frame).set_index("feature")
    assert pd.isna(screen.loc["VE", "corr_health"])
    assert "vs SoHQ" not in screen.loc["VE", "signal"]
    compared = compare_cells(pd.concat([
        frame.assign(cell_id="a"), frame.assign(cell_id="b"),
    ], ignore_index=True), late_frac=.5)
    assert not compared.empty
    assert compared.spread_trend_r.isna().all()