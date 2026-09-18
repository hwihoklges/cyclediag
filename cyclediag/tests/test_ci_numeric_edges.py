"""Undefined statistics must stay unknown without emitting runtime warnings."""

import warnings

import numpy as np
import pandas as pd
import pytest

from cyclediag.features.cycle_indicators_plots import fit_sohq_from_rest_v_end
from cyclediag.models.indicator_scoring import _summarize_indicators


@pytest.mark.parametrize("constant_target", [False, True])
def test_proxy_constant_series_has_undefined_correlation(constant_target):
    frame = pd.DataFrame({
        "EoC_restV_end": [4.0] * 4,
        "EoD_restV_end": [3.0] * 4,
        "SoHQ": [90.0] * 4 if constant_target else [100.0, 95.0, 90.0, 85.0],
    })
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        fit, scored = fit_sohq_from_rest_v_end(frame)
    assert np.isnan(fit.pearson_r)
    assert np.isfinite(scored["SoHQ_hat"]).all()


def test_proxy_varying_series_retains_correlation():
    frame = pd.DataFrame({
        "EoC_restV_end": [4.0, 4.1, 4.2, 4.3],
        "EoD_restV_end": [3.0] * 4,
        "SoHQ": [85.0, 90.0, 95.0, 100.0],
    })
    fit, _ = fit_sohq_from_rest_v_end(frame)
    assert fit.pearson_r == pytest.approx(1.0)


@pytest.mark.parametrize("early_z", [np.nan, 1.0])
def test_summary_missing_late_scores_remain_unknown(early_z):
    frame = pd.DataFrame({"cycle": [1, 2, 3], "VE": [0.9, 0.8, 0.7]})
    z = pd.DataFrame({"VE": [early_z, np.nan, np.nan]})
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        summary = _summarize_indicators(frame, ["VE"], z)
    assert np.isnan(summary.loc[0, "late_median_abs_z"])
    if np.isnan(early_z):
        assert summary.loc[0, "status"] == "unknown"
        assert np.isnan(summary.loc[0, "indicator_score"])
    else:
        assert summary.loc[0, "median_abs_z"] == early_z