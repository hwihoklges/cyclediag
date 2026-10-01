"""C5: pulse quality must use an explicit dimensional current and contiguous clocks."""

import pandas as pd
import pytest

from cyclediag.features.quality import cycle_quality_metrics


def quality(current, time=None, step=None, **kwargs):
    data = {"voltage": [3.5] * len(current), "current": current}
    if time is not None:
        data["step_time"] = time
    if step is not None:
        data["StepNo"] = step
    return cycle_quality_metrics(pd.DataFrame(data), **kwargs)


def test_explicit_amperes_and_capacity_times_rate_agree_without_fixed_ampere_threshold():
    values = [0, -0.6, -0.6, -0.6]
    a = quality(values, [0, 1, 1.5, 2], pulse_current_threshold_a=0.5)
    c = quality(values, [0, 1, 1.5, 2], nominal_capacity_ah=0.5, pulse_c_rate=1)
    assert a["pulse_sample_count_1s"] == c["pulse_sample_count_1s"] == 3
    assert a["pulse_threshold_source"] == "explicit_A"
    assert c["pulse_threshold_source"] == "nominal_Ah_times_C_rate"
    assert a["pulse_current_stability"] == 0
    assert quality(values, [0, 1, 1.5, 2])["pulse_sample_count_1s"] is None
    assert quality([0, -30, -30], [0, 1, 2])["pulse_current_stability"] is None


@pytest.mark.parametrize("kwargs", [
    {"pulse_current_threshold_a": -1},
    {"pulse_current_threshold_a": float("nan")},
    {"nominal_capacity_ah": 0, "pulse_c_rate": 1},
])
def test_invalid_threshold_fails_closed(kwargs):
    with pytest.raises(ValueError):
        quality([-1, -1], [0, 1], **kwargs)


def test_first_contiguous_pulse_only_and_step_boundary():
    q = quality([0, -1, -1.1, 0, -5, -5, 0, 0],
                [0, 0, .5, 1, 0, .5, 0, 100],
                [1, 2, 2, 2, 3, 3, 4, 5], pulse_current_threshold_a=.5,
                tau_relax_est=1)
    assert q["pulse_sample_count_1s"] == 2
    assert q["pulse_current_stability"] == pytest.approx(.05 / 1.05)
    assert q["rest_sufficiency"] == 0


def test_missing_time_nan_and_clock_reset_cannot_join_pulses_or_rest():
    assert quality([1, 1, 1], pulse_current_threshold_a=.5)["pulse_sample_count_1s"] is None
    q = quality([-1, -1, -1, -1], [0, .2, float("nan"), .3], pulse_current_threshold_a=.5)
    assert q["pulse_sample_count_1s"] == 2
    assert q["pulse_current_stability"] == 0
    reset = quality([-1, -1, -1, -1], [0, .2, 0, .3], pulse_current_threshold_a=.5)
    assert reset["pulse_sample_count_1s"] == 2
    assert reset["pulse_current_stability"] == 0
    direction = quality([-1, -1, 1, 1], [0, .2, .4, .6], pulse_current_threshold_a=.5)
    assert direction["pulse_sample_count_1s"] == 2
