"""Tests for peak_evolution module."""

from __future__ import annotations

import numpy as np
import pandas as pd

from cyclediag.features.peak_evolution import (
    PeakEvolutionConfig,
    build_evolution_map,
    evaluate_synthetic_mot,
    extract_tracks,
    make_synthetic_map,
    preflight_passes,
    run_preflight_checks,
    track_peaks,
    viterbi_ridge,
)
from cyclediag.features.dqdv_peaks import DqdvPeakConfig


def _synthetic_raw(n_cycles: int = 30) -> pd.DataFrame:
    rows = []
    for cyc in range(1, n_cycles + 1):
        q = np.linspace(0, 70, 120)
        v = 4.2 - 1.2 * (q / 70) + 0.05 * np.sin(6 * np.pi * q / 70 + cyc * 0.02)
        for qi, vi in zip(q, v):
            rows.append({
                "cycle": cyc,
                "voltage": vi,
                "current": -38.0,
                "discharge_capacity": qi,
                "step_type": "discharge",
            })
    return pd.DataFrame(rows)


def test_viterbi_single_ridge_synthetic():
    emap = make_synthetic_map(n_cycles=60, n_grid=150, seed=1)
    lam = 0.5
    path, score, conf = viterbi_ridge(emap.M, lam)
    assert len(path) == emap.n_cycles
    assert score > 0
    assert conf.min() >= 0


def test_extract_multiple_tracks():
    emap = make_synthetic_map(n_cycles=80, n_grid=200, seed=2)
    tracks = extract_tracks(emap, lam=0.3, max_tracks=4, snr_stop=0.75)
    assert 1 <= len(tracks) <= 4
    assert all("path" in t for t in tracks)


def test_build_evolution_map_from_raw():
    raw = _synthetic_raw(15)
    cfg = PeakEvolutionConfig(
        quality_gate=False,
        normalize="none",
        n_grid=200,
        use_roi_extract=False,
        dqdv_config=DqdvPeakConfig(n_interp=200, sg_window=7, merge_v_sep_v=0.003),
    )
    emap = build_evolution_map(raw, config=cfg, rate="0.5C", step_df=None)
    assert emap.M.shape[0] == len(emap.cycles)
    assert emap.n_grid == cfg.n_grid
    assert np.isfinite(emap.M).any()


def test_track_peaks_pipeline_synthetic():
    emap = make_synthetic_map(n_cycles=50, n_grid=120, seed=3)
    cfg = PeakEvolutionConfig(lam=0.4, max_tracks=3, use_roi_extract=False, normalize="none")
    result = track_peaks(emap, config=cfg)
    assert not result.tracks.empty
    assert not result.trajectories.empty
    assert "track_id" in result.tracks.columns
    # initial confirmation should not be labeled "born"
    if not result.events.empty:
        assert set(result.events["event_type"]).issubset(
            {"initial", "born", "died", "merged", "split", "resumed", "ambiguous"}
        )


def test_synthetic_mot_metrics():
    emap = make_synthetic_map(n_cycles=40, n_grid=100, seed=4)
    lam = 0.5
    path, _, _ = viterbi_ridge(emap.M, lam)
    true_paths = [path]
    pred_paths = [path]
    metrics = evaluate_synthetic_mot(
        true_paths, pred_paths,
        true_events={(10, "born")},
        pred_events={(10, "born")},
    )
    assert metrics["id_switch_rate"] == 0.0
    assert metrics["mostly_tracked_ratio"] >= 0.9


def test_preflight_runs_on_synthetic():
    raw = _synthetic_raw(10)
    cfg = PeakEvolutionConfig(
        dqdv_config=DqdvPeakConfig(n_interp=200, sg_window=7),
        n_grid=200,
    )
    table = run_preflight_checks(raw, step_df=None, config=cfg)
    checks = table["check_id"].astype(str)
    assert checks.is_unique
    assert set(checks) == {"1", "2", "3", "4", "4b", "5"}
    items = table.assign(check_id=checks).set_index("check_id")["item"]
    assert items["4"] == "RPT_C3_resolved_peaks_full_dqdv"
    assert items["4b"] == "RPT_C3_resolved_peaks_interior_dvdq"
    assert "status" in table.columns


def test_assign_mode_evolution_raises():
    from cyclediag.features.peak_assign import assign_peaks_for_leg
    import pytest

    v = np.linspace(3.0, 4.1, 50)
    q = np.linspace(0, 60, 50)
    with pytest.raises(ValueError, match="multi-cycle"):
        assign_peaks_for_leg(v, q, "discharge", assign_mode="evolution")
