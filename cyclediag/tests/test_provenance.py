"""Bounded reproducibility contracts; synthetic inputs, no LFS fixtures."""
from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import pytest

from cyclediag import api
from cyclediag.features.lges_extract import LgesExtractConfig
from cyclediag.io.cycler_csv import ColumnMap, normalize_cycler_dataframe
from cyclediag.provenance import (
    build_provenance, canonical_json, code_metadata, config_hash, dataframe_hash,
    file_hash, read_feature_csv, save_features_csv,
)


def raw():
    return pd.DataFrame({"CycleIndex": [1, 1, 2, 2], "Voltage(V)": [3., 4., 3., 4.],
                         "ChargeCapacity": [0., 1., 0., 1.], "Current(A)": [1.] * 4,
                         "TotalTime_sec": [0., 1., 2., 3.], "StepType": ["charge"] * 4})


def test_canonical_config_and_defaults():
    cfg = LgesExtractConfig()
    before = deepcopy(cfg)
    assert config_hash(cfg) == config_hash(asdict(cfg))
    assert config_hash({"b": 2, "a": 1}) == config_hash({"a": 1, "b": 2})
    assert config_hash(ColumnMap()) != config_hash(ColumnMap(units={"capacity": "mAh"}))
    assert cfg == before
    with pytest.raises(TypeError):
        canonical_json(object())


def test_input_units_schema_and_selectors(tmp_path, monkeypatch):
    monkeypatch.setenv("CYCLEDIAG_SECRET", "do-not-serialize-me")
    df = normalize_cycler_dataframe(raw())
    cfg = LgesExtractConfig(cell_id="sensitive-cell")
    def meta(source=df, **kwargs):
        return build_provenance(source, df, config=cfg, column_map=ColumnMap(), **kwargs)
    a = meta(cycles=[1])
    assert a["config_sha256"] != meta(cycles=[2])["config_sha256"]
    assert a["config_sha256"] != meta(cycles=[1], selectors={"cv_only": True})["config_sha256"]
    changed = df.copy()
    changed.loc[0, "voltage"] += .1
    assert dataframe_hash(df) != dataframe_hash(changed)
    changed = df.copy(deep=True)
    changed.attrs = deepcopy(df.attrs)
    changed.attrs["unit_schema"] = "other"
    assert dataframe_hash(df) != dataframe_hash(changed)
    assert dataframe_hash(df) != dataframe_hash(df.astype({"cycle": float}))
    path = tmp_path / "sensitive-cell.csv"
    raw().to_csv(path, index=False)
    b = meta(path)
    assert b["input"] == {"kind": "file_bytes", "sha256": file_hash(path)}
    path.write_bytes(path.read_bytes() + b"\n")
    assert b["input"]["sha256"] != meta(path)["input"]["sha256"]
    text = canonical_json(a)
    assert "sensitive-cell" not in text and "do-not-serialize-me" not in text
    assert str(tmp_path) not in canonical_json(b)


def test_api_does_not_mutate_config_and_explicit_result(tmp_path, monkeypatch):
    path = tmp_path / "cell.csv"
    raw().to_csv(path, index=False)
    cfg = LgesExtractConfig(with_diagnosis=False, enrich_assb=False)
    before = deepcopy(cfg)
    observed = []
    def extract(df, **kwargs):
        observed.append(kwargs["cycles"])
        kwargs["config"].rest_labels = "mutated-by-downstream"
        return pd.DataFrame({"cycle": [1], "scientific_warnings": [["incomplete"]]})
    monkeypatch.setattr(api, "extract_lges_features_table", extract)
    feats = api.extract_features(path, cycles=iter([1]), config=cfg)
    assert cfg == before and observed == [[1]]
    assert feats.attrs["unit_warnings"]
    monkeypatch.setattr(api, "diagnose_dataframe", lambda features, **kw: {"features": features})
    yes = api.diagnose_csv(path, config=cfg)
    no = api.diagnose_csv(path, config=cfg, with_screen=False)
    assert yes["provenance"]["config_sha256"] != no["provenance"]["config_sha256"]
    assert yes["provenance"]["extraction_config_sha256"] == yes["features"].attrs["provenance"]["config_sha256"]


def test_checked_roundtrip_and_warning_preservation(tmp_path):
    frame = pd.DataFrame({"cycle": [1, 2], "value": [np.nan, 2.],
                          "scientific_warnings": [["missing_time"], []],
                          "feature_units": [{"value": "Ah"}, {"value": "Ah"}]})
    frame.attrs = {"unit_schema": "canonical_A_Ah_V_s_v1", "unit_warnings": ["assumed"],
                   "scientific_warnings": ["not_validated"], "secret": "must-not-export",
                   "provenance": {"schema_version": 1}}
    original = deepcopy(frame)
    path = tmp_path / "features.csv"
    sidecar = save_features_csv(frame, path)
    restored = read_feature_csv(path)
    assert restored.scientific_warnings.tolist() == frame.scientific_warnings.tolist()
    assert restored.feature_units.tolist() == frame.feature_units.tolist()
    assert restored.attrs["scientific_warnings"] == ["not_validated"]
    assert restored.attrs["unit_warnings"] == ["assumed"]
    assert "secret" not in sidecar.read_text()
    pd.testing.assert_frame_equal(frame, original)
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="integrity"):
        read_feature_csv(path)


def test_export_source_protection_and_empty(tmp_path):
    path = tmp_path / "raw.csv"
    path.write_text("original")
    frame = pd.DataFrame({"file": [str(path)]})
    with pytest.raises(ValueError, match="source"):
        save_features_csv(frame, path)
    side = path.with_name(path.name + ".metadata.json")
    side.write_text("original-side")
    with pytest.raises(ValueError, match="source"):
        save_features_csv(pd.DataFrame(), path, source_path=side)
    assert path.read_text() == "original"
    empty = tmp_path / "empty.csv"
    save_features_csv(pd.DataFrame(), empty)
    assert read_feature_csv(empty).empty


def test_structured_nulls_and_git_timeout(tmp_path, monkeypatch):
    path = tmp_path / "warnings.csv"
    save_features_csv(pd.DataFrame({"warnings": [["missing"], None, []]}), path)
    assert read_feature_csv(path).warnings.tolist() == [["missing"], None, []]
    def timeout(*args, **kwargs):
        assert kwargs["shell"] is False and kwargs["timeout"] == 3
        raise subprocess.TimeoutExpired("git", 3)
    monkeypatch.setattr(subprocess, "run", timeout)
    assert code_metadata(tmp_path)["commit"] is None


def test_real_checkout_code_ownership():
    meta = code_metadata()
    assert meta["commit"] is not None
    assert isinstance(meta["dirty"], bool)


def test_unrelated_parent_git_not_attributed(tmp_path):
    subprocess.run(["git", "init", str(tmp_path)], check=True, capture_output=True)
    package = tmp_path / "installed" / "cyclediag"
    package.mkdir(parents=True)
    for name in ("__init__.py", "api.py"):
        (package / name).write_text("# installed code")
    meta = code_metadata(package)
    assert meta["commit"] is None and meta["dirty"] is None
    original = meta["source_sha256"]
    (package / "api.py").write_text("# changed code")
    assert original != code_metadata(package)["source_sha256"]


@pytest.mark.parametrize("feature_set,cv", [("vp_v1_basic", False), ("vp_lges_cycle_v1", False),
                                           ("vp_lges_cycle_v2", False), ("vp_v1_basic", True)])
def test_cli_extract_sidecars(tmp_path, feature_set, cv):
    from cyclediag.__main__ import main
    source = tmp_path / "raw.csv"
    target = tmp_path / "features.csv"
    raw().to_csv(source, index=False)
    args = ["extract", "--input", str(source), "--out", str(target), "--feature-set", feature_set]
    assert main(args + (["--cv-only"] if cv else [])) == 0
    result = read_feature_csv(target)
    assert result.attrs["provenance"]["input"]["sha256"] == file_hash(source)
    assert "unit_warnings" in result.attrs