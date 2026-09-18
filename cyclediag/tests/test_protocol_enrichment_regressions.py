import numpy as np
import pandas as pd
import pytest

from cyclediag.io.cycler_csv import normalize_cycler_dataframe
from cyclediag.io.cycle_protocol import _dchg_ah, build_protocol_exclusion
from cyclediag.features.enrich_assb import _set_cell, enrich_feature_table
from cyclediag.diagnosis.constraints import scientific_validity


@pytest.mark.parametrize("header,value", [("DischargeCapacity(mAh)", 300000.), ("DischargeCapacity(Ah)", 300.), ("discharge_capacity", 300.)])
def test_large_canonical_capacity_stays_ah(header, value):
    frame = normalize_cycler_dataframe(pd.DataFrame({
        "CycleIndex": [1] * 4, header: [value] * 4,
        "Current(A)": [38., 0., -38., 0.],
    }))
    assert _dchg_ah(frame) == 300.
    protocol = build_protocol_exclusion(frame)
    assert protocol.flags.iloc[0]["dchg_Ah"] == 300.
    assert protocol.flags.iloc[0]["protocol_kind"] == "routine"
    assert protocol.excluded == set()
    assert _dchg_ah(pd.DataFrame({"discharge_capacity": [300.]})) == 300.


def test_structured_cells_preserve_nondefault_index_and_row_independence():
    out = pd.DataFrame({"cycle": [1, 2, 3], "payload": [np.nan] * 3}, index=[10, 20, 30])
    _set_cell(out, out.cycle < 3, "payload", {"warnings": ["a"]})
    _set_cell(out, out.cycle == 3, "items", [1, 2])
    assert out.at[10, "payload"] == {"warnings": ["a"]}
    out.at[10, "payload"]["warnings"].append("b")
    assert out.at[20, "payload"] == {"warnings": ["a"]}
    assert out.at[30, "items"] == [1, 2]


@pytest.mark.parametrize("cycles", [[1], [1, 2, 3]])
def test_enrichment_preserves_quality_and_ecm_validity(cycles):
    raw = pd.concat([pd.DataFrame({
        "cycle": c, "current": [0., -100., -100., 0.],
        "voltage": [4., 3.9, 3.8, 3.95], "step_time": [0., 1., 2., 3.],
        "discharge_capacity": [0., 0.1, 0.2, 0.2],
    }) for c in cycles], ignore_index=True)
    features = pd.DataFrame({"cycle": cycles, "dchgCapa": [70.] * len(cycles)}, index=[10 * c for c in cycles])
    features["scientific_validity"] = [{"existing": c} for c in cycles]
    out, _ = enrich_feature_table(features, raw, expected_pulse_current=100.)
    for i, (_, row) in enumerate(out.iterrows()):
        quality = row["quality_scientific_validity"]
        key = "ecm_scientific_validity" + (f"_soc{[80, 50, 20][i]}" if len(cycles) == 3 else "")
        ecm = row[key]
        assert isinstance(quality, dict) and quality["warnings"]
        assert isinstance(ecm, dict) and ecm["warnings"]
        assert row["scientific_validity"] == {"existing": cycles[i]}
        validity = scientific_validity(row.to_dict(), {})
        assert set(quality["warnings"] + ecm["warnings"]) <= set(validity["warnings"])
    records = out.to_dict(orient="records")
    assert isinstance(records[0]["quality_scientific_validity"], dict)