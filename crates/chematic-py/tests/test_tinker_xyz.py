"""Tinker's declared data survives public Python editing and serialization."""
import copy
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / "validation/tinker_xyz/cases.json").read_text(encoding="utf-8"))

@pytest.mark.parametrize("case", PACKET["cases"], ids=lambda c: c["id"])
def test_record_roundtrip(case):
    record = chematic.parse_tinker_xyz(case["input"].encode())
    assert record == case["expected"]
    assert chematic.parse_tinker_xyz(chematic.write_tinker_xyz(record)) == record

@pytest.mark.parametrize("case", PACKET["refusals"], ids=lambda c: c["id"])
def test_refusal(case):
    with pytest.raises(ValueError):
        chematic.parse_tinker_xyz(case["input"])

def test_edit_and_reject_unrepresentable_metadata():
    record = copy.deepcopy(PACKET["cases"][0]["expected"])
    record["atoms"][0].update(name="OW", atom_type=90001, position_angstrom=[.12345678901234567, 0, 0])
    record["cell"] = dict(lengths_angstrom=[10, 11, 12], angles_degrees=[70, 80, 90])
    assert chematic.parse_tinker_xyz(chematic.write_tinker_xyz(record)) == record
    record["atoms"][0]["isotope"] = 18
    with pytest.raises(ValueError):
        chematic.write_tinker_xyz(record)
    del record["atoms"][0]["isotope"]
    record["atoms"][0]["position_angstrom"][0] = float("nan")
    with pytest.raises(ValueError):
        chematic.write_tinker_xyz(record)
