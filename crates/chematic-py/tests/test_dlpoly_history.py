"""The same original DL-POLY HISTORY cases exercised through the built Python extension."""
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / "validation/dlpoly_history/cases.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", PACKET["valid"], ids=lambda c: c["id"])
def test_dlpoly_history_roundtrip(case):
    for source in (case["input"], case["input"].encode("ascii")):
        record = chematic.parse_dlpoly_history(source)
        assert record == case["expected"]
        assert chematic.parse_dlpoly_history(chematic.write_dlpoly_history(record)) == record


@pytest.mark.parametrize("case", PACKET["invalid"], ids=lambda c: c["id"])
def test_dlpoly_history_refusal(case):
    with pytest.raises(ValueError):
        chematic.parse_dlpoly_history(case["input"])


def test_dlpoly_history_edits_and_unknown_fields():
    record = chematic.parse_dlpoly_history(PACKET["valid"][0]["input"])
    record["frames"][0]["atoms"][0]["charge"] = -0.98765
    assert chematic.parse_dlpoly_history(chematic.write_dlpoly_history(record)) == record
    record["frames"][0]["atoms"][0]["charge"] = 0.123456789012345
    with pytest.raises(ValueError):
        chematic.write_dlpoly_history(record)
    record["frames"][0]["atoms"][0]["charge"] = float("nan")
    with pytest.raises(ValueError):
        chematic.write_dlpoly_history(record)
    record["frames"][0]["atoms"][0]["charge"] = 0.0
    record["frames"][0]["atoms"][0]["unrecognized"] = 1
    with pytest.raises(ValueError):
        chematic.write_dlpoly_history(record)
