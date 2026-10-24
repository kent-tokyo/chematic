"""The same original BGF cases exercised through the built Python extension."""
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / "validation/bgf/cases.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", PACKET["valid"], ids=lambda c: c["id"])
def test_bgf_roundtrip(case):
    for source in (case["input"], case["input"].encode("ascii")):
        record = chematic.parse_bgf(source)
        assert record == case["expected"]
        assert chematic.parse_bgf(chematic.write_bgf(record)) == record


@pytest.mark.parametrize("case", PACKET["invalid"], ids=lambda c: c["id"])
def test_bgf_refusal(case):
    with pytest.raises(ValueError):
        chematic.parse_bgf(case["input"])


def test_bgf_edits_and_unknown_fields():
    record = chematic.parse_bgf(PACKET["valid"][0]["input"])
    record["atoms"][0]["partial_charge"] = -0.98765
    assert chematic.parse_bgf(chematic.write_bgf(record)) == record
    record["atoms"][0]["partial_charge"] = 0.123456
    with pytest.raises(ValueError):
        chematic.write_bgf(record)
    record["atoms"][0]["partial_charge"] = float("nan")
    with pytest.raises(ValueError):
        chematic.write_bgf(record)
    record["atoms"][0]["partial_charge"] = 0.0
    record["atoms"][0]["unrecognized"] = 1
    with pytest.raises(ValueError):
        chematic.write_bgf(record)
