"""Retain parsed binary64 numeric values through the public mmCIF writer."""
import json
from pathlib import Path
import struct

import chematic
import pytest

CASES = json.loads((Path(__file__).parents[3] /
                    "validation/mmcif_numeric_preservation/cases.json").read_text())["cases"]


def bits(values):
    return [struct.pack(">d", float(value)).hex() for value in values]


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_mmcif_numeric_record_round_trip(case):
    before = chematic.parse_mmcif(case["input"])
    output = chematic.write_mmcif(before["atoms"], before["cell"], before["space_group"], case["id"])
    after = chematic.parse_mmcif(output)
    assert before["atoms"] == after["atoms"] == case["expected"]["atoms"]
    assert before["space_group"] == after["space_group"]
    assert after["cell"] == case["expected"]["cell"]
    for atom, expected in zip(after["atoms"], case["expected"]["atoms"], strict=True):
        assert bits(atom[key] for key in ("x", "y", "z", "occupancy", "b_iso")) == bits(
            expected[key] for key in ("x", "y", "z", "occupancy", "b_iso"))
    if after["cell"] is not None:
        keys = ("a", "b", "c", "alpha", "beta", "gamma")
        assert bits(after["cell"][k] for k in keys) == bits(case["expected"]["cell"][k] for k in keys)
