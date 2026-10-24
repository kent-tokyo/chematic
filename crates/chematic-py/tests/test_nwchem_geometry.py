"""The same original NWChem geometry cases exercised through the built Python extension."""
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3] / "validation/nwchem_geometry/cases.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", PACKET["valid"], ids=lambda c: c["id"])
def test_nwchem_geometry_roundtrip(case):
    for source in (case["input"], case["input"].encode("ascii")):
        record = chematic.parse_nwchem_geometry(source)
        assert record == case["expected"]
        assert chematic.parse_nwchem_geometry(chematic.write_nwchem_geometry(record)) == record


@pytest.mark.parametrize("case", PACKET["invalid"], ids=lambda c: c["id"])
def test_nwchem_geometry_refusal(case):
    with pytest.raises(ValueError):
        chematic.parse_nwchem_geometry(case["input"])


def test_nwchem_geometry_edits_and_unknown_fields():
    record = chematic.parse_nwchem_geometry(PACKET["valid"][0]["input"])
    record["centers"][0]["nuclear_charge"] = -0.98765
    assert chematic.parse_nwchem_geometry(chematic.write_nwchem_geometry(record)) == record
    record["centers"][0]["mass_amu"] = -1.0
    with pytest.raises(ValueError):
        chematic.write_nwchem_geometry(record)
    record["centers"][0]["mass_amu"] = 1.0
    record["centers"][0]["nuclear_charge"] = float("nan")
    with pytest.raises(ValueError):
        chematic.write_nwchem_geometry(record)
    record["centers"][0]["nuclear_charge"] = 0.0
    record["centers"][0]["unrecognized"] = 1
    with pytest.raises(ValueError):
        chematic.write_nwchem_geometry(record)
