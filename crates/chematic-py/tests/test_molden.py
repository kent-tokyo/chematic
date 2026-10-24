"""Molden geometry/vibration data and explicit units survive public API edits."""
import copy
import json
from pathlib import Path
import pytest
import chematic
PACKET=json.loads((Path(__file__).resolve().parents[3]/"validation/molden_vibrations/cases.json").read_text(encoding="utf-8"))

@pytest.mark.parametrize("case",PACKET["cases"],ids=lambda c:c["id"])
def test_roundtrip(case):
    record=chematic.parse_molden(case["input"].encode())
    assert record==case["expected"]
    assert chematic.parse_molden(chematic.write_molden(record))==record

@pytest.mark.parametrize("case",PACKET["refusals"],ids=lambda c:c["id"])
def test_refusal(case):
    with pytest.raises(ValueError):chematic.parse_molden(case["input"])

def test_edits_and_absent_raman_are_retained():
    record=copy.deepcopy(PACKET["cases"][2]["expected"])
    record["frequencies"][0]=-123.45678901234567
    record["normal_modes"][0]["displacements_bohr"][0]=[1e-130,.12345678901234567,-3]
    record["intensities"][0]["raman"]=0.
    assert chematic.parse_molden(chematic.write_molden(record))==record
    record["orbitals"]=[]
    with pytest.raises(ValueError):chematic.write_molden(record)
    del record["orbitals"]
    record["intensities"][0]["ir"]=float("nan")
    with pytest.raises(ValueError):chematic.write_molden(record)
