"""Public CONFIG APIs retain simulation vectors without guessing units."""
import copy
import json
from pathlib import Path

import pytest
import chematic

PACKET = json.loads((Path(__file__).resolve().parents[3]/"validation/dlpoly_config/cases.json").read_text(encoding="utf-8"))

@pytest.mark.parametrize("case",PACKET["cases"],ids=lambda c:c["id"])
def test_roundtrip(case):
    record=chematic.parse_dlpoly_config(case["input"].encode())
    assert record==case["expected"]
    assert chematic.parse_dlpoly_config(chematic.write_dlpoly_config(record))==record

@pytest.mark.parametrize("case",PACKET["refusals"],ids=lambda c:c["id"])
def test_refusal(case):
    with pytest.raises(ValueError):chematic.parse_dlpoly_config(case["input"])

def test_edited_vectors_and_opaque_context():
    record=copy.deepcopy(PACKET["cases"][2]["expected"])
    record["atoms"][0]["velocity"]=[1.2345678901234567,0,0]
    record["atoms"][0]["force"]=[-123.45678901234567,0,0]
    record["atoms"][0]["header_extra"]=["opaque","-0.1250"]
    assert chematic.parse_dlpoly_config(chematic.write_dlpoly_config(record))==record
    record["units"]="angstrom"
    with pytest.raises(ValueError):chematic.write_dlpoly_config(record)
    del record["units"]
    record["atoms"][0]["force"]=None
    with pytest.raises(ValueError):chematic.write_dlpoly_config(record)
    record["atoms"][0]["force"]=[float("nan"),0,0]
    with pytest.raises(ValueError):chematic.write_dlpoly_config(record)
