"""Source retention and bounded CJSON projection; development fixtures only."""
import json
import struct
from pathlib import Path
import pytest
import chematic
SOURCE = (Path(__file__).resolve().parents[3] / "validation/cjson_preservation/charged.cjson").read_text(encoding="utf-8")
PARTIAL = (Path(__file__).resolve().parents[3] / "validation/cjson_preservation/partial-charges.cjson").read_text(encoding="utf-8")

@pytest.mark.parametrize("root,legacy", [("ヘリウム β", None), (None, "legacy"), ("same", "same")])
def test_root_and_legacy_names_survive_both_conversion_directions(root, legacy):
    # Original synthetic title-only fixture under the project license.
    value = {"atoms": {"elements": {"number": [2]}, "coords": {"3d": [0,0,0]}}}
    if root is not None: value["name"] = root
    if legacy is not None: value["properties"] = {"name": legacy}
    source = json.dumps(value)
    expected = root if root is not None else legacy
    record = json.loads(chematic.parse_cjson_record_json(source))
    assert record["name"] == expected and record["name_conflict"] is False
    mol2 = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")["output"]
    after = json.loads(chematic.convert_with_report(mol2, "mol2", "cjson", policy="allow_loss")["output"])
    assert after["name"] == after["properties"]["name"] == expected

def test_conflicting_names_are_retained_without_selecting_a_mol2_title():
    source = json.dumps({"atoms": {"elements": {"number": [2]}, "coords": {"3d": [0,0,0]}},
                         "name": "root", "properties": {"name": "legacy"}})
    record = json.loads(chematic.parse_cjson_record_json(source))
    assert record["name"] is None and record["name_conflict"] is True
    assert chematic.roundtrip_cjson_record(source) == source
    report = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")
    assert report["status"] == "refused" and report["output"] is None
    assert report["candidate_atom_mapping"] == []

def test_ionic_mol2_formal_charges_are_separate_from_partial_charges():
    # Original synthetic NaCl representation under the project license.
    source = json.dumps({"atoms": {"elements": {"number": [11, 17]},
        "formalCharges": [1, -1], "coords": {"3d": [0, 0, 0, 3, 0, 0]}},
        "partialCharges": {"USER_CHARGES": [0.8, -0.8]}})
    report = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")
    assert report["candidate_graph_verified"] is True
    fields = {o["field"]: o["status"] for o in report["observations"]}
    assert fields["formal_charges"] == fields["partial_charges"] == "preserved"
    after = json.loads(chematic.parse_mol2_record_json(report["output"]))
    assert after["formal_charges"] == [1, -1]
    assert [a["partial_charge"] for a in after["atoms"]] == [0.8, -0.8]
    assert chematic.convert_with_report(source, "cjson", "mol2")["output"] is None

@pytest.mark.parametrize("source", [SOURCE, SOURCE.replace("\n", "\r\n")])
def test_exact_record_and_inspection(source):
    assert chematic.roundtrip_cjson_record(source) == source
    record = json.loads(chematic.parse_cjson_record_json(source))
    assert record["schema"] == "chematic.cjson-record.v1"
    assert record["source"] == source
    assert record["formal_charges"] == [0, 1, -1, 0]
    assert record["isotopes"] == [13, 15, None, 2]
    assert record["opaque_fields"] == ["/atoms/coords/2d", "/atoms/elements/symbol", "/atoms/partialCharges", "/vendor~1a~0b"]
    assert record["name"] == "isotope and charge fixture"
    assert record["name_conflict"] is False
    assert struct.pack(">d", record["coordinates"][0][0]) == struct.pack(">d", -1.4175225513503917)
    assert struct.pack(">d", record["coordinates"][0][1]) == struct.pack(">d", -0.0)
    assert struct.pack(">d", record["coordinates"][0][2]) == bytes.fromhex("0000000000000001")
    for policy in ["strict", "allow_loss", "preserve_source"]:
        report = chematic.convert_with_report(source, "cjson", "cjson", policy=policy)
        assert report["status"] == "completed"
        assert report["output"] == source
        assert report["verification_scope"] == "exact_source_bytes"
        assert report["candidate_graph_verified"] is False

@pytest.mark.parametrize("target", ["smiles", "smi"])
def test_checked_projection(target):
    allowed = chematic.convert_with_report(SOURCE, "cjson", target, policy="allow_loss")
    assert allowed["status"] == "completed_with_loss"
    assert allowed["candidate_graph_verified"] is True
    assert len(allowed["candidate_atom_mapping"]) == 4
    fields = {o["field"]: o["status"] for o in allowed["observations"]}
    assert fields["coordinates"] == fields["opaque_fields"] == "dropped"
    assert "13C" in allowed["output"] and "15N" in allowed["output"] and "2H" in allowed["output"]
    strict = chematic.convert_with_report(SOURCE, "cjson", target)
    assert strict["status"] == "refused"
    assert strict["output"] is None

@pytest.mark.parametrize("source", [
    SOURCE.replace("[0, 1, -1, 0]", "[]"),
    SOURCE.replace("[13, 15, 0, 2]", "null"),
    SOURCE.replace('"order": [1, 1]', '"order": [2.2, 1]'),
    SOURCE.replace("[0, 1, 1, 3]", "[0, 1, 1, 0]"),
    SOURCE.replace('"stereo": "opaque"', '"stereo": "opaque", "stereo": "other"'),
    SOURCE.replace('"chemicalJson": 1', '"chemicalJson": 0'),
])
def test_ambiguous_input_is_rejected(source):
    for call in [chematic.roundtrip_cjson_record, chematic.parse_cjson_record_json,
                 lambda text: chematic.convert_with_report(text, "cjson", "cjson")]:
        with pytest.raises(ValueError):
            call(source)

@pytest.mark.parametrize("n", [-0.0, -1.4175225513503917, 5e-324, 1.7976931348623157e308])
def test_cjson_to_mol2_exact_partial_and_coordinate_values(n):
    value = json.loads(PARTIAL)
    value["partialCharges"]["USER_CHARGES"][0] = n
    value["atoms"]["coords"]["3d"][0] = n
    source = json.dumps(value)
    record = json.loads(chematic.parse_cjson_record_json(source))
    assert record["opaque_fields"] == ["/properties/vendor"]
    assert record["name"] == "acetamide charges"
    assert struct.pack(">d", record["partial_charge_datasets"]["USER_CHARGES"][0]) == struct.pack(">d", n)
    allowed = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")
    assert allowed["status"] == "completed_with_loss"
    assert allowed["candidate_graph_verified"] is True
    after = json.loads(chematic.parse_mol2_record_json(allowed["output"]))
    assert struct.pack(">d", after["atoms"][0]["partial_charge"]) == struct.pack(">d", n)
    assert struct.pack(">d", after["coordinates"][0][0]) == struct.pack(">d", n)
    assert after["name"] == record["name"]
    assert after["charge_type"] == "USER_CHARGES"
    fields = {o["field"]: o["status"] for o in allowed["observations"]}
    for field in ["coordinates", "partial_charges", "charge_method", "molecule_name"]:
        assert fields[field] == "preserved"
    strict = chematic.convert_with_report(source, "cjson", "mol2")
    assert strict["status"] == "refused" and strict["output"] is None

@pytest.mark.parametrize("mode", ["multiple", "missing-coordinates", "name-newline", "method-whitespace"])
def test_valid_cjson_refused_without_guessing_mol2_representation(mode):
    value = json.loads(PARTIAL)
    if mode == "multiple": value["partialCharges"]["SECOND"] = [0] * 4
    elif mode == "missing-coordinates": del value["atoms"]["coords"]
    elif mode == "name-newline": value["properties"]["name"] = "a\n@<TRIPOS>ATOM"
    else: value["partialCharges"] = {"two methods": [0] * 4}
    source = json.dumps(value)
    assert chematic.roundtrip_cjson_record(source) == source
    report = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")
    assert report["status"] == "refused" and report["output"] is None
    assert report["candidate_atom_mapping"] == []

@pytest.mark.parametrize("invalid", [None, [], {"method": []}, {"method": [0, 0, 0]}, {"method": [0, 0, 0, "NaN"]}])
def test_malformed_partial_dataset_rejected(invalid):
    value = json.loads(PARTIAL)
    value["partialCharges"] = invalid
    source = json.dumps(value)
    for call in [chematic.roundtrip_cjson_record, chematic.parse_cjson_record_json,
                 lambda text: chematic.convert_with_report(text, "cjson", "mol2", policy="allow_loss")]:
        with pytest.raises(ValueError): call(source)
