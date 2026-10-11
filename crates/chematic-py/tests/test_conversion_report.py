"""Checked conversion: hand-derived chemistry and shared envelope contract."""
import json
from pathlib import Path
import pytest
import chematic

SOURCE = (Path(__file__).parents[3] / "validation/mol2_preservation/amide.mol2").read_text()

@pytest.mark.parametrize("policy,status", [("strict", "refused"), ("allow_loss", "completed_with_loss")])
def test_checked_cross_format(policy, status):
    report = chematic.convert_with_report(SOURCE, "mol2", "smiles", policy=policy)
    assert report == json.loads(chematic.convert_with_report_json(SOURCE, "mol2", "smiles", policy))
    assert report["schema"] == "chematic.conversion-report.v1"
    assert report["status"] == status
    assert report["candidate_graph_verified"] is True
    assert len(report["candidate_atom_mapping"]) == 4
    assert len(report["candidate_bond_mapping"]) == 3
    fields = {o["field"]: o["status"] for o in report["observations"]}
    assert fields["coordinates"] == fields["partial_charges"] == fields["opaque_sections"] == "dropped"
    assert fields["formal_charge_interpretation"] == fields["stereochemistry"] == "unverified"
    if policy == "strict":
        assert report["output"] is None
    else:
        mol = chematic.from_smiles(report["output"])
        graph = json.loads(mol.to_moljson())
        assert len(graph["atoms"]) == 4
        assert len(graph["bonds"]) == 3
        assert all(atom["charge"] == 0 for atom in graph["atoms"])

@pytest.mark.parametrize("policy", ["strict", "allow_loss", "preserve_source"])
def test_source_retention_has_separate_scope(policy):
    source = SOURCE.replace("\n", "\r\n")
    report = chematic.convert_with_report(source, ".MOL2", "mol2", policy=policy)
    assert report["output"] == source
    assert report["status"] == "completed"
    assert report["verification_scope"] == "exact_source_bytes"
    assert report["candidate_graph_verified"] is False
    assert any(o["status"] == "opaque_retained" for o in report["observations"])

@pytest.mark.parametrize("text,source,target,policy", [
    ("invalid", "mol2", "mol2", "strict"),
    (SOURCE, "sdf", "smiles", "strict"),
    (SOURCE, "mol2", "cml", "allow_loss"),
    (SOURCE, "mol2", "smiles", "preserve_source"),
    (SOURCE, "mol2", "smiles", "typo"),
    (SOURCE.replace("ACM -0.70", "ACM NaN"), "mol2", "mol2", "preserve_source"),
])
def test_errors_are_distinct_from_policy_refusals(text, source, target, policy):
    with pytest.raises(ValueError):
        chematic.convert_with_report(text, source, target, policy=policy)

def test_default_is_strict():
    assert chematic.convert_with_report(SOURCE, "mol2", "smiles")["status"] == "refused"
    assert json.loads(chematic.convert_with_report_json(SOURCE, "mol2", "smiles"))["output"] is None

@pytest.mark.parametrize("policy", ["strict", "allow_loss"])
def test_cjson_keeps_3d_geometry_and_projected_charge_with_loss_report(policy):
    source = SOURCE.replace("0.123456789 0 0 C.3", "0.123456789 2.5 -3.75 N.4")
    report = chematic.convert_with_report(source, "mol2", "cjson", policy=policy)
    fields = {o["field"]: o["status"] for o in report["observations"]}
    assert fields["coordinates"] == fields["projected_formal_charges"] == "preserved"
    assert fields["partial_charges"] == fields["charge_method"] == fields["molecule_name"] == "preserved"
    assert fields["formal_charge_interpretation"] == "unverified"
    assert "implicit_hydrogens" not in fields
    assert report["candidate_graph_verified"] is True
    if policy == "strict":
        assert report["status"] == "refused"
        assert report["output"] is None
    else:
        mol, coords = chematic.from_cjson(report["output"])
        assert coords[0] == [0.123456789, 2.5, -3.75]
        assert json.loads(mol.to_cjson(coords))["atoms"]["formalCharges"] == [1, 0, 0, 0]

@pytest.mark.parametrize("value", ["-1.4175225513503917", "-0.0", "5e-324", "1.7976931348623157e308"])
def test_cjson_coordinate_roundtrip_matches_python_float_bits(value):
    import struct
    source = SOURCE.replace("0.123456789", value)
    report = chematic.convert_with_report(source, "mol2", "cjson", policy="allow_loss")
    mol, coords = chematic.from_cjson(report["output"])
    assert struct.pack(">d", coords[0][0]) == struct.pack(">d", float(value))
    # Exercise the existing CJSON API too: the precision fix is in shared Rust.
    _, reparsed_coords = chematic.from_cjson(mol.to_cjson(coords))
    assert struct.pack(">d", reparsed_coords[0][0]) == struct.pack(">d", float(value))
    assert any(o["field"] == "coordinates" and o["status"] == "preserved" for o in report["observations"])

@pytest.mark.parametrize("token", ["nc", "du", "un"])
def test_cjson_changed_bonds_stay_unverified(token):
    source = SOURCE.replace("12 7 3 am", f"12 7 3 {token}")
    report = chematic.convert_with_report(source, "mol2", "cjson", policy="allow_loss")
    assert report["candidate_graph_verified"] is False
    assert report["candidate_atom_mapping"] == []
    assert any(o["field"] == "projected_graph" and o["status"] == "unverified" for o in report["observations"])

def test_cjson_rejects_preserve_source_policy():
    with pytest.raises(ValueError):
        chematic.convert_with_report(SOURCE, "mol2", "cjson", policy="preserve_source")

@pytest.mark.parametrize("literal", ["-0.0", "-1.4175225513503917", "5e-324", "1.7976931348623157e308"])
def test_cjson_preserves_declared_partial_charge_bits_and_name(literal):
    import struct
    source = SOURCE.replace("ACM 0.70", f"ACM {literal}")
    report = chematic.convert_with_report(source, "mol2", "cjson", policy="allow_loss")
    output = json.loads(report["output"])
    assert output["properties"]["name"] == "acetamide"
    assert struct.pack(">d", output["partialCharges"]["USER_CHARGES"][0]) == struct.pack(">d", float(literal))
    fields = {o["field"]: o["status"] for o in report["observations"]}
    assert fields["partial_charges"] == fields["charge_method"] == "preserved"

@pytest.mark.parametrize("source", [SOURCE.replace("ACM 0.70", "ACM"), SOURCE.replace("USER_CHARGES", "NO_CHARGES")])
def test_cjson_incomplete_or_undeclared_partial_charges_are_not_filled(source):
    report = chematic.convert_with_report(source, "mol2", "cjson", policy="allow_loss")
    assert "partialCharges" not in json.loads(report["output"])
    assert any(o["field"] == "partial_charges" and o["status"] == "dropped" for o in report["observations"])
