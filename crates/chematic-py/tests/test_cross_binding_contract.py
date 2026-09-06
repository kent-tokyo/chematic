"""Consume the shared Rust/Python/Node/WASM contract fixture."""

import json
from pathlib import Path

import pytest

import chematic


_FIXTURE_PATH = Path(__file__).parents[3] / "validation" / "cross_binding_contract.json"
_DOCUMENT = json.loads(_FIXTURE_PATH.read_text())


def test_shared_fixture_schema_is_stable():
    assert _DOCUMENT["schema_version"] == 1
    assert len(_DOCUMENT["fixtures"]) == 4
    assert _DOCUMENT["descriptor_contract"]["schema_version"] == 1
    assert _DOCUMENT["descriptor_contract"]["fields"]["tpsa"]["unit"] == "A2"
    assert _DOCUMENT["fingerprint_contract"]["schema_version"] == 1
    assert _DOCUMENT["fingerprint_contract"]["operations"]["ecfp4"]["bytes"] == 256
    assert _DOCUMENT["fingerprint_contract"]["operations"]["maccs"]["bytes"] == 21
    assert _DOCUMENT["batch_canonicalization_contract"]["schema_version"] == 1


def test_python_binding_matches_shared_batch_canonicalization_contract():
    contract = _DOCUMENT["batch_canonicalization_contract"]
    actual = json.loads(chematic.canonicalize_smiles_batch_json(contract["inputs"]))
    assert actual["schema_version"] == 1
    assert actual["operation"] == "canonicalize_smiles"
    assert actual["status"] == "complete"
    assert actual["record_count"] == len(contract["expected"])
    observed = []
    for record in actual["records"]:
        result = {
            "input_index": record["input_index"],
            "input": record["input"],
            "status": record["status"],
        }
        if record["status"] == "accepted":
            result["canonical_smiles"] = record["canonical_smiles"]
        observed.append(result)
    assert observed == contract["expected"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["descriptor_contract"]["fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_descriptor_contract(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    assert mol.mw == pytest.approx(fixture["molecular_weight"], abs=1e-6)
    assert mol.tpsa == pytest.approx(fixture["tpsa"], abs=1e-6)
    assert mol.hbd == fixture["hbd"]
    assert mol.hba == fixture["hba"]
    assert mol.heavy_atoms == fixture["heavy_atoms"]


@pytest.mark.parametrize("fixture", _DOCUMENT["fixtures"], ids=lambda item: item["id"])
def test_python_binding_matches_shared_fixture(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    assert mol.smiles == fixture["canonical_smiles"]
    assert mol.heavy_atoms == fixture["heavy_atoms"]


@pytest.mark.parametrize(
    "fixture",
    _DOCUMENT["fingerprint_contract"]["fixtures"],
    ids=lambda item: item["id"],
)
def test_python_binding_matches_shared_fingerprint_shape(fixture):
    mol = chematic.from_smiles(fixture["smiles"])
    assert len(mol.ecfp4()) == 256
    assert len(mol.maccs()) == 21
    assert any(mol.ecfp4())
    assert any(mol.maccs())


@pytest.mark.parametrize(
    "case",
    _DOCUMENT["adversarial"],
    ids=lambda item: item["id"],
)
def test_python_binding_rejects_shared_adversarial_fixture(case):
    parsers = {
        "smiles": chematic.from_smiles,
        "mol": chematic.from_mol_block,
        "mol_v3000": chematic.from_mol_v3000,
        "mol2": chematic.from_mol2,
        "cml": chematic.from_cml,
        "cjson": chematic.from_cjson,
        "moljson": chematic.from_moljson,
        "cdxml": chematic.from_cdxml,
    }
    with pytest.raises((ValueError, RuntimeError, TypeError)):
        parsers[case["format"]](case["input"])


def test_python_cdxml_document_edit_preserves_multi_page_presentation():
    cdxml = '<CDXML>\n<page id="p1">\n<arrow id="a1"/>\n</page>\n<page id="p2">\n<text id="t1"/>\n</page>\n</CDXML>'
    edited = chematic.edit_cdxml_document_json(
        cdxml,
        json.dumps({"kind": "set_page_attribute", "page_id": "p2", "key": "title", "value": "Page 2"}),
    )
    assert 'title="Page 2"' in edited
    assert '<arrow id="a1"/>' in edited


def test_python_semantic_markush_contract_expands_with_mapping():
    model = {
        "schema": "chematic.semantic.v1",
        "atom_ids": ["a1", "a2"],
        "bond_ids": [],
        "r_groups": [{
            "id": "r1",
            "attachment_atoms": ["a2"],
            "alternatives": ["[*]O"],
            "selected_alternative": None,
        }],
        "polymer_units": [],
        "extensions": {},
    }
    selected = chematic.semantic_apply_json_command(
        json.dumps(model), json.dumps({"group_id": "r1", "alternative": 0})
    )
    expanded = json.loads(chematic.semantic_expand_json("CC", selected))
    assert expanded["schema"] == "chematic.semantic-expanded.v1"
    assert expanded["source_to_expanded"]["r1"] == [2]


def test_python_rxn_document_contract_is_loss_aware():
    document = {
        "id": "rxn-contract",
        "steps": [{
            "id": "step-1",
            "components": [
                {"id": "reactant-1", "role": "reactant", "smiles": "CC", "coefficient": 1, "origin": "authored"},
                {"id": "product-1", "role": "product", "smiles": "CC", "coefficient": 1, "origin": "authored"},
            ],
            "conditions": [],
            "provenance": [],
            "origin": "authored",
        }],
        "provenance": [],
    }
    rxn = chematic.to_rxn_document_json(json.dumps(document))
    decoded = json.loads(chematic.from_rxn_document_json(rxn))
    assert [c["role"] for c in decoded["steps"][0]["components"]] == ["reactant", "product"]

    document["steps"][0]["conditions"] = [{"key": "temperature", "value": "25 C"}]
    with pytest.raises(ValueError):
        chematic.to_rxn_document_json(json.dumps(document))
