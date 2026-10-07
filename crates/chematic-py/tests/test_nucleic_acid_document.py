import json
from pathlib import Path

import chematic


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = json.loads((ROOT / "validation/nucleic_acid_document_contract.json").read_text())


def _at_path(value, path):
    for segment in path:
        value = value[segment]
    return value


def test_nucleic_acid_valid_documents_roundtrip_and_edit():
    for case in FIXTURE["valid_cases"]:
        envelope = json.loads(chematic.nucleic_acid_validate_json(json.dumps(case["document"])))
        assert envelope["ok"] is True, case["id"]
        assert envelope["document"] == case["document"], case["id"]

        edited = json.loads(
            chematic.nucleic_acid_apply_json_command(
                json.dumps(case["document"]), json.dumps(case["command"])
            )
        )
        assert edited["ok"] is True, case["id"]
        assert _at_path(edited["document"], case["edited_path"]) == case["edited_value"]


def test_nucleic_acid_invalid_documents_have_typed_categories():
    for case in FIXTURE["invalid_cases"]:
        envelope = json.loads(chematic.nucleic_acid_validate_json(json.dumps(case["document"])))
        assert envelope["ok"] is False, case["id"]
        assert envelope["error"]["code"] == case["expected_code"], case["id"]
        assert envelope["error"]["path"].startswith("/")


def test_nucleic_acid_binding_enforces_collection_limits():
    strands = []
    atom_ids = []
    for index in range(65):
        atom_id = f"a{index}"
        atom_ids.append(atom_id)
        strands.append({
            "id": f"s{index}",
            "kind": "dna",
            "residues": [{
                "id": f"r{index}",
                "base": "A",
                "sugar": "deoxyribose",
                "atom_refs": [atom_id],
                "annotations": {},
            }],
            "annotations": {},
        })
    document = {
        "schema": "chematic.nucleic-acid.v1",
        "atom_ids": atom_ids,
        "strands": strands,
        "linkages": [],
        "annotations": {},
    }
    envelope = json.loads(chematic.nucleic_acid_validate_json(json.dumps(document)))
    assert envelope["ok"] is False
    assert envelope["error"]["code"] == "resource_limit"
    assert envelope["error"]["path"] == "/strands"
