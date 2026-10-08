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


SEQUENCES = json.loads((ROOT / "validation/nucleic_acid_edit_sequences.json").read_text())


def _ownership(document):
    owner = {}
    for strand in document["strands"]:
        for residue in strand["residues"]:
            for atom in residue["atom_refs"]:
                assert atom not in owner, f"atom {atom} owned twice"
                owner[atom] = residue["id"]
    assert set(document["atom_ids"]) <= set(owner)
    return owner


def test_nucleic_acid_edit_sequences_reload_and_keep_atom_ownership():
    starts = {case["id"]: case["document"] for case in FIXTURE["valid_cases"]}
    for sequence in SEQUENCES["sequences"]:
        document = starts[sequence["start"]]
        owners = _ownership(document)
        for k, step in enumerate(sequence["steps"]):
            envelope = json.loads(
                chematic.nucleic_acid_apply_json_command(json.dumps(document), json.dumps(step["command"]))
            )
            assert envelope["ok"] is step["ok"], (sequence["id"], k, envelope)
            if step["ok"]:
                assert envelope["document"] == step["document"], (sequence["id"], k)
                document = envelope["document"]
                reread = json.loads(chematic.nucleic_acid_validate_json(json.dumps(document)))
                assert reread["ok"] is True and reread["document"] == document, (sequence["id"], k)
                assert _ownership(document) == owners, (sequence["id"], k)
            else:
                assert envelope["error"]["code"] == step["code"], (sequence["id"], k)
                assert envelope["error"]["path"] == step["path"], (sequence["id"], k)
