#!/usr/bin/env python3
"""Write (or check) the nucleic-acid edit-sequence fixture.

Each sequence starts from a valid document of
``validation/nucleic_acid_document_contract.json`` and applies commands one
after another; every step records the envelope chematic returns (the edited
document, or the typed error with the document left as it was). The Rust,
Python and WASM tests replay the sequences and require the same envelopes, a
document that validates again after serialization, and unchanged atom
ownership (every atom id owned by exactly one residue, as in the start).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import chematic

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "validation/nucleic_acid_document_contract.json"
OUTPUT = ROOT / "validation/nucleic_acid_edit_sequences.json"


def residue(rid):
    return {"kind": "residue", "id": rid}


SEQUENCES = {
    "short_dna": [
        {"kind": "set_annotation", "target": residue("d2"), "key": "selected", "value": True},
        {"kind": "set_annotation", "target": {"kind": "linkage", "id": "dl1"}, "key": "note", "value": "checked"},
        {"kind": "set_residue_identity", "residue_id": "d1", "base": "G", "sugar": "deoxyribose"},
        {"kind": "set_annotation", "target": residue("missing"), "key": "selected", "value": True},
        {"kind": "remove_annotation", "target": residue("d2"), "key": "selected"},
        {"kind": "set_annotation", "target": {"kind": "document"}, "key": "revision", "value": 2},
    ],
    "short_rna": [
        {"kind": "set_annotation", "target": {"kind": "strand", "id": "rna-1"}, "key": "label", "value": "RNA sample"},
        {"kind": "set_residue_identity", "residue_id": "r1", "base": "U", "sugar": "ribose"},
        {"kind": "set_residue_identity", "residue_id": "r1", "base": "Q", "sugar": "ribose"},
        {"kind": "remove_annotation", "target": {"kind": "strand", "id": "rna-1"}, "key": "label"},
    ],
    "explicit_modification": [
        {"kind": "set_residue_identity", "residue_id": "m1", "base": "C", "modification": "5mC", "sugar": "deoxyribose"},
        {"kind": "set_residue_identity", "residue_id": "m1", "base": "C", "modification": "not-a-modification", "sugar": "deoxyribose"},
        {"kind": "set_residue_identity", "residue_id": "m1", "base": "C", "sugar": "deoxyribose"},
        {"kind": "set_annotation", "target": residue("m1"), "key": "", "value": 1},
    ],
}


def build() -> dict:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    starts = {case["id"]: case["document"] for case in contract["valid_cases"]}
    sequences = []
    for sid, commands in SEQUENCES.items():
        document = starts[sid]
        steps = []
        for command in commands:
            envelope = json.loads(
                chematic.nucleic_acid_apply_json_command(json.dumps(document), json.dumps(command))
            )
            if envelope["ok"]:
                document = envelope["document"]
                steps.append({"command": command, "ok": True, "document": document})
            else:
                steps.append({"command": command, "ok": False,
                              "code": envelope["error"]["code"], "path": envelope["error"]["path"]})
        sequences.append({"id": sid, "start": sid, "steps": steps})
    return {
        "schema_version": 1,
        "contract": "validation/nucleic_acid_document_contract.json",
        "generated_with": f"chematic {chematic.__version__}",
        "sequences": sequences,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    data = build()
    text = json.dumps(data, indent=1, sort_keys=False) + "\n"
    if args.check:
        old = json.loads(OUTPUT.read_text(encoding="utf-8"))
        old.pop("generated_with", None)
        data.pop("generated_with", None)
        same = old == data
        print("edit sequences match" if same else "edit sequences differ")
        return 0 if same else 1
    OUTPUT.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
