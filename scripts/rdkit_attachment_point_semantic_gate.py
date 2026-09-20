#!/usr/bin/env python3
"""Check a bounded RDKit -> chematic -> RDKit CXSMILES attachment contract.

RDKit creates each source CXSMILES payload. CheMatic's CLI reports the
supported metadata it retained and emits a CXSMILES payload that RDKit reopens.
The gate is deliberately narrow: it proves atom-map keyed labels, ordinary
topology, and degree-one wildcard attachment identity. It does not claim MDL
ATTCHPT collapse, V3000 ENDPTS/ATTACH semantics, or generic CX query editing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from rdkit import Chem


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = ROOT / "target" / "debug" / "chematic"

# Each label is installed through RDKit before it writes the CXSMILES source.
# Expected attachment points are keyed by atom-map number, not emitted atom
# order, so the test catches a label accidentally moving during serialization.
CASES = {
    "attachment_one": {
        "smiles": "[*:11][C:12]",
        "labels": {11: "_AP1"},
        "attachment_points": [{"atom_map": 11, "identifier": 1}],
    },
    "attachment_u32_max": {
        "smiles": "[*:11][C:12]",
        "labels": {11: "_AP4294967295"},
        "attachment_points": [{"atom_map": 11, "identifier": 4294967295}],
    },
    "two_attachment_points": {
        "smiles": "[C:1]([*:11])[*:12]",
        "labels": {11: "_AP1", 12: "_AP2"},
        "attachment_points": [
            {"atom_map": 11, "identifier": 1},
            {"atom_map": 12, "identifier": 2},
        ],
    },
    "non_dummy_label": {
        "smiles": "[C:1][*:11]",
        "labels": {1: "_AP7"},
        "attachment_points": [],
    },
    "zero_is_not_attachment": {
        "smiles": "[*:11][C:12]",
        "labels": {11: "_AP0"},
        "attachment_points": [],
    },
    "signed_is_not_attachment": {
        "smiles": "[*:11][C:12]",
        "labels": {11: "_AP-1"},
        "attachment_points": [],
    },
    "text_is_not_attachment": {
        "smiles": "[*:11][C:12]",
        "labels": {11: "_AP1x"},
        "attachment_points": [],
    },
    "overflow_is_not_attachment": {
        "smiles": "[*:11][C:12]",
        "labels": {11: "_AP4294967296"},
        "attachment_points": [],
    },
}


def cli_provenance(cli: Path) -> dict[str, str]:
    with cli.open("rb") as handle:
        return {
            "cli_name": cli.name,
            "cli_sha256": hashlib.file_digest(handle, "sha256").hexdigest(),
        }


def labeled_rdkit_molecule(case: dict) -> Chem.Mol:
    molecule = Chem.MolFromSmiles(case["smiles"])
    if molecule is None:
        raise RuntimeError(f"RDKit rejected fixture: {case['smiles']}")
    for atom_map, label in case["labels"].items():
        matches = [atom for atom in molecule.GetAtoms() if atom.GetAtomMapNum() == atom_map]
        if len(matches) != 1:
            raise RuntimeError(f"fixture atom map {atom_map} is not unique")
        matches[0].SetProp("atomLabel", label)
    return molecule


def signature(molecule: Chem.Mol) -> dict:
    atom_by_map = {}
    for atom in molecule.GetAtoms():
        atom_map = atom.GetAtomMapNum()
        if not atom_map:
            raise RuntimeError("fixture/output omitted required atom map")
        atom_by_map[str(atom_map)] = {
            "atomic_number": atom.GetAtomicNum(),
            "formal_charge": atom.GetFormalCharge(),
            "isotope": atom.GetIsotope(),
            "label": atom.GetProp("atomLabel") if atom.HasProp("atomLabel") else None,
        }
    bonds = sorted(
        (
            min(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()),
            max(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()),
            str(bond.GetBondType()),
            str(bond.GetStereo()),
        )
        for bond in molecule.GetBonds()
    )
    return {"atoms_by_map": atom_by_map, "bonds_by_map": bonds}


def run_cli(cli: Path, source: str) -> tuple[int, dict | None, str]:
    run = subprocess.run(
        [str(cli), "cxsmiles", source],
        capture_output=True,
        text=True,
    )
    if run.returncode:
        return run.returncode, None, run.stderr.strip()
    try:
        return run.returncode, json.loads(run.stdout), ""
    except json.JSONDecodeError as error:
        return run.returncode, None, f"invalid CLI JSON: {error}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cli", type=Path, default=DEFAULT_CLI)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cli = args.cli.resolve()
    if not cli.is_file():
        parser.error(f"chematic CLI not found: {cli}")

    rows = []
    for case_id, case in CASES.items():
        original = labeled_rdkit_molecule(case)
        source = Chem.MolToCXSmiles(original, canonical=False)
        exit_code, record, error = run_cli(cli, source)
        row: dict[str, object] = {
            "id": case_id,
            "source_smiles": case["smiles"],
            "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "cli_exit_code": exit_code,
        }
        if record is None:
            row["error"] = error
            rows.append(row)
            continue

        written = record["cxsmiles"]
        reopened = Chem.MolFromSmiles(written)
        if reopened is None:
            row["error"] = "RDKit rejected CheMatic CXSMILES output"
            row["chematic_cxsmiles"] = written
            rows.append(row)
            continue
        expected_points = sorted(case["attachment_points"], key=lambda point: point["atom_map"])
        observed_points = sorted(
            (
                {"atom_map": point["atom_map"], "identifier": point["identifier"]}
                for point in record["marked_attachment_points"]
            ),
            key=lambda point: point["atom_map"],
        )
        before = signature(original)
        after = signature(reopened)
        row.update(
            {
                "chematic_cxsmiles": written,
                "output_sha256": hashlib.sha256(written.encode()).hexdigest(),
                "before": before,
                "after": after,
                "rdkit_semantic_equal": before == after,
                "expected_attachment_points": expected_points,
                "observed_attachment_points": observed_points,
                "attachment_identity_equal": expected_points == observed_points,
            }
        )
        rows.append(row)

    failures = [
        row["id"]
        for row in rows
        if row.get("cli_exit_code") != 0
        or not row.get("rdkit_semantic_equal", False)
        or not row.get("attachment_identity_equal", False)
    ]
    result = {
        "schema_version": 1,
        "profile": "rdkit_cxsmiles_attachment_point_semantic_roundtrip_v1",
        "rdkit_version": Chem.rdBase.rdkitVersion,
        **cli_provenance(cli),
        "cases": len(rows),
        "semantic_equal": len(rows) - len(failures),
        "failures": failures,
        "scope": "RDKit-written CXSMILES atom-map keyed labels, ordinary topology, and degree-one wildcard attachment identity",
        "not_claimed": [
            "MDL ATTCHPT collapse or position semantics",
            "V3000 ENDPTS/ATTACH or haptic chemistry",
            "generic CX query editing or unknown-field preservation",
            "cross-binding parity beyond the CLI route",
        ],
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("cases", "semantic_equal", "failures")}))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
