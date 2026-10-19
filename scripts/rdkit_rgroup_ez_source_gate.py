#!/usr/bin/env python3
"""Compare source-built E/Z R-group decompositions with pinned RDKit defaults."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdRGroupDecomposition

ROOT = Path(__file__).resolve().parents[1]
RDKIT_VERSION = "2026.09.1"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(smiles: str, *, isomeric: bool) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=isomeric)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture", type=Path,
        default=ROOT / "validation/rdkit/rgroup_ez_holdout.json",
    )
    parser.add_argument("--source-jsonl", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"expected RDKit {RDKIT_VERSION}; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    fixture = json.loads(args.fixture.read_text(encoding="utf-8"))
    if fixture["schema_version"] != 1 or fixture["core"] != "C=C":
        raise ValueError("unsupported holdout fixture")
    substituents = fixture["substituents"]
    if substituents != ["F", "Cl", "Br", "I", "C", "CO"]:
        raise ValueError("unexpected substituent holdout")
    expected_pairs = list(itertools.combinations(substituents, 2))
    source_rows = [
        json.loads(line) for line in args.source_jsonl.read_text(encoding="utf-8").splitlines()
    ]
    if len(source_rows) != len(expected_pairs):
        raise ValueError("source row count differs from fixture")

    cases = []
    for (left, right), source in zip(expected_pairs, source_rows, strict=True):
        inputs = [f"{left}/C=C/{right}", f"{left}/C=C\\{right}"]
        if [source["left"], source["right"]] != [left, right] or source["inputs"] != inputs:
            raise ValueError(f"source pair or input mismatch: {left}, {right}")
        if len(source["rows"]) != 2:
            raise ValueError(f"source decomposition missing input: {left}, {right}")

        input_isomers = [canonical(value, isomeric=True) for value in inputs]
        input_constitutions = [canonical(value, isomeric=False) for value in inputs]
        input_stereo_only = (
            input_isomers[0] != input_isomers[1]
            and input_constitutions[0] == input_constitutions[1]
        )

        cores = [row["Core"] for row in source["rows"]]
        source_core_isomers = [canonical(value, isomeric=True) for value in cores]
        source_core_constitutions = [canonical(value, isomeric=False) for value in cores]
        source_groups_same = source["rows"][0]["R"] == source["rows"][1]["R"]
        source_preserves_geometry = (
            source_core_isomers[0] != source_core_isomers[1]
            and source_core_constitutions[0] == source_core_constitutions[1]
            and source_groups_same
        )

        rdkit_rows, unmatched = rdRGroupDecomposition.RGroupDecompose(
            [Chem.MolFromSmarts(fixture["core"])],
            [Chem.MolFromSmiles(value) for value in inputs],
            asSmiles=True,
            asRows=True,
        )
        if unmatched or len(rdkit_rows) != 2:
            raise ValueError(f"RDKit decomposition missing input: {left}, {right}")
        rdkit_cores = [row["Core"] for row in rdkit_rows]
        rdkit_core_isomers = [canonical(value, isomeric=True) for value in rdkit_cores]
        rdkit_core_constitutions = [canonical(value, isomeric=False) for value in rdkit_cores]
        rdkit_preserves_geometry = (
            rdkit_core_isomers[0] != rdkit_core_isomers[1]
            and rdkit_core_constitutions[0] == rdkit_core_constitutions[1]
            and all(
                rdkit_rows[0][key] == rdkit_rows[1][key]
                for key in rdkit_rows[0].keys() - {"Core"}
            )
        )
        cases.append({
            "substituents": [left, right],
            "inputs": inputs,
            "input_isomers": input_isomers,
            "input_stereo_only": input_stereo_only,
            "chematic_rows": source["rows"],
            "chematic_core_isomers": source_core_isomers,
            "chematic_preserves_geometry": source_preserves_geometry,
            "rdkit_default_rows": rdkit_rows,
            "rdkit_default_rows_identical": rdkit_rows[0] == rdkit_rows[1],
            "rdkit_preserves_geometry": rdkit_preserves_geometry,
        })

    paths = [
        "crates/chematic-chem/src/rgroup.rs",
        "crates/chematic-chem/src/rgroup_tests.rs",
        "crates/chematic-chem/examples/rgroup_ez_dump.rs",
        "scripts/rdkit_rgroup_ez_source_gate.py",
        "validation/rdkit/rgroup_ez_holdout.json",
    ]
    report = {
        "schema_version": 1,
        "profile": "source_rgroup_alkene_core_ez_holdout",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "rdkit_version": rdBase.rdkitVersion,
        "source_file_sha256": {path: digest(ROOT / path) for path in paths},
        "source_dump_sha256": digest(args.source_jsonl),
        "criterion": fixture["criterion"],
        "counts": {
            "total_pairs": len(cases),
            "stereo_only_inputs": sum(row["input_stereo_only"] for row in cases),
            "chematic_preserves_geometry": sum(
                row["chematic_preserves_geometry"] for row in cases
            ),
            "rdkit_default_preserves_geometry": sum(
                row["rdkit_preserves_geometry"] for row in cases
            ),
            "rdkit_default_rows_identical": sum(
                row["rdkit_default_rows_identical"] for row in cases
            ),
        },
        "cases": cases,
        "limits": [
            "15 systematic pairs from six substituents around one C=C core",
            "RDKit RGroupDecompose default arguments with asSmiles=True and asRows=True",
            "This measures retention of specified alkene geometry, not overall R-group accuracy.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report["counts"], sort_keys=True))
    return 0 if all(
        row["input_stereo_only"] and row["chematic_preserves_geometry"]
        for row in cases
    ) else 1


if __name__ == "__main__":
    raise SystemExit(main())
