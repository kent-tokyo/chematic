#!/usr/bin/env python3
"""Compare automatic R-group columns and multi-substituent rows with RDKit."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdRGroupDecomposition

ROOT = Path(__file__).resolve().parents[1]
CASES = [
    ("benzene_positions", "c1ccccc1", [
        "c1ccccc1", "Cc1ccccc1", "Cc1ccccc1C", "Cc1cccc(C)c1", "Cc1ccc(C)cc1",
    ]),
    ("shared_column", "CC", ["CC", "CCC", "CC(C)C", "CC(O)C"]),
    ("heteroatom_sites", "CO", ["CO", "CCO", "COC", "CC(O)C"]),
    ("asymmetric_ring", "c1ccncc1", ["c1ccncc1", "Cc1ccncc1", "Cc1ccc(C)nc1"]),
]


def canonical(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"unreadable R-group output: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    cases = []
    for name, core, smiles in CASES:
        rdkit_rows, unmatched = rdRGroupDecomposition.RGroupDecompose(
            [Chem.MolFromSmarts(core)],
            [Chem.MolFromSmiles(value) for value in smiles],
            asSmiles=True,
            asRows=True,
        )
        chematic_rows = chematic.rgroup_decompose_labeled(
            core, [chematic.from_smiles(value) for value in smiles]
        )
        if unmatched or len(rdkit_rows) != len(smiles) or len(chematic_rows) != len(smiles):
            raise ValueError(f"incomplete decomposition for {name}: {unmatched}")
        for value, rdkit_row, chematic_row in zip(smiles, rdkit_rows, chematic_rows):
            exact = (
                chematic_row is not None
                and rdkit_row.keys() == (chematic_row.keys() - {"mol_idx"})
                and all(canonical(group) == canonical(chematic_row[key])
                        for key, group in rdkit_row.items())
            )
            cases.append({
                "set": name, "core": core, "molecule": value,
                "rdkit_row": rdkit_row, "chematic_row": chematic_row,
                "graph_and_labels_match": exact,
            })

    extension = Path(importlib.import_module("chematic.chematic").__file__)
    report = {
        "profile": "source_automatic_rgroup_columns",
        "versions": {"chematic": chematic.__version__, "rdkit": rdBase.rdkitVersion},
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "source_sha256": {
            "crates/chematic-chem/src/rgroup.rs": hashlib.sha256(
                (ROOT / "crates/chematic-chem/src/rgroup.rs").read_bytes()
            ).hexdigest(),
        },
        "agreement": {"matched": sum(row["graph_and_labels_match"] for row in cases),
                      "total": len(cases)},
        "cases": cases,
        "limits": [
            "Four selected small-molecule batches, not a representative R-group corpus",
            "Graph and R-label comparison after canonicalizing with RDKit, not raw SMILES text",
            "One SMARTS core per batch and default RDKit RGroupDecompose parameters",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"R-group automatic columns: {report['agreement']['matched']}/{report['agreement']['total']} exact")
    if report["agreement"]["matched"] != report["agreement"]["total"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
