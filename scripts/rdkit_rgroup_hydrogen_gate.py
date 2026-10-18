#!/usr/bin/env python3
"""Check labelled R-group hydrogen-column semantics against RDKit."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path

import chematic
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdRGroupDecomposition

CASES = [
    ("[c:1]1ccccc1", "c1ccccc1", "Cc1ccccc1"),
    ("[C:1]1CCCCC1", "C1CCCCC1", "CC1CCCCC1"),
    ("[O:1]C", "CO", "COC"),
    ("[N:1]C", "CN", "CNC"),
]


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unreadable R-group output: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1; found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    results = []
    for core, no_group, group in CASES:
        rdkit_rows, unmatched = rdRGroupDecomposition.RGroupDecompose(
            [Chem.MolFromSmarts(core)],
            [Chem.MolFromSmiles(no_group), Chem.MolFromSmiles(group)],
            asSmiles=True,
            asRows=True,
        )
        source_rows = chematic.rgroup_decompose_labeled(
            core, [chematic.from_smiles(no_group), chematic.from_smiles(group)]
        )
        single = chematic.rgroup_decompose_labeled(core, [chematic.from_smiles(no_group)])
        exact = (
            not unmatched
            and all(row is not None for row in source_rows)
            and all(
                canonical(rdkit_rows[index][key]) == canonical(source_rows[index][key])
                for index in range(2)
                for key in ("Core", "R1")
            )
            and single[0] is not None
            and "R1" not in single[0]
        )
        results.append({"core": core, "molecules": [no_group, group],
                        "rdkit_rows": rdkit_rows, "chematic_rows": source_rows,
                        "single_unsubstituted_row": single[0], "graph_and_labels_match": exact})

    extension = Path(importlib.import_module("chematic.chematic").__file__)
    report = {
        "profile": "source_labelled_rgroup_hydrogen_columns",
        "versions": {"chematic": chematic.__version__, "rdkit": rdBase.rdkitVersion},
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "agreement": {"matched": sum(row["graph_and_labels_match"] for row in results),
                      "total": len(results)},
        "cases": results,
        "limits": ["one explicitly mapped core atom", "two molecules per core",
                   "graph and label comparison, not raw SMILES string comparison"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(f"R-group: {report['agreement']['matched']}/{report['agreement']['total']} exact")
    if not all(row["graph_and_labels_match"] for row in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
