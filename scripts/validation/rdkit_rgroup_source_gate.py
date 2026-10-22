#!/usr/bin/env python3
"""Compare a source-built labelled R-group dump with RDKit 2026.09.1."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdRGroupDecomposition

ROOT = Path(__file__).resolve().parents[2]
CASES = [
    ("c1ccccc1", ["c1ccccc1", "Cc1ccccc1", "Clc1ccccc1"]),
    ("c1ccccc1", ["Cc1ccccc1", "Cc1ccc(Cl)cc1", "Cc1cc(Cl)ccc1"]),
    ("c1ccccc1", ["Cc1c(Cl)cccc1", "Cc1ccc(Cl)cc1", "Cc1cc(Cl)ccc1"]),
    ("c1ccccc1", ["Cc1ccccc1", "CCc1ccccc1", "Oc1ccccc1"]),
    ("c1ccncc1", ["c1ccncc1", "Cc1ccncc1", "Clc1ccncc1"]),
    ("c1ccncc1", ["Cc1ccncc1", "Cc1ccc(Cl)nc1", "Cc1cc(Cl)ncc1"]),
    ("C1CCCCC1", ["C1CCCCC1", "CC1CCCCC1", "ClC1CCCCC1"]),
    ("C1CCCCC1", ["CC1CCCCC1", "CC1CCC(Cl)CC1", "CC1CC(Cl)CCC1"]),
    ("c1cc([*:1])ccc1[*:2]", ["Cc1ccccc1", "Cc1ccc(CC)cc1", "Cc1cc(CC)ccc1"]),
    ("C([*:1])([*:2])O", ["CO", "CCO", "CC(O)N"]),
    ("C([*:1])([*:2])([*:3])O", ["CO", "CCO", "CC(O)N"]),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"unreadable R-group fragment: {smiles}")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")

    binary = args.binary.resolve(strict=True)
    payload = "".join(
        json.dumps({"core": core, "molecules": molecules}) + "\n"
        for core, molecules in CASES
    )
    process = subprocess.run(
        [str(binary)], input=payload, text=True, capture_output=True, check=True
    )
    source_rows = [json.loads(line) for line in process.stdout.splitlines()]
    if len(source_rows) != len(CASES):
        raise ValueError(f"source emitted {len(source_rows)} rows for {len(CASES)} cases")

    cases = []
    for (core, molecules), source in zip(CASES, source_rows, strict=True):
        oracle, unmatched = rdRGroupDecomposition.RGroupDecompose(
            [Chem.MolFromSmarts(core)],
            [Chem.MolFromSmiles(value) for value in molecules],
            asSmiles=True,
            asRows=True,
        )
        exact = not unmatched and len(oracle) == len(source)
        if exact:
            for expected, actual in zip(oracle, source, strict=True):
                if actual is None or expected.keys() != actual.keys():
                    exact = False
                    break
                if any(canonical(value) != canonical(actual[key])
                       for key, value in expected.items()):
                    exact = False
                    break
        cases.append({
            "core": core,
            "molecules": molecules,
            "rdkit_rows": oracle,
            "unmatched": unmatched,
            "source_rows": source,
            "exact_graph_and_labels": exact,
        })

    report = {
        "profile": "source_rdkit_2026_09_1_labelled_rgroup_symmetry",
        "rdkit_version": rdBase.rdkitVersion,
        "rdkit_rgroup_module_sha256": sha256(Path(rdRGroupDecomposition.__file__)),
        "python_version": sys.version.split()[0],
        "source_binary_sha256": sha256(binary),
        "source_file_sha256": {
            path: sha256(ROOT / path)
            for path in (
                "crates/chematic-chem/src/rgroup.rs",
                "crates/chematic-chem/src/rgroup_tests.rs",
                "crates/chematic-chem/examples/rgroup_gate_probe.rs",
                "scripts/validation/rdkit_rgroup_source_gate.py",
            )
        },
        "agreement": {
            "exact": sum(row["exact_graph_and_labels"] for row in cases),
            "total": len(cases),
        },
        "cases": cases,
        "limits": [
            "Eleven exposed small batches, three molecules per batch.",
            "RDKit graph and label agreement is a compatibility measure, not general R-group accuracy.",
            "The CheMatic arm is a source-built Rust example, not a published Python wheel.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"R-group: {report['agreement']['exact']}/{report['agreement']['total']} exact")
    return 0 if report["agreement"]["exact"] == len(CASES) else 1


if __name__ == "__main__":
    raise SystemExit(main())
