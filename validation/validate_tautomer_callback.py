#!/usr/bin/env python3
"""Compare real Rust checkpoint/partial-result behavior with an installed RDKit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import chematic as ch
from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize as standard

CASES = [
    "",
    "CCO",
    "[Na+].[Cl-]",
    "CC=O",
    "OC=C",
    "CC(=O)C",
    "CC(=O)CC=O",
    "CC(=O)CC(=O)C",
    "Oc1ccccn1",
    "O=c1cccc[nH]1",
    "NC(=N)N",
    "NC(=N)NC(=N)N",
    "c1ccccc1O",
    "Cc1ccc(O)cc1",
    "CC(=O)/C=C/c1ccc(O)cc1",
    "C[C@@H](C=O)C(=O)O",
    "O=P(O)(O)O",
    "OP(O)O",
    "C[N+](=O)[O-]",
    "C[N+](O)=O",
    "ON=C(C)C",
    "CCN=O",
    "N=c1cc[nH]cc1",
    "O=C1CCCCC1",
]


def run_case(smiles, stop, max_tautomers, max_transforms):
    mol = ch.from_smiles(smiles) if smiles else ch.RWMol().GetMol()
    trace = []
    transforms = []

    def callback(progress):
        trace.append(list(progress.keys))
        transforms.append(progress.transforms)
        return stop is None or len(trace) < stop

    got = ch.enumerate_tautomers_checked(
        mol,
        max_tautomers=max_tautomers,
        max_transforms=max_transforms,
        callback=callback,
    )
    oracle_trace = []

    class Callback(standard.TautomerEnumeratorCallback):
        def __call__(self, original, result):
            oracle_trace.append(list(result.smiles))
            return stop is None or len(oracle_trace) < stop

    enum = standard.TautomerEnumerator()
    enum.SetMaxTautomers(max_tautomers)
    enum.SetMaxTransforms(max_transforms)
    enum.SetCallback(Callback())
    result = enum.Enumerate(Chem.MolFromSmiles(smiles))
    expected = {
        "keys": list(result.smiles),
        "smiles": [Chem.MolToSmiles(mol) for mol in result],
        "status": str(result.status),
        "trace": oracle_trace,
    }
    actual = {
        "keys": list(got.keys),
        "smiles": list(got.smiles),
        "status": got.status,
        "trace": trace,
    }
    return {
        "input": smiles,
        "stop_on_call": stop,
        "max_tautomers": max_tautomers,
        "max_transforms": max_transforms,
        "match": actual == expected,
        "actual": actual,
        "expected": expected,
        "transforms": transforms,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for smiles in CASES:
        for stop in [None, 1, 2, 3, 5, 10]:
            records.append(run_case(smiles, stop, 1000, 1000))
        for max_tautomers, max_transforms in [
            (1, 1000),
            (1000, 1),
            (1, 1),
            (2, 3),
            (3, 10),
            (10, 3),
        ]:
            records.append(run_case(smiles, None, max_tautomers, max_transforms))
    extension = Path(ch.chematic.__file__)
    report = {
        "rdkit_version": rdBase.rdkitVersion,
        "chematic_version": ch.__version__,
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
        "cases": len(records),
        "matched": sum(r["match"] for r in records),
        "checkpoint_calls": sum(len(r["actual"]["trace"]) for r in records),
        "records": records,
        "scope": "Concrete molecule profile, immutable key snapshots, matched-transform checkpoints; no wall-clock preemption or SetCallback API compatibility claim.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "records"}, indent=2))
    for record in records:
        if not record["match"]:
            print(
                "MISMATCH",
                record["input"],
                record["stop_on_call"],
                record["max_tautomers"],
                record["max_transforms"],
            )
    return 0 if report["matched"] == len(records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
