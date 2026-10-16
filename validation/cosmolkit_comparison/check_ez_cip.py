#!/usr/bin/env python3
"""Compare source-built E/Z CIP labels with RDKit and COSMolKit by bond atoms.

Build ``cargo build -p chematic-chem --example ez_cip_dump --offline`` first.
Run with RDKit 2025.09.3 and COSMolKit 0.5.0 installed::

    python check_ez_cip.py --corpus ../benchmark_corpora/rdkit-js-browser-10k-v1.smi \\
        --chematic ../../target/debug/examples/ez_cip_dump
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import defaultdict
from pathlib import Path

import cosmolkit
from rdkit import Chem, RDLogger


def score(expected: dict, actual: dict) -> dict[str, int]:
    return {
        "matched": sum(actual.get(key) == value for key, value in expected.items()),
        "missing": sum(key not in actual for key in expected),
        "wrong": sum(key in actual and actual[key] != value for key, value in expected.items()),
        "extra": sum(key not in expected for key in actual),
    }


def display_labels(labels: dict[tuple[int, int], str]) -> dict[str, str]:
    return {f"{begin}-{end}": code for (begin, end), code in labels.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--chematic", type=Path, required=True)
    args = parser.parse_args()
    if Chem.rdBase.rdkitVersion != "2025.09.3" or cosmolkit.__version__ != "0.5.0":
        raise SystemExit("requires RDKit 2025.09.3 and COSMolKit 0.5.0")
    RDLogger.DisableLog("rdApp.*")
    corpus, binary = args.corpus.resolve(), args.chematic.resolve()
    completed = subprocess.run(
        [str(binary), str(corpus)], text=True, capture_output=True, check=True
    )
    chematic: dict[tuple[int, str], dict[tuple[int, int], str]] = defaultdict(dict)
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[1] == "ERROR":
            raise RuntimeError(f"CheMatic parse failed on row {fields[0]}")
        if len(fields) != 5:
            raise RuntimeError(f"invalid CheMatic output: {line}")
        index, mode, begin, end, code = fields
        if mode not in ("LegacyFast", "Accurate") or code not in ("E", "Z"):
            raise RuntimeError(f"invalid CheMatic label: {line}")
        key = tuple(sorted((int(begin), int(end))))
        chematic[int(index), mode][key] = code

    totals = {
        name: {"matched": 0, "missing": 0, "wrong": 0, "extra": 0}
        for name in ("chematic_default", "chematic_accurate", "cosmolkit")
    }
    eligible_rows = rows_with_ez = expected_bonds = 0
    mismatches: list[dict] = []
    for index, line in enumerate(corpus.read_text().splitlines()):
        tokens = line.split()
        if not tokens:
            continue
        source = tokens[0]
        rdkit = Chem.MolFromSmiles(source)
        if rdkit is None:
            continue
        Chem.rdCIPLabeler.AssignCIPLabels(rdkit)
        expected = {
            tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))): bond.GetProp(
                "_CIPCode"
            )
            for bond in rdkit.GetBonds()
            if bond.HasProp("_CIPCode") and bond.GetProp("_CIPCode") in ("E", "Z")
        }
        eligible_rows += 1
        rows_with_ez += bool(expected)
        expected_bonds += len(expected)
        other = cosmolkit.Molecule.from_smiles(source).with_cip_labels()
        cosmolkit_bonds = {
            tuple(sorted((bond.begin(), bond.end()))): code
            for bond in other.bonds()
            if (code := other.bond_property_string(bond.id(), "_CIPCode")) in ("E", "Z")
        }
        for name, actual in (
            ("chematic_default", chematic.get((index, "LegacyFast"), {})),
            ("chematic_accurate", chematic.get((index, "Accurate"), {})),
            ("cosmolkit", cosmolkit_bonds),
        ):
            counts = score(expected, actual)
            for key, count in counts.items():
                totals[name][key] += count
            if any(counts[key] for key in ("missing", "wrong", "extra")) and len(mismatches) < 20:
                mismatches.append({
                    "index": index,
                    "engine": name,
                    "expected": display_labels(expected),
                    "actual": display_labels(actual),
                })

    print(json.dumps({
        "corpus_sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(),
        "chematic_binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "cosmolkit_version": cosmolkit.__version__,
        "eligible_rows": eligible_rows,
        "rows_with_ez": rows_with_ez,
        "expected_bonds": expected_bonds,
        "totals": totals,
        "mismatches": mismatches,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
