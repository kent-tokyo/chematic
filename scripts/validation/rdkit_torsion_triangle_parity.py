#!/usr/bin/env python3
"""Compare RDKit/CheMatic torsion codes and bits on substituted cyclopropanes.

RDKit 2026.09.1 is the fixed oracle. The enumerated small-ring structures
and three rooted SMILES spellings are deterministic; this is a development
stress set, not evidence of universal fingerprint agreement.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdFingerprintGenerator


def structures() -> list[str]:
    base = Chem.MolFromSmiles("C1CC1")
    substituents = (None, 6, 7, 8, 9, 17, 35)
    smiles = set()
    for atomic_numbers in itertools.product(substituents, repeat=3):
        rw = Chem.RWMol(base)
        for ring_atom, atomic_number in enumerate(atomic_numbers):
            if atomic_number is not None:
                substituent = rw.AddAtom(Chem.Atom(atomic_number))
                rw.AddBond(ring_atom, substituent, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        for root in range(3):
            smiles.add(Chem.MolToSmiles(mol, rootedAtAtom=root, canonical=False))
    return sorted(smiles)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.09.1":
        parser.error(f"expected RDKit 2026.09.1, found {rdBase.rdkitVersion}")
    RDLogger.DisableLog("rdApp.*")
    selected = structures()
    run = subprocess.run(
        [str(args.binary.resolve())], input="\n".join(selected) + "\n",
        capture_output=True, text=True, check=True,
    )
    lines = run.stdout.splitlines()
    if len(lines) != len(selected):
        raise ValueError("CheMatic output count differs from input count")
    generator = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
    mismatches = []
    for smiles, line in zip(selected, lines, strict=True):
        source, sparse_text, bits_text = line.split("\t")
        if source != smiles:
            raise ValueError("CheMatic output order differs from input order")
        mol = Chem.MolFromSmiles(smiles)
        expected_sparse = generator.GetSparseCountFingerprint(mol).GetNonzeroElements()
        expected_bits = list(generator.GetFingerprint(mol).GetOnBits())
        actual_sparse = dict(ast.literal_eval(sparse_text))
        actual_bits = ast.literal_eval(bits_text)
        if actual_sparse != expected_sparse or actual_bits != expected_bits:
            mismatches.append({
                "smiles": smiles,
                "sparse_equal": actual_sparse == expected_sparse,
                "bits_equal": actual_bits == expected_bits,
            })
    script = Path(__file__)
    report = {
        "schema_version": 1,
        "rdkit_version": rdBase.rdkitVersion,
        "selection": "All 7^3 choices from unsubstituted/C/N/O/F/Cl/Br at each cyclopropane carbon, each rooted at atoms 0, 1, and 2; unique SMILES sorted",
        "rows": len(selected),
        "sparse_and_bit_exact": len(selected) - len(mismatches),
        "mismatches": mismatches,
        "hashes": {
            "binary_sha256": hashlib.sha256(args.binary.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
            "selected_smiles_sha256": hashlib.sha256(("\n".join(selected) + "\n").encode()).hexdigest(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"{report['sparse_and_bit_exact']}/{report['rows']} sparse and bit exact")
    if mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
