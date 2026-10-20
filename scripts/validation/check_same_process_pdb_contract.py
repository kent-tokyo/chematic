#!/usr/bin/env python3
"""Same-process schematic/RDKit PDB semantic contract."""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.same_process_report import make_report, write_report

ROOT = Path(__file__).resolve().parents[2]
MALFORMED = {
    "bad_atom": "ATOM      1  BAD\n",
    "bad_coordinate": "ATOM      1  CA  ALA A   1       nope   0.000   0.000\n",
}


def schematic(text: str) -> tuple[int, tuple[int, int, bool] | None]:
    import chematic

    try:
        molecule, coords = chematic.from_pdb(text)
        finite = all(math.isfinite(value) for row in coords for value in row)
        return 1, (len(coords), molecule.heavy_atoms, finite)
    except Exception:
        return 0, None


def rdkit(text: str) -> tuple[int, tuple[int, int, bool] | None]:
    from rdkit import Chem

    molecule = Chem.MolFromPDBBlock(text, sanitize=True, removeHs=False)
    if molecule is None:
        return 0, None
    if molecule.GetNumConformers() == 0:
        return 0, None
    conformer = molecule.GetConformer()
    finite = all(math.isfinite(conformer.GetAtomPosition(i).x) and math.isfinite(conformer.GetAtomPosition(i).y) and math.isfinite(conformer.GetAtomPosition(i).z) for i in range(molecule.GetNumAtoms()))
    return 1, (molecule.GetNumAtoms(), molecule.GetNumHeavyAtoms(), finite)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", type=Path, default=Path("benchmarks/fixtures/minimal.pdb"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.pdb if args.pdb.is_absolute() else ROOT / args.pdb
    payload = path.read_bytes()
    text = payload.decode("utf-8")
    started = time.perf_counter()
    schematic_rows = [schematic(text) for _ in range(args.repeats)]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [rdkit(text) for _ in range(args.repeats)]
    rdkit_seconds = time.perf_counter() - started
    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    mismatches = sum(left[1] != right[1] for left, right in zip(schematic_rows, rdkit_rows))
    errors = []
    if schematic_records != args.repeats or rdkit_records != args.repeats:
        errors.append("valid-record count mismatch")
    if mismatches:
        errors.append(f"PDB signature mismatch in {mismatches} repetitions")
    malformed = {}
    malformed_rejection_mismatches = []
    for case_id, case in MALFORMED.items():
        left, right = schematic(case), rdkit(case)
        malformed[case_id] = {"chematic": {"records": left[0], "failures": 1 - left[0]}, "rdkit": {"records": right[0], "failures": 1 - right[0]}}
        if left[0] != right[0]:
            malformed_rejection_mismatches.append(case_id)
    report = make_report(
        ROOT, path, payload, args.repeats,
        "same_process_pdb_semantic_contract",
        rows={
            "chematic": {"records": schematic_records, "failures": args.repeats - schematic_records, "seconds": round(schematic_seconds, 6), "boundary": "schematic.from_pdb in the current process"},
            "rdkit": {"records": rdkit_records, "failures": args.repeats - rdkit_records, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolFromPDBBlock in the current process"},
        },
        comparison={"expected_records": args.repeats, "signature": "total atoms, heavy atoms, finite coordinates", "signature_mismatch_repetitions": mismatches, "malformed_cases": malformed, "malformed_rejection_mismatches": malformed_rejection_mismatches, "malformed_contract": "reported only; schematic PDB parser is intentionally lenient", "bond_parity": "not claimed", "speed_is_non_ranking_context": True},
        errors=errors,
    )
    return write_report(ROOT, args.output, report)


if __name__ == "__main__":
    raise SystemExit(main())
