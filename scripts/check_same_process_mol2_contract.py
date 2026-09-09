#!/usr/bin/env python3
"""Same-process schematic/RDKit MOL2 semantic contract.

MOL2 charge typing is intentionally not compared as canonical SMILES: the
two readers may infer missing charges differently. The gate compares accepted
record counts plus heavy-atom and ring-count signatures and records that
charge/SMILES parity is outside this bounded contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from benchmark_version import workspace_version

ROOT = Path(__file__).resolve().parents[1]
MALFORMED = {
    "missing_sections": "@<TRIPOS>MOLECULE\nname\n1 1 0 0 0\n",
    "bad_atom": "@<TRIPOS>MOLECULE\nname\n1 0 0 0 0\n@<TRIPOS>ATOM\n1 Xx 0 0 0\n",
}


def schematic(block: str) -> tuple[int, tuple[int, int] | None]:
    import chematic

    try:
        molecule = chematic.from_mol2(block)
        return 1, (molecule.heavy_atoms, molecule.ring_count)
    except Exception:
        return 0, None


def rdkit(block: str) -> tuple[int, tuple[int, int] | None]:
    from rdkit import Chem

    molecule = Chem.MolFromMol2Block(block, sanitize=True, removeHs=False)
    if molecule is None:
        return 0, None
    return 1, (molecule.GetNumHeavyAtoms(), molecule.GetRingInfo().NumRings())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mol2", type=Path, default=Path("benchmarks/fixtures/ethanol.mol2"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.mol2 if args.mol2.is_absolute() else ROOT / args.mol2
    payload = path.read_bytes()
    block = payload.decode("utf-8")
    started = time.perf_counter()
    schematic_rows = [schematic(block) for _ in range(args.repeats)]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [rdkit(block) for _ in range(args.repeats)]
    rdkit_seconds = time.perf_counter() - started
    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    mismatches = sum(left[1] != right[1] for left, right in zip(schematic_rows, rdkit_rows))
    errors = []
    if schematic_records != args.repeats or rdkit_records != args.repeats:
        errors.append("valid-record count mismatch")
    if mismatches:
        errors.append(f"structural signature mismatch in {mismatches} repetitions")
    malformed = {}
    for case_id, case in MALFORMED.items():
        left, right = schematic(case), rdkit(case)
        malformed[case_id] = {"chematic": {"records": left[0], "failures": 1 - left[0]}, "rdkit": {"records": right[0], "failures": 1 - right[0]}}
        if left[0] != 0 or right[0] != 0:
            errors.append(f"malformed case {case_id} was accepted")
    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified" if not errors else "failed",
        "gate": "same_process_mol2_semantic_contract",
        "fixture": {"path": str(path.relative_to(ROOT)), "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()},
        "repeats": args.repeats,
        "rows": {
            "chematic": {"records": schematic_records, "failures": args.repeats - schematic_records, "seconds": round(schematic_seconds, 6), "boundary": "schematic.from_mol2 in the current process"},
            "rdkit": {"records": rdkit_records, "failures": args.repeats - rdkit_records, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolFromMol2Block in the current process"},
        },
        "comparison": {"expected_records": args.repeats, "structural_signature_mismatch_repetitions": mismatches, "malformed_cases": malformed, "charge_and_smiles_parity": "not claimed", "speed_is_non_ranking_context": True},
        "tool_versions": {},
        "errors": errors,
    }
    from rdkit import rdBase

    report["tool_versions"]["rdkit"] = rdBase.rdkitVersion
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
