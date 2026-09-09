#!/usr/bin/env python3
"""Same-process schematic/RDKit V3000 MOL semantic contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from benchmark_version import workspace_version

ROOT = Path(__file__).resolve().parents[1]
MALFORMED = {
    "bad_counts": "bad\n  chematic\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS not-numbers\nM  V30 END CTAB\nM  END\n",
    "bad_atom": "bad\n  chematic\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 Xx 0 0 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n",
}


def parse_schematic(block: str) -> tuple[int, str | None, tuple[int, int] | None]:
    import chematic

    try:
        molecule = chematic.from_mol_v3000(block)
        return 1, molecule.smiles, (molecule.heavy_atoms, molecule.ring_count)
    except Exception:
        return 0, None, None


def parse_rdkit(block: str) -> tuple[int, str | None, tuple[int, int] | None]:
    from rdkit import Chem

    molecule = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    if molecule is None:
        return 0, None, None
    return 1, Chem.MolToSmiles(molecule, canonical=True), (molecule.GetNumHeavyAtoms(), molecule.GetRingInfo().NumRings())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mol", type=Path, default=Path("benchmarks/fixtures/ethanol.v3000"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.mol if args.mol.is_absolute() else ROOT / args.mol
    payload = path.read_bytes()
    block = payload.decode("utf-8")
    started = time.perf_counter()
    schematic_rows = [parse_schematic(block) for _ in range(args.repeats)]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [parse_rdkit(block) for _ in range(args.repeats)]
    rdkit_seconds = time.perf_counter() - started
    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    mismatches = sum(left[2] != right[2] for left, right in zip(schematic_rows, rdkit_rows))
    errors = []
    if schematic_records != args.repeats:
        errors.append("schematic record count mismatch")
    if rdkit_records != args.repeats:
        errors.append("RDKit record count mismatch")
    if mismatches:
        errors.append(f"canonical SMILES mismatch in {mismatches} repetitions")
    malformed = {}
    for case_id, malformed_block in MALFORMED.items():
        left = parse_schematic(malformed_block)
        right = parse_rdkit(malformed_block)
        malformed[case_id] = {"chematic": {"records": left[0], "failures": 1 - left[0]}, "rdkit": {"records": right[0], "failures": 1 - right[0]}}
        if left[0] != 0 or right[0] != 0:
            errors.append(f"malformed case {case_id} was accepted")
    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified" if not errors else "failed",
        "gate": "same_process_v3000_mol_semantic_contract",
        "fixture": {"path": str(path.relative_to(ROOT)), "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()},
        "repeats": args.repeats,
        "rows": {
            "chematic": {"records": schematic_records, "failures": args.repeats - schematic_records, "seconds": round(schematic_seconds, 6), "boundary": "schematic.from_mol_v3000 in the current process"},
            "rdkit": {"records": rdkit_records, "failures": args.repeats - rdkit_records, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolFromMolBlock in the current process"},
        },
        "comparison": {"expected_records": args.repeats, "canonical_smiles_mismatch_repetitions": mismatches, "malformed_cases": malformed, "speed_is_non_ranking_context": True},
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
