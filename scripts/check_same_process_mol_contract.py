#!/usr/bin/env python3
"""Compare schematic and RDKit V2000 MOL parsing in one Python process.

The input is an SDF fixture split at record delimiters. This is a semantic
contract gate for materialized MOL blocks; timings are recorded only as
context and must not be used as a cross-engine speed ranking.
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
    "bad_atom_line": "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\nnot an atom line\nM  END\n",
    "bad_coordinate": "broken\n  chematic\n\n  1  0  0  0  0  0            999 V2000\n  nope 0.0 0.0 C  0  0  0  0  0  0  0  0  0  0\nM  END\n",
}


def blocks(payload: str) -> list[str]:
    return [block.lstrip("\r\n") for block in payload.split("$$$$") if block.strip()]


def parse_schematic(block: str) -> tuple[int, list[str]]:
    import chematic

    try:
        molecule = chematic.from_mol_block(block)
    except Exception:
        return 0, []
    return 1, [molecule.smiles]


def parse_rdkit(block: str) -> tuple[int, list[str]]:
    from rdkit import Chem

    molecule = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    if molecule is None:
        return 0, []
    return 1, [Chem.MolToSmiles(molecule, canonical=True)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.sdf if args.sdf.is_absolute() else ROOT / args.sdf
    payload = path.read_bytes()
    text = payload.decode("utf-8")
    input_blocks = blocks(text)
    expected_records = len(input_blocks) * args.repeats

    started = time.perf_counter()
    schematic_rows = [parse_schematic(block) for _ in range(args.repeats) for block in input_blocks]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [parse_rdkit(block) for _ in range(args.repeats) for block in input_blocks]
    rdkit_seconds = time.perf_counter() - started
    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    schematic_failures = expected_records - schematic_records
    rdkit_failures = expected_records - rdkit_records
    mismatches = sum(
        1 for left, right in zip(schematic_rows, rdkit_rows) if left[1] != right[1]
    )
    errors = []
    if (schematic_records, schematic_failures) != (expected_records, 0):
        errors.append("schematic valid-record/failure count mismatch")
    if (rdkit_records, rdkit_failures) != (expected_records, 0):
        errors.append("RDKit valid-record/failure count mismatch")
    if mismatches:
        errors.append(f"canonical SMILES mismatch in {mismatches} block repetitions")

    malformed = {}
    for case_id, case in MALFORMED.items():
        left = parse_schematic(case)
        right = parse_rdkit(case)
        malformed[case_id] = {
            "chematic": {"records": left[0], "failures": 1 - left[0]},
            "rdkit": {"records": right[0], "failures": 1 - right[0]},
        }
        if left[0] != 0 or right[0] != 0:
            errors.append(f"malformed case {case_id} was accepted")

    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified" if not errors else "failed",
        "gate": "same_process_v2000_mol_semantic_contract",
        "fixture": {"path": str(path.relative_to(ROOT)), "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()},
        "repeats": args.repeats,
        "input_blocks_per_repeat": len(input_blocks),
        "rows": {
            "chematic": {"records": schematic_records, "failures": schematic_failures, "seconds": round(schematic_seconds, 6), "boundary": "schematic.from_mol_block in the current process"},
            "rdkit": {"records": rdkit_records, "failures": rdkit_failures, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolFromMolBlock in the current process"},
        },
        "comparison": {"expected_records": expected_records, "canonical_smiles_mismatch_block_repetitions": mismatches, "malformed_cases": malformed, "speed_is_non_ranking_context": True},
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
