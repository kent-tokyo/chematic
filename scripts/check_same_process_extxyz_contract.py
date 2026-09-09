#!/usr/bin/env python3
"""Same-process schematic/RDKit Extended XYZ boundary contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from benchmark_version import workspace_version

ROOT = Path(__file__).resolve().parents[1]
MALFORMED = {
    "missing_atom": "2\nProperties=species:S:1:pos:R:3\nC 0 0 0\n",
    "bad_coordinate": "1\nProperties=species:S:1:pos:R:3\nC nope 0 0\n",
}


def frames(text: str) -> list[str]:
    lines = text.splitlines(keepends=True)
    result = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        count = int(lines[i].strip())
        end = i + 2 + count
        result.append("".join(lines[i:end]))
        i = end
    return result


def rdkit_xyz(frame: str) -> tuple[int, tuple[int, bool] | None]:
    from rdkit import Chem

    lines = frame.splitlines()
    count = int(lines[0])
    atoms = lines[2 : 2 + count]
    xyz = "\n".join([str(count), lines[1], *[" ".join(row.split()[:4]) for row in atoms]]) + "\n"
    molecule = Chem.MolFromXYZBlock(xyz)
    if molecule is None:
        return 0, None
    conformer = molecule.GetConformer()
    finite = all(math.isfinite(conformer.GetAtomPosition(i).x) and math.isfinite(conformer.GetAtomPosition(i).y) and math.isfinite(conformer.GetAtomPosition(i).z) for i in range(molecule.GetNumAtoms()))
    return 1, (molecule.GetNumAtoms(), finite)


def schematic_extxyz(frame: str) -> tuple[int, tuple[int, bool, bool] | None]:
    import chematic

    try:
        result = chematic.from_extxyz(frame)
        coords = result["coords"]
        finite = all(math.isfinite(value) for row in coords for value in row)
        return 1, (len(coords), finite, "energy" in result["info"])
    except Exception:
        return 0, None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--extxyz", type=Path, default=Path("benchmarks/fixtures/streaming.extxyz"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.extxyz if args.extxyz.is_absolute() else ROOT / args.extxyz
    payload = path.read_bytes()
    input_frames = frames(payload.decode("utf-8"))
    expected = len(input_frames) * args.repeats
    started = time.perf_counter()
    schematic_rows = [schematic_extxyz(frame) for _ in range(args.repeats) for frame in input_frames]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [rdkit_xyz(frame) for _ in range(args.repeats) for frame in input_frames]
    rdkit_seconds = time.perf_counter() - started
    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    mismatches = sum(left[1][:2] != right[1] for left, right in zip(schematic_rows, rdkit_rows))
    errors = []
    if schematic_records != expected or rdkit_records != expected:
        errors.append("valid-frame count mismatch")
    if mismatches:
        errors.append(f"frame signature mismatch in {mismatches} repetitions")
    malformed = {}
    for case_id, case in MALFORMED.items():
        left = schematic_extxyz(case)
        try:
            right = rdkit_xyz(case)
        except Exception:
            right = (0, None)
        malformed[case_id] = {"chematic": {"records": left[0], "failures": 1 - left[0]}, "rdkit": {"records": right[0], "failures": 1 - right[0]}}
        if left[0] != 0 or right[0] != 0:
            errors.append(f"malformed case {case_id} was accepted")
    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "status": "local-verified" if not errors else "failed",
        "gate": "same_process_extxyz_semantic_contract",
        "fixture": {"path": str(path.relative_to(ROOT)), "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()},
        "repeats": args.repeats,
        "input_frames_per_repeat": len(input_frames),
        "rows": {
            "chematic": {"records": schematic_records, "failures": expected - schematic_records, "seconds": round(schematic_seconds, 6), "boundary": "schematic.from_extxyz in the current process"},
            "rdkit": {"records": rdkit_records, "failures": expected - rdkit_records, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolFromXYZBlock over the coordinate projection in the current process"},
        },
        "comparison": {"expected_records": expected, "frame_signature_mismatch_repetitions": mismatches, "malformed_cases": malformed, "extended_properties": "schematic-only; RDKit coordinate projection does not claim property parity", "speed_is_non_ranking_context": True},
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
