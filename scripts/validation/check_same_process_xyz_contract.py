#!/usr/bin/env python3
"""Same-process schematic/RDKit XYZ frame contract.

XYZ readers do not share a bond model here, so this gate compares frame
acceptance, total/ heavy-atom counts, and finite coordinate values only.
"""

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
    "missing_atom": "2\nmissing second atom\nC 0 0 0\n",
    "bad_coordinate": "1\nnon-numeric coordinate\nC nope 0 0\n",
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


def schematic(frame: str) -> tuple[int, tuple[int, int] | None, bool]:
    import chematic

    try:
        molecule, coords = chematic.from_xyz(frame)
        finite = all(math.isfinite(value) for row in coords for value in row)
        return 1, (len(coords), molecule.heavy_atoms), finite
    except Exception:
        return 0, None, False


def rdkit(frame: str) -> tuple[int, tuple[int, int] | None, bool]:
    from rdkit import Chem

    molecule = Chem.MolFromXYZBlock(frame)
    if molecule is None:
        return 0, None, False
    conformer = molecule.GetConformer()
    finite = all(math.isfinite(conformer.GetAtomPosition(i).x) and math.isfinite(conformer.GetAtomPosition(i).y) and math.isfinite(conformer.GetAtomPosition(i).z) for i in range(molecule.GetNumAtoms()))
    return 1, (molecule.GetNumAtoms(), molecule.GetNumHeavyAtoms()), finite


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xyz", type=Path, default=Path("benchmarks/fixtures/streaming.xyz"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    path = args.xyz if args.xyz.is_absolute() else ROOT / args.xyz
    payload = path.read_bytes()
    input_frames = frames(payload.decode("utf-8"))
    expected = len(input_frames) * args.repeats
    started = time.perf_counter()
    schematic_rows = [schematic(frame) for _ in range(args.repeats) for frame in input_frames]
    schematic_seconds = time.perf_counter() - started
    started = time.perf_counter()
    rdkit_rows = [rdkit(frame) for _ in range(args.repeats) for frame in input_frames]
    rdkit_seconds = time.perf_counter() - started
    schematic_records = sum(row[0] for row in schematic_rows)
    rdkit_records = sum(row[0] for row in rdkit_rows)
    signature_mismatches = sum(left[1] != right[1] or left[2] != right[2] for left, right in zip(schematic_rows, rdkit_rows))
    errors = []
    if schematic_records != expected or rdkit_records != expected:
        errors.append("valid-frame count mismatch")
    if signature_mismatches:
        errors.append(f"frame signature mismatch in {signature_mismatches} repetitions")
    malformed = {}
    for case_id, case in MALFORMED.items():
        left, right = schematic(case), rdkit(case)
        malformed[case_id] = {"chematic": {"records": left[0], "failures": 1 - left[0]}, "rdkit": {"records": right[0], "failures": 1 - right[0]}}
        if left[0] != 0 or right[0] != 0:
            errors.append(f"malformed case {case_id} was accepted")
    report = make_report(
        ROOT, path, payload, args.repeats,
        "same_process_xyz_semantic_contract",
        rows={
            "chematic": {"records": schematic_records, "failures": expected - schematic_records, "seconds": round(schematic_seconds, 6), "boundary": "schematic.from_xyz in the current process"},
            "rdkit": {"records": rdkit_records, "failures": expected - rdkit_records, "seconds": round(rdkit_seconds, 6), "boundary": "RDKit MolFromXYZBlock in the current process"},
        },
        comparison={"expected_records": expected, "frame_signature_mismatch_repetitions": signature_mismatches, "malformed_cases": malformed, "bond_parity": "not claimed", "speed_is_non_ranking_context": True},
        errors=errors,
        extra={"input_frames_per_repeat": len(input_frames)},
    )
    return write_report(ROOT, args.output, report)


if __name__ == "__main__":
    raise SystemExit(main())
