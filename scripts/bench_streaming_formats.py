#!/usr/bin/env python3
"""Benchmark RDKit's SDF/MOL/XYZ parser lanes on P1 fixtures.

This is intentionally separate from the Rust file-backed runner: RDKit's
Python block API and its file-backed supplier are distinct operations. The
report therefore records the comparison boundary for each mode.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path

from rdkit import Chem, RDLogger


# The fixture intentionally contains a mixed 2D/3D coordinate edge; suppress
# its non-fatal diagnostic so benchmark output stays machine-readable.
RDLogger.DisableLog("rdApp.warning")


def sdf_blocks(text: str) -> list[str]:
    blocks = [block.lstrip("\r\n") for block in text.split("$$$$") if block.strip()]
    # RDKit's block constructor accepts the MOL portion, while SDMolSupplier
    # additionally handles data fields. Keep the parser comparison focused on
    # the same graph block used by chematic's fast streaming path.
    return [block.split("M  END", 1)[0] + "M  END\n" for block in blocks]


def xyz_frames(text: str) -> list[str]:
    lines = text.splitlines()
    frames: list[str] = []
    offset = 0
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        count = int(lines[offset].strip())
        end = offset + count + 2
        frames.append("\n".join(lines[offset:end]) + "\n")
        offset = end
    return frames


def mol2_blocks(text: str) -> list[str]:
    return ["@<TRIPOS>MOLECULE\n" + block.strip() + "\n" for block in text.split("@<TRIPOS>MOLECULE") if block.strip()]


def measure(label: str, blocks: list[str], repeats: int, source_bytes: int | None = None) -> dict[str, object]:
    started = time.perf_counter()
    records = 0
    for _ in range(repeats):
        for block in blocks:
            if label == "xyz":
                molecule = Chem.MolFromXYZBlock(block)
            elif label == "mol2":
                molecule = Chem.MolFromMol2Block(block, sanitize=False, cleanupSubstructures=False)
            else:
                molecule = Chem.MolFromMolBlock(block)
            if molecule is not None:
                records += 1
    elapsed = time.perf_counter() - started
    total_bytes = (source_bytes if source_bytes is not None else sum(len(block.encode()) for block in blocks)) * repeats
    return {
        "engine": "rdkit",
        "format": label,
        "repeats": repeats,
        "records": records,
        "failures": repeats * len(blocks) - records,
        "input_bytes": total_bytes,
        "seconds": round(elapsed, 6),
        "records_per_second": round(records / elapsed, 2),
        "bytes_per_second": round(total_bytes / elapsed, 2),
        "comparison_boundary": "Python block API, not file-backed BufRead streaming",
    }


def measure_rdkit_file_backed_sdf(path: Path, repeats: int) -> dict[str, object]:
    started = time.perf_counter()
    records = 0
    failures = 0
    for _ in range(repeats):
        with path.open("rb") as handle:
            supplier = Chem.ForwardSDMolSupplier(
                handle,
                sanitize=False,
                removeHs=False,
                strictParsing=True,
            )
            for molecule in supplier:
                if molecule is None:
                    failures += 1
                else:
                    records += 1
    elapsed = time.perf_counter() - started
    total_bytes = path.stat().st_size * repeats
    return {
        "engine": "rdkit",
        "format": "sdf",
        "mode": "file_backed_supplier",
        "repeats": repeats,
        "records": records,
        "failures": failures,
        "input_bytes": total_bytes,
        "seconds": round(elapsed, 6),
        "records_per_second": round(records / elapsed, 2),
        "bytes_per_second": round(total_bytes / elapsed, 2),
        "comparison_boundary": "RDKit ForwardSDMolSupplier in one Python process",
    }


def measure_rdkit_file_backed_xyz(path: Path, repeats: int) -> dict[str, object]:
    started = time.perf_counter()
    records = 0
    failures = 0
    payload = path.read_text()
    frames = xyz_frames(payload)
    for _ in range(repeats):
        for frame in frames:
            if Chem.MolFromXYZBlock(frame) is None:
                failures += 1
            else:
                records += 1
    elapsed = time.perf_counter() - started
    total_bytes = path.stat().st_size * repeats
    return {
        "engine": "rdkit",
        "format": "xyz",
        "mode": "file_backed_block_parser",
        "repeats": repeats,
        "records": records,
        "failures": failures,
        "input_bytes": total_bytes,
        "seconds": round(elapsed, 6),
        "records_per_second": round(records / elapsed, 2),
        "bytes_per_second": round(total_bytes / elapsed, 2),
        "comparison_boundary": "RDKit MolFromXYZBlock after Python frame splitting, not the Rust BufRead parser",
    }


def measure_openbabel_sdf(path: Path, repeats: int, executable: str) -> dict[str, object]:
    started = time.perf_counter()
    records = 0
    failures = 0
    converted_pattern = re.compile(r"^(\d+) molecules converted$")
    for _ in range(repeats):
        completed = subprocess.run(
            [executable, "-isdf", str(path), "-osmi", "-O", "/dev/null"],
            check=False,
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        match = next(
            (converted_pattern.match(line.strip()) for line in completed.stderr.splitlines()),
            None,
        )
        converted = int(match.group(1)) if match else 0
        records += converted
        if completed.returncode != 0 or not match:
            failures += 1
    elapsed = time.perf_counter() - started
    total_bytes = path.stat().st_size * repeats
    return {
        "engine": "openbabel",
        "format": "sdf",
        "mode": "file_backed_cli",
        "repeats": repeats,
        "records": records,
        "failures": failures,
        "input_bytes": total_bytes,
        "seconds": round(elapsed, 6),
        "records_per_second": round(records / elapsed, 2),
        "bytes_per_second": round(total_bytes / elapsed, 2),
        "comparison_boundary": "Open Babel CLI process per repetition; startup and conversion included",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--xyz", type=Path, default=Path("benchmarks/fixtures/streaming.xyz"))
    parser.add_argument("--repeats", type=int, default=200)
    parser.add_argument(
        "--mode",
        choices=("block", "mol-block", "v3000-block", "mol2-block", "file-backed", "xyz-file-backed"),
        default="block",
        help="RDKit block constructors, file-backed SDF suppliers, or XYZ blocks over a file",
    )
    parser.add_argument("--openbabel", help="also measure Open Babel's file-backed SDF CLI")
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    if args.mode in ("file-backed", "xyz-file-backed"):
        if args.mode == "xyz-file-backed":
            print(json.dumps([measure_rdkit_file_backed_xyz(args.xyz, args.repeats)], indent=2))
            return
        results: list[dict[str, object]] = [
            measure_rdkit_file_backed_sdf(args.sdf, args.repeats),
        ]
        if args.openbabel:
            results.append(measure_openbabel_sdf(args.sdf, args.repeats, args.openbabel))
        print(json.dumps(results, indent=2))
        return
    sdf = sdf_blocks(args.sdf.read_text())
    xyz = xyz_frames(args.xyz.read_text())
    if args.mode in ("mol-block", "v3000-block", "mol2-block"):
        if args.mode == "v3000-block":
            sdf = [args.sdf.read_text()]
            label = "mol"
            source = args.sdf.stat().st_size
        elif args.mode == "mol2-block":
            sdf = mol2_blocks(args.sdf.read_text())
            label = "mol2"
            source = args.sdf.stat().st_size
        else:
            label = "mol"
            source = args.sdf.stat().st_size
        print(json.dumps([measure(label, sdf, args.repeats, source)], indent=2))
        return
    print(json.dumps([
        measure("sdf", sdf, args.repeats),
        measure("mol", sdf, args.repeats),
        measure("xyz", xyz, args.repeats),
    ], indent=2))


if __name__ == "__main__":
    main()
