#!/usr/bin/env python3
"""Benchmark RDKit's equivalent SDF/MOL/XYZ block parsing on P1 fixtures.

This is intentionally separate from the Rust file-backed runner: RDKit's
Python API parses a supplied block rather than exposing the same BufRead
iterator contract. The report therefore records the comparison boundary.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from rdkit import Chem


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


def measure(label: str, blocks: list[str], repeats: int) -> dict[str, object]:
    started = time.perf_counter()
    records = 0
    for _ in range(repeats):
        for block in blocks:
            if (Chem.MolFromXYZBlock(block) if label == "xyz" else Chem.MolFromMolBlock(block)) is not None:
                records += 1
    elapsed = time.perf_counter() - started
    total_bytes = sum(len(block.encode()) for block in blocks) * repeats
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--xyz", type=Path, default=Path("benchmarks/fixtures/streaming.xyz"))
    parser.add_argument("--repeats", type=int, default=200)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    sdf = sdf_blocks(args.sdf.read_text())
    xyz = xyz_frames(args.xyz.read_text())
    print(json.dumps([
        measure("sdf", sdf, args.repeats),
        measure("mol", sdf, args.repeats),
        measure("xyz", xyz, args.repeats),
    ], indent=2))


if __name__ == "__main__":
    main()
