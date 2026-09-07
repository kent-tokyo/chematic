#!/usr/bin/env python3
"""Collect one reproducible chematic streaming benchmark matrix.

The matrix is deliberately a Rust-runner report, not a cross-engine ranking.
It keeps file-backed ``BufRead`` formats separate from the materialized
one-shot formats and records compressed versus uncompressed input bytes.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORMATS = {
    "sdf": (Path("benchmarks/fixtures/streaming.sdf"), 2),
    "mol": (Path("benchmarks/fixtures/streaming.sdf"), 2),
    "xyz": (Path("benchmarks/fixtures/streaming.xyz"), 2),
    "v3000": (Path("benchmarks/fixtures/ethanol.v3000"), 1),
    "mol2": (Path("benchmarks/fixtures/ethanol.mol2"), 1),
    "cml": (Path("benchmarks/fixtures/ethanol.cml"), 1),
    "cdxml": (Path("benchmarks/fixtures/ethanol.cdxml"), 1),
    "mmcif": (Path("benchmarks/fixtures/minimal.mmcif"), 1),
    "pdb": (Path("benchmarks/fixtures/minimal.pdb"), 1),
}


def run(binary: str, fmt: str, path: Path, repeats: int, gzip_input: bool) -> dict[str, object]:
    command = [binary, "--format", fmt, "--path", str(path), "--repeats", str(repeats)]
    if gzip_input:
        command.append("--gzip")
    completed = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
    return json.loads(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary",
        default="target/debug/examples/streaming_benchmark",
        help="built streaming_benchmark binary",
    )
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--gzip", action="store_true", help="add gzip rows for every format")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="chematic-streaming-matrix-") as directory:
        temp = Path(directory)
        for fmt, (relative_fixture, records_per_pass) in FORMATS.items():
            fixture = ROOT / relative_fixture
            payload = fixture.read_bytes()
            digest = hashlib.sha256(payload).hexdigest()
            cases = [(False, fixture)]
            if args.gzip:
                compressed = temp / f"{fmt}.gz"
                with gzip.open(compressed, "wb") as output:
                    output.write(payload)
                cases.append((True, compressed))
            for compressed, path in cases:
                row = run(args.binary, fmt, path, args.repeats, compressed)
                expected_records = records_per_pass * args.repeats
                if row.get("records") != expected_records or row.get("failures") != 0:
                    errors.append(f"{fmt} {'gzip' if compressed else 'plain'} count mismatch: {row}")
                rows.append(
                    {
                        "format": fmt,
                        "compression": "gzip" if compressed else "none",
                        "fixture": {
                            "path": str(relative_fixture),
                            "bytes": len(payload),
                            "sha256": digest,
                            "compressed_bytes": path.stat().st_size if compressed else None,
                        },
                        "expected_records": expected_records,
                        "runner": row,
                    }
                )

    report = {
        "schema_version": 1,
        "target_version": "1.0.9",
        "repeats": args.repeats,
        "formats": list(FORMATS),
        "rows": rows,
        "comparison_boundary": (
            "chematic Rust runner only; file-backed BufRead for SDF/MOL/XYZ, "
            "materialized one-shot parsing for V3000/MOL2/CML/CDXML/mmCIF/PDB; "
            "gzip rows measure compressed input while parser limits apply after decompression"
        ),
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    if errors:
        print("streaming benchmark matrix failures:", *errors, sep="\n")
        return 1
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
