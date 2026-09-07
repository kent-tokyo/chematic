#!/usr/bin/env python3
"""Check record/failure agreement for one identical file-backed SDF input.

This is a contract check, not a throughput ranking. The three engines use
different parser and process boundaries; the report keeps those boundaries
explicit while requiring the same valid-record and failure counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_json(command: list[str]) -> object:
    completed = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
    return json.loads(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument(
        "--binary",
        nargs="+",
        default=["cargo", "run", "-p", "chematic-mol", "--example", "streaming_benchmark", "--offline", "--"],
        help="chematic runner command before its format arguments",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")
    sdf = args.sdf if args.sdf.is_absolute() else ROOT / args.sdf
    payload = sdf.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    expected_records = 2 * args.repeats
    chematic = run_json(
        [
            *args.binary,
            "--format",
            "sdf",
            "--path",
            str(sdf),
            "--repeats",
            str(args.repeats),
        ]
    )
    comparator = run_json(
        [
            "python3",
            "scripts/bench_streaming_formats.py",
            "--mode",
            "file-backed",
            "--sdf",
            str(sdf),
            "--repeats",
            str(args.repeats),
            "--openbabel",
            args.openbabel,
        ]
    )
    if not isinstance(comparator, list) or len(comparator) != 2:
        raise SystemExit("comparison runner did not return RDKit and Open Babel rows")
    rdkit, openbabel = comparator
    rows = {"chematic": chematic, "rdkit": rdkit, "openbabel": openbabel}
    errors: list[str] = []
    for engine, row in rows.items():
        if not isinstance(row, dict):
            errors.append(f"{engine} row is not an object")
            continue
        if row.get("records") != expected_records or row.get("failures") != 0:
            errors.append(f"{engine} count mismatch: {row}")
        if row.get("input_bytes") != len(payload) * args.repeats:
            errors.append(f"{engine} input byte mismatch: {row.get('input_bytes')}")
    report = {
        "schema_version": 1,
        "target_version": "1.0.9",
        "fixture": {"path": str(sdf.relative_to(ROOT)), "bytes": len(payload), "sha256": digest},
        "repeats": args.repeats,
        "expected_records": expected_records,
        "rows": rows,
        "comparison_boundary": {
            "chematic": "Rust SdfFileReader over file-backed BufRead",
            "rdkit": "RDKit ForwardSDMolSupplier in one Python process",
            "openbabel": "Open Babel CLI conversion per repetition, including process startup",
        },
    }
    if errors:
        print("streaming cross-engine contract failures:", *errors, sep="\n", flush=True)
        return 1
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
