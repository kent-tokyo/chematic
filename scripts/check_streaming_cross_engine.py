#!/usr/bin/env python3
"""Check record/failure agreement for one identical file-backed SDF or XYZ input.

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
    parser.add_argument("--format", choices=("sdf", "mol", "xyz"), default="sdf")
    parser.add_argument("--sdf", type=Path, default=Path("benchmarks/fixtures/streaming.sdf"))
    parser.add_argument("--xyz", type=Path, default=Path("benchmarks/fixtures/streaming.xyz"))
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--openbabel", default="obabel", help="Open Babel executable; used only for SDF")
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
    path = (args.sdf if args.format in ("sdf", "mol") else args.xyz)
    path = path if path.is_absolute() else ROOT / path
    payload = path.read_bytes()
    expected_records = 2 * args.repeats
    digest = hashlib.sha256(payload).hexdigest()
    chematic = run_json(
        [
            *args.binary,
            "--format",
            args.format,
            "--path",
            str(path),
            "--repeats",
            str(args.repeats),
        ]
    )
    comparator_command = [
        "python3",
        "scripts/bench_streaming_formats.py",
        "--mode",
        "file-backed" if args.format == "sdf" else ("mol-block" if args.format == "mol" else "xyz-file-backed"),
        "--sdf" if args.format in ("sdf", "mol") else "--xyz",
        str(path),
        "--repeats",
        str(args.repeats),
    ]
    if args.format == "sdf":
        comparator_command.extend(["--openbabel", args.openbabel])
    comparator = run_json(comparator_command)
    expected_comparators = 2 if args.format == "sdf" else 1
    if not isinstance(comparator, list) or len(comparator) != expected_comparators:
        raise SystemExit("comparison runner returned an unexpected row count")
    rows = {"chematic": chematic, "rdkit": comparator[0]}
    if args.format == "sdf":
        rows["openbabel"] = comparator[1]
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
        "format": args.format,
        "fixture": {"path": str(path.relative_to(ROOT)), "bytes": len(payload), "sha256": digest},
        "repeats": args.repeats,
        "expected_records": expected_records,
        "rows": rows,
        "comparison_boundary": {
            "chematic": f"Rust {args.format.upper()} file-backed BufRead reader",
            "rdkit": "RDKit Python block parser over blocks split from the identical file",
            **({"openbabel": "Open Babel CLI conversion per repetition, including process startup"} if args.format == "sdf" else {}),
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
