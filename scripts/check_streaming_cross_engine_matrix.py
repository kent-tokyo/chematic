#!/usr/bin/env python3
"""Check same-input record/failure accounting across installed engines.

This is deliberately a contract matrix, not a normalized speed benchmark.
Each row records the parser/process boundary because RDKit's Python block
constructors and Open Babel's per-repetition CLI are not equivalent to the
Rust runner's file-backed or materialized paths.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import time
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FORMATS = {
    "sdf": (Path("benchmarks/fixtures/streaming.sdf"), 2, "rdkit_sdf", "sdf"),
    "mol": (Path("benchmarks/fixtures/streaming.sdf"), 2, "rdkit_mol", "mol"),
    "xyz": (Path("benchmarks/fixtures/streaming.xyz"), 2, "rdkit_xyz", "xyz"),
    "extxyz": (Path("benchmarks/fixtures/streaming.extxyz"), 2, None, "xyz"),
    "v3000": (Path("benchmarks/fixtures/ethanol.v3000"), 1, "rdkit_v3000", "mol"),
    "mol2": (Path("benchmarks/fixtures/ethanol.mol2"), 1, "rdkit_mol2", "mol2"),
    "cml": (Path("benchmarks/fixtures/ethanol.cml"), 1, None, "cml"),
    "cdxml": (Path("benchmarks/fixtures/ethanol.cdxml"), 1, None, "cdxml"),
    "mmcif": (Path("benchmarks/fixtures/minimal.mmcif"), 1, None, "cif"),
    "pdb": (Path("benchmarks/fixtures/minimal.pdb"), 1, None, "pdb"),
}
CONVERTED = re.compile(r"^(\d+) molecule(?:s)? converted$")


def run_json(command: list[str]) -> dict[str, object]:
    completed = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
    return json.loads(completed.stdout)


def version(executable: str, flag: str = "-V") -> str:
    completed = subprocess.run([executable, flag], check=True, text=True, capture_output=True)
    return (completed.stdout or completed.stderr).strip().splitlines()[0]


def rdkit_version() -> str:
    completed = subprocess.run(
        ["python3", "-c", "import rdkit; print(rdkit.__version__)"],
        check=True,
        text=True,
        capture_output=True,
    )
    return completed.stdout.strip()


def openbabel_row(executable: str, input_format: str, path: Path, repeats: int) -> dict[str, object]:
    started = time.perf_counter()
    records = 0
    failures = 0
    for _ in range(repeats):
        completed = subprocess.run(
            [executable, f"-i{input_format}", str(path), "-osmi", "-O", "/dev/null"],
            check=False,
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        match = next(
            (CONVERTED.match(line.strip()) for line in completed.stderr.splitlines()),
            None,
        )
        if completed.returncode != 0 or match is None:
            failures += 1
        else:
            records += int(match.group(1))
    elapsed = time.perf_counter() - started
    return {
        "engine": "openbabel",
        "records": records,
        "failures": failures,
        "input_bytes": path.stat().st_size * repeats,
        "seconds": round(elapsed, 6),
        "records_per_second": round(records / elapsed, 2),
        "boundary": "Open Babel CLI conversion per repetition, including process startup",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="target/debug/examples/streaming_benchmark")
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument("--repeats", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.repeats <= 0:
        raise SystemExit("--repeats must be positive")

    rows: list[dict[str, object]] = []
    errors: list[str] = []
    for fmt, (relative_fixture, records_per_pass, rdkit_mode, openbabel_format) in FORMATS.items():
        fixture = ROOT / relative_fixture
        payload = fixture.read_bytes()
        expected = records_per_pass * args.repeats
        digest = hashlib.sha256(payload).hexdigest()
        chematic = run_json(
            [
                args.binary,
                "--format",
                fmt,
                "--path",
                str(fixture),
                "--repeats",
                str(args.repeats),
            ]
        )
        engines: list[dict[str, object]] = [
            {
                "engine": "chematic",
                "records": chematic["records"],
                "failures": chematic["failures"],
                "input_bytes": chematic["input_bytes"],
                "seconds": chematic["seconds"],
                "records_per_second": chematic["records_per_second"],
                "boundary": (
                    "Rust file-backed BufRead reader"
                    if chematic.get("execution_mode") == "file_backed_bufread"
                    else "Rust materialized one-shot parser"
                ),
            }
        ]
        if rdkit_mode:
            mode_path = "--xyz" if rdkit_mode == "rdkit_xyz" else "--sdf"
            mode_name = {
                "rdkit_sdf": "file-backed",
                "rdkit_mol": "mol-block",
                "rdkit_xyz": "xyz-file-backed",
                "rdkit_v3000": "v3000-block",
                "rdkit_mol2": "mol2-block",
            }[rdkit_mode]
            rdkit_rows = run_json(
                [
                    "python3",
                    "scripts/bench_streaming_formats.py",
                    "--mode",
                    mode_name,
                    mode_path,
                    str(fixture),
                    "--repeats",
                    str(args.repeats),
                ]
            )
            rdkit = rdkit_rows[0]
            engines.append(
                {
                    "engine": "rdkit",
                    "records": rdkit["records"],
                    "failures": rdkit["failures"],
                    "input_bytes": rdkit["input_bytes"],
                    "seconds": rdkit["seconds"],
                    "records_per_second": rdkit["records_per_second"],
                    "boundary": rdkit["comparison_boundary"],
                }
            )
        engines.append(openbabel_row(args.openbabel, openbabel_format, fixture, args.repeats))
        for engine in engines:
            if engine["records"] != expected or engine["failures"] != 0:
                errors.append(f"{fmt}/{engine['engine']} expected {expected}: {engine}")
        rows.append(
            {
                "format": fmt,
                "fixture": {
                    "path": str(relative_fixture),
                    "bytes": len(payload),
                    "sha256": digest,
                },
                "expected_records": expected,
                "engines": engines,
            }
        )

    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "repeats": args.repeats,
        "formats": list(FORMATS),
        "tool_versions": {"rdkit": rdkit_version(), "openbabel": version(args.openbabel)},
        "rows": rows,
        "comparison_boundary": "Record/failure contract only; parser/process boundaries are intentionally not normalized into a speed claim.",
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    if errors:
        print("streaming cross-engine matrix failures:", *errors, sep="\n")
        return 1
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
