#!/usr/bin/env python3
"""Paired same-process parse/write/round-trip timings against Open Babel."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Any

try:
    from benchmark_version import workspace_version
except ModuleNotFoundError:  # imported as scripts.* by pytest
    from scripts.benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = {
    "v3000": ROOT / "benchmarks/fixtures/ethanol.v3000",
    "mol2": ROOT / "benchmarks/fixtures/ethanol.mol2",
    "cml": ROOT / "benchmarks/fixtures/ethanol.cml",
    "cdxml": ROOT / "benchmarks/fixtures/ethanol.cdxml",
}
OPERATIONS = ("parse", "write", "roundtrip")


def ci95_lower(values: list[float]) -> float:
    if len(values) < 2:
        return values[0]
    return statistics.mean(values) - 2.086 * statistics.stdev(values) / math.sqrt(len(values))


def run_row(binary: Path, fmt: str, fixture: Path, operation: str, repeats: int) -> dict[str, Any]:
    completed = subprocess.run(
        [
            str(binary),
            "--format",
            fmt,
            "--path",
            str(fixture),
            "--benchmark-operation",
            operation,
            "--repeats",
            str(repeats),
        ],
        cwd=ROOT,
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"benchmark failed: {completed.stderr.strip()}")
    row = json.loads(completed.stdout)
    if row.get("records") != repeats:
        raise RuntimeError(f"unexpected record count: {row}")
    return row


def source_provenance() -> dict[str, Any]:
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout.strip()
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=normal"],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    ).stdout
    return {"commit": revision, "dirty": bool(status.strip())}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chematic-harness",
        type=Path,
        default=ROOT / "target/release/examples/file_io_semantic_probe",
    )
    parser.add_argument(
        "--openbabel-harness",
        type=Path,
        default=ROOT / "target/release/openbabel_file_io_harness",
    )
    parser.add_argument("--blocks", type=int, default=21)
    parser.add_argument("--repeats", type=int, default=500)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.blocks < 21:
        parser.error("--blocks must be at least 21")
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    for harness in (args.chematic_harness, args.openbabel_harness):
        if not harness.exists():
            parser.error(f"harness does not exist: {harness}")

    lanes: dict[str, Any] = {}
    errors: list[str] = []
    for fmt, fixture in FIXTURES.items():
        for operation in OPERATIONS:
            name = f"{fmt}_{operation}"
            rows = []
            for block in range(args.blocks):
                order = (
                    (("chematic", args.chematic_harness), ("openbabel", args.openbabel_harness))
                    if block % 2 == 0
                    else (("openbabel", args.openbabel_harness), ("chematic", args.chematic_harness))
                )
                values: dict[str, dict[str, Any]] = {}
                for engine, harness in order:
                    values[engine] = run_row(harness, fmt, fixture, operation, args.repeats)
                chematic_ns = int(values["chematic"]["elapsed_ns"])
                openbabel_ns = int(values["openbabel"]["elapsed_ns"])
                rows.append(
                    {
                        "block": block + 1,
                        "order": [engine for engine, _ in order],
                        "chematic_ns": chematic_ns,
                        "openbabel_ns": openbabel_ns,
                        "speedup": openbabel_ns / chematic_ns,
                    }
                )
            speedups = [float(row["speedup"]) for row in rows]
            lower = ci95_lower(speedups)
            every = all(value > 1.0 for value in speedups)
            passed = every and lower > 1.0
            if not passed:
                errors.append(f"{name}: every={every}, ci95_lower={lower:.4f}")
            lanes[name] = {
                "format": fmt,
                "operation": operation,
                "fixture": str(fixture.relative_to(ROOT)),
                "blocks": rows,
                "median_chematic_ns": statistics.median(
                    int(row["chematic_ns"]) for row in rows
                ),
                "median_openbabel_ns": statistics.median(
                    int(row["openbabel_ns"]) for row in rows
                ),
                "median_speedup": statistics.median(speedups),
                "paired_speedup_ci95_lower_bound": lower,
                "every_block_faster": every,
                "status": "faster" if passed else "not_proven_faster",
            }

    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "source": source_provenance(),
        "gate": "openbabel_file_io_paired_same_process_v1",
        "status": "local-verified" if not errors else "diagnostic",
        "blocks": args.blocks,
        "repeats_per_process": args.repeats,
        "alternating_process_order": True,
        "timer_boundary": "internal steady-clock loop; process startup excluded",
        "semantic_prerequisite": "scripts/check_openbabel_file_io_semantics.py",
        "lanes": lanes,
        "claim_boundary": (
            "tiny checked-in fixtures; parser/write/round-trip hot loops only; "
            "not large-file throughput, peak memory, or published-package evidence"
        ),
        "errors": errors,
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
