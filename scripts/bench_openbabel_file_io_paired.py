#!/usr/bin/env python3
"""Paired CLI round-trip timings for semantically verified format fixtures."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import tempfile
import time
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


def command_for_chematic(probe: Path, fmt: str, source: Path, output: Path) -> list[str]:
    return [
        str(probe),
        "--format",
        fmt,
        "--path",
        str(source),
        "--rewrite-output",
        str(output),
    ]


def command_for_openbabel(executable: str, fmt: str, source: Path, output: Path) -> list[str]:
    if fmt == "v3000":
        return [executable, "-imol", str(source), "-omol", "-x3", "-O", str(output)]
    return [executable, f"-i{fmt}", str(source), f"-o{fmt}", "-O", str(output)]


def run_iterations(command: list[str], iterations: int) -> int:
    started = time.perf_counter_ns()
    for _ in range(iterations):
        completed = subprocess.run(
            command,
            cwd=ROOT,
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"benchmark command failed: {' '.join(command)}")
    return time.perf_counter_ns() - started


def ci95_lower(values: list[float]) -> float:
    if len(values) < 2:
        return values[0]
    # n >= 21 is enforced; t(0.975, 20) is conservative for larger samples.
    margin = 2.086 * statistics.stdev(values) / math.sqrt(len(values))
    return statistics.mean(values) - margin


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


def executable_version(executable: str) -> str:
    completed = subprocess.run(
        [executable, "-V"], check=True, text=True, capture_output=True
    )
    output = (completed.stdout or completed.stderr).strip()
    return output.splitlines()[0] if output else "unknown"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--openbabel", default="obabel")
    parser.add_argument(
        "--probe",
        type=Path,
        default=ROOT / "target/release/examples/file_io_semantic_probe",
    )
    parser.add_argument("--blocks", type=int, default=21)
    parser.add_argument("--iterations-per-block", type=int, default=10)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.blocks < 21:
        parser.error("--blocks must be at least 21")
    if args.iterations_per_block < 1:
        parser.error("--iterations-per-block must be positive")
    if not args.probe.exists():
        parser.error(f"semantic probe does not exist: {args.probe}")

    formats: dict[str, Any] = {}
    errors: list[str] = []
    with tempfile.TemporaryDirectory(prefix="chematic-openbabel-speed-") as temp:
        temp_dir = Path(temp)
        for fmt, source in FIXTURES.items():
            chematic_command = command_for_chematic(
                args.probe, fmt, source, temp_dir / f"chematic.{fmt}"
            )
            openbabel_command = command_for_openbabel(
                args.openbabel, fmt, source, temp_dir / f"openbabel.{fmt}"
            )
            # Warm both paths once before entering the paired sequence.
            run_iterations(chematic_command, 1)
            run_iterations(openbabel_command, 1)
            rows: list[dict[str, Any]] = []
            for block in range(args.blocks):
                order = (
                    (("chematic", chematic_command), ("openbabel", openbabel_command))
                    if block % 2 == 0
                    else (("openbabel", openbabel_command), ("chematic", chematic_command))
                )
                elapsed: dict[str, int] = {}
                for engine, command in order:
                    elapsed[engine] = run_iterations(command, args.iterations_per_block)
                rows.append(
                    {
                        "block": block + 1,
                        "order": [engine for engine, _ in order],
                        "chematic_ns": elapsed["chematic"],
                        "openbabel_ns": elapsed["openbabel"],
                        "speedup": elapsed["openbabel"] / elapsed["chematic"],
                    }
                )
            speedups = [float(row["speedup"]) for row in rows]
            lower = ci95_lower(speedups)
            every_block_faster = all(speedup > 1.0 for speedup in speedups)
            passed = every_block_faster and lower > 1.0
            if not passed:
                errors.append(
                    f"{fmt}: every_block_faster={every_block_faster}, ci95_lower={lower:.4f}"
                )
            formats[fmt] = {
                "fixture": str(source.relative_to(ROOT)),
                "blocks": rows,
                "median_chematic_ns": statistics.median(
                    int(row["chematic_ns"]) for row in rows
                ),
                "median_openbabel_ns": statistics.median(
                    int(row["openbabel_ns"]) for row in rows
                ),
                "mean_speedup": statistics.mean(speedups),
                "median_speedup": statistics.median(speedups),
                "paired_speedup_ci95_lower_bound": lower,
                "every_block_faster": every_block_faster,
                "status": "faster" if passed else "not_proven_faster",
            }

    report = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "source": source_provenance(),
        "gate": "openbabel_file_io_paired_cli_roundtrip_v1",
        "status": "local-verified" if not errors else "failed",
        "lane": "fresh CLI process, parse plus same-format write, process startup included",
        "blocks": args.blocks,
        "iterations_per_block": args.iterations_per_block,
        "alternating_order": True,
        "semantic_prerequisite": "scripts/check_openbabel_file_io_semantics.py",
        "formats": formats,
        "tool_versions": {"openbabel": executable_version(args.openbabel)},
        "claim_boundary": (
            "tiny checked-in fixture CLI round trips only; not parser-only, writer-only, "
            "large-file throughput, peak memory, or total-format superiority"
        ),
        "errors": errors,
    }
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        target = args.output if args.output.is_absolute() else ROOT / args.output
        target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
