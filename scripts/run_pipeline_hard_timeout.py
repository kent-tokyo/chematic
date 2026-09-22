#!/usr/bin/env python3
"""Run one pipeline arm per subprocess with a real wall-clock timeout.

PipelineV2's cooperative timeout is useful for normal callers, but a single
stage can run past a checkpoint.  This runner is deliberately an external
benchmark boundary: every manifest row is retained, including rows whose
child process is killed by the hard timeout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFESTS = {
    "A": ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    "B": ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command: list[str]) -> str | None:
    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.stdout.strip() if completed.returncode == 0 else None


def run_one(binary: Path, tier: str, index: int, arm: str, timeout_s: float) -> dict:
    manifest = json.loads(MANIFESTS[tier].read_text(encoding="utf-8"))
    molecule = manifest["molecules"][index]
    command = [
        str(binary),
        "--tier",
        tier,
        "--start",
        str(index),
        "--count",
        "1",
        "--only-arm",
        arm,
    ]
    try:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "tier": tier,
            "name": molecule["name"],
            "smiles": molecule["smiles"],
            "primary_category": molecule.get("primary_category", "unknown"),
            "arm": arm,
            "status": "timeout",
            "external_timeout_s": timeout_s,
            "failure_cause": "external_hard_timeout",
        }
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    if completed.returncode != 0 or len(lines) != 1:
        return {
            "tier": tier,
            "name": molecule["name"],
            "smiles": molecule["smiles"],
            "primary_category": molecule.get("primary_category", "unknown"),
            "arm": arm,
            "status": "runner_failure",
            "returncode": completed.returncode,
            "stdout_lines": len(lines),
            "stderr_tail": completed.stderr[-500:],
        }
    try:
        row = json.loads(lines[0])
    except json.JSONDecodeError as exc:
        return {
            "tier": tier,
            "name": molecule["name"],
            "smiles": molecule["smiles"],
            "primary_category": molecule.get("primary_category", "unknown"),
            "arm": arm,
            "status": "runner_failure",
            "failure_cause": f"invalid_child_json:{exc}",
        }
    row["external_timeout_s"] = timeout_s
    row["original_index"] = index
    return row


def build_metadata(binary: Path, arm: str, timeout_s: float, tiers: list[str]) -> dict:
    revision = command_output(["git", "rev-parse", "HEAD"])
    dirty = command_output(["git", "status", "--short"])
    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_revision": revision,
        "source_dirty": dirty != "",
        "source_dirty_paths": dirty.splitlines() if dirty else [],
        "arm": arm,
        "external_timeout_s": timeout_s,
        "tiers": tiers,
        "binary": {
            "path": str(binary.resolve()),
            "sha256": sha256_file(binary),
        },
        "manifests": {
            tier: {
                "path": str(MANIFESTS[tier].relative_to(ROOT)),
                "sha256": sha256_file(MANIFESTS[tier]),
                "row_count": len(
                    json.loads(MANIFESTS[tier].read_text(encoding="utf-8"))["molecules"]
                ),
            }
            for tier in tiers
        },
        "runtime": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
            "rustc": command_output(["rustc", "--version"]),
        },
        "environment": {
            "SCHEMATIC_MMFF94_MAX_ITERATIONS": os.environ.get(
                "SCHEMATIC_MMFF94_MAX_ITERATIONS"
            ),
        },
        "command": {
            "binary": str(binary),
            "arm": arm,
            "timeout_s": timeout_s,
            "tiers": tiers,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/debug/examples/pipeline_v2_vs_rdkit_dump")
    parser.add_argument("--arm", default="chematic_pipeline_v2_mmff94_strict")
    parser.add_argument("--timeout-s", type=float, default=30.0)
    parser.add_argument("--tier", choices=sorted(MANIFESTS), action="append")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--metadata-output",
        type=Path,
        help="optional provenance sidecar; JSONL row schema remains unchanged",
    )
    args = parser.parse_args()
    if args.timeout_s <= 0:
        parser.error("--timeout-s must be positive")
    if not args.binary.is_file():
        parser.error(f"pipeline binary does not exist: {args.binary}")
    tiers = args.tier or list(MANIFESTS)
    if args.metadata_output:
        metadata = build_metadata(args.binary, args.arm, args.timeout_s, tiers)
        args.metadata_output.write_text(
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    with args.output.open("w", encoding="utf-8") as handle:
        for tier in tiers:
            molecules = json.loads(MANIFESTS[tier].read_text(encoding="utf-8"))["molecules"]
            for index in range(len(molecules)):
                row = run_one(args.binary, tier, index, args.arm, args.timeout_s)
                row["original_index"] = index
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                handle.flush()
                print(f"{tier} {index + 1}/{len(molecules)} {row['name']} {row['status']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
