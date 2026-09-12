#!/usr/bin/env python3
"""Run one pipeline arm per subprocess with a real wall-clock timeout.

PipelineV2's cooperative timeout is useful for normal callers, but a single
stage can run past a checkpoint.  This runner is deliberately an external
benchmark boundary: every manifest row is retained, including rows whose
child process is killed by the hard timeout.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFESTS = {
    "A": ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_a.json",
    "B": ROOT / "validation/manifests/pipeline_v2_vs_rdkit_etkdgv3_tier_b.json",
}


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
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/debug/examples/pipeline_v2_vs_rdkit_dump")
    parser.add_argument("--arm", default="chematic_pipeline_v2_mmff94_strict")
    parser.add_argument("--timeout-s", type=float, default=30.0)
    parser.add_argument("--tier", choices=sorted(MANIFESTS), action="append")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.timeout_s <= 0:
        parser.error("--timeout-s must be positive")
    if not args.binary.is_file():
        parser.error(f"pipeline binary does not exist: {args.binary}")
    tiers = args.tier or list(MANIFESTS)
    with args.output.open("w", encoding="utf-8") as handle:
        for tier in tiers:
            molecules = json.loads(MANIFESTS[tier].read_text(encoding="utf-8"))["molecules"]
            for index in range(len(molecules)):
                row = run_one(args.binary, tier, index, args.arm, args.timeout_s)
                handle.write(json.dumps(row, sort_keys=True) + "\n")
                handle.flush()
                print(f"{tier} {index + 1}/{len(molecules)} {row['name']} {row['status']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
