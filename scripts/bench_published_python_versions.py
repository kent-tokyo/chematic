#!/usr/bin/env python3
"""Counterbalanced paired timing of two installed published Python wheels.

Each arm runs in its own isolated interpreter and reports output digest plus
process peak RSS. A block contains both arms, alternating order every block.
This compares versions, not RDKit; use bench_python_op_matrix_vs_rdkit.py's
--paired mode for the separate oracle comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

OPERATIONS = ("parse_smiles", "hba", "rdkit_ecfp4", "ring_count", "has_substructure")
MODES = ("parse_inclusive", "prepared", "precomputed")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def corpus_rows(path: Path, limit: int) -> list[str]:
    return [line.split()[0] for line in path.read_text(encoding="utf-8").splitlines() if line.strip()][:limit]


def evaluate(mol, operation: str):
    if operation == "hba":
        return str(mol.rdkit_hba)
    if operation == "rdkit_ecfp4":
        return bytes(mol.rdkit_ecfp4()).hex()
    if operation == "ring_count":
        return str(mol.ring_count)
    if operation == "has_substructure":
        return "1" if mol.has_substructure("c1ccccc1") else "0"
    raise ValueError(f"unsupported operation: {operation}")


def worker(args: argparse.Namespace) -> int:
    import chematic

    if chematic.__version__ != args.expected_version:
        raise RuntimeError(f"imported {chematic.__version__}, expected {args.expected_version}")
    rows = corpus_rows(args.corpus, args.limit)
    if not rows:
        raise RuntimeError("empty corpus")
    inputs = rows * args.multiplier
    if args.mode == "parse_inclusive" and args.operation != "parse_smiles":
        def call_all():
            return [evaluate(chematic.from_smiles(smiles), args.operation) for smiles in inputs]
    elif args.operation == "parse_smiles":
        if args.mode != "parse_inclusive":
            raise RuntimeError("parse_smiles only has parse_inclusive mode")
        def call_all():
            # Time parsing alone. Serializing back to SMILES would turn this
            # into a parse+canonical-writer lane with very different cost.
            return ["1" if chematic.from_smiles(smiles) is not None else "0" for smiles in inputs]
    else:
        molecules = [chematic.from_smiles(smiles) for smiles in inputs]
        if args.mode == "precomputed":
            [evaluate(mol, args.operation) for mol in molecules]
        def call_all():
            return [evaluate(mol, args.operation) for mol in molecules]
    # Exclude one warmup from the measured operation. Fresh processes prevent
    # version A's caches from leaking into B, but repeated prepared calls are
    # labeled as precomputed rather than cold.
    call_all()
    if args.mode == "prepared" and args.operation != "parse_smiles":
        # Rebuild inputs after warmup so this lane times first-use perception.
        molecules = [chematic.from_smiles(smiles) for smiles in inputs]
    start = time.perf_counter_ns()
    values = call_all()
    elapsed = time.perf_counter_ns() - start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = peak if sys.platform == "darwin" else peak * 1024
    digest = hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()
    print(json.dumps({"version": chematic.__version__, "operation": args.operation,
                      "mode": args.mode, "input_count": len(inputs), "elapsed_ns": elapsed,
                      "output_sha256": digest, "peak_process_rss_bytes": peak_bytes,
                      "runtime": platform.platform()}))
    return 0


def ci(ratios: list[float]) -> dict:
    rng = random.Random(0)
    n = len(ratios)
    medians = sorted(statistics.median(ratios[rng.randrange(n)] for _ in range(n)) for _ in range(4000))
    return {"median_b_over_a": statistics.median(ratios),
            "bootstrap_95pct": [medians[99], medians[3899]], "bootstrap_samples": 4000, "seed": 0}


def parent(args: argparse.Namespace) -> int:
    if args.blocks < 20:
        raise ValueError("at least 20 paired blocks are required")
    if args.multiplier < 1 or args.limit < 1:
        raise ValueError("limit and multiplier must be positive")
    if args.operation == "parse_smiles" and args.mode != "parse_inclusive":
        raise ValueError("parse_smiles only has parse_inclusive mode")
    artifacts = {"a": {"version": args.version_a, "python": str(args.python_a),
                       "wheel_sha256": sha256(args.wheel_a)},
                 "b": {"version": args.version_b, "python": str(args.python_b),
                       "wheel_sha256": sha256(args.wheel_b)}}
    pairs = []
    for block in range(args.blocks):
        order = ("a", "b") if block % 2 == 0 else ("b", "a")
        values = {}
        for arm in order:
            env = artifacts[arm]
            command = [env["python"], __file__, "--worker", "--corpus", str(args.corpus),
                       "--limit", str(args.limit), "--multiplier", str(args.multiplier),
                       "--operation", args.operation, "--mode", args.mode,
                       "--expected-version", env["version"]]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            values[arm] = json.loads(completed.stdout)
        pairs.append({"block": block, "order": list(order), **values})
        print(f"block {block + 1}/{args.blocks}: {values['a']['elapsed_ns']} vs {values['b']['elapsed_ns']} ns", flush=True)
    ratios = [pair["a"]["elapsed_ns"] / pair["b"]["elapsed_ns"] for pair in pairs]
    digests = {arm: sorted({pair[arm]["output_sha256"] for pair in pairs}) for arm in ("a", "b")}
    report = {"schema": "paired-published-python-versions/v1", "artifacts": artifacts,
              "corpus": {"sha256": sha256(args.corpus), "limit": args.limit, "multiplier": args.multiplier},
              "operation": args.operation, "mode": args.mode,
              "protocol": {"blocks": args.blocks, "order": "AB/BA alternating, fresh process per arm",
                           "warmup": "one operation call per arm outside timer", "memory": "whole-process peak RSS"},
              "paired": ci(ratios), "output_sha256_sets": digests,
              "outputs_equal_across_versions": len(digests["a"]) == len(digests["b"]) == 1 and digests["a"] == digests["b"],
              "pairs": pairs}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"paired": report["paired"], "outputs_equal": report["outputs_equal_across_versions"]}))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--python-a", type=Path)
    parser.add_argument("--python-b", type=Path)
    parser.add_argument("--version-a")
    parser.add_argument("--version-b")
    parser.add_argument("--wheel-a", type=Path)
    parser.add_argument("--wheel-b", type=Path)
    parser.add_argument("--expected-version")
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--multiplier", type=int, default=10)
    parser.add_argument("--operation", choices=OPERATIONS, required=True)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker:
        return worker(args)
    for key in ("python_a", "python_b", "version_a", "version_b", "wheel_a", "wheel_b", "output"):
        if getattr(args, key) is None:
            parser.error(f"--{key.replace('_', '-')} is required")
    return parent(args)


if __name__ == "__main__":
    raise SystemExit(main())
