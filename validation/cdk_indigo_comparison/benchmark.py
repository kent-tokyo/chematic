#!/usr/bin/env python3
"""Alternating-order persistent-process benchmark for equivalent operations."""

from __future__ import annotations

import argparse
import json
import random
import shlex
import statistics
import subprocess
from pathlib import Path


def adapter_arg(value: str) -> tuple[str, str]:
    engine, separator, command = value.partition("=")
    if not separator or not engine or not command:
        raise argparse.ArgumentTypeError("adapter must use ENGINE=COMMAND")
    return engine, command


def ratio_ci(numerator: list[int], denominator: list[int], seed: int) -> dict:
    ratios = [left / right for left, right in zip(numerator, denominator)]
    rng = random.Random(seed)
    boot = [statistics.median(ratios[rng.randrange(len(ratios))] for _ in ratios) for _ in range(10_000)]
    boot.sort()
    return {
        "median": statistics.median(ratios),
        "ci95": [boot[249], boot[9749]],
        "winning_blocks": sum(value > 1.0 for value in ratios),
        "blocks": len(ratios),
    }


class Adapter:
    def __init__(self, name: str, command: str):
        self.name = name
        self.process = subprocess.Popen(
            shlex.split(command) + ["--server"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        self.metadata = self.call({"command": "metadata"})
        if self.metadata.get("engine") != name:
            raise RuntimeError(f"adapter {name} identified itself as {self.metadata!r}")

    def call(self, request: dict) -> dict:
        assert self.process.stdin is not None and self.process.stdout is not None
        self.process.stdin.write(json.dumps(request, separators=(",", ":")) + "\n")
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            stderr = self.process.stderr.read() if self.process.stderr else ""
            raise RuntimeError(f"adapter {self.name} exited early: {stderr}")
        response = json.loads(line)
        if "error" in response:
            raise RuntimeError(f"adapter {self.name}: {response['error']}")
        return response

    def close(self) -> None:
        if self.process.stdin:
            self.process.stdin.close()
        self.process.wait(timeout=10)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adapter", action="append", required=True, type=adapter_arg)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--scorecard", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blocks", type=int, default=21)
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument("--operation", action="append")
    args = parser.parse_args()
    if args.blocks < 3 or args.iterations < 1:
        parser.error("--blocks must be at least 3 and --iterations must be positive")
    commands = dict(args.adapter)
    if len(commands) != len(args.adapter) or "chematic" not in commands:
        parser.error("adapter names must be unique and include chematic")
    rows = [json.loads(line) for line in args.corpus.read_text().splitlines() if line.strip()]
    smiles = [row["smiles"] for row in rows]
    scorecard = json.loads(args.scorecard.read_text())
    operations = args.operation or ["parse", "formula", "tpsa_milli", "hba", "hbd", "rings"]
    adapters = {name: Adapter(name, command) for name, command in commands.items()}
    samples = {
        lane: {operation: {name: [] for name in adapters} for operation in operations}
        for lane in ("pipeline", "prepared")
    }
    try:
        for lane in samples:
            for operation in operations:
                if operation == "parse" and lane == "prepared":
                    continue
                request = {
                    "command": "benchmark",
                    "operation": operation,
                    "lane": lane,
                    "smiles": smiles,
                    "iterations": args.iterations,
                }
                for adapter in adapters.values():
                    adapter.call(request)  # warm-up
                names = list(adapters)
                for block in range(args.blocks):
                    order = names[block % len(names):] + names[:block % len(names)]
                    for name in order:
                        result = adapters[name].call(request)
                        samples[lane][operation][name].append(result)
    finally:
        for adapter in adapters.values():
            adapter.close()

    comparisons = {}
    reference_engine = scorecard["reference_engine"]
    corpus_records = scorecard["corpus_records"]
    for lane, lane_values in samples.items():
        comparisons[lane] = {}
        for operation, by_engine in lane_values.items():
            if not by_engine["chematic"]:
                continue
            operation_score = scorecard.get("operations", {}).get(operation, {})
            comparisons[lane][operation] = {}
            for engine, values in by_engine.items():
                if engine == "chematic":
                    continue
                if engine == reference_engine:
                    counts = {
                        "match": corpus_records,
                        "mismatch": 0,
                        "uncomparable": 0,
                    }
                else:
                    counts = operation_score.get("against_reference", {}).get(engine, {})
                chematic_counts = operation_score.get("against_reference", {}).get("chematic", {})
                rankable = (
                    counts.get("match", 0) == corpus_records
                    and chematic_counts.get("match", 0) == corpus_records
                    and counts.get("mismatch", 0) == 0
                    and counts.get("uncomparable", 0) == 0
                    and chematic_counts.get("mismatch", 0) == 0
                    and chematic_counts.get("uncomparable", 0) == 0
                    and all(row["errors"] == 0 for row in values + by_engine["chematic"])
                )
                entry = {"rankable": rankable}
                if rankable:
                    entry["competitor_over_chematic"] = ratio_ci(
                        [row["elapsed_ns"] for row in values],
                        [row["elapsed_ns"] for row in by_engine["chematic"]],
                        seed=17,
                    )
                comparisons[lane][operation][engine] = entry
    report = {
        "schema_version": 1,
        "blocks": args.blocks,
        "iterations": args.iterations,
        "corpus": str(args.corpus),
        "engines": {name: adapter.metadata for name, adapter in adapters.items()},
        "samples": samples,
        "comparisons": comparisons,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
