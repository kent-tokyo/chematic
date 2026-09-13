#!/usr/bin/env python3
"""Evaluate an adjudicated A5 paired result packet.

The input is JSONL.  Each row represents one predeclared, independently
adjudicated case and must contain ``id``, ``cluster``, ``adjudicated`` and
boolean ``chematic_correct``/``rdkit_correct`` fields.  Engine refusals and
timeouts remain rows in the packet and should be represented by a false
correctness flag plus an optional ``chematic_status``/``rdkit_status``.

Unresolved rows are retained in the accounting but are not silently counted
as correct.  The report includes both the conservative all-row result and the
adjudicated-valid result, so a small or incomplete packet cannot look like an
equivalence proof.
"""

from __future__ import annotations

import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path


MARGIN_PP = 0.1
BOOTSTRAP_REPETITIONS = 10_000
ALLOWED_STATUSES = {"ok", "refused", "timeout", "unsupported", "error"}


def fail(message: str) -> "NoReturn":
    raise SystemExit(f"A5 paired evaluation invalid: {message}")


def load_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    seen: set[str] = set()
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError as exc:
            fail(f"line {line_number} is not JSON: {exc}")
        if not isinstance(row, dict):
            fail(f"line {line_number} must be an object")
        case_id = row.get("id")
        cluster = row.get("cluster")
        if not isinstance(case_id, str) or not case_id or case_id in seen:
            fail(f"line {line_number} has a missing or duplicate id")
        if not isinstance(cluster, str) or not cluster:
            fail(f"line {line_number} has no cluster")
        if not isinstance(row.get("adjudicated"), bool):
            fail(f"line {line_number} adjudicated must be boolean")
        for engine in ("chematic", "rdkit"):
            key = f"{engine}_correct"
            if not isinstance(row.get(key), bool):
                fail(f"line {line_number} {key} must be boolean")
            status = row.get(f"{engine}_status", "ok")
            if status not in ALLOWED_STATUSES:
                fail(f"line {line_number} has invalid {engine}_status")
            if status != "ok" and row[key]:
                fail(f"line {line_number} marks non-ok {engine} status as correct")
        seen.add(case_id)
        rows.append(row)
    if not rows:
        fail("packet is empty")
    return rows


def category(row: dict) -> str:
    if not row["adjudicated"]:
        return "unresolved"
    schematic = row["chematic_correct"]
    rdkit = row["rdkit_correct"]
    if schematic and rdkit:
        return "both_correct"
    if schematic:
        return "chematic_only_correct"
    if rdkit:
        return "rdkit_only_correct"
    return "both_incorrect"


def score(rows: list[dict]) -> dict:
    valid = [row for row in rows if row["adjudicated"]]
    return {
        "rows": len(rows),
        "adjudicated_rows": len(valid),
        "schematic_correct": sum(row["chematic_correct"] for row in valid),
        "rdkit_correct": sum(row["rdkit_correct"] for row in valid),
        "schematic_accuracy": (sum(row["chematic_correct"] for row in valid) / len(valid)) if valid else None,
        "rdkit_accuracy": (sum(row["rdkit_correct"] for row in valid) / len(valid)) if valid else None,
    }


def cluster_bootstrap(rows: list[dict], repetitions: int, seed: int) -> dict:
    valid = [row for row in rows if row["adjudicated"]]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in valid:
        grouped[row["cluster"]].append(row)
    clusters = list(grouped)
    if not clusters:
        return {"clusters": 0, "repetitions": 0, "difference_pp": None, "ci95_pp": None}
    rng = random.Random(seed)
    effects: list[float] = []
    for _ in range(repetitions):
        sampled = [rng.choice(clusters) for _ in clusters]
        sampled_rows = [row for cluster in sampled for row in grouped[cluster]]
        n = len(sampled_rows)
        effects.append(100.0 * sum(row["chematic_correct"] - row["rdkit_correct"] for row in sampled_rows) / n)
    effects.sort()
    low = effects[int(0.025 * (len(effects) - 1))]
    high = effects[int(0.975 * (len(effects) - 1))]
    observed = 100.0 * sum(row["chematic_correct"] - row["rdkit_correct"] for row in valid) / len(valid)
    return {
        "clusters": len(clusters),
        "repetitions": repetitions,
        "seed": seed,
        "difference_pp": observed,
        "ci95_pp": [low, high],
    }


def verdict(bootstrap: dict, valid_rows: int) -> str:
    interval = bootstrap.get("ci95_pp")
    if valid_rows == 0 or bootstrap.get("clusters", 0) < 2 or interval is None:
        return "insufficient_evidence"
    low, high = interval
    if low > 0:
        return "schematic_superior"
    if low >= -MARGIN_PP and high <= MARGIN_PP:
        return "equivalent"
    return "not_equivalent_or_inconclusive"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", type=int, default=20260913)
    parser.add_argument("--repetitions", type=int, default=BOOTSTRAP_REPETITIONS)
    args = parser.parse_args()
    if args.repetitions < 100:
        fail("bootstrap repetitions must be at least 100")
    rows = load_rows(args.input)
    categories = Counter(category(row) for row in rows)
    summary = score(rows)
    bootstrap = cluster_bootstrap(rows, args.repetitions, args.seed)
    report = {
        "schema_version": 1,
        "protocol": "a5-paired-adjudication-v1",
        "input": str(args.input),
        "category_counts": dict(sorted(categories.items())),
        "all_rows": {
            "rows": len(rows),
            "schematic_correct": sum(row["adjudicated"] and row["chematic_correct"] for row in rows),
            "rdkit_correct": sum(row["adjudicated"] and row["rdkit_correct"] for row in rows),
        },
        "adjudicated_valid": summary,
        "paired_bootstrap": bootstrap,
        "equivalence_margin_pp": MARGIN_PP,
        "verdict": verdict(bootstrap, summary["adjudicated_rows"]),
        "unresolved_policy": "retained_in_all_rows; not counted as correct; valid-rate sensitivity required",
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
