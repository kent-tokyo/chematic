#!/usr/bin/env python3
"""Summarize published-wheel 3D timing, coverage, and common-judge quality."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path


ARM_PAIRS = {
    "raw_embedding": (
        "chematic_pipeline_v2_no_ff",
        "rdkit_etkdgv3_raw",
    ),
    "embedding_plus_uff": (
        "chematic_pipeline_v2_uff_only",
        "rdkit_etkdgv3_uff",
    ),
    "embedding_plus_mmff94": (
        "chematic_pipeline_v2_mmff94_strict",
        "rdkit_etkdgv3_mmff94",
    ),
    "best_of_10_plus_uff": (
        "chematic_pipeline_v2_uff_best_of_10",
        "rdkit_etkdgv3_best_of_n",
    ),
}

T_975 = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
    11: 2.201,
    12: 2.179,
    13: 2.160,
    14: 2.145,
    15: 2.131,
    16: 2.120,
    17: 2.110,
    18: 2.101,
    19: 2.093,
    20: 2.086,
    25: 2.060,
    30: 2.042,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_jsonl(path: Path) -> list[dict[str, object]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def timing_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    elapsed = [float(row["elapsed_ms"]) for row in rows]
    successful = [
        float(row["elapsed_ms"]) for row in rows if row.get("status") == "success"
    ]
    return {
        "rows": len(rows),
        "status_counts": dict(
            sorted(Counter(str(row.get("status")) for row in rows).items())
        ),
        "all_rows": {
            "p50_ms": percentile(elapsed, 0.5),
            "p95_ms": percentile(elapsed, 0.95),
            "total_ms": sum(elapsed),
        },
        "successful_rows": {
            "count": len(successful),
            "p50_ms": percentile(successful, 0.5),
            "p95_ms": percentile(successful, 0.95),
        },
    }


def critical_value(degrees_of_freedom: int) -> float:
    for maximum_df in sorted(T_975):
        if degrees_of_freedom <= maximum_df:
            return T_975[maximum_df]
    return 1.960


def paired_speedup(
    chematic_rows: list[dict[str, object]], rdkit_rows: list[dict[str, object]]
) -> dict[str, object]:
    def successful(rows: list[dict[str, object]]) -> dict[tuple[str, str], float]:
        return {
            (str(row["tier"]), str(row["name"])): float(row["elapsed_ms"])
            for row in rows
            if row.get("status") == "success" and float(row["elapsed_ms"]) > 0
        }

    left = successful(chematic_rows)
    right = successful(rdkit_rows)
    keys = sorted(left.keys() & right.keys())
    log_ratios = [math.log(right[key] / left[key]) for key in keys]
    if len(log_ratios) < 2:
        return {"common_successes": len(keys), "status": "insufficient_pairs"}
    mean_log = statistics.fmean(log_ratios)
    standard_error = statistics.stdev(log_ratios) / math.sqrt(len(log_ratios))
    lower = math.exp(mean_log - critical_value(len(log_ratios) - 1) * standard_error)
    ratios = [math.exp(value) for value in log_ratios]
    return {
        "common_successes": len(keys),
        "geometric_mean_rdkit_over_chematic": math.exp(mean_log),
        "lower_95_rdkit_over_chematic": lower,
        "median_rdkit_over_chematic": percentile(ratios, 0.5),
        "minimum_rdkit_over_chematic": min(ratios),
        "chematic_faster_gate": lower > 1.0,
    }


def quality_summary(
    scored_rows: list[dict[str, object]], engine: str, arm: str, denominator: int
) -> dict[str, object]:
    rows = [
        row
        for row in scored_rows
        if row.get("engine") == engine and row.get("arm") == arm
    ]
    sound = sum(row.get("independently_sound") is True for row in rows)
    stereo_clean = sum(
        isinstance(row.get("stereo"), dict) and row["stereo"].get("violated") == 0  # type: ignore[index]
        for row in rows
    )
    usable = sum(
        row.get("independently_sound") is True
        and isinstance(row.get("stereo"), dict)
        and row["stereo"].get("violated") == 0  # type: ignore[index]
        for row in rows
    )
    return {
        "scored_successes": len(rows),
        "independently_sound": sound,
        "stereo_clean": stereo_clean,
        "usable": usable,
        "usable_over_all_inputs": usable / denominator if denominator else None,
    }


def rows_for_arm(rows: list[dict[str, object]], arm: str) -> list[dict[str, object]]:
    return [row for row in rows if row.get("arm") == arm]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chematic", type=Path, required=True)
    parser.add_argument("--rdkit", type=Path, required=True)
    parser.add_argument("--scored", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--pair",
        nargs=3,
        action="append",
        metavar=("LABEL", "CHEMATIC_ARM", "RDKIT_ARM"),
        help="summarize an explicit arm pair; repeat for multiple pairs",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    chematic_rows = load_jsonl(args.chematic)
    rdkit_rows = load_jsonl(args.rdkit)
    scored_rows = load_jsonl(args.scored)
    result: dict[str, object] = {
        "schema_version": 1,
        "benchmark": "public-package-3d-vs-rdkit-summary",
        "inputs": {
            "chematic": {"path": str(args.chematic), "sha256": sha256(args.chematic)},
            "rdkit": {"path": str(args.rdkit), "sha256": sha256(args.rdkit)},
            "scored": {"path": str(args.scored), "sha256": sha256(args.scored)},
        },
        "arm_pairs": {},
        "boundary": "Speed is paired only over common successful molecules. Quality uses the same external geometry/stereo scorer for both engines and retains every failure in the all-input denominator. A lane passes speed only when the paired 95% lower bound exceeds 1.0; this does not by itself establish quality non-inferiority.",
    }
    arm_results: dict[str, object] = {}
    arm_pairs = (
        {
            label: (chematic_arm, rdkit_arm)
            for label, chematic_arm, rdkit_arm in args.pair
        }
        if args.pair
        else ARM_PAIRS
    )
    for label, (chematic_arm, rdkit_arm) in arm_pairs.items():
        left = rows_for_arm(chematic_rows, chematic_arm)
        right = rows_for_arm(rdkit_rows, rdkit_arm)
        if len(left) != len(right):
            raise SystemExit(
                f"row-count mismatch for {label}: chematic={len(left)} rdkit={len(right)}"
            )
        denominator = len(left)
        arm_results[label] = {
            "arms": {"chematic": chematic_arm, "rdkit": rdkit_arm},
            "timing": {
                "chematic": timing_summary(left),
                "rdkit": timing_summary(right),
                "paired": paired_speedup(left, right),
            },
            "quality": {
                "chematic": quality_summary(
                    scored_rows, "chematic", chematic_arm, denominator
                ),
                "rdkit": quality_summary(scored_rows, "rdkit", rdkit_arm, denominator),
            },
        }
    result["arm_pairs"] = arm_results
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(arm_results, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
