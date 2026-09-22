#!/usr/bin/env python3
"""Validate a paired Parse + compatible-Morgan browser measurement.

The isolated browser runner already verifies the paired packed-fingerprint
digest and set-bit count. This checker supplies the separate statistical
acceptance condition: the one-sided lower bound of the paired log-speedup
(`RDKit time / chematic time`) must exceed one.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path


# Two-sided 95% Student-t critical values for the small, fixed CI sample
# sizes used by the workflow. Using a table avoids a SciPy dependency in CI.
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


def critical_value(degrees_of_freedom: int) -> float:
    """Return a conservative 95% t critical value for supported sample sizes."""
    for maximum_df in sorted(T_975):
        if degrees_of_freedom <= maximum_df:
            return T_975[maximum_df]
    return 1.960


def operation_mean(run: dict[str, object], operation: str) -> float:
    try:
        value = run["operations"][operation]["mean_ms"]  # type: ignore[index]
    except (KeyError, TypeError) as error:
        raise ValueError(
            f"missing raw_runs.*.operations.{operation}.mean_ms"
        ) from error
    if not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"invalid {operation} mean: {value!r}")
    return float(value)


def paired_summary(
    document: dict[str, object], operation: str = "parse_fp"
) -> dict[str, float | int | str]:
    if document.get("schema_version") != 2:
        raise ValueError("expected benchmark schema_version 2")
    raw = document.get("raw_runs")
    if not isinstance(raw, dict):
        raise ValueError("missing raw_runs")
    chematic = raw.get("chematic")
    rdkit = raw.get("rdkit")
    if not isinstance(chematic, list) or not isinstance(rdkit, list):
        raise ValueError("raw_runs must contain chematic and rdkit lists")
    if len(chematic) != len(rdkit) or len(chematic) < 2:
        raise ValueError("need at least two paired chematic/RDKit repetitions")

    log_speedups = [
        math.log(operation_mean(right, operation) / operation_mean(left, operation))
        for left, right in zip(chematic, rdkit, strict=True)
    ]
    mean_log = statistics.fmean(log_speedups)
    standard_error = statistics.stdev(log_speedups) / math.sqrt(len(log_speedups))
    lower = math.exp(mean_log - critical_value(len(log_speedups) - 1) * standard_error)
    return {
        "operation": operation,
        "pairs": len(log_speedups),
        "geometric_mean_speedup": math.exp(mean_log),
        "lower_95_speedup": lower,
        "minimum_observed_speedup": math.exp(min(log_speedups)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument(
        "--operation", choices=("parse_fp", "prepared_fp"), default="parse_fp"
    )
    parser.add_argument("--min-lower-speedup", type=float, default=1.0)
    args = parser.parse_args()
    if not math.isfinite(args.min_lower_speedup) or args.min_lower_speedup <= 0:
        parser.error("--min-lower-speedup must be finite and positive")

    try:
        document = json.loads(args.input.read_text(encoding="utf-8"))
        summary = paired_summary(document, args.operation)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"Parse + Morgan speed gate invalid: {error}", file=sys.stderr)
        return 2

    print(json.dumps(summary, indent=2))
    if summary["lower_95_speedup"] <= args.min_lower_speedup:
        print(
            "Parse + Morgan speed gate failed: paired 95% lower bound is not above "
            f"{args.min_lower_speedup}",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
