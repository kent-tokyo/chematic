#!/usr/bin/env python3
"""Summarize recorded attempt counts and elapsed times of a timed holdout."""
import argparse
import json
import math
import statistics
from pathlib import Path


def distribution(values):
    if not values:
        return None
    values = sorted(values)
    # Nearest-rank percentiles retain an observed value, including for small sets.
    return {
        "min": values[0],
        "median": statistics.median(values),
        "p95": values[math.ceil(0.95 * len(values)) - 1],
        "max": values[-1],
    }


def summarize(data):
    seconds = data["seconds"]
    cap = data["max_candidates"]
    rows = data["results"]
    if not rows or not math.isfinite(seconds) or seconds <= 0 or cap < 1:
        raise ValueError("invalid timing budget, candidate cap, or empty results")
    output = {
        "rows": len(rows),
        "platform": data["platform"],
        "seconds": seconds,
        "max_candidates": cap,
        "timing": data["timing"],
        "errors": sum(bool(row.get("error") or row.get("engine_errors")) for row in rows),
        "engine_errors": {engine: sum(engine in row.get("engine_errors", {}) for row in rows)
                          for engine in ("chematic", "rdkit")},
        "engines": {},
    }
    for engine in ("chematic", "rdkit"):
        records = [row[engine] for row in rows if engine in row]
        for record in records:
            attempted, accepted, elapsed = (record[k] for k in ("attempted", "accepted", "elapsed"))
            if (not isinstance(attempted, int) or not isinstance(accepted, int)
                    or not 0 <= accepted <= attempted <= cap
                    or not math.isfinite(elapsed) or elapsed < 0):
                raise ValueError(f"invalid {engine} timing record")
        output["engines"][engine] = {
            "records": len(records),
            "missing_records": len(rows) - len(records),
            "attempted": distribution([r["attempted"] for r in records]),
            "accepted": distribution([r["accepted"] for r in records]),
            "elapsed_seconds": distribution([r["elapsed"] for r in records]),
            "elapsed_over_budget": sum(r["elapsed"] > seconds for r in records),
            "reached_attempt_cap": sum(r["attempted"] == cap for r in records),
        }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(json.loads(args.input.read_text()))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
