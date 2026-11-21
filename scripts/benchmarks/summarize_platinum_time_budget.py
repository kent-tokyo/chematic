#!/usr/bin/env python3
"""Summarize paired equal-time recovery and candidate counts."""

import argparse
import json
import random
import statistics
from pathlib import Path


def recovered(row, engine):
    record = row.get(engine, {})
    rmsd = record.get("rmsd", [])
    stereo = record.get("stereo", [])
    if row.get("error") or len(rmsd) != len(stereo):
        return False
    return any(
        distance <= 1.0 and not status["atoms"] and not status["bonds"]
        for distance, status in zip(rmsd, stereo)
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = json.loads(args.input.read_text())
    rows = report["results"]
    if not rows:
        raise ValueError("no rows")
    pairs = [(int(recovered(r, "chematic")), int(recovered(r, "rdkit"))) for r in rows]
    rng = random.Random(42)
    samples = []
    for _ in range(10000):
        draw = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        samples.append(sum(a - b for a, b in draw) / len(draw))
    samples.sort()
    summary = {
        "rows": len(rows),
        "seconds": report["seconds"],
        "errors": sum(bool(row.get("error")) for row in rows),
        "recovery_difference": sum(a - b for a, b in pairs) / len(pairs),
        "paired_bootstrap_95_percentile_interval": [samples[249], samples[9749]],
    }
    for engine in ("chematic", "rdkit"):
        records = [row[engine] for row in rows if engine in row]
        summary[engine] = {
            "recovered": sum(recovered(row, engine) for row in rows),
            "median_accepted": statistics.median(record["accepted"] for record in records) if records else None,
            "reached_candidate_cap": sum(record["accepted"] == report["max_candidates"] for record in records),
            "stereo_violations": sum(bool(v["atoms"] or v["bonds"]) for record in records for v in record["stereo"]),
        }
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
