#!/usr/bin/env python3
"""Summarize paired equal-time recovery and candidate counts."""

import argparse
import json
import random
import statistics
from pathlib import Path


def recovered(row, engine, require_geometry=False, require_converged=False):
    record = row.get(engine, {})
    rmsd = record.get("rmsd", [])
    stereo = record.get("stereo", [])
    if row.get("error") or len(rmsd) != len(stereo):
        return False
    geometry = record.get("geometry", [])
    if require_geometry:
        reference = row.get("geometry_reference", {})
        if len(geometry) != len(rmsd) or not reference:
            return False
        if (
            reference["bonds_supported"] != reference["bonds_total"]
            or reference["angles_supported"] != reference["angles_total"]
        ):
            return False
    else:
        geometry = [None] * len(rmsd)
    statuses = record.get("statuses", [])
    if require_converged and len(statuses) != len(rmsd):
        return False
    if not require_converged:
        statuses = [0] * len(rmsd)
    return any(
        distance <= 1.0
        and optimizer_status == 0
        and not status["atoms"]
        and not status["bonds"]
        and (
            screen is None
            or (
                screen["finite"]
                and screen["bond_violations"] == 0
                and screen["angle_violations"] == 0
                and screen["heavy_atom_clashes"] == 0
            )
        )
        for distance, status, screen, optimizer_status in zip(
            rmsd, stereo, geometry, statuses
        )
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
    geometry = report.get("geometry_checked", False)
    pairs = [
        (int(recovered(r, "chematic", geometry)), int(recovered(r, "rdkit", geometry)))
        for r in rows
    ]
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
    if all(
        row.get("error") or "statuses" in row.get(engine, {})
        for row in rows
        for engine in ("chematic", "rdkit")
    ):
        converged_pairs = [
            (
                int(recovered(r, "chematic", geometry, True)),
                int(recovered(r, "rdkit", geometry, True)),
            )
            for r in rows
        ]
        draws = []
        for _ in range(10000):
            draw = [
                converged_pairs[rng.randrange(len(converged_pairs))]
                for _ in converged_pairs
            ]
            draws.append(sum(a - b for a, b in draw) / len(draw))
        draws.sort()
        summary["converged_recovery"] = {
            "chematic": sum(a for a, _ in converged_pairs),
            "rdkit": sum(b for _, b in converged_pairs),
            "difference": sum(a - b for a, b in converged_pairs) / len(converged_pairs),
            "paired_bootstrap_95_percentile_interval": [draws[249], draws[9749]],
        }
    for engine in ("chematic", "rdkit"):
        records = [row[engine] for row in rows if engine in row]
        summary[engine] = {
            "recovered": sum(recovered(row, engine, geometry) for row in rows),
            "median_accepted": (
                statistics.median(record["accepted"] for record in records)
                if records
                else None
            ),
            "reached_candidate_cap": sum(
                record["accepted"] == report["max_candidates"] for record in records
            ),
            "stereo_violations": sum(
                bool(v["atoms"] or v["bonds"])
                for record in records
                for v in record.get("stereo", [])
            ),
            "converged": (
                sum(status == 0 for record in records for status in record["statuses"])
                if all("statuses" in record for record in records)
                else None
            ),
        }
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
