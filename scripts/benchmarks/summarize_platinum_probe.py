#!/usr/bin/env python3
"""Summarize paired Platinum development recovery, counting failures in the denominator."""

import argparse
import json
import math
import random
from collections import Counter
from pathlib import Path


def recovered(row, engine, stage, lowest_energy=False):
    prefix = engine + ("_raw" if stage == "raw" else "")
    distances = row.get(prefix + "_rmsd", [])
    stereo = row.get(prefix + "_stereo", [])
    if row.get("status") != "ok" or row.get(engine + "_error") or len(distances) != len(stereo):
        return False
    if lowest_energy:
        energies = row.get(engine + "_energies", [])
        if len(energies) != len(distances):
            return False
        finite = [i for i, energy in enumerate(energies) if math.isfinite(energy)]
        if not finite:
            return False
        index = min(finite, key=lambda i: energies[i])
        distances, stereo = [distances[index]], [stereo[index]]
    return any(math.isfinite(rms) and rms <= 1.0 and not error["atoms"] and not error["bonds"] for rms, error in zip(distances, stereo))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10000)
    args = parser.parse_args()
    rows = json.loads(args.input.read_text())["results"]
    if not rows:
        raise ValueError("no rows")
    if args.bootstrap < 1:
        raise ValueError("bootstrap must be positive")
    rng = random.Random(42)
    summary = {"rows": len(rows), "statuses": dict(Counter(row.get("status", "missing") for row in rows)), "bootstrap_seed": 42, "bootstrap_samples": args.bootstrap, "stages": {}}
    for stage, lowest in (("raw", False), ("minimized", False), ("minimized_lowest_energy", True)):
        pairs = [(int(recovered(r, "chematic", stage, lowest)), int(recovered(r, "rdkit", stage, lowest))) for r in rows]
        diff = sum(a - b for a, b in pairs) / len(pairs)
        # Paired resampling must draw one shared index for both engines.
        samples = []
        for _ in range(args.bootstrap):
            draw = [pairs[rng.randrange(len(pairs))] for _ in pairs]
            samples.append(sum(a - b for a, b in draw) / len(draw))
        samples.sort()
        summary["stages"][stage] = {"chematic_recovered": sum(a for a, _ in pairs), "rdkit_recovered": sum(b for _, b in pairs), "difference": diff, "paired_bootstrap_95_percentile_interval": [samples[int(0.025 * (len(samples)-1))], samples[int(0.975 * (len(samples)-1))]]}
    for engine in ("chematic", "rdkit"):
        summary[engine] = {}
        for stage in ("raw", "minimized"):
            key = engine + ("_raw" if stage == "raw" else "") + "_stereo"
            violations = [v for row in rows for v in row.get(key, [])]
            summary[engine][stage] = {"conformers_checked": len(violations), "stereo_violations": sum(bool(v["atoms"] or v["bonds"]) for v in violations)}
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
