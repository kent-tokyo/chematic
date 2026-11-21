#!/usr/bin/env python3
"""Report paired scaffold-cluster intervals for the full time-budget holdout."""
import argparse
import json
import random
from collections import defaultdict
from pathlib import Path

from summarize_platinum_time_budget import recovered


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("split", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=10000)
    args = p.parse_args()
    data = json.loads(args.input.read_text())
    split = json.loads(args.split.read_text())
    if data["source_sha256"] != split["source_sha256"]:
        raise ValueError("source hash mismatch")
    expected = {r["index"]: r["scaffold_sha256"] for r in split["rows"] if r.get("split") == "holdout"}
    rows = data["results"]
    indices = [row["index"] for row in rows]
    if set(indices) != set(expected) or len(indices) != len(expected) or args.bootstrap < 1:
        raise ValueError("holdout is incomplete, duplicated, or bootstrap count is invalid")
    geometry = data.get("geometry_checked", False)
    result = {"rows": len(rows), "seconds": data["seconds"], "bootstrap_samples": args.bootstrap, "bootstrap_seed": 42, "stages": {}}
    rng = random.Random(42)
    for label, converged in (("recovery", False), ("converged_recovery", True)):
        clusters = defaultdict(list)
        for row in rows:
            delta = int(recovered(row, "chematic", geometry, converged)) - int(recovered(row, "rdkit", geometry, converged))
            clusters[expected[row["index"]]].append(delta)
        totals = [(sum(values), len(values)) for values in clusters.values()]
        samples = []
        for _ in range(args.bootstrap):
            draw = [totals[rng.randrange(len(totals))] for _ in totals]
            samples.append(sum(value for value, _ in draw) / sum(size for _, size in draw))
        samples.sort()
        lo = int(0.025 * (args.bootstrap - 1))
        hi = int(0.975 * (args.bootstrap - 1))
        result["stages"][label] = {"difference": sum(v for v, _ in totals) / len(rows), "scaffold_clusters": len(totals), "paired_scaffold_cluster_95_interval": [samples[lo], samples[hi]]}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
