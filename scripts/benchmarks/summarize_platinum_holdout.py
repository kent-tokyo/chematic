#!/usr/bin/env python3
"""Summarize paired recovery and scaffold-cluster uncertainty for a fixed probe."""
import argparse
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from summarize_platinum_probe import recovered


def converged_recovery(row, engine, geometry):
    rmsd = row.get(f"{engine}_rmsd", [])
    statuses = row.get(f"{engine}_statuses", [])
    if len(rmsd) != len(statuses):
        return False
    filtered = dict(row)
    filtered[f"{engine}_rmsd"] = [value if status == 0 else 1e100 for value, status in zip(rmsd, statuses)]
    return recovered(filtered, engine, "minimized", require_geometry=geometry)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("split", type=Path)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=10000)
    args = p.parse_args()
    data = json.loads(args.input.read_text())
    split = json.loads(args.split.read_text())
    rows = data["results"]
    if data["split_sha256"] != hashlib.sha256(args.split.read_bytes()).hexdigest():
        raise ValueError("split hash mismatch")
    expected = {r["index"] for r in split["rows"] if r.get("split") == "holdout"}
    actual = [r["index"] for r in rows]
    if set(actual) != expected or len(actual) != len(expected):
        raise ValueError("holdout is incomplete or contains duplicate records")
    if not rows or args.bootstrap < 1:
        raise ValueError("rows and bootstrap samples must be positive")
    scaffold = {r["index"]: r["scaffold_sha256"] for r in split["rows"] if "scaffold_sha256" in r}
    geometry = data.get("geometry_checked", False)
    result = {"rows": len(rows), "statuses": dict(Counter(r.get("status", "missing") for r in rows)), "bootstrap_samples": args.bootstrap, "bootstrap_seed": 42, "stages": {}}
    rng = random.Random(42)
    for label, fn in (("recovery", lambda r,e: recovered(r,e,"minimized",require_geometry=geometry)), ("converged_recovery", lambda r,e: converged_recovery(r,e,geometry))):
        pairs = [(int(fn(r,"chematic")), int(fn(r,"rdkit"))) for r in rows]
        clusters = defaultdict(list)
        for row, pair in zip(rows, pairs):
            clusters[scaffold[row["index"]]].append(pair[0] - pair[1])
        cluster_totals = [(sum(v),len(v)) for v in clusters.values()]
        n = len(pairs)
        individual = []
        clustered = []
        for _ in range(args.bootstrap):
            sample = [pairs[rng.randrange(n)] for _ in range(n)]
            individual.append(sum(a-b for a,b in sample) / n)
            draw = [cluster_totals[rng.randrange(len(cluster_totals))] for _ in cluster_totals]
            clustered.append(sum(x for x,_ in draw) / sum(y for _,y in draw))
        individual.sort(); clustered.sort()
        lo = int(0.025 * (args.bootstrap - 1)); hi = int(0.975 * (args.bootstrap - 1))
        result["stages"][label] = {"chematic_recovered":sum(a for a,_ in pairs), "rdkit_recovered":sum(b for _,b in pairs), "difference":sum(a-b for a,b in pairs)/n, "paired_molecule_95_interval":[individual[lo],individual[hi]], "scaffold_clusters":len(clusters), "paired_scaffold_cluster_95_interval":[clustered[lo],clustered[hi]]}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
