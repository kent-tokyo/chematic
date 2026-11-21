#!/usr/bin/env python3
"""Select complete scaffold groups from a precomputed development/holdout split."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--development-scaffolds", type=int, default=200)
    parser.add_argument("--holdout-scaffolds", type=int, default=2000)
    args = parser.parse_args()
    data = json.loads(args.split.read_text(encoding="utf-8"))
    for partition, limit in (("development", args.development_scaffolds), ("holdout", args.holdout_scaffolds)):
        if limit < 1:
            parser.error("scaffold counts must be positive")
        scaffolds = {r["scaffold_sha256"] for r in data["rows"] if r.get("split") == partition}
        if len(scaffolds) < limit:
            raise ValueError(f"{partition}: only {len(scaffolds)} scaffolds available")
        key = lambda scaffold: (hashlib.sha256((partition + ":" + scaffold).encode()).digest(), scaffold)
        selected = set(sorted(scaffolds, key=key)[:limit])
        for row in data["rows"]:
            if row.get("split") == partition and row["scaffold_sha256"] not in selected:
                row["split"] = "unused_" + partition
    data["counts"] = dict(Counter(r.get("split", r.get("status", "unknown")) for r in data["rows"]))
    data["rule"] += "; within each nonexcluded partition retain complete scaffold groups with lowest SHA256(partition + ':' + scaffold hash)"
    data["selection"] = {"development_scaffolds": args.development_scaffolds, "holdout_scaffolds": args.holdout_scaffolds, "input_split_sha256": hashlib.sha256(args.split.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
