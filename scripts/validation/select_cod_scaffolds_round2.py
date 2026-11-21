#!/usr/bin/env python3
"""Reserve untouched COD scaffold groups after the first candidate's evaluation."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("split", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    raw = args.split.read_bytes()
    data = json.loads(raw)
    for row in data["rows"]:
        if row.get("split") in ("development", "holdout"):
            row["split"] = "excluded_previously_evaluated"
    for old, new, limit in (("unused_development", "development", 200), ("unused_holdout", "holdout", 2000)):
        scaffolds = {r["scaffold_sha256"] for r in data["rows"] if r.get("split") == old}
        if len(scaffolds) < limit:
            raise ValueError(f"{old} has fewer than {limit} scaffold groups")
        key = lambda scaffold: (hashlib.sha256(("round2:" + new + ":" + scaffold).encode()).digest(), scaffold)
        selected = set(sorted(scaffolds, key=key)[:limit])
        for row in data["rows"]:
            if row.get("split") == old and row["scaffold_sha256"] in selected:
                row["split"] = new
    data["counts"] = dict(Counter(r.get("split", r.get("status", "unknown")) for r in data["rows"]))
    data["rule"] += "; second round excludes first-round development and holdout, then selects unused scaffold groups by SHA256(round2:partition:scaffold hash)"
    data["second_round_selection"] = {"input_split_sha256": hashlib.sha256(raw).hexdigest(), "development_scaffolds": 200, "holdout_scaffolds": 2000}
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
