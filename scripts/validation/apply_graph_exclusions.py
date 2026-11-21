#!/usr/bin/env python3
"""Exclude entire eligible scaffold groups flagged by a conservative graph audit."""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("split", type=Path)
    p.add_argument("audit", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    split_bytes = args.split.read_bytes()
    data = json.loads(split_bytes)
    audit_bytes = args.audit.read_bytes()
    audit = json.loads(audit_bytes)
    if audit["split_sha256"] != hashlib.sha256(split_bytes).hexdigest() or audit["source_sha256"] != data["source_sha256"]:
        raise ValueError("audit does not match split")
    indexed = {r["index"]: r for r in data["rows"]}
    scaffolds = {indexed[r["index"]]["scaffold_sha256"] for r in audit["overlaps"]}
    excluded_indices = []
    for row in data["rows"]:
        if row.get("split") in ("development", "holdout") and row["scaffold_sha256"] in scaffolds:
            row["split"] = "excluded_invalid_exclusion_graph"
            excluded_indices.append(row["index"])
    data["counts"] = dict(Counter(r.get("split", r.get("status", "unknown")) for r in data["rows"]))
    data["conservative_graph_exclusion"] = {"audit_sha256": hashlib.sha256(audit_bytes).hexdigest(), "excluded_indices": excluded_indices, "excluded_scaffolds": len(scaffolds)}
    args.output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(data["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
