#!/usr/bin/env python3
"""Partition an already-audited common-domain COD holdout by scaffold hash."""
import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(data):
    return (json.dumps(data, indent=2) + "\n").encode("utf-8")


def shard(split_bytes, audit_bytes, shard_count):
    if shard_count < 1:
        raise ValueError("shard count must be positive")
    split, audit = json.loads(split_bytes), json.loads(audit_bytes)
    parent_hash = digest(split_bytes)
    if (
        audit["split_sha256"] != parent_hash
        or audit["source_sha256"] != split["source_sha256"]
        or audit["invalid_exclusion_count"] != len(split["invalid_exclusions"])
        or audit["overlaps"]
    ):
        raise ValueError("parent invalid-exclusion audit is not clean or does not match")
    selected = [r for r in split["rows"] if r.get("split") == "holdout"]
    if not selected or len({r["index"] for r in selected}) != len(selected):
        raise ValueError("holdout must be nonempty and have unique record indices")
    assignments = {}
    for row in selected:
        key = row["scaffold_sha256"]
        if len(key) != 64 or any(c not in "0123456789abcdef" for c in key):
            raise ValueError("invalid scaffold SHA-256")
        assignments[row["index"]] = int(key, 16) % shard_count
    outputs = []
    for part in range(shard_count):
        item = deepcopy(split)
        for row in item["rows"]:
            if row.get("split") == "holdout" and assignments[row["index"]] != part:
                row["split"] = "heldout_other_shard"
        holdout = [r for r in item["rows"] if r.get("split") == "holdout"]
        if not holdout:
            raise ValueError(f"shard {part} is empty")
        item["counts"] = dict(Counter(r.get("split", r.get("status", "unknown")) for r in item["rows"]))
        if "eligibility" in item:
            item["eligibility"]["parent_eligible_holdout_records"] = item["eligibility"]["eligible_holdout_records"]
            item["eligibility"]["eligible_holdout_records"] = len(holdout)
        item["sharding"] = {
            "rule": "integer scaffold SHA-256 modulo shard count",
            "parent_split_sha256": parent_hash,
            "shard_index": part,
            "shard_count": shard_count,
            "holdout_records": len(holdout),
        }
        item_bytes = json_bytes(item)
        part_audit = deepcopy(audit)
        part_audit["split_sha256"] = digest(item_bytes)
        part_audit["derived_from_audit_sha256"] = digest(audit_bytes)
        part_audit["eligible_scaffolds"] = len({r["scaffold_sha256"] for r in item["rows"] if r.get("split") in ("development", "holdout")})
        part_audit["rule"] += " Filtered to a deterministic scaffold shard."
        outputs.append((item_bytes, json_bytes(part_audit)))
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("split", type=Path)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--shards", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outputs = shard(args.split.read_bytes(), args.audit.read_bytes(), args.shards)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for i, (split_bytes, audit_bytes) in enumerate(outputs):
        (args.output_dir / f"split-{i}.json").write_bytes(split_bytes)
        (args.output_dir / f"audit-{i}.json").write_bytes(audit_bytes)
    print(json.dumps({"shard_counts": [json.loads(x)["sharding"]["holdout_records"] for x, _ in outputs]}))


if __name__ == "__main__":
    main()
