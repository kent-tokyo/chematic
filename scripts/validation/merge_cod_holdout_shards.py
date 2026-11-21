#!/usr/bin/env python3
"""Merge complete, independently audited COD result shards without changing rows."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(data):
    return hashlib.sha256(data).hexdigest()


def merge(parent_bytes, pairs, kind):
    if kind not in ("fixed", "time"):
        raise ValueError("kind must be fixed or time")
    parent = json.loads(parent_bytes)
    parent_sha = digest(parent_bytes)
    expected = {r["index"]: r for r in parent["rows"] if r.get("split") == "holdout"}
    if not expected or len(expected) != sum(r.get("split") == "holdout" for r in parent["rows"]):
        raise ValueError("parent holdout is empty or has duplicate indices")
    combined, seen, parts, provenance = [], set(), set(), []
    common, template = None, None
    count = len(pairs)
    for split_bytes, data_bytes in pairs:
        split, data = json.loads(split_bytes), json.loads(data_bytes)
        shard = split["sharding"]
        part = shard["shard_index"]
        if shard["parent_split_sha256"] != parent_sha or shard["shard_count"] != count or not 0 <= part < count or part in parts:
            raise ValueError("invalid, duplicate, or unrelated shard")
        parts.add(part)
        expected_part = {i for i, row in expected.items() if int(row["scaffold_sha256"], 16) % count == part}
        shard_rows = [r for r in split["rows"] if r.get("split") == "holdout"]
        observed_part = {r["index"] for r in shard_rows}
        if (observed_part != expected_part or len(shard_rows) != len(expected_part)
                or shard["holdout_records"] != len(expected_part)
                or split["source_sha256"] != parent["source_sha256"]
                or any(r["name"] != expected[r["index"]]["name"]
                       or r["scaffold_sha256"] != expected[r["index"]]["scaffold_sha256"]
                       for r in shard_rows)):
            raise ValueError("shard membership differs from deterministic parent assignment")
        rows = data["results"]
        indices = [r["index"] for r in rows]
        if set(indices) != expected_part or len(indices) != len(expected_part) or seen.intersection(indices):
            raise ValueError("shard results are incomplete, duplicated, or misassigned")
        if any(r["name"] != expected[r["index"]]["name"] for r in rows):
            raise ValueError("result record identity mismatch")
        if data["source_sha256"] != parent["source_sha256"]:
            raise ValueError("reference source mismatch")
        if kind == "fixed" and data["split_sha256"] != digest(split_bytes):
            raise ValueError("fixed result split hash mismatch")
        metadata = {k: v for k, v in data.items() if k not in ("platform", "split_sha256", "results")}
        if common is None:
            common, template = metadata, dict(data)
        elif metadata != common:
            raise ValueError("shard run settings or software identity differ")
        provenance.append({"shard_index": part, "platform": data["platform"], "split_sha256": digest(split_bytes), "result_sha256": digest(data_bytes), "records": len(rows)})
        combined.extend(rows)
        seen.update(indices)
    if seen != set(expected) or parts != set(range(count)):
        raise ValueError("merged holdout is incomplete")
    template["results"] = sorted(combined, key=lambda r: r["index"])
    template["platform"] = "multiple runners; see shard_provenance"
    if kind == "fixed":
        template["split_sha256"] = parent_sha
    template["shard_provenance"] = sorted(provenance, key=lambda x: x["shard_index"])
    return template


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent_split", type=Path)
    parser.add_argument("--kind", choices=("fixed", "time"), required=True)
    parser.add_argument("--shards-dir", type=Path, required=True)
    parser.add_argument("--shards", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.shards < 1:
        parser.error("shards must be positive")
    pairs = [((args.shards_dir / f"split-{i}.json").read_bytes(), (args.shards_dir / f"{args.kind}-{i}.json").read_bytes()) for i in range(args.shards)]
    result = merge(args.parent_split.read_bytes(), pairs, args.kind)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"kind": args.kind, "records": len(result["results"]), "shards": len(pairs)}))


if __name__ == "__main__":
    main()
