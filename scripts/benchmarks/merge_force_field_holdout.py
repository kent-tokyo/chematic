#!/usr/bin/env python3
"""Merge complete force-field shards against the sealed holdout manifest."""
import argparse
import hashlib
import json
from pathlib import Path


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def merge(manifest_path, shard_dirs):
    manifest_bytes = manifest_path.read_bytes()
    manifest = json.loads(manifest_bytes)
    expected = {s["index"]: s for s in manifest["shards"]}
    if len(shard_dirs) != manifest["shard_count"] or len(expected) != len(shard_dirs):
        raise ValueError("missing or duplicate manifest shards")
    results, platforms, provenance = [], set(), []
    versions = None
    seen = set()
    for directory in shard_dirs:
        metadata = json.loads((directory / "metadata.json").read_bytes())
        index = metadata["index"]
        if index not in expected or index in seen:
            raise ValueError("unknown or duplicate result shard")
        seen.add(index)
        shard = expected[index]
        if metadata["wheel_sha256"] != manifest["wheel_sha256"]:
            raise ValueError("wheel mismatch")
        if metadata["probe_sha256"] != "4bbc4dcd44c1632263f97e341ed8c955efe1f9d290cb4589a1d81b0d5cdd3e18":
            raise ValueError("force-field probe mismatch")
        if metadata["manifest_sha256"] != sha256(manifest_bytes):
            raise ValueError("manifest mismatch")
        input_bytes = (directory / "input.json").read_bytes()
        result_bytes = (directory / "result.json").read_bytes()
        input_data, result = json.loads(input_bytes), json.loads(result_bytes)
        if (input_data["source_sha256"] != manifest["source_sha256"]
                or input_data["split_sha256"] != shard["split_sha256"]
                or input_data["partition"] != "holdout"
                or result["development_sha256"] != sha256(input_bytes)):
            raise ValueError("shard input provenance mismatch")
        ids = [(r["index"], r["name"]) for r in input_data["results"]]
        want = [(r["index"], r["name"]) for r in shard["records"]]
        actual = [(r["index"], r["name"]) for r in result["results"]]
        if len(ids) != shard["count"] or len(set(ids)) != len(ids) or set(ids) != set(want) or ids != actual:
            raise ValueError("shard records are missing, duplicated or reordered")
        current_versions = (result["rdkit_version"], result["chematic_version"], result["seed"])
        if current_versions[0] != "2026.09.2" or current_versions[2] != 42:
            raise ValueError("unexpected engine or seed")
        if versions is not None and versions != current_versions:
            raise ValueError("inconsistent engine versions")
        versions = current_versions
        platforms.add(result["platform"])
        provenance.append({"index": index, "input_sha256": sha256(input_bytes),
                           "result_sha256": sha256(result_bytes)})
        results.extend(result["results"])
    if seen != set(expected):
        raise ValueError("missing result shard")
    return {"source_sha256": manifest["source_sha256"],
            "split_sha256": manifest["original_split_sha256"],
            "wheel_sha256": manifest["wheel_sha256"],
            "manifest_sha256": sha256(manifest_bytes),
            "rdkit_version": versions[0], "chematic_version": versions[1],
            "seed": versions[2], "platform": sorted(platforms),
            "shards": sorted(provenance, key=lambda r: r["index"]),
            "results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("shard_dirs", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = merge(args.manifest, args.shard_dirs)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"rows": len(output["results"]), "shards": len(output["shards"])}))


if __name__ == "__main__":
    main()
