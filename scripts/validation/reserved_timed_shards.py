#!/usr/bin/env python3
"""Prepare and merge timed holdout shards without filtering failed inputs."""

import argparse
import hashlib
import json
import re
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(sdf, split_path, out, shard_count):
    from rdkit import Chem, RDLogger, rdBase
    RDLogger.DisableLog("rdApp.*")
    split = json.loads(split_path.read_text())
    if sha(sdf) != split["source_sha256"]:
        raise ValueError("source hash differs from reserved split")
    if shard_count < 1:
        raise ValueError("shard count must be positive")
    supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
    selected = [row for row in split["rows"] if row.get("split") == "holdout"]
    if not selected or len({row["index"] for row in selected}) != len(selected):
        raise ValueError("holdout indices are empty or duplicated")
    out.mkdir(parents=True, exist_ok=True)
    shards = [[] for _ in range(shard_count)]
    for row in selected:
        item = {"index": row["index"], "name": row["name"]}
        try:
            reference = supplier[row["index"]]
            if reference is None or reference.GetProp("_Name") != row["name"]:
                raise ValueError("reference parse or identity mismatch")
            item["smiles"] = Chem.MolToSmiles(Chem.RemoveHs(reference), isomericSmiles=True)
        except Exception as error:
            item["input_error"] = repr(error)
        shards[int(row["scaffold_sha256"], 16) % shard_count].append(item)
    manifest = {"source_sha256": split["source_sha256"], "split_sha256": sha(split_path),
                "rdkit_version": rdBase.rdkitVersion, "shard_count": shard_count, "shards": []}
    for i, rows in enumerate(shards):
        path = out / f"baseline-{i:03d}.json"
        path.write_text(json.dumps({"source_sha256": split["source_sha256"], "results": rows}, indent=2) + "\n")
        manifest["shards"].append({"index": i, "count": len(rows), "sha256": sha(path),
                                   "records": [{"index": r["index"], "name": r["name"]} for r in rows]})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def merge(manifest_path, result_paths, out, wheel_sha256, source_commit, *, allow_mixed_platforms=False):
    manifest = json.loads(manifest_path.read_text())
    if len(result_paths) != manifest["shard_count"]:
        raise ValueError("wrong number of shards")
    records = {item["index"]: item for shard in manifest["shards"] for item in shard["records"]}
    if sum(s["count"] for s in manifest["shards"]) != len(records):
        raise ValueError("duplicate records in manifest")
    results, template, signatures, platforms = [], None, set(), {}
    for path in result_paths:
        data = json.loads(path.read_text())
        if (data.get("source_sha256") != manifest["source_sha256"]
                or data.get("wheel_sha256") != wheel_sha256
                or data.get("source_commit") != source_commit
                or data.get("rdkit") != manifest["rdkit_version"]
                or data.get("profile") != "mixed-fallback"
                or data.get("stereo_policy") != "strict"
                or data.get("hard_timeout") is not True
                or data.get("geometry_checked") is not True
                or data.get("seconds") != 0.2 or data.get("max_candidates") != 20
                or data.get("repair_iters") != 400):
            raise ValueError(f"candidate or evaluation settings differ in {path}")
        signature = (data.get("script_sha256"), data.get("stereo_verifier_sha256"),
                     data.get("platform"), data.get("chematic"))
        if None in signature:
            raise ValueError("missing scorer, platform or package provenance")
        benchmark_dir = Path(__file__).resolve().parents[1] / "benchmarks"
        if (signature[0] != sha(benchmark_dir / "platinum_time_budget.py")
                or signature[1] != sha(benchmark_dir / "platinum_stereo_strict.py")):
            raise ValueError("scorer or stereo verifier differs from checkout")
        signatures.add(signature if not allow_mixed_platforms else (signature[0], signature[1], signature[3]))
        platform = data["platform"]
        platforms[platform] = platforms.get(platform, 0) + len(data["results"])
        template = data if template is None else template
        results.extend([{**row, "evaluation_platform": platform} for row in data["results"]])
    if len(signatures) != 1 or len(results) != len(records):
        raise ValueError("scorer mismatch or incomplete results")
    if len(platforms) > 1:
        families = set()
        for platform in platforms:
            match = re.search(r"-(x86_64|aarch64|arm64)(?:-|$)", platform)
            if not match or platform.split("-", 1)[0] not in ("Linux", "macOS"):
                raise ValueError("unrecognized mixed platform")
            families.add((platform.split("-", 1)[0], match.group(1)))
        if len(families) != 1:
            raise ValueError("mixed operating systems or architectures")
    observed = {row["index"]: row for row in results}
    if len(observed) != len(results) or set(observed) != set(records):
        raise ValueError("missing or duplicated result identities")
    if any(observed[i].get("name") != records[i]["name"] for i in records):
        raise ValueError("result name differs from manifest")
    template["results"] = sorted(results, key=lambda row: row["index"])
    template["platforms"] = [{"platform": p, "records": n} for p, n in sorted(platforms.items())]
    if len(platforms) > 1:
        template["platform"] = "mixed; see platforms and each result evaluation_platform"
    template["split_sha256"] = manifest["split_sha256"]
    template["shard_manifest_sha256"] = sha(manifest_path)
    out.write_text(json.dumps(template, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("sdf", type=Path)
    prep.add_argument("split", type=Path)
    prep.add_argument("output", type=Path)
    prep.add_argument("--shards", type=int, default=8)
    combine = sub.add_parser("merge")
    combine.add_argument("manifest", type=Path)
    combine.add_argument("output", type=Path)
    combine.add_argument("results", nargs="+", type=Path)
    combine.add_argument("--wheel-sha256", required=True)
    combine.add_argument("--source-commit", required=True)
    combine.add_argument("--allow-mixed-platforms", action="store_true",
                         help="allow same-OS, same-architecture runners with differing platform strings; record each result platform")
    args = parser.parse_args()
    if args.command == "prepare":
        manifest = prepare(args.sdf, args.split, args.output, args.shards)
        print(json.dumps([s["count"] for s in manifest["shards"]]))
    else:
        merge(args.manifest, args.results, args.output, args.wheel_sha256, args.source_commit,
              allow_mixed_platforms=args.allow_mixed_platforms)


if __name__ == "__main__":
    main()
