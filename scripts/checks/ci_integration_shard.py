#!/usr/bin/env python3
"""Run every Rust integration test across a small, deterministic CI matrix.

Cargo metadata discovers integration-test targets, so a new test file does not
need a parallel workflow edit. Packages are assigned as units to reuse builds.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from typing import Any


def assign_packages(metadata: dict[str, Any], shard_count: int) -> list[list[str]]:
    if shard_count < 1:
        raise ValueError("shard count must be positive")
    members = set(metadata["workspace_members"])
    packages = []
    for package in metadata["packages"]:
        if package["id"] not in members:
            continue
        count = sum("test" in target["kind"] for target in package["targets"])
        if count:
            packages.append((package["name"], count))
    packages.sort(key=lambda entry: (-entry[1], entry[0]))
    if len(packages) < shard_count:
        raise ValueError("more shards than packages with integration tests")

    shards: list[list[str]] = [[] for _ in range(shard_count)]
    weights = [0] * shard_count
    for name, count in packages:
        shard = min(range(shard_count), key=lambda index: (weights[index], index))
        shards[shard].append(name)
        weights[shard] += count
    return shards


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", type=int, default=4)
    parser.add_argument("--shard", type=int)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    if args.run and args.shard is None:
        parser.error("--run requires --shard")
    metadata = json.loads(subprocess.check_output(
        ["cargo", "metadata", "--no-deps", "--format-version", "1"], text=True
    ))
    shards = assign_packages(metadata, args.shards)
    if args.shard is not None and not 0 <= args.shard < args.shards:
        parser.error("--shard is outside the matrix")
    for index, packages in enumerate(shards):
        print(f"shard {index}: {', '.join(packages)}", flush=True)
    if not args.run:
        return
    packages = shards[args.shard]
    command = ["cargo", "test", "--locked", "--tests", "--quiet"]
    for package in packages:
        command.extend(["-p", package])
    subprocess.run(command, check=True)
    if "chematic-crystal" in packages:
        subprocess.run([
            "cargo", "test", "--locked", "-p", "chematic-crystal",
            "--test", "serde_roundtrip", "--features", "serde", "--quiet",
        ], check=True)


if __name__ == "__main__":
    main()
