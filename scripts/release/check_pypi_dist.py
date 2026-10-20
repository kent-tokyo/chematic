#!/usr/bin/env python3
"""Check release wheel completeness and PyPI hashes for the exact build."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.request
from pathlib import Path

PYTHONS = ("39", "310", "311", "312", "313")
PLATFORMS = ("linux-x86_64", "linux-aarch64", "macos-x86_64", "macos-arm64", "windows-amd64")
ATTESTATION_SUFFIX = ".publish.attestation"


def platform_for(tag: str) -> str:
    if tag.startswith("manylinux_") and tag.endswith("_x86_64"):
        return "linux-x86_64"
    if tag.startswith("manylinux_") and tag.endswith("_aarch64"):
        return "linux-aarch64"
    if tag.startswith("macosx_") and tag.endswith("_x86_64"):
        return "macos-x86_64"
    if tag.startswith("macosx_") and tag.endswith("_arm64"):
        return "macos-arm64"
    if tag == "win_amd64":
        return "windows-amd64"
    raise ValueError(f"unexpected wheel platform: {tag}")


def check_layout(names: set[str], version: str) -> None:
    sdist = f"chematic-{version}.tar.gz"
    if sdist not in names:
        raise ValueError(f"missing sdist: {sdist}")
    pattern = re.compile(rf"chematic-{re.escape(version)}-cp(39|310|311|312|313)-cp\1-(.+)\.whl\Z")
    found = set()
    for name in names - {sdist}:
        match = pattern.fullmatch(name)
        if not match:
            raise ValueError(f"unexpected distribution: {name}")
        key = (match.group(1), platform_for(match.group(2)))
        if key in found:
            raise ValueError(f"duplicate wheel target: {key}")
        found.add(key)
    expected = {(python, platform) for python in PYTHONS for platform in PLATFORMS}
    if found != expected:
        raise ValueError(f"wheel matrix missing={sorted(expected - found)} extra={sorted(found - expected)}")


def local_hashes(directory: Path) -> dict[str, str]:
    files = {path.name: path for path in directory.iterdir() if path.is_file()}
    distributions = {
        name: path for name, path in files.items() if name.endswith((".whl", ".tar.gz"))
    }
    attestations = {name for name in files if name.endswith(ATTESTATION_SUFFIX)}
    unexpected = set(files) - set(distributions) - attestations
    if unexpected:
        raise ValueError(f"unexpected release files: {sorted(unexpected)}")
    orphaned = {
        name for name in attestations if name.removesuffix(ATTESTATION_SUFFIX) not in distributions
    }
    if orphaned:
        raise ValueError(f"orphaned publish attestations: {sorted(orphaned)}")
    return {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in distributions.items()
    }


def published_hashes(version: str) -> dict[str, str]:
    url = (
        f"https://pypi.org/pypi/chematic/{version}/json"
        f"?release-verification={time.time_ns()}"
    )
    request = urllib.request.Request(url, headers={"Cache-Control": "no-cache"})
    with urllib.request.urlopen(request, timeout=20) as response:
        payload = json.load(response)
    return {entry["filename"]: entry["digests"]["sha256"] for entry in payload["urls"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dist", type=Path, required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--verify-published", action="store_true")
    args = parser.parse_args()
    local = local_hashes(args.dist)
    check_layout(set(local), args.version)
    print(f"Complete PyPI distribution matrix: {len(local)} files")
    if not args.verify_published:
        return
    for attempt in range(20):
        try:
            published = published_hashes(args.version)
        except (OSError, ValueError):
            published = {}
        if published == local:
            print(f"PyPI hashes match all {len(local)} built distributions")
            return
        if set(published) == set(local):
            raise SystemExit("PyPI file names match, but one or more SHA-256 hashes differ")
        if attempt < 19:
            time.sleep(15)
    raise SystemExit(f"PyPI file set differs: missing={sorted(set(local) - set(published))} extra={sorted(set(published) - set(local))}")


if __name__ == "__main__":
    main()
