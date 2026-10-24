#!/usr/bin/env python3
"""Fail closed if frozen trust holdouts or their post-candidate policy drift."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "validation/trust_holdouts"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_ids(base: Path) -> set[str]:
    result: set[str] = set()
    for path in base.glob("*.jsonl"):
        for line in path.read_text().splitlines():
            row = json.loads(line)
            for value in row.values():
                if isinstance(value, dict) and "source_id" in value:
                    result.add(value["source_id"])
                elif isinstance(value, list):
                    result.update(
                        item["source_id"] for item in value
                        if isinstance(item, dict) and "source_id" in item
                    )
            if "source_id" in row:
                result.add(row["source_id"])
    return result


def check_version(version: str) -> set[str]:
    base = BASE / version
    manifest = json.loads((base / "manifest.json").read_text())
    assert manifest["schema"] == f"chematic.trust-holdout-manifest.{version}"
    candidate = manifest["candidate_frozen_commit"]
    assert len(candidate) == 40
    subprocess.run(["git", "cat-file", "-e", f"{candidate}^{{commit}}"], cwd=ROOT, check=True)
    candidate_time = datetime.fromisoformat(
        subprocess.check_output(
            ["git", "show", "-s", "--format=%cI", candidate], cwd=ROOT, text=True
        ).strip()
    )
    for source in manifest["sources"].values():
        acquired_time = datetime.fromisoformat(source["acquired_at"].replace("Z", "+00:00"))
        assert acquired_time > candidate_time
    assert all(source["acquired_after_candidate_freeze"] for source in manifest["sources"].values())
    assert manifest["sources"]["nci"]["decompressed_sha256"] == (
        "031736b438a9596d89cf11802fa97934c3878ebbf021924c417a57e02c7930cd"
    )
    forbidden = {"expected", "expected_output", "result", "oracle_output"}
    for name, contract in manifest["lanes"].items():
        path = base / name
        assert sha(path) == contract["sha256"], name
        rows = [json.loads(line) for line in path.read_text().splitlines() if line]
        assert len(rows) == contract["rows"], name
        assert contract["contains_expected_outputs"] is False
        assert all(not (set(row) & forbidden) for row in rows), name
    return source_ids(base)


def main() -> None:
    v1_ids = check_version("v1")
    versions = ["v1"]
    if (BASE / "v2" / "manifest.json").is_file():
        v2_ids = check_version("v2")
        assert not (v1_ids & v2_ids), "v2 must not reuse v1 source rows"
        versions.append("v2")
    print(f"Frozen trust holdouts {', '.join(versions)}: provenance, hashes, rows, and policy passed")


if __name__ == "__main__":
    main()
