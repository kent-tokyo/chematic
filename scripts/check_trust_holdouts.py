#!/usr/bin/env python3
"""Fail closed if frozen trust holdouts or their post-candidate policy drift."""

from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "validation/trust_holdouts/v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads((BASE / "manifest.json").read_text())
    assert manifest["schema"] == "chematic.trust-holdout-manifest.v1"
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
        path = BASE / name
        assert sha(path) == contract["sha256"], name
        rows = [json.loads(line) for line in path.read_text().splitlines() if line]
        assert len(rows) == contract["rows"], name
        assert contract["contains_expected_outputs"] is False
        assert all(not (set(row) & forbidden) for row in rows), name
    print("Frozen trust holdouts: post-candidate provenance, hashes, rows, and policy passed")


if __name__ == "__main__":
    main()
