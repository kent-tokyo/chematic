#!/usr/bin/env python3
"""Fail closed if the frozen trust holdouts or their policy drift."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "validation/trust_holdouts/v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    manifest = json.loads((BASE / "manifest.json").read_text())
    assert manifest["schema"] == "chematic.trust-holdout-manifest.v1"
    source = ROOT / manifest["source"]["path"]
    assert sha(source) == manifest["source"]["sha256"]
    forbidden = {"expected", "expected_output", "result", "oracle_output"}
    for name, contract in manifest["lanes"].items():
        path = BASE / name
        assert sha(path) == contract["sha256"], name
        rows = [json.loads(line) for line in path.read_text().splitlines() if line]
        assert len(rows) == contract["rows"], name
        assert contract["contains_expected_outputs"] is False
        assert all(not (set(row) & forbidden) for row in rows), name
    print("Frozen trust holdouts: hashes, row counts, and outcome-blind policy passed")


if __name__ == "__main__":
    main()
