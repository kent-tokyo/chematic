#!/usr/bin/env python3
"""Fail-closed validator for the checked-in streaming contract matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_FORMATS = ["sdf", "mol", "xyz", "extxyz", "v3000", "mol2", "cml", "cdxml", "mmcif", "pdb"]


def fail(message: str) -> None:
    raise SystemExit(f"streaming matrix invalid: {message}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "manifest",
        type=Path,
        nargs="?",
        help="matrix JSON (defaults to the current workspace-versioned result)",
    )
    args = parser.parse_args()
    manifest = args.manifest or Path(
        "validation/results"
    ) / f"cross-engine-matrix-v{workspace_version(ROOT)}.json"
    path = manifest if manifest.is_absolute() else ROOT / manifest
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read JSON: {error}")
    if data.get("schema_version") != 1 or data.get("target_version") != workspace_version(ROOT):
        fail("schema or target version is stale")
    repeats = data.get("repeats")
    if not isinstance(repeats, int) or repeats <= 0:
        fail("repeats must be a positive integer")
    if data.get("formats") != EXPECTED_FORMATS:
        fail(f"formats must be {EXPECTED_FORMATS}")
    rows = data.get("rows")
    if not isinstance(rows, list) or len(rows) != len(EXPECTED_FORMATS):
        fail("rows must contain exactly one entry per format")
    seen: set[str] = set()
    for row in rows:
        fmt = row.get("format") if isinstance(row, dict) else None
        if fmt not in EXPECTED_FORMATS or fmt in seen:
            fail(f"invalid or duplicate format row: {fmt!r}")
        seen.add(fmt)
        fixture = row.get("fixture")
        if not isinstance(fixture, dict):
            fail(f"{fmt}: missing fixture")
        fixture_path = fixture.get("path")
        fixture_file = ROOT / fixture_path if isinstance(fixture_path, str) else None
        if fixture_file is None or not fixture_file.is_file():
            fail(f"{fmt}: fixture does not exist: {fixture_path!r}")
        payload = fixture_file.read_bytes()
        if fixture.get("bytes") != len(payload) or fixture.get("sha256") != hashlib.sha256(payload).hexdigest():
            fail(f"{fmt}: fixture hash/size mismatch")
        expected = row.get("expected_records")
        if not isinstance(expected, int) or expected <= 0:
            fail(f"{fmt}: invalid expected_records")
        engines = row.get("engines")
        if not isinstance(engines, list) or not engines:
            fail(f"{fmt}: missing engines")
        names: set[str] = set()
        for engine in engines:
            if not isinstance(engine, dict):
                fail(f"{fmt}: malformed engine row")
            name = engine.get("engine")
            if not isinstance(name, str) or name in names:
                fail(f"{fmt}: duplicate/missing engine name")
            names.add(name)
            if engine.get("records") != expected or engine.get("failures") != 0:
                fail(f"{fmt}/{name}: record/failure contract mismatch")
            if not isinstance(engine.get("boundary"), str) or not engine["boundary"].strip():
                fail(f"{fmt}/{name}: missing comparison boundary")
        if "chematic" not in names:
            fail(f"{fmt}: chematic result is required")
    if seen != set(EXPECTED_FORMATS):
        fail("format rows are incomplete")
    print(f"streaming matrix OK: {len(rows)} formats, {repeats} repeats")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
