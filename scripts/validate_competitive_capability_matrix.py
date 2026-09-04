#!/usr/bin/env python3
"""Validate the operation-level capability inventory without network access."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
MATRIX = ROOT / "validation" / "competitive_capability_matrix.json"
ALLOWED_STATES = {"supported", "partial", "unsupported", "not_measured"}
REQUIRED_ENGINES = {"chematic", "rdkit", "openbabel", "cdk", "sdfrust", "kekule"}


def fail(message: str) -> None:
    print(f"competitive capability matrix: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    try:
        data = json.loads(MATRIX.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read JSON: {exc}")
    if data.get("schema_version") != 1:
        fail("schema_version must be 1")
    version = data.get("target_version", "")
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        fail("target_version must be semantic version")
    if "cosmolkit" not in data.get("excluded_comparators", []):
        fail("cosmolkit must remain explicitly excluded from this matrix")
    engines = data.get("engines", [])
    engine_ids = [engine.get("id") for engine in engines]
    if set(engine_ids) != REQUIRED_ENGINES or len(engine_ids) != len(set(engine_ids)):
        fail("engine inventory must contain exactly the six in-scope engines")
    operations = data.get("operations", [])
    if not operations:
        fail("at least one operation is required")
    seen = set()
    for operation in operations:
        operation_id = operation.get("id")
        if not operation_id or operation_id in seen:
            fail(f"operation ids must be non-empty and unique: {operation_id}")
        seen.add(operation_id)
        statuses = operation.get("status")
        if set(statuses or {}) != REQUIRED_ENGINES:
            fail(f"{operation_id} must state every in-scope engine")
        invalid = set(statuses.values()) - ALLOWED_STATES
        if invalid:
            fail(f"{operation_id} has invalid states: {sorted(invalid)}")
    print(f"Competitive capability matrix OK: {len(operations)} operations, {len(engines)} engines")


if __name__ == "__main__":
    main()
