#!/usr/bin/env python3
"""Recompute descriptor scorecard counters from schema-v2 raw rows.

The summary fields in a diagnostic artifact are claims.  This validator
rebuilds their row/field accounting from ``raw_rows`` so a dropped row,
silently changed status, or stale summary cannot become a positive result.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path


REQUIRED_FIELDS = (
    "molecular_weight",
    "hba",
    "hbd",
    "tpsa",
    "logp",
    "molar_refractivity",
    "fsp3",
    "aromatic_ring_count",
)


def validate(path: Path) -> list[str]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return [f"cannot read {path}: {error}"]
    errors: list[str] = []
    raw_rows = document.get("raw_rows")
    if not isinstance(raw_rows, list):
        return [f"{path}: raw_rows must be a list"]
    rows = document.get("rows")
    if not isinstance(rows, int) or isinstance(rows, bool) or rows != len(raw_rows):
        errors.append(f"{path}: rows does not equal raw_rows length")
    raw_payload = json.dumps(raw_rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    if document.get("raw_rows_sha256") != hashlib.sha256(raw_payload).hexdigest():
        errors.append(f"{path}: raw_rows_sha256 mismatch")

    ids = [row.get("index") for row in raw_rows if isinstance(row, dict)]
    if len(ids) != len(set(ids)) or set(ids) != set(range(len(raw_rows))):
        errors.append(f"{path}: raw row indices are missing or duplicated")

    recomputed = {
        field: {"parsed": 0, "matches": 0, "strict_matches": 0, "mismatches": 0}
        for field in REQUIRED_FIELDS
    }
    for row in raw_rows:
        if not isinstance(row, dict) or row.get("status") != "ok":
            continue
        fields = row.get("fields")
        if not isinstance(fields, dict) or set(fields) != set(REQUIRED_FIELDS):
            errors.append(f"{path}: row {row.get('index')} fields do not match the contract")
            continue
        for name in REQUIRED_FIELDS:
            value = fields[name]
            if not isinstance(value, dict) or value.get("status") != "ok":
                continue
            delta = value.get("absolute_error")
            if not isinstance(delta, (int, float)) or isinstance(delta, bool) or not math.isfinite(float(delta)):
                errors.append(f"{path}: row {row.get('index')} field {name} has invalid absolute_error")
                continue
            stats = recomputed[name]
            stats["parsed"] += 1
            if value.get("matched") is True:
                stats["matches"] += 1
            if value.get("strict_passed") is True:
                stats["strict_matches"] += 1
            else:
                stats["mismatches"] += 1

    summary_fields = document.get("fields")
    if not isinstance(summary_fields, dict) or set(summary_fields) != set(REQUIRED_FIELDS):
        errors.append(f"{path}: summary fields do not match the contract")
    else:
        for name, expected in recomputed.items():
            actual = summary_fields[name]
            if not isinstance(actual, dict):
                errors.append(f"{path}: summary for {name} is not an object")
                continue
            for key, expected_value in expected.items():
                if actual.get(key) != expected_value:
                    errors.append(
                        f"{path}: {name}.{key}={actual.get(key)!r}, expected {expected_value!r}"
                    )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("diagnostics", type=Path, nargs="+")
    args = parser.parse_args()
    errors = [error for path in args.diagnostics for error in validate(path)]
    if errors:
        print("descriptor raw accounting: BLOCKED", file=sys.stderr)
        print("\n".join(f"- {error}" for error in errors), file=sys.stderr)
        return 1
    print(f"descriptor raw accounting: PASS ({len(args.diagnostics)} artifacts)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
