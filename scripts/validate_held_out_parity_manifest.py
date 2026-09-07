#!/usr/bin/env python3
"""Validate the v1.0.8 held-out parity coverage ledger without dependencies.

The ledger is deliberately conservative: ``not_measured`` is a valid state,
but it cannot be mistaken for a completed parity report. A measured operation
must point to an existing report and name its comparison boundary.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "validation" / "held_out_parity_manifest.json"
REQUIRED_OPERATIONS = {
    "morgan_ecfp",
    "maccs",
    "topological",
    "torsion",
    "descriptors",
    "standardization",
}
STATUSES = {"measured", "not_measured"}
BINDINGS = {"rust", "python", "node", "wasm"}


def validate(document: dict, root: Path = ROOT) -> list[str]:
    errors: list[str] = []
    if document.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if document.get("target_version") != "1.0.8":
        errors.append("target_version must be 1.0.8")
    operations = document.get("operations")
    if not isinstance(operations, dict):
        return ["operations must be an object"]
    missing = REQUIRED_OPERATIONS - set(operations)
    extra = set(operations) - REQUIRED_OPERATIONS
    if missing:
        errors.append(f"missing operations: {sorted(missing)}")
    if extra:
        errors.append(f"unknown operations: {sorted(extra)}")

    for name, record in operations.items():
        if not isinstance(record, dict):
            errors.append(f"operation {name!r} must be an object")
            continue
        status = record.get("status")
        if status not in STATUSES:
            errors.append(f"operation {name!r} has invalid status {status!r}")
        bindings = record.get("bindings")
        if not isinstance(bindings, list) or not bindings or any(b not in BINDINGS for b in bindings):
            errors.append(f"operation {name!r} has invalid bindings")
        boundary = record.get("boundary")
        if not isinstance(boundary, str) or not boundary.strip():
            errors.append(f"operation {name!r} is missing boundary")
        if status == "measured":
            report = record.get("report")
            if not isinstance(report, str) or not report:
                errors.append(f"measured operation {name!r} is missing report")
            elif not (root / report).is_file():
                errors.append(f"measured operation {name!r} report does not exist: {report}")
        elif "report" in record:
            errors.append(f"not_measured operation {name!r} must not declare report")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    try:
        document = json.loads(args.manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"manifest read failure: {exc}", file=sys.stderr)
        return 1
    errors = validate(document)
    if errors:
        print("held-out parity manifest failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    measured = sum(r["status"] == "measured" for r in document["operations"].values())
    print(f"held-out parity manifest OK: {measured} measured, {len(REQUIRED_OPERATIONS) - measured} not_measured")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
