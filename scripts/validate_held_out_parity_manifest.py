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


def validate_report(name: str, report_path: Path, target_version: str) -> list[str]:
    """Validate arithmetic invariants in a measured report.

    The reports intentionally use operation-specific field names, so this is
    a small explicit schema gate rather than a permissive duck-typing check.
    It catches stale or hand-edited headline counts without declaring parity
    successful merely because a JSON file exists.
    """
    errors: list[str] = []
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"operation {name!r} report is not valid JSON: {exc}"]
    if not isinstance(report, dict):
        return [f"operation {name!r} report must be a JSON object"]
    report_version = report.get("target_version")
    if report_version is not None and report_version != target_version:
        errors.append(
            f"operation {name!r} report target_version {report_version!r} "
            f"does not match {target_version!r}"
        )

    if name == "morgan_ecfp":
        total = report.get("total_inputs")
        success = report.get("success_count")
        error_count = report.get("error_count")
        exact = report.get("exact_match_among_success")
        if not all(isinstance(value, int) for value in (total, success, error_count, exact)):
            return [f"operation {name!r} headline counts must be integers"]
        if total != success + error_count:
            errors.append(f"operation {name!r} total_inputs != success_count + error_count")
        if exact != success:
            errors.append(f"operation {name!r} exact_match must equal success_count")
        if report.get("gate_passed") is not True or report.get("gate_failures") != []:
            errors.append(f"operation {name!r} parity gate is not passed")
        return errors

    if name == "descriptors":
        parsed = report.get("parsed")
        if not isinstance(parsed, int):
            return [f"operation {name!r} parsed must be an integer"]
        fields = report.get("fields")
        if not isinstance(fields, dict):
            return [f"operation {name!r} fields must be an object"]
        for field, values in fields.items():
            if not isinstance(values, dict) or not all(
                isinstance(values.get(key), int) for key in ("matches", "mismatches")
            ) or values["matches"] + values["mismatches"] != parsed:
                errors.append(f"operation {name!r} field {field!r} counts do not equal parsed")
        return errors

    rows = report.get("rows")
    if not isinstance(rows, int):
        return [f"operation {name!r} rows must be an integer"]
    if name == "standardization":
        valid = report.get("valid")
        error_count = report.get("errors")
        matches = report.get("matches")
        mismatches = report.get("mismatches")
        if not all(isinstance(value, int) for value in (valid, error_count, matches, mismatches)):
            return [f"operation {name!r} standardization counts must be integers"]
        if valid + error_count != rows:
            errors.append(f"operation {name!r} valid + errors != rows")
        if matches + mismatches != valid:
            errors.append(f"operation {name!r} matches + mismatches != valid")
        return errors

    schematic_failures = report.get("chematic_failures")
    rdkit_failures = report.get("rdkit_failures")
    exact_matches = report.get("exact_matches")
    exact_mismatches = report.get("exact_mismatches")
    if not all(
        isinstance(value, int)
        for value in (schematic_failures, rdkit_failures, exact_matches, exact_mismatches)
    ):
        return [f"operation {name!r} comparison counts must be integers"]
    common_success = rows - schematic_failures - rdkit_failures
    if exact_matches + exact_mismatches != common_success:
        errors.append(
            f"operation {name!r} exact matches + mismatches != common successful rows"
        )
    return errors


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
            else:
                errors.extend(validate_report(name, root / report, document["target_version"]))
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
