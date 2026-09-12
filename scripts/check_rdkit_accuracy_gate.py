#!/usr/bin/env python3
"""Fail-closed A0 validator for the RDKit accuracy comparison contract.

This is deliberately independent of RDKit and the chematic runtime. It checks
whether an evidence record is complete enough to be considered; it does not
produce chemical results and cannot turn shared-success agreement into full
coverage.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "manifests" / "rdkit_accuracy_v1.json"


def load(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON value must be an object")
    return value


def finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def nonnegative_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def positive_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


def validate_holdout(holdout: dict, required: list[str]) -> list[str]:
    errors: list[str] = []
    rows = holdout.get("rows")
    fields = holdout.get("fields")
    if not positive_int(rows):
        return ["holdout must contain a positive row count"]
    if not isinstance(fields, dict) or set(fields) != set(required):
        errors.append("holdout must contain exactly the manifest descriptor fields")
        return errors
    for name in required:
        record = fields[name]
        if not isinstance(record, dict):
            errors.append(f"holdout {name}: field record must be an object")
            continue
        if record.get("rows") != rows:
            errors.append(f"holdout {name}: row count is not {rows}")
        checks = record.get("checks")
        if not isinstance(checks, list) or len(checks) != rows:
            errors.append(f"holdout {name}: checks do not cover every row")
            continue
        if record.get("failed") != 0 or record.get("strict_passed") != rows:
            errors.append(f"holdout {name}: strict failures are present")
        seen_ids: set[str] = set()
        for item in checks:
            if (
                not isinstance(item, dict)
                or not isinstance(item.get("id"), str)
                or not item["id"]
                or item["id"] in seen_ids
            ):
                errors.append(f"holdout {name}: missing or duplicate row id")
                continue
            seen_ids.add(item["id"])
            for key in ("chematic", "rdkit", "absolute_error"):
                if item.get(key) is not None and not finite(item[key]):
                    errors.append(f"holdout {name}: {key} is non-finite")
    return errors


def validate_raw_rows(diagnostics: dict, required: list[str], tolerances: dict) -> list[str]:
    """Recompute the aggregate counters from schema-v2 raw rows."""
    raw_rows = diagnostics.get("raw_rows")
    rows = diagnostics.get("rows")
    errors: list[str] = []
    if not isinstance(raw_rows, list) or not positive_int(rows) or len(raw_rows) != rows:
        return ["schema-v2 raw_rows must cover every diagnostics row"]
    payload = json.dumps(raw_rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    if diagnostics.get("raw_rows_sha256") != hashlib.sha256(payload).hexdigest():
        errors.append("raw_rows_sha256 does not match raw_rows")
    seen: set[int] = set()
    parsed_rows = invalid_rows = unsupported = 0
    field_counts = {
        name: {"parsed": 0, "matches": 0, "strict_matches": 0, "mismatches": 0}
        for name in required
    }
    for row in raw_rows:
        if not isinstance(row, dict) or not isinstance(row.get("index"), int) or isinstance(row.get("index"), bool):
            errors.append("raw row index must be an integer")
            continue
        index = row["index"]
        if index in seen:
            errors.append(f"duplicate raw row index: {index}")
        seen.add(index)
        if not isinstance(row.get("smiles"), str) or not row["smiles"]:
            errors.append(f"raw row {index}: smiles is missing")
        status = row.get("status")
        if status == "invalid_input":
            invalid_rows += 1
            continue
        if status != "ok":
            errors.append(f"raw row {index}: unknown status {status!r}")
            continue
        parsed_rows += 1
        fields = row.get("fields")
        if not isinstance(fields, dict) or set(fields) != set(required):
            errors.append(f"raw row {index}: fields do not match the manifest")
            continue
        for name in required:
            item = fields[name]
            if not isinstance(item, dict):
                errors.append(f"raw row {index} {name}: field record is not an object")
                continue
            if item.get("status") == "unsupported":
                unsupported += 1
                continue
            if item.get("status") != "ok":
                errors.append(f"raw row {index} {name}: unknown field status")
                continue
            if any(not finite(item.get(key)) for key in ("chematic", "rdkit", "absolute_error")):
                errors.append(f"raw row {index} {name}: non-finite value")
                continue
            delta = abs(float(item["chematic"]) - float(item["rdkit"]))
            if abs(delta - float(item["absolute_error"])) > 1e-12:
                errors.append(f"raw row {index} {name}: absolute_error is inconsistent")
            expected_match = delta <= diagnostics["fields"][name]["tolerance"]
            expected_strict = delta <= diagnostics["fields"][name]["strict_tolerance"]
            if item.get("matched") is not expected_match or item.get("strict_passed") is not expected_strict:
                errors.append(f"raw row {index} {name}: pass flags are inconsistent")
            field_counts[name]["parsed"] += 1
            field_counts[name]["matches"] += int(expected_match)
            field_counts[name]["strict_matches"] += int(expected_strict)
            field_counts[name]["mismatches"] += int(not expected_strict)
    if seen != set(range(rows)):
        errors.append("raw row indices are not a complete 0..rows-1 sequence")
    if parsed_rows != diagnostics.get("parsed") or invalid_rows != diagnostics.get("parse_failures"):
        errors.append("raw row status accounting differs from diagnostics")
    if unsupported != diagnostics.get("unsupported_values"):
        errors.append("raw unsupported-value accounting differs from diagnostics")
    for name in required:
        actual = diagnostics["fields"].get(name, {})
        if any(actual.get(key) != value for key, value in field_counts[name].items()):
            errors.append(f"{name}: aggregate counters differ from raw rows")
    return errors


def validate(manifest: dict, diagnostics: dict, holdout: dict) -> list[str]:
    errors: list[str] = []
    fields = diagnostics.get("fields")
    required = manifest["descriptor_contract"]["required_fields"]
    if not isinstance(fields, dict):
        return ["diagnostics.fields must be an object"]
    if diagnostics.get("schema_version") not in (1, 2):
        errors.append("diagnostics schema_version must be 1 or 2")
    if not isinstance(diagnostics.get("profile"), str) or not diagnostics["profile"]:
        errors.append("diagnostics profile is missing")
    if not isinstance(diagnostics.get("chematic_version"), str) or not diagnostics["chematic_version"]:
        errors.append("chematic version is missing")
    if diagnostics.get("corpus") != manifest["descriptor_contract"]["corpus"]:
        errors.append("diagnostics corpus does not match manifest")
    if set(fields) != set(required):
        errors.append("diagnostics fields must exactly match the manifest descriptor fields")
    missing = [name for name in required if name not in fields]
    if missing:
        errors.append(f"missing required descriptor fields: {', '.join(missing)}")
    rows = diagnostics.get("rows")
    parsed = diagnostics.get("parsed")
    failures = diagnostics.get("parse_failures")
    if not all(nonnegative_int(value) for value in (rows, parsed, failures)):
        errors.append("diagnostics rows/parsed/parse_failures must be non-negative integers")
    elif parsed + failures != rows:
        errors.append(f"parser accounting mismatch: parsed={parsed} failures={failures} rows={rows}")
    elif positive_int(rows) and parsed == 0:
        errors.append("diagnostics contains no successfully parsed rows")
    if not diagnostics.get("rdkit_version"):
        errors.append("RDKit version is missing")
    if not diagnostics.get("corpus"):
        errors.append("corpus identity is missing")
    if diagnostics.get("unsupported_values") != 0:
        errors.append("diagnostics contains unsupported descriptor values")

    tolerances = manifest["descriptor_contract"]["strict_tolerances"]
    for name in required:
        record = fields.get(name)
        if not isinstance(record, dict):
            errors.append(f"{name}: evidence record is missing or not an object")
            continue
        for key in ("parsed", "matches", "strict_matches", "mismatches"):
            if not nonnegative_int(record.get(key)):
                errors.append(f"{name}: {key} must be a non-negative integer")
        if isinstance(rows, int) and record.get("parsed") != rows:
            errors.append(f"{name}: coverage is incomplete ({record.get('parsed')}/{rows})")
        if record.get("matches", 0) > record.get("parsed", 0):
            errors.append(f"{name}: matches exceed parsed rows")
        if record.get("mismatches") != record.get("parsed", -1) - record.get("matches", -2):
            errors.append(f"{name}: matches/mismatches accounting is inconsistent")
        if record.get("strict_matches") != record.get("parsed", -1):
            errors.append(f"{name}: strict mismatch or incomplete coverage")
        if record.get("strict_tolerance") != tolerances[name]:
            errors.append(f"{name}: strict_tolerance does not match manifest")
        for key in ("mae", "median_abs_error", "p95_abs_error", "max_abs_error"):
            if not finite(record.get(key)):
                errors.append(f"{name}: {key} is absent or non-finite")
        field_errors = record.get("errors")
        if not isinstance(field_errors, list):
            errors.append(f"{name}: errors must be an array")
            field_errors = []
        for error in field_errors:
            if not isinstance(error, dict) or not error.get("smiles") or not error.get("cause"):
                errors.append(f"{name}: every error must include smiles and cause")

    if diagnostics.get("schema_version") == 2:
        errors.extend(validate_raw_rows(diagnostics, required, tolerances))
    errors.extend(validate_holdout(holdout, required))
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("diagnostics", type=Path)
    parser.add_argument("holdout", type=Path)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    try:
        errors = validate(load(args.manifest), load(args.diagnostics), load(args.holdout))
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"A0 accuracy gate input error: {exc}", file=sys.stderr)
        return 2
    if errors:
        print("A0 RDKit accuracy gate: BLOCKED")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("A0 RDKit accuracy gate: PASS (evidence contract complete)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
