#!/usr/bin/env python3
"""Validate comparison JSONL records using the checked-in contract.

This intentionally has no third-party dependency so the comparison gate can run
before RDKit or an external competitor is installed.
"""

import argparse
import hashlib
import json
from pathlib import Path

CORPUS = Path(__file__).with_name("smoke_corpus.jsonl")
MANIFEST = Path(__file__).with_name("corpus_manifest.json")
STATUSES = {"ok", "parse_error", "unsupported", "error"}


def corpus() -> tuple[str, dict[str, str]]:
    rows = [json.loads(line) for line in CORPUS.read_text().splitlines() if line.strip()]
    manifest = json.loads(MANIFEST.read_text())
    actual_hash = hashlib.sha256(CORPUS.read_bytes()).hexdigest()
    expected_rows = manifest["records"]
    if manifest["corpus"] != CORPUS.name:
        raise ValueError("corpus manifest names a different corpus file")
    if manifest["sha256"] != actual_hash:
        raise ValueError("corpus manifest sha256 does not match smoke_corpus.jsonl")
    if rows != expected_rows:
        raise ValueError("smoke_corpus.jsonl rows do not match corpus_manifest.json")
    return actual_hash, {row["id"]: row["smiles"] for row in rows}


def validate(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        expected_hash, expected_smiles = corpus()
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return [f"corpus manifest invalid: {exc}"]
    seen: set[str] = set()
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError as exc:
            errors.append(f"line {number}: invalid JSON: {exc.msg}")
            continue
        required = {"schema_version", "engine", "engine_version", "corpus_sha256",
                    "id", "smiles", "status", "operations"}
        missing = required - record.keys()
        if missing:
            errors.append(f"line {number}: missing fields: {sorted(missing)}")
            continue
        if record["schema_version"] != 1:
            errors.append(f"line {number}: unsupported schema_version {record['schema_version']!r}")
        if record["corpus_sha256"] != expected_hash:
            errors.append(f"line {number}: corpus_sha256 does not match smoke_corpus.jsonl")
        if record["status"] not in STATUSES:
            errors.append(f"line {number}: invalid record status {record['status']!r}")
        ident = record["id"]
        if ident in seen:
            errors.append(f"line {number}: duplicate id {ident!r}")
        seen.add(ident)
        if ident in expected_smiles and record["smiles"] != expected_smiles[ident]:
            errors.append(f"line {number}: smiles does not match corpus record {ident!r}")
        if not isinstance(record["operations"], dict):
            errors.append(f"line {number}: operations must be an object")
            continue
        for operation, result in record["operations"].items():
            if not isinstance(result, dict) or result.get("status") not in STATUSES:
                errors.append(f"line {number}: invalid status for operation {operation!r}")
    expected_ids = set(expected_smiles)
    missing_ids = expected_ids - seen
    extra_ids = seen - expected_ids
    if missing_ids:
        errors.append(f"missing corpus records: {sorted(missing_ids)}")
    if extra_ids:
        errors.append(f"unknown corpus records: {sorted(extra_ids)}")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    errors = validate(args.results)
    if errors:
        print(json.dumps({"valid": False, "errors": errors}, indent=2, sort_keys=True))
        raise SystemExit(1)
    print(json.dumps({"valid": True, "records": len([line for line in args.results.read_text().splitlines() if line.strip()])}, sort_keys=True))


if __name__ == "__main__":
    main()
