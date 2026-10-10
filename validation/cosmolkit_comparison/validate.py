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


def corpus(
    corpus_path: Path | None = None, manifest_path: Path | None = None
) -> tuple[str, dict[str, str]]:
    corpus_path = CORPUS if corpus_path is None else corpus_path
    manifest_path = MANIFEST if manifest_path is None else manifest_path
    rows = [json.loads(line) for line in corpus_path.read_text().splitlines() if line.strip()]
    manifest = json.loads(manifest_path.read_text())
    actual_hash = hashlib.sha256(corpus_path.read_bytes()).hexdigest()
    expected_rows = manifest["records"]
    if manifest["corpus"] != corpus_path.name:
        raise ValueError("corpus manifest names a different corpus file")
    if manifest["sha256"] != actual_hash:
        raise ValueError(f"corpus manifest sha256 does not match {corpus_path.name}")
    if rows != expected_rows:
        raise ValueError(f"{corpus_path.name} rows do not match {manifest_path.name}")
    return actual_hash, {row["id"]: row["smiles"] for row in rows}


def validate(
    path: Path,
    expected_engine: str | None = None,
    corpus_path: Path | None = None,
    manifest_path: Path | None = None,
) -> list[str]:
    errors: list[str] = []
    corpus_path = CORPUS if corpus_path is None else corpus_path
    manifest_path = MANIFEST if manifest_path is None else manifest_path
    try:
        expected_hash, expected_smiles = corpus(corpus_path, manifest_path)
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
        if expected_engine is not None and record["engine"] != expected_engine:
            errors.append(f"line {number}: engine does not match requested {expected_engine!r}")
        if record["corpus_sha256"] != expected_hash:
            errors.append(f"line {number}: corpus_sha256 does not match {corpus_path.name}")
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
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    args = parser.parse_args()
    errors = validate(args.results, corpus_path=args.corpus, manifest_path=args.manifest)
    if errors:
        print(json.dumps({"valid": False, "errors": errors}, indent=2, sort_keys=True))
        raise SystemExit(1)
    print(json.dumps({"valid": True, "records": len([line for line in args.results.read_text().splitlines() if line.strip()])}, sort_keys=True))


if __name__ == "__main__":
    main()
