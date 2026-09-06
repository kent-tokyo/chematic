#!/usr/bin/env python3
"""Validate a versioned operation scorecard without third-party packages.

The validator is intentionally conservative: unsupported, failed, missing, or
not-measured rows cannot be used as a positive claim. A scorecard may still
contain those rows, but any entry under ``claims`` must point to an ``ok`` row.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "release-metadata" / "v1.0.8.json"
BAD_STATUSES = {"unsupported", "failure", "failed", "not_measured", "missing"}
GOOD_STATUSES = {"ok", "supported", "match"}


def validate(document: dict, expected_version: str) -> list[str]:
    errors: list[str] = []
    if document.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    actual_version = document.get("target_version") or document.get("release", {}).get("version")
    if actual_version != expected_version:
        errors.append(f"target_version must be {expected_version}, got {actual_version!r}")
    corpus = document.get("corpus_sha256")
    if not isinstance(corpus, str) or not re.fullmatch(r"[0-9a-f]{64}", corpus):
        errors.append("corpus_sha256 must be a lowercase SHA-256 digest")

    engines = document.get("engines")
    if not isinstance(engines, dict) or not engines:
        errors.append("engines must be a non-empty mapping")
    else:
        for engine, metadata in engines.items():
            if not isinstance(metadata, dict) or not metadata.get("engine_version"):
                errors.append(f"engine {engine!r} is missing engine_version")
            commits = metadata.get("source_commits")
            if not isinstance(commits, list) or not commits or any(
                commit is not None and not re.fullmatch(r"[0-9a-f]{7,40}", str(commit))
                for commit in commits
            ):
                errors.append(f"engine {engine!r} has invalid source_commits")

    operations = document.get("operations")
    if not isinstance(operations, dict) or not operations:
        errors.append("operations must be a non-empty mapping")
    else:
        for operation, record in operations.items():
            if not isinstance(record, dict):
                errors.append(f"operation {operation!r} is not an object")
                continue
            status_counts = record.get("status_counts")
            if not isinstance(status_counts, dict) or not status_counts:
                errors.append(f"operation {operation!r} is missing status_counts")
                continue
            for engine, counts in status_counts.items():
                if not isinstance(counts, dict) or any(
                    not isinstance(count, int) or count < 0 for count in counts.values()
                ):
                    errors.append(f"operation {operation!r}/{engine!r} has invalid counts")
                unknown = set(counts) - GOOD_STATUSES - BAD_STATUSES - {"mismatch", "uncomparable"}
                if unknown:
                    errors.append(f"operation {operation!r}/{engine!r} has unknown statuses: {sorted(unknown)}")

    for index, claim in enumerate(document.get("claims", [])):
        if not isinstance(claim, dict):
            errors.append(f"claim {index} is not an object")
            continue
        status = claim.get("status")
        if status not in GOOD_STATUSES:
            errors.append(f"claim {index} uses non-positive status {status!r}")
        operation = claim.get("operation")
        if not isinstance(operations, dict) or operation not in operations:
            errors.append(f"claim {index} references unknown operation {operation!r}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scorecard", type=Path)
    parser.add_argument("--version")
    args = parser.parse_args()
    try:
        document = json.loads(args.scorecard.read_text(encoding="utf-8"))
        metadata = json.loads(DEFAULT_METADATA.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"scorecard read failure: {exc}", file=sys.stderr)
        return 1
    expected = args.version or metadata.get("release", {}).get("version")
    errors = validate(document, expected)
    if errors:
        print("scorecard validation failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"scorecard OK: {args.scorecard} ({expected})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
