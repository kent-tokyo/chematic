#!/usr/bin/env python3
"""Validate the current streaming failure-taxonomy evidence."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FORMATS = ("sdf", "mol", "xyz", "extxyz", "v3000", "mol2", "cml", "cdxml", "mmcif", "pdb")
EXPECTED_SUPPLEMENTAL = {
    "sdf": 0,
    "mol": 0,
    "xyz": 0,
    "extxyz": 6,
    "v3000": 0,
    "mol2": 2,
    "cml": 5,
    "cdxml": 3,
    "mmcif": 4,
    "pdb": 0,
}


def fail(message: str) -> int:
    print(f"streaming failure taxonomy invalid: {message}", file=sys.stderr)
    return 1


def main() -> int:
    version = workspace_version(ROOT)
    path = ROOT / "validation" / "results" / f"streaming-failure-taxonomy-v{version}.json"
    corpus = ROOT / "validation" / "streaming_format_safety_cases.json"
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
        corpus_document = json.loads(corpus.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return fail(f"cannot read report or corpus: {error}")

    if report.get("schema_version") != 1 or report.get("target_version") != version:
        return fail("schema or target version is stale")
    if report.get("gate") != "streaming_failure_taxonomy" or report.get("status") != "local-verified":
        return fail("report is not a verified taxonomy result")
    if report.get("errors") != []:
        return fail("report contains recorded runner errors")
    boundary = report.get("boundary", "")
    if "diagnostic variant evidence only" not in boundary:
        return fail("diagnostic/non-safety boundary is missing")

    if corpus_document.get("schema_version") != 1 or tuple(corpus_document.get("cases", {})) != FORMATS:
        return fail("base corpus format order changed")
    corpus_meta = report.get("corpus")
    if not isinstance(corpus_meta, dict) or corpus_meta.get("cases_per_format") != 12:
        return fail("base corpus metadata is stale")
    if corpus_meta.get("formats") != len(FORMATS) or corpus_meta.get("total_cases") != 140:
        return fail("total corpus metadata is stale")
    if corpus_meta.get("sha256") != hashlib.sha256(corpus.read_bytes()).hexdigest():
        return fail("base corpus digest does not match")

    rows = report.get("rows")
    if not isinstance(rows, dict) or set(rows) != set(FORMATS):
        return fail("taxonomy rows are missing or incomplete")
    for fmt in FORMATS:
        row = rows[fmt]
        if not isinstance(row, dict):
            return fail(f"{fmt}: row is not an object")
        if row.get("base_cases") != 12 or row.get("supplemental_cases") != EXPECTED_SUPPLEMENTAL[fmt]:
            return fail(f"{fmt}: case counts are stale")
        expected_cases = 12 + EXPECTED_SUPPLEMENTAL[fmt]
        kinds = row.get("failure_kinds")
        if not isinstance(kinds, dict) or not kinds or sum(kinds.values()) != expected_cases:
            return fail(f"{fmt}: failure-kind counts do not cover {expected_cases} cases")
        if any(not isinstance(name, str) or not name or not isinstance(count, int) or count <= 0 for name, count in kinds.items()):
            return fail(f"{fmt}: failure-kind entries must be non-empty positive integers")
        if row.get("cases") != expected_cases:
            return fail(f"{fmt}: total case count is stale")

    print("streaming failure taxonomy OK: 140 cases, 10 formats, current corpus digest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
