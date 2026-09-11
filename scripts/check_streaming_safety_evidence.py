#!/usr/bin/env python3
"""Validate the checked-in streaming safety summary for the current candidate."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_FORMATS = ("sdf", "mol", "xyz", "extxyz", "v3000", "mol2", "cml", "cdxml", "mmcif", "pdb")
EXPECTED_MALFORMED_PER_FORMAT = 80
EXPECTED_GENERATED_PER_FORMAT = 48


def fail(message: str) -> None:
    raise SystemExit(f"streaming safety evidence invalid: {message}")


def main() -> int:
    version = workspace_version(ROOT)
    candidates = sorted((ROOT / "benchmarks").glob(f"*-streaming-safety-v{version}.json"))
    if len(candidates) != 1:
        fail(
            f"expected exactly one current streaming safety report for v{version}, "
            f"found {len(candidates)}"
        )
    report_path = candidates[0]
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read report: {error}")
    if report.get("schema_version") != 1:
        fail("schema_version must be 1")
    if report.get("target_version") != version:
        fail("target version is stale")
    if report.get("status") != "local-current-candidate":
        fail("report is not for the current local candidate")
    environment = report.get("environment")
    if not isinstance(environment, dict) or tuple(environment.get("formats", ())) != EXPECTED_FORMATS:
        fail("format list is missing or out of order")
    results = report.get("results")
    if not isinstance(results, dict):
        fail("results are missing")
    expected_formats = len(EXPECTED_FORMATS)
    expected_malformed = expected_formats * EXPECTED_MALFORMED_PER_FORMAT
    expected_generated = expected_formats * EXPECTED_GENERATED_PER_FORMAT
    if results.get("malformed_negative_cases") != expected_malformed:
        fail(f"malformed case count must be {expected_malformed}")
    if results.get("malformed_unique_payloads") != expected_malformed:
        fail("malformed payloads must all be unique")
    if results.get("malformed_duplicate_reuses") != 0:
        fail("malformed duplicate reuses must be zero")
    if results.get("generated_parser_entry_cases") != expected_generated:
        fail(f"generated parser-entry count must be {expected_generated}")
    if results.get("generated_parser_entry_cases_per_format") != EXPECTED_GENERATED_PER_FORMAT:
        fail("generated parser-entry per-format count is stale")
    if results.get("generated_parser_entry_families_per_format") != 8:
        fail("generated parser-entry family count is stale")
    if results.get("generated_parser_entry_cases_per_family") != 6:
        fail("generated parser-entry cases-per-family count is stale")
    if results.get("generated_parser_entry_duplicate_reuses") != 0:
        fail("generated parser-entry duplicate reuses must be zero")
    per_format = results.get("malformed_cases_per_format")
    if per_format != {fmt: EXPECTED_MALFORMED_PER_FORMAT for fmt in EXPECTED_FORMATS}:
        fail("per-format malformed counts do not match the current runner wave")
    if results.get("oversized_input_cases") != expected_formats or results.get("gzip_cases") != expected_formats * 2:
        fail("oversized/gzip counts are stale")
    if results.get("failures") != 0:
        fail("recorded failures must be zero")
    print(
        f"streaming safety evidence OK: {expected_malformed} malformed, "
        f"{expected_generated} generated parser-entry, {expected_formats} formats"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
