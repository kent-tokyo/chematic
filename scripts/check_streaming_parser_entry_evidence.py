#!/usr/bin/env python3
"""Validate the checked-in parser-entry coverage evidence.

The streaming runner is the authoritative executable gate. This validator
ensures that its checked-in summary still describes the current category
declarations and workspace version, rather than silently becoming stale after
the generated wave changes.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from benchmark_version import workspace_version
from check_streaming_format_limits import GENERATED_PARSER_ENTRY_CATEGORIES


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    version = workspace_version(ROOT)
    path = ROOT / "validation" / "results" / f"streaming-parser-entry-categories-v{version}.json"
    kinds_path = ROOT / "validation" / "results" / f"streaming-parser-entry-failure-kinds-v{version}.json"
    errors: list[str] = []
    try:
        evidence = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"streaming parser-entry evidence read failure: {exc}", file=sys.stderr)
        return 1
    try:
        kinds_evidence = json.loads(kinds_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"streaming parser-entry failure-kind evidence read failure: {exc}", file=sys.stderr)
        return 1

    expected_cases = sum(len(categories) * 6 for categories in GENERATED_PARSER_ENTRY_CATEGORIES.values())
    if evidence.get("schema_version") != 1:
        errors.append("schema_version must be 1")
    if evidence.get("target_version") != version:
        errors.append(f"target_version must be {version}")
    if evidence.get("result") != "pass":
        errors.append("result must be pass")
    if evidence.get("generated_cases") != expected_cases:
        errors.append(f"generated_cases must be {expected_cases}")
    if evidence.get("unique_generated_cases") != expected_cases:
        errors.append(f"unique_generated_cases must be {expected_cases}")
    family_count = len(next(iter(GENERATED_PARSER_ENTRY_CATEGORIES.values())))
    if evidence.get("formats") != len(GENERATED_PARSER_ENTRY_CATEGORIES):
        errors.append("formats count does not match runner declarations")
    if evidence.get("cases_per_format") != family_count * 6:
        errors.append("cases_per_format must equal families_per_format * cases_per_family")
    if evidence.get("families_per_format") != family_count:
        errors.append(f"families_per_format must be {family_count}")
    if evidence.get("cases_per_family") != 6:
        errors.append("cases_per_family must be 6")
    if evidence.get("duplicate_reuses") != 0:
        errors.append("duplicate_reuses must be 0")

    categories = evidence.get("categories")
    if categories != {fmt: list(values) for fmt, values in GENERATED_PARSER_ENTRY_CATEGORIES.items()}:
        errors.append("evidence categories do not match the runner declarations")

    if errors:
        print("Streaming parser-entry evidence failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    if kinds_evidence.get("target_version") != version:
        errors.append(f"failure-kind target_version must be {version}")
    if kinds_evidence.get("result") != "pass":
        errors.append("failure-kind result must be pass")
    if kinds_evidence.get("generated_cases") != expected_cases:
        errors.append(f"failure-kind generated_cases must be {expected_cases}")
    if kinds_evidence.get("failure_cases") != expected_cases:
        errors.append(f"failure_cases must be {expected_cases}")
    if kinds_evidence.get("empty_categories") != 0:
        errors.append("empty_categories must be 0")
    kind_categories = kinds_evidence.get("category_failure_kinds", {})
    expected_categories = {
        fmt: list(values) for fmt, values in GENERATED_PARSER_ENTRY_CATEGORIES.items()
    }
    if set(kind_categories) != set(expected_categories):
        errors.append("failure-kind evidence format set does not match runner declarations")
    for fmt, categories in expected_categories.items():
        if set(kind_categories.get(fmt, {})) != set(categories):
            errors.append(f"failure-kind categories for {fmt} do not match runner declarations")
        for category in categories:
            if not kind_categories.get(fmt, {}).get(category):
                errors.append(f"failure-kind category is empty: {fmt}/{category}")
    if errors:
        print("Streaming parser-entry evidence failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        "Streaming parser-entry evidence OK: "
        f"{expected_cases} generated cases, {len(GENERATED_PARSER_ENTRY_CATEGORIES)} formats, "
        f"{family_count} families per format"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
