#!/usr/bin/env python3
"""Validate the checked-in competitive benchmark protocol offline."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "validation" / "competitive_benchmark_manifest.json"


def fail(message: str) -> None:
    print(f"competitive benchmark manifest: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    try:
        data = json.loads(MANIFEST.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read JSON: {exc}")

    if data.get("schema_version") != 1:
        fail("schema_version must be 1")
    if data.get("status") != "prepared":
        fail("status must remain prepared until an actual result is recorded")
    if not re.fullmatch(r"\d+\.\d+\.\d+", data.get("target_version", "")):
        fail("target_version must be a semantic version")

    required_metadata = set(data.get("required_metadata", []))
    if len(required_metadata) < 10:
        fail("required_metadata is too small")
    if len(data.get("engines", [])) < 2:
        fail("at least chematic and one comparison engine are required")

    for corpus in data.get("corpora", []):
        path = ROOT / corpus.get("path", "")
        if not path.exists():
            fail(f"corpus path does not exist: {corpus.get('path')}")
        if not corpus.get("id") or not corpus.get("role"):
            fail("every corpus needs id and role")

    for operation in data.get("operations", []):
        if not operation.get("id") or not operation.get("metric"):
            fail("every operation needs id and metric")
        runner = ROOT / operation.get("runner", "")
        if not runner.exists():
            fail(f"runner does not exist: {operation.get('runner')}")

    rules = data.get("fairness_rules", [])
    if len(rules) < 6 or any(not isinstance(rule, str) or not rule for rule in rules):
        fail("fairness_rules must contain at least six non-empty rules")

    print(
        "Competitive benchmark manifest OK: "
        f"{len(data['engines'])} engines, {len(data['corpora'])} corpora, "
        f"{len(data['operations'])} operations; preparation only"
    )


if __name__ == "__main__":
    main()
