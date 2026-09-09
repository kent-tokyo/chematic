#!/usr/bin/env python3
"""Validate the reproducible local #337 determinism evidence bundle."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "validation" / "results" / f"mmff94-issue337-determinism-v{workspace_version(ROOT)}.json"


def fail(message: str) -> None:
    raise SystemExit(f"#337 determinism evidence invalid: {message}")


def main() -> int:
    try:
        report = json.loads(REPORT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read report: {error}")
    expected = {
        "schema_version": 1,
        "target_version": workspace_version(ROOT),
        "issue": 337,
        "result": "pass",
        "fixtures": 6,
        "seeded_relabelings_per_fixture": 256,
        "seeded_relabeling_checks": 1536,
        "rdkit_representative_family_parity": "open",
    }
    for key, value in expected.items():
        if report.get(key) != value:
            fail(f"{key} must be {value!r}")
    if report.get("observed_macrocycle_ring_sizes") != [28, 30, 32, 32, 31, 29]:
        fail("macrocycle ring-size boundary is stale")
    if not isinstance(report.get("boundary"), str) or not report["boundary"].strip():
        fail("boundary disclosure is missing")
    print("#337 determinism evidence OK: 6 fixtures, 1,536 seeded checks; parity remains open")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
