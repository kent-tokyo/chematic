#!/usr/bin/env python3
"""Validate the checked-in triclinic periodic-neighbor evidence envelope."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    raise SystemExit(f"triclinic neighbor evidence invalid: {message}")


def main() -> int:
    version = workspace_version(ROOT)
    path = ROOT / "validation" / "results" / f"periodic-neighbor-triclinic-cell-list-v{version}.json"
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read report: {error}")

    if report.get("schema_version") != 1:
        fail("schema_version must be 1")
    if report.get("target_version") != version:
        fail("target version is stale")
    if report.get("status") != "local-verified":
        fail("report is not local-verified")
    if "seeded_triclinic_sweep_matches_brute_force_oracle" not in str(report.get("fixture", "")):
        fail("seeded sweep fixture is missing")
    if report.get("command") != "cargo test -p chematic-crystal --offline --test neighbor --quiet":
        fail("reproduction command is missing or changed")

    candidate = report.get("candidate")
    if not isinstance(candidate, dict):
        fail("candidate section is missing")
    if candidate.get("fixtures") != 16 or candidate.get("seeded_fixtures") != 12:
        fail("fixture counts must be 16 total and 12 seeded")
    if candidate.get("seed") != "0x9e3779b97f4a7c15":
        fail("seed is missing or changed")
    cutoffs = candidate.get("cutoffs_angstrom")
    expected = [1.5, 1.75, 2.0, 2.75, 3.25, 4.0, 4.25, 5.5, 5.75, 6.0]
    if cutoffs != expected:
        fail("cutoff sweep is missing or changed")
    if "independent fractional-delta brute-force oracle" not in str(candidate.get("result", "")):
        fail("independent oracle claim is missing")
    print("Triclinic neighbor evidence OK: 16 fixtures, 12 seeded, 10 cutoffs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
