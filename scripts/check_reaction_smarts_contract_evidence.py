#!/usr/bin/env python3
"""Validate the bounded reaction-SMARTS contract evidence and fixture."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "validation" / "cross_binding_contract.json"
REPORT = ROOT / "validation" / "results" / f"reaction-smarts-bounded-contract-v{workspace_version(ROOT)}.json"


def fail(message: str) -> None:
    raise SystemExit(f"reaction SMARTS contract evidence invalid: {message}")


def main() -> int:
    try:
        fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
        report = json.loads(REPORT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read fixture or report: {error}")
    cases = fixture.get("reaction_smarts_contract", {}).get("cases")
    if not isinstance(cases, list) or len(cases) != 38:
        fail("fixture must contain exactly 38 cases")
    if report.get("schema_version") != 1 or report.get("target_version") != workspace_version(ROOT):
        fail("schema or target version is stale")
    if report.get("fixture") != "validation/cross_binding_contract.json::reaction_smarts_contract":
        fail("report fixture reference is incorrect")
    if report.get("cases") != len(cases) or report.get("result") != "pass":
        fail("report case count or result is invalid")
    surfaces = report.get("surfaces")
    if not isinstance(surfaces, dict):
        fail("surface results are missing")
    if surfaces.get("rust", {}).get("contract_cases") != 38:
        fail("Rust contract case count is invalid")
    if surfaces.get("python", {}).get("contract_cases") != 38:
        fail("Python contract case count is invalid")
    if surfaces.get("wasm", {}).get("contract_cases") != 38:
        fail("WASM contract case count is invalid")
    if not any(">[Pd]>" in case["smarts"] for case in cases):
        fail("agent-section case is missing")
    if not any(case["matches"] and ":1]>>[O:1]" in case["smarts"] and ":1]>>[O:1]" in case["reaction"] for case in cases):
        fail("mapped positive case is missing")
    if not any(not case["matches"] and ":1]>>[O:1]" in case["smarts"] and ":2]" in case["reaction"] for case in cases):
        fail("mapped negative case is missing")
    boundary = report.get("boundary")
    if not isinstance(boundary, str) or not boundary.strip():
        fail("boundary disclosure is missing")
    print("reaction SMARTS contract evidence OK: 38 cases; Rust/Python/WASM map and agent boundaries present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
