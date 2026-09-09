#!/usr/bin/env python3
"""Validate the shared reaction-application evidence against its fixture."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "validation" / "cross_binding_contract.json"
REPORT = ROOT / "validation" / "results" / f"reaction-application-breadth-v{workspace_version(ROOT)}.json"


def fail(message: str) -> None:
    raise SystemExit(f"reaction application evidence invalid: {message}")


def main() -> int:
    try:
        fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
        report = json.loads(REPORT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read fixture or report: {error}")

    contract = fixture.get("reaction_application_contract", {})
    additional = contract.get("additional_cases")
    negative = contract.get("negative_cases")
    if not isinstance(additional, list) or not isinstance(negative, list):
        fail("additional_cases and negative_cases must be arrays")
    positive_count = len(additional) + 1
    if report.get("schema_version") != 1:
        fail("schema_version must be 1")
    if report.get("target_version") != workspace_version(ROOT):
        fail("target version is stale")
    if report.get("fixture") != "validation/cross_binding_contract.json":
        fail("fixture reference is incorrect")
    if report.get("cases") != positive_count:
        fail("positive case count does not match fixture")
    if report.get("negative_cases") != len(negative):
        fail("negative case count does not match fixture")
    if report.get("status") != "local-verified":
        fail("report is not local-verified")
    bindings = report.get("bindings")
    if not isinstance(bindings, dict) or any(bindings.get(name) != "passed" for name in ("rust", "python", "node_wasm", "wasm_source")):
        fail("all binding results must be passed")
    boundary = report.get("negative_boundary")
    if not isinstance(boundary, str) or not boundary.strip():
        fail("negative boundary disclosure is missing")
    print(
        f"reaction application evidence OK: {positive_count} positive, "
        f"{len(negative)} negative cases; Rust/Python/Node-WASM/WASM passed"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
