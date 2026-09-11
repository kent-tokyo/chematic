#!/usr/bin/env python3
"""Validate the shared reaction-application evidence against its fixture."""

from __future__ import annotations

import json
import re
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "validation" / "cross_binding_contract.json"
REPORT = ROOT / "validation" / "results" / f"reaction-application-breadth-v{workspace_version(ROOT)}.json"


def current_workspace_version() -> str:
    cargo = (ROOT / "Cargo.toml").read_text(encoding="utf-8")
    match = re.search(r'(?ms)^\[workspace\.package\]\s*$.*?^version\s*=\s*"([^"]+)"\s*$', cargo)
    if not match:
        fail("Cargo.toml has no workspace package version")
    return match.group(1)


QUALITY_REPORT = ROOT / "validation" / "results" / f"reaction-quality-boundary-v{current_workspace_version()}.json"


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
    # The v1.0.12 fixture added two deliberately separate quality-boundary
    # cases.  The older breadth report remains a valid 15-case historical
    # measurement; do not make it stale merely because the shared fixture was
    # extended. Validate those two cases through their dedicated current
    # report below.
    quality_cases = [case for case in additional if case.get("quality_boundary")]
    breadth_cases = [case for case in additional if not case.get("quality_boundary")]
    positive_count = len(breadth_cases) + 1
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

    try:
        quality_report = json.loads(QUALITY_REPORT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read quality-boundary report: {error}")
    if quality_report.get("target_version") != current_workspace_version():
        fail("quality-boundary report target version is stale")
    if quality_report.get("fixture") != "validation/cross_binding_contract.json#reaction_application_contract":
        fail("quality-boundary fixture reference is incorrect")
    quality = quality_report.get("cases")
    if not isinstance(quality, dict) or set(quality) != {case.get("quality_boundary") for case in quality_cases}:
        fail("quality-boundary case set does not match fixture")
    if any(quality[name].get("status") != "pass" for name in quality):
        fail("quality-boundary cases must be passed")
    quality_bindings = quality_report.get("bindings", {})
    if quality_bindings.get("rust", {}).get("status") != "pass":
        fail("quality-boundary Rust result must be passed")
    if quality_bindings.get("node_wasm", {}).get("status") != "pass":
        fail("quality-boundary Node/WASM result must be passed")
    python_status = quality_bindings.get("python", {}).get("status")
    if python_status not in {"pass", "environment_blocked"}:
        fail("quality-boundary Python result must be pass or an explicit environment boundary")
    print(
        f"reaction application evidence OK: {positive_count} historical positive, "
        f"{len(quality_cases)} supplemental quality cases, {len(negative)} negative; "
        "Rust/Python/Node-WASM/WASM boundaries explicit"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
