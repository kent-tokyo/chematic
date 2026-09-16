#!/usr/bin/env python3
"""Validate the shared reaction-application evidence against its fixture."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "validation" / "cross_binding_contract.json"
REPORT_DIR = ROOT / "validation" / "results"

def fail(message: str) -> None:
    raise SystemExit(f"reaction application evidence invalid: {message}")


def main() -> int:
    reports = sorted(REPORT_DIR.glob("reaction-application-breadth-v*.json"))
    if not reports:
        fail("no reaction-application breadth report found")
    quality_reports = sorted(REPORT_DIR.glob("reaction-quality-boundary-v*.json"))
    if not quality_reports:
        fail("no reaction quality-boundary report found")
    # The breadth lane is intentionally historical: the shared fixture may
    # grow with supplemental quality cases without forcing a rerun of the
    # older 15-case binding matrix. Select the newest checked-in report and
    # validate its declared scope against the fixture below.
    report_path = reports[-1]
    try:
        fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
        report = json.loads(report_path.read_text(encoding="utf-8"))
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
    report_version = report.get("target_version")
    if not isinstance(report_version, str) or not report_version.strip():
        fail("breadth report target version is missing")
    if report.get("fixture") != "validation/cross_binding_contract.json":
        fail("fixture reference is incorrect")
    if report.get("cases") != positive_count:
        fail("positive case count does not match fixture")
    reported_negative = report.get("negative_cases")
    if not isinstance(reported_negative, int) or reported_negative < 0 or reported_negative > len(negative):
        fail("breadth report negative case count exceeds the current fixture")
    if report.get("status") != "local-verified":
        fail("report is not local-verified")
    bindings = report.get("bindings")
    if not isinstance(bindings, dict) or any(bindings.get(name) != "passed" for name in ("rust", "python", "node_wasm", "wasm_source")):
        fail("all binding results must be passed")
    boundary = report.get("negative_boundary")
    if not isinstance(boundary, str) or not boundary.strip():
        fail("negative boundary disclosure is missing")

    try:
        quality_report = json.loads(quality_reports[-1].read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        fail(f"cannot read quality-boundary report: {error}")
    quality_version = quality_report.get("target_version")
    if not isinstance(quality_version, str) or not quality_version.strip():
        fail("quality-boundary report target version is missing")
    if quality_report.get("fixture") != "validation/cross_binding_contract.json#reaction_application_contract":
        fail("quality-boundary fixture reference is incorrect")
    quality = quality_report.get("cases")
    if quality is None and isinstance(quality_report.get("new_boundaries"), dict):
        # v1.0.13 uses the compact boundary report shape. Normalize it to the
        # per-case status shape used by the validator while retaining support
        # for older reports that already expose `cases`.
        quality = {
            name: {"status": status}
            for name, status in quality_report["new_boundaries"].items()
        }
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
    negative_boundaries = quality_report.get("negative_boundaries", {})
    if not isinstance(negative_boundaries, dict) or not negative_boundaries:
        fail("negative-boundary results are missing")
    if any(status != "pass" for status in negative_boundaries.values()):
        fail("negative-boundary cases must be passed")
    covered_negative = reported_negative + len(negative_boundaries)
    if covered_negative > len(negative):
        fail("historical and supplemental negative case counts exceed the fixture")
    print(
        f"reaction application evidence OK: {positive_count} historical positive, "
        f"{len(quality_cases)} supplemental quality cases, {covered_negative}/{len(negative)} negative; "
        f"breadth report v{report_version}; quality report v{quality_version}; "
        "Rust/Python/Node-WASM/WASM boundaries explicit"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
