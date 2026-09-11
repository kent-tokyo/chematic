#!/usr/bin/env python3
"""Validate the local Criterion calibration contract and workflow boundary."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "criterion-gate-calibration.json"
WORKFLOW = ROOT / ".github" / "workflows" / "bench-pr-gate.yml"


def fail(message: str) -> int:
    print(f"Criterion calibration contract invalid: {message}", file=sys.stderr)
    return 1


def main() -> int:
    try:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        workflow = WORKFLOW.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        return fail(f"cannot read manifest or workflow: {error}")

    if manifest.get("schema_version") != 1:
        return fail("schema_version must be 1")
    thresholds = manifest.get("thresholds")
    if thresholds != {
        "stage1_route_ratio": 1.04,
        "stage2_fail_ratio": 1.04,
        "stage1_blocks": 3,
        "stage2_blocks": 10,
    }:
        return fail("threshold or block contract changed")
    cases = manifest.get("cases")
    if not isinstance(cases, list) or {case.get("id") for case in cases} != {
        "plus5", "plus10", "plus6", "clean_noop", "stage2_build_noise",
        "stage2_real_effect", "stage2_contaminated",
    }:
        return fail("calibration case set is incomplete")
    by_id = {case["id"]: case for case in cases}
    for case_id in ("plus5", "plus10", "plus6"):
        case = by_id[case_id]
        if case.get("stage") != 1 or case.get("expected_route") != "route" or case.get("expected_gate") != "eligible_for_stage2":
            return fail(f"{case_id}: stage-1 routing contract changed")
    if by_id["clean_noop"].get("expected_route") != "no-route":
        return fail("clean_noop must remain no-route")
    for case_id, expected_gate in (("stage2_build_noise", "inconclusive"), ("stage2_real_effect", "fail"), ("stage2_contaminated", "inconclusive")):
        case = by_id[case_id]
        if case.get("stage") != 2 or case.get("sign_test") != "fail" or case.get("expected_gate") != expected_gate:
            return fail(f"{case_id}: stage-2 outcome contract changed")
    if by_id["stage2_contaminated"].get("environment_contaminated") is not True:
        return fail("contaminated case must be explicitly marked")
    if by_id["stage2_real_effect"].get("environment_contaminated") is not False:
        return fail("real-effect case must be uncontaminated")

    required_fragments = (
        "STAGE1_ROUTE_THRESHOLD=1.04",
        "STAGE2_FAIL_THRESHOLD=1.04",
        "if [ \"$environment_contaminated\" -eq 1 ]; then",
        "not blocking.",
        "any_fail=1",
    )
    for fragment in required_fragments:
        if fragment not in workflow:
            return fail(f"workflow boundary is missing: {fragment}")
    if re.search(r'if \[ "\$environment_contaminated" \] -eq 1; then[\s\S]{0,700}any_fail=1', workflow):
        return fail("workflow can set any_fail inside contaminated branch")

    print("Criterion calibration contract OK: 7 cases, threshold 1.04, contaminated runs remain non-blocking")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
