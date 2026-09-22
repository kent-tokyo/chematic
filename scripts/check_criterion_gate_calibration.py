#!/usr/bin/env python3
"""Validate the local Criterion calibration contract and workflow boundary."""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "criterion-gate-calibration.json"
HOSTED_EVIDENCE = ROOT / "validation" / "results" / "criterion-gate-hosted-calibration-2026-09-22.json"
WORKFLOW = ROOT / ".github" / "workflows" / "bench-pr-gate.yml"


def fail(message: str) -> int:
    print(f"Criterion calibration contract invalid: {message}", file=sys.stderr)
    return 1


def main() -> int:
    try:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        hosted = json.loads(HOSTED_EVIDENCE.read_text(encoding="utf-8"))
        workflow = WORKFLOW.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        return fail(f"cannot read manifest or workflow: {error}")

    if manifest.get("schema_version") != 2:
        return fail("schema_version must be 2")
    thresholds = manifest.get("thresholds")
    if thresholds != {
        "stage1_route_ratio": 1.04,
        "stage2_fail_ratio": 1.04,
        "null_contamination_ratio": 1.04,
        "stage1_blocks": 3,
        "stage2_blocks": 10,
        "null_blocks": 10,
    }:
        return fail("threshold or block contract changed")
    cases = manifest.get("cases")
    expected_ids = {
        "plus5", "plus10", "plus6", "clean_noop", "stage2_build_noise",
        "stage2_real_effect", "stage2_contaminated",
        "null_small_build_noise", "null_one_sided_contamination",
    }
    if not isinstance(cases, list) or len(cases) != len(expected_ids) or {case.get("id") for case in cases} != expected_ids:
        return fail("calibration case set is incomplete")
    by_id = {case["id"]: case for case in cases}
    if any(not isinstance(case.get("median_ratio"), (int, float)) or not math.isfinite(case["median_ratio"]) or case["median_ratio"] <= 0 for case in cases):
        return fail("every calibration case needs a finite positive median_ratio")
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
    if by_id["null_small_build_noise"].get("expected_environment_contaminated") is not False:
        return fail("small null build noise must remain clean")
    if by_id["null_one_sided_contamination"].get("expected_environment_contaminated") is not True:
        return fail("material one-sided null bias must mark contamination")

    if hosted.get("schema_version") != 1 or hosted.get("issue") != 70:
        return fail("hosted evidence schema/issue is invalid")
    if hosted.get("thresholds") != {
        "stage1_route_ratio": 1.04,
        "stage2_fail_ratio": 1.04,
        "null_contamination_ratio": 1.04,
        "stage2_blocks": 10,
        "null_blocks": 10,
    }:
        return fail("hosted evidence thresholds changed")

    plus10 = hosted.get("plus10", {})
    if not (
        plus10.get("accepted") is True
        and plus10.get("stage2_verdict") == "fail"
        and plus10.get("stage2_baseline_wins") == 10
        and plus10.get("stage2_candidate_wins") == 0
        and plus10.get("stage2_median_ratio", 0) >= 1.04
        and plus10.get("null_verdict") == "inconclusive"
        and plus10.get("workflow_conclusion") == "failure"
    ):
        return fail("hosted +10% calibration is incomplete")

    plus5 = hosted.get("plus5", {})
    plus5_runs = plus5.get("runs")
    if not isinstance(plus5_runs, list) or len(plus5_runs) < 10:
        return fail("hosted +5% calibration needs at least 10 runs")
    detected = 0
    for run in plus5_runs:
        ratios = (
            run.get("stage1_median_ratio"),
            run.get("stage2_median_ratio"),
            run.get("null_median_ratio"),
        )
        if any(not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0 for value in ratios):
            return fail("hosted +5% run has an invalid ratio")
        if run.get("null_verdict") != "inconclusive":
            return fail("hosted +5% run has a contaminated null control")
        if run.get("detected") is True:
            detected += 1
            if not (
                run.get("stage2_verdict") == "fail"
                and run.get("stage2_baseline_wins") == 10
                and run.get("stage2_candidate_wins") == 0
                and run["stage2_median_ratio"] >= 1.04
            ):
                return fail("hosted +5% detected run lacks the full gate signal")
    if not (
        detected == plus5.get("detected_runs")
        and len(plus5_runs) == plus5.get("total_runs")
        and detected > len(plus5_runs) / 2
        and plus5.get("accepted") is True
    ):
        return fail("hosted +5% calibration did not meet strict-majority acceptance")

    contaminated = hosted.get("one_sided_contamination", {})
    if not (
        contaminated.get("accepted") is True
        and contaminated.get("environment_contaminated") is True
        and contaminated.get("null_verdict") == "fail"
        and contaminated.get("null_median_ratio", 0) >= 1.04
        and contaminated.get("stage2_verdict") == "fail"
        and contaminated.get("stage2_median_ratio", 0) >= 1.04
        and contaminated.get("reported_outcome") == "environment-inconclusive"
        and contaminated.get("workflow_conclusion") == "success"
    ):
        return fail("hosted one-sided-contamination calibration is incomplete")
    if "must never be merged" not in hosted.get("temporary_injection_policy", ""):
        return fail("temporary calibration branch policy is missing")
    if hosted.get("overall_accepted") is not True:
        return fail("hosted calibration is not accepted")

    required_fragments = (
        "STAGE1_ROUTE_THRESHOLD=1.04",
        "STAGE2_FAIL_THRESHOLD=1.04",
        "NULL_CONTAMINATION_THRESHOLD=1.04",
        "NULL_CONTROL_BLOCKS=10",
        "contamination-check",
        "if [ \"$environment_contaminated\" -eq 1 ]; then",
        "not blocking.",
        "any_fail=1",
    )
    for fragment in required_fragments:
        if fragment not in workflow:
            return fail(f"workflow boundary is missing: {fragment}")
    if re.search(r'if \[ "\$environment_contaminated" \] -eq 1; then[\s\S]{0,700}any_fail=1', workflow):
        return fail("workflow can set any_fail inside contaminated branch")

    print(
        "Criterion calibration contract OK: 9 local cases; hosted +10 blocked, "
        f"+5 detected {detected}/{len(plus5_runs)}, contamination downgraded"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
