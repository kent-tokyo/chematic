#!/usr/bin/env python3
"""Validate the bounded Issue #503 residual evidence.

This gate protects the distinction between a reproducible residual and a
completed canonical-invariance proof.  It intentionally validates evidence
shape and fail-closed conclusions; it does not promote an unstable canonical
winner.
"""

from __future__ import annotations

import json
import sys
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_audit_1024_2026-09-11.json"
EXPERIMENT = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_close_side_experiment_2026-09-11.json"
CURRENT_SUMMARY = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-12.json"


def fail(message: str) -> int:
    print(f"E/Z residual evidence invalid: {message}", file=sys.stderr)
    return 1


def read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("top-level value is not an object")
    return value


def validate_current_summary(summary: dict) -> str | None:
    if summary.get("relabelings_per_molecule") != 1024:
        return "current summary is not the required 1024-relabeling run"
    topology = summary.get("topology")
    axis1 = summary.get("axis1_summary")
    provenance = summary.get("source_provenance")
    if not isinstance(topology, dict) or not isinstance(axis1, dict):
        return "current summary is missing topology or axis-1 summary"
    if topology.get("n_coupled_components") != 28:
        return "current coupled-component count changed"
    if topology.get("coupled_component_sizes") != [2]:
        return "current component-size boundary changed"
    if topology.get("coupled_component_shapes") != ["path"]:
        return "current component-shape boundary changed"
    if axis1.get("n_divergent") != 4:
        return "current expected residual count is not 4/28"
    if axis1.get("n_cross_correspondence_failures_total") != 0:
        return "current cross-correspondence failures are present"
    verdict = summary.get("verdict")
    if not isinstance(verdict, dict) or verdict.get("verdict") != "NEEDS-RESEARCH, confirmed residuals found":
        return "current summary does not preserve the fail-closed research verdict"
    if not isinstance(provenance, dict) or not isinstance(provenance.get("source_commit"), str):
        return "current source commit provenance is missing"
    if not isinstance(provenance.get("corpus_sha256"), str) or len(provenance["corpus_sha256"]) != 64:
        return "current corpus hash provenance is missing"
    if provenance.get("worktree_dirty") is not True:
        return "current run must record the dirty worktree boundary"
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-summary", type=Path, default=CURRENT_SUMMARY)
    args = parser.parse_args()
    try:
        audit = read(AUDIT)
        experiment = read(EXPERIMENT)
        current_summary = read(args.current_summary)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        return fail(f"cannot read evidence: {error}")

    current_error = validate_current_summary(current_summary)
    if current_error:
        return fail(current_error)

    if audit.get("schema_version") != 1 or audit.get("issue") != 503:
        return fail("audit schema or issue number changed")
    if audit.get("relabelings_per_molecule") != 1024:
        return fail("audit is not the required 1024-relabeling run")
    if audit.get("coupled_components") != 28:
        return fail("coupled-component count changed")
    if audit.get("component_sizes") != [2] or audit.get("component_shapes") != ["path"]:
        return fail("component shape boundary changed")
    if audit.get("cross_correspondence_failures") != 0:
        return fail("cross-correspondence failures are present")
    if audit.get("divergent_components") != 4:
        return fail("expected stable residual count is not 4/28")
    inputs = audit.get("divergent_inputs")
    if not isinstance(inputs, list) or len(inputs) != 4 or not all(isinstance(item, str) and item for item in inputs):
        return fail("the four residual inputs are not recorded")
    interpretation = audit.get("interpretation", "")
    if "fail-closed" not in interpretation or "No canonical winner" not in interpretation:
        return fail("audit does not record the fail-closed/no-winner boundary")

    if experiment.get("schema_version") != 1 or experiment.get("issue") != 503:
        return fail("experiment schema or issue number changed")
    if experiment.get("status") != "rejected":
        return fail("rejected experiment was promoted")
    baseline = experiment.get("baseline")
    result = experiment.get("experiment_result")
    if not isinstance(baseline, dict) or not isinstance(result, dict):
        return fail("baseline or experiment result is missing")
    if baseline.get("coupled_components") != 28 or baseline.get("divergent_components") != 4:
        return fail("baseline no longer records 4/28")
    if result.get("coupled_components") != 28 or result.get("divergent_components") != 7:
        return fail("rejected close-side experiment no longer records worsened 7/28")
    if result.get("cross_correspondence_failures") != 0:
        return fail("rejected experiment introduced correspondence failures")
    conclusion = experiment.get("conclusion", "")
    if "must not be promoted" not in conclusion or "remains in production" not in conclusion:
        return fail("rejected conclusion does not preserve the production guard")

    print("E/Z residual evidence OK: current 1024-relabeling summary confirms 4/28 residual; historical 7/28 experiment remains rejected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
