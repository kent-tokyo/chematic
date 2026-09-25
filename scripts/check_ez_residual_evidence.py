#!/usr/bin/env python3
"""Validate the bounded Issue #149/#503 canonical E/Z evidence.

The historical files preserve the measured 4/28 residual and a rejected
7/28 experiment. The current files prove that the adopted complete-slot
planner resolves all 28 corpus components under the pinned 1,024-relabeling
gate. The current run was measured on v1.0.23 (after the Issue #632 canonical
E/Z change); the 2026-09-22 run at `2fe0cdd7` predates #632 and is kept as
history. This remains bounded evidence, not a claim about every possible E/Z
spelling.
"""

from __future__ import annotations

import json
import sys
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_audit_1024_2026-09-11.json"
EXPERIMENT = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_close_side_experiment_2026-09-11.json"
CURRENT_SUMMARY = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-25.json"
CURRENT_ROWS = ROOT / "validation" / "results" / "ez_shared_carrier_coupling_mechanism_audit_1024_2026-09-25.jsonl"


def fail(message: str) -> int:
    print(f"E/Z residual evidence invalid: {message}", file=sys.stderr)
    return 1


def read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("top-level value is not an object")
    return value


def read_jsonl(path: Path) -> list[dict]:
    rows = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        try:
            value = json.loads(line)
        except json.JSONDecodeError as error:
            raise ValueError(f"{path}:{line_number}: {error}") from error
        if not isinstance(value, dict):
            raise ValueError(f"{path}:{line_number}: row is not an object")
        rows.append(value)
    return rows


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
    if axis1.get("n_molecules_tested") != 28:
        return "current axis-1 molecule count changed"
    if axis1.get("n_divergent") != 0:
        return "current 1,024-relabeling gate is not 0/28"
    if axis1.get("n_cross_correspondence_failures_total") != 0:
        return "current cross-correspondence failures are present"
    verdict = summary.get("verdict")
    if not isinstance(verdict, dict) or verdict.get("verdict") != "PASS, no sampled residuals":
        return "current summary does not record the bounded pass verdict"
    if not isinstance(provenance, dict) or not isinstance(provenance.get("source_commit"), str):
        return "current source commit provenance is missing"
    if not isinstance(provenance.get("corpus_sha256"), str) or len(provenance["corpus_sha256"]) != 64:
        return "current corpus hash provenance is missing"
    if provenance.get("worktree_dirty") is not False:
        return "current run was not measured from a clean worktree"
    if provenance.get("tracked_diff_sha256") != "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855":
        return "current run does not record an empty tracked diff"
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-summary", type=Path, default=CURRENT_SUMMARY)
    args = parser.parse_args()
    try:
        audit = read(AUDIT)
        experiment = read(EXPERIMENT)
        current_summary = read(args.current_summary)
        current_rows = read_jsonl(CURRENT_ROWS)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        return fail(f"cannot read evidence: {error}")

    current_error = validate_current_summary(current_summary)
    if current_error:
        return fail(current_error)
    if len(current_rows) != 28:
        return fail("current detail file does not contain 28 components")
    if any(row.get("axis1", {}).get("divergent") is not False for row in current_rows):
        return fail("current detail file contains a divergent component")
    if any(
        row.get("axis1", {}).get("n_cross_correspondence_failures") != 0
        for row in current_rows
    ):
        return fail("current detail file contains correspondence failures")

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

    print(
        "E/Z evidence OK: current 1024-relabeling gate resolves 28/28; "
        "historical 4/28 and rejected 7/28 evidence remain preserved"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
