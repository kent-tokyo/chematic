#!/usr/bin/env python3
"""Account for every source-wheel A6 row without treating diagnostics as quality parity."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

ARM = "chematic_pipeline_v2_mmff94_strict_stereo_safe"
REASONS = {
    "gradient_converged",
    "iteration_limit",
    "constraint_rejected_fallback",
}


def summarize(rows: list[dict]) -> dict:
    if len(rows) != 265 or {row.get("row_index") for row in rows} != set(range(265)):
        raise ValueError("expected exactly one result for each fixed A+B row 0..264")
    statuses = Counter(row.get("status") for row in rows)
    reasons: Counter[str] = Counter()
    failures = []
    sound = 0
    for row in rows:
        if row.get("arm") != ARM:
            raise ValueError(f"unexpected arm on row {row['row_index']}")
        if row["status"] != "success":
            failures.append(
                {
                    "row_index": row["row_index"],
                    "status": row["status"],
                    "cause": row.get("failure_cause"),
                }
            )
            continue
        ff = row.get("force_field")
        if not isinstance(ff, dict):
            raise ValueError(f"missing force-field report on row {row['row_index']}")
        reason = ff.get("mmff94_termination")
        if reason not in REASONS:
            raise ValueError(f"missing/unknown MMFF94 stop reason on row {row['row_index']}")
        if ff.get("converged") is not (reason == "gradient_converged"):
            raise ValueError(f"convergence and stop reason disagree on row {row['row_index']}")
        reasons[reason] += 1
        sound += row.get("final_validation", {}).get("sound") is True
    return {
        "rows": len(rows),
        "statuses": dict(sorted(statuses.items())),
        "success_terminations": dict(sorted(reasons.items())),
        "successful_rows_with_internal_sound_geometry": sound,
        "failures": failures,
        "success_retention_gate": statuses["success"] == 265 and sound == 265,
        "full_a6_quality_gate": "not_run",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--build-profile", required=True, choices=("dev", "release"))
    args = parser.parse_args()
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    configuration = metadata["configuration"]
    if (
        metadata["package"]["kind"] != "source_candidate"
        or metadata["result"]["termination"] != "completed"
        or metadata["result"]["selected_molecules"] != 265
        or metadata["result"]["emitted_rows"] != 265
        or configuration["tiers"] != "AB"
        or configuration["arms"] != [ARM]
        or configuration["random_seed"] != 20260801
        or configuration["max_attempts"] != 8
        or configuration["ring_torsion_policy"] != "diagnostic_only"
        or configuration["force_field_max_iterations"] != 200
    ):
        raise ValueError("not the fixed 265-row source-wheel MMFF94 lane")
    digest = hashlib.sha256(args.rows.read_bytes()).hexdigest()
    if digest != metadata["result"]["output_sha256"]:
        raise ValueError("runner output SHA-256 differs from the raw rows")
    rows = [
        json.loads(line)
        for line in args.rows.read_text(encoding="utf-8").splitlines()
        if line
    ]
    report = {
        "schema": "a6-source-mmff94-termination/v1",
        "build_profile": args.build_profile,
        "wheel_sha256": metadata["package"]["wheel"]["sha256"],
        "source_revision": metadata["package"]["source_revision"],
        "rows_sha256": digest,
        **summarize(rows),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: report[key]
                for key in ("statuses", "success_terminations", "success_retention_gate")
            },
            sort_keys=True,
        )
    )
    if not report["success_retention_gate"] and os.environ.get("GITHUB_ACTIONS") == "true":
        print("::warning::A6 source wheel lost successful fixed-265 rows; adoption remains blocked")
    return 0  # Diagnostic until release-profile rows are adjudicated; never claim A6 passed.


if __name__ == "__main__":
    raise SystemExit(main())
