#!/usr/bin/env python3
"""Validate the checked-in same-process cross-engine contract bundle.

The individual runners intentionally keep their format-specific signatures and
boundaries.  This gate only verifies that the complete bundle is present,
current, reproducible in structure, and has not silently turned a timing
context into a ranking claim.  PDB's lenient-parser mismatch is an expected
documented boundary, so it is checked for explicit reporting rather than
treated as parity.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from benchmark_version import workspace_version


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_GATES = {
    "sdf": "same_process_sdf_semantic_contract",
    "v2000": "same_process_v2000_mol_semantic_contract",
    "v3000": "same_process_v3000_mol_semantic_contract",
    "mol2": "same_process_mol2_semantic_contract",
    "xyz": "same_process_xyz_semantic_contract",
    "extxyz": "same_process_extxyz_semantic_contract",
    "pdb": "same_process_pdb_semantic_contract",
    "cdxml": "same_process_cdxml_semantic_contract",
}


def main() -> int:
    errors: list[str] = []
    version = workspace_version(ROOT)
    for fmt, gate in EXPECTED_GATES.items():
        relative = f"benchmarks/2026-09-10-same-process-{fmt}-contract-v{version}.json"
        if fmt in {"v2000", "v3000"}:
            relative = f"benchmarks/2026-09-10-same-process-{fmt}-mol-contract-v{version}.json"
        path = ROOT / relative
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{fmt}: cannot read report: {exc}")
            continue
        if report.get("schema_version") != 1:
            errors.append(f"{fmt}: schema_version must be 1")
        if report.get("target_version") != version:
            errors.append(f"{fmt}: target_version is {report.get('target_version')!r}, expected {version}")
        if report.get("gate") != gate:
            errors.append(f"{fmt}: unexpected gate {report.get('gate')!r}")
        if report.get("status") != "local-verified":
            errors.append(f"{fmt}: report status is not local-verified")
        if report.get("repeats") != 20:
            errors.append(f"{fmt}: expected 20 repetitions")
        comparison = report.get("comparison")
        if not isinstance(comparison, dict):
            errors.append(f"{fmt}: comparison object is missing")
            continue
        if comparison.get("speed_is_non_ranking_context") is not True:
            errors.append(f"{fmt}: timing boundary is not explicitly non-ranking")
        expected_records = comparison.get("expected_records")
        if not isinstance(expected_records, int) or expected_records <= 0:
            errors.append(f"{fmt}: expected_records must be a positive integer")
        rows = report.get("rows")
        if not isinstance(rows, dict):
            errors.append(f"{fmt}: timing rows are missing")
        else:
            for lane in ("chematic", "rdkit"):
                row = rows.get(lane)
                if not isinstance(row, dict):
                    errors.append(f"{fmt}: {lane} timing row is missing")
                    continue
                seconds = row.get("seconds")
                if not isinstance(seconds, (int, float)) or not math.isfinite(seconds) or seconds <= 0:
                    errors.append(f"{fmt}: {lane} timing must be finite and positive")
                records = row.get("records")
                if records != expected_records:
                    errors.append(f"{fmt}: {lane} records must equal expected_records ({expected_records})")
                failures = row.get("failures", 0)
                if not isinstance(failures, int) or failures < 0:
                    errors.append(f"{fmt}: {lane} failures must be a non-negative integer")
                elif isinstance(records, int) and records + failures != expected_records:
                    errors.append(f"{fmt}: {lane} records plus failures must equal expected_records")
                boundary = str(row.get("boundary", ""))
                if "current process" not in boundary:
                    errors.append(f"{fmt}: {lane} timing boundary must name the current process")
        # Every lane must report malformed cases and make its acceptance
        # boundary inspectable.  PDB intentionally records mismatches because
        # schematic's default parser is lenient; that is not hidden here.
        malformed = comparison.get("malformed_cases")
        if not isinstance(malformed, dict) or not malformed:
            errors.append(f"{fmt}: malformed case report is missing")
        elif any(
            not isinstance(lanes, dict)
            or any(
                not isinstance(result, dict)
                or not isinstance(result.get("records"), int)
                or not isinstance(result.get("failures"), int)
                or result.get("records", 0) < 0
                or result.get("failures", 0) < 0
                for result in lanes.values()
            )
            for lanes in malformed.values()
        ):
            errors.append(f"{fmt}: malformed case records/failures must be integers")
        if fmt == "pdb":
            if "malformed_rejection_mismatches" not in comparison:
                errors.append("pdb: lenient-parser mismatch list is missing")
            if "lenient" not in str(comparison.get("malformed_contract", "")).lower():
                errors.append("pdb: lenient parser boundary is not documented")
        elif fmt == "cdxml":
            if "malformed_rejection_mismatches" not in comparison:
                errors.append("cdxml: malformed rejection mismatch list is missing")
        else:
            if comparison.get("malformed_rejection_mismatches", []) not in ([], None):
                errors.append(f"{fmt}: unexpected malformed rejection mismatch")
        tool_versions = report.get("tool_versions")
        if not isinstance(tool_versions, dict) or not str(tool_versions.get("rdkit", "")).strip():
            errors.append(f"{fmt}: RDKit version is missing")

    if errors:
        print("Same-process contract bundle failures:", file=sys.stderr)
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(
        "Same-process contract bundle OK: "
        f"{len(EXPECTED_GATES)} formats, 20 repetitions each, "
        "PDB leniency and CDXML parser boundaries explicitly bounded"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
