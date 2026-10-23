#!/usr/bin/env python3
"""Fail closed on a comparable pair of SMARTS/RDKit evidence summaries.

This gate deliberately validates measurement identity before discussing an
improvement. It does not treat RDKit parse errors or residual cells as success:
they remain separately counted in the emitted packet.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path


ROLES = {"baseline", "candidate"}
REQUIRED_PROVENANCE = {
    "source_commit",
    "source_tree_sha256",
    "source_diff_sha256",
    "binary_sha256",
    "corpus_sha256",
    "query_sha256",
}
REQUIRED_SUMMARY = {
    "schema_version",
    "role",
    "provenance",
    "rdkit_version",
    "rdkit_pinned_source",
    "comparison_config",
    "n_rows_in_dump",
    "n_molecules_compared",
    "n_alignment_checked",
    "n_alignment_failures",
    "total_cells",
    "bucket_counts",
}
DISALLOWED_SUCCESS_BUCKETS = {
    "parity_regresses",
    "parity_worsens_agreement",
    "chematic_parse_error",
    "chematic_parity_error",
}


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON must be an object")
    return value


def valid_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def validate_arm(arm: dict, role: str) -> list[str]:
    errors: list[str] = []
    missing = sorted(REQUIRED_SUMMARY - arm.keys())
    if missing:
        return [f"{role}: missing fields: {', '.join(missing)}"]
    if arm["schema_version"] != 1:
        errors.append(f"{role}: schema_version must be 1")
    if arm["role"] != role:
        errors.append(f"{role}: role must be {role!r}")
    provenance = arm["provenance"]
    if not isinstance(provenance, dict) or set(provenance) != REQUIRED_PROVENANCE:
        errors.append(f"{role}: provenance must contain exactly the required keys")
    elif not all(valid_sha256(provenance[key]) for key in REQUIRED_PROVENANCE - {"source_commit"}):
        errors.append(f"{role}: provenance digests must be lowercase SHA-256")
    elif not isinstance(provenance["source_commit"], str) or len(provenance["source_commit"]) != 40:
        errors.append(f"{role}: provenance.source_commit must be a 40-character SHA")
    if not isinstance(arm["bucket_counts"], dict) or sum(arm["bucket_counts"].values()) != arm["total_cells"]:
        errors.append(f"{role}: bucket counts must sum to total_cells")
    if arm["n_rows_in_dump"] != arm["n_molecules_compared"]:
        errors.append(f"{role}: dump row accounting is incomplete")
    if arm["n_alignment_checked"] != arm["n_molecules_compared"] or arm["n_alignment_failures"] != 0:
        errors.append(f"{role}: atom alignment is incomplete")
    return errors


def validate_pair(baseline: dict, candidate: dict) -> list[str]:
    errors = validate_arm(baseline, "baseline") + validate_arm(candidate, "candidate")
    bp = baseline.get("provenance", {})
    cp = candidate.get("provenance", {})
    for key in ("corpus_sha256", "query_sha256"):
        if bp.get(key) != cp.get(key):
            errors.append(f"baseline/candidate {key} differs")
    for key in ("rdkit_version", "rdkit_pinned_source", "comparison_config", "n_rows_in_dump", "total_cells"):
        if baseline.get(key) != candidate.get(key):
            errors.append(f"baseline/candidate {key} differs")
    for key in ("source_commit", "source_tree_sha256", "binary_sha256"):
        if bp.get(key) == cp.get(key):
            errors.append(f"baseline and candidate reuse the same {key}")
    return errors


def atomic_write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary_path = Path(handle.name)
    os.replace(temporary_path, path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--require-non-regression",
        action="store_true",
        help="fail unless the candidate is comparable and introduces no residual or error-bucket regression",
    )
    args = parser.parse_args()
    try:
        baseline = load(args.baseline)
        candidate = load(args.candidate)
        errors = validate_pair(baseline, candidate)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"SMARTS baseline/candidate packet error: {exc}", file=sys.stderr)
        return 2

    def bucket(arm: dict, name: str) -> int:
        return int(arm["bucket_counts"].get(name, 0))

    adoption_errors: list[str] = []
    if not errors:
        baseline_residuals = bucket(baseline, "both_disagree_same_as_default")
        candidate_residuals = bucket(candidate, "both_disagree_same_as_default")
        if candidate_residuals > baseline_residuals:
            adoption_errors.append(
                f"candidate residual cells regressed: {candidate_residuals} > {baseline_residuals}"
            )
        for name in sorted(DISALLOWED_SUCCESS_BUCKETS):
            if bucket(candidate, name) != 0:
                adoption_errors.append(f"candidate has {bucket(candidate, name)} {name} cells")

    result = {
        "schema_version": 1,
        "profile": "rdkit_smarts_baseline_candidate_v1",
        "baseline": str(args.baseline),
        "candidate": str(args.candidate),
        "comparison": {
            "residual_cells": {"baseline": bucket(baseline, "both_disagree_same_as_default"), "candidate": bucket(candidate, "both_disagree_same_as_default")},
            "rdkit_parse_error_cells": {"baseline": bucket(baseline, "rdkit_smarts_parse_error"), "candidate": bucket(candidate, "rdkit_smarts_parse_error")},
            "candidate_disallowed_buckets": {name: bucket(candidate, name) for name in sorted(DISALLOWED_SUCCESS_BUCKETS)},
        },
        "errors": errors,
        "comparison_valid": not errors,
        "adoption_allowed": not errors and not adoption_errors,
        "adoption_errors": adoption_errors,
        "gate_passed": not errors,
    }
    atomic_write(args.output, result)
    print("SMARTS baseline/candidate comparison: " + ("PASS" if not errors else "BLOCKED"))
    if errors:
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("SMARTS candidate adoption: " + ("ALLOWED" if not adoption_errors else "BLOCKED"))
    if adoption_errors:
        print("\n".join(f"- {error}" for error in adoption_errors))
        if args.require_non_regression:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
