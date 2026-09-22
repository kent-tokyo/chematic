#!/usr/bin/env python3
"""Validate the checked-in A0 core-eight sealed acceptance summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUMMARY = ROOT / "validation/results/a0-core-eight-sealed-acceptance-20260922.json"
FIELDS = {
    "molecular_weight": 1e-6,
    "hba": 0.0,
    "hbd": 0.0,
    "tpsa": 1e-6,
    "logp": 1e-6,
    "molar_refractivity": 1e-6,
    "fsp3": 1e-6,
    "aromatic_ring_count": 0.0,
}


def validate(summary: object) -> list[str]:
    if not isinstance(summary, dict):
        return ["summary must be an object"]
    errors: list[str] = []
    rows = 8_000
    if summary.get("profile") != "rdkit_core_eight_frozen_candidate_evaluation_v1":
        errors.append("unexpected evaluation profile")
    if summary.get("status") != "accepted" or summary.get("decision", {}).get("accepted") is not True:
        errors.append("A0 decision is not accepted")
    candidate = summary.get("candidate", {})
    if not isinstance(candidate.get("commit"), str) or len(candidate["commit"]) != 40:
        errors.append("candidate commit is missing or not a full SHA")
    if not isinstance(candidate.get("tag"), str) or not candidate["tag"].startswith("trust-eval-candidate-"):
        errors.append("annotated candidate tag is missing")
    oracle = summary.get("oracle", {})
    if oracle != {"engine": "RDKit", "version": "2025.09.3", "profile": "rdkit_compat_v2"}:
        errors.append("oracle contract changed")
    cohort = summary.get("cohort", {})
    if cohort.get("split") != "sealed_holdout" or cohort.get("rows") != rows:
        errors.append("sealed cohort must contain exactly 8,000 rows")
    for key in ("sha256", "manifest_sha256", "source_sha256"):
        if not isinstance(cohort.get(key), str) or len(cohort[key]) != 64:
            errors.append(f"cohort {key} is missing")
    accounting = summary.get("accounting", {})
    if accounting != {
        "rows": rows,
        "parsed": rows,
        "candidate_parse_failures": 0,
        "oracle_parse_failures": 0,
        "all_rows_accounted": True,
    }:
        errors.append("row accounting is incomplete")
    if summary.get("strict_tolerances") != FIELDS:
        errors.append("strict tolerance contract changed")
    observed = summary.get("fields", {})
    if set(observed) != set(FIELDS):
        errors.append("declared field set changed")
    else:
        for field in FIELDS:
            result = observed[field]
            if (
                result.get("strict_matches") != rows
                or result.get("mismatches") != 0
                or result.get("unsupported") != 0
            ):
                errors.append(f"{field} is not strict-green on all rows")
    provenance = summary.get("provenance", {})
    for key in (
        "candidate_binary_sha256",
        "evaluator_script_sha256",
        "unused_data_attestation_sha256",
        "raw_result_sha256",
    ):
        if not isinstance(provenance.get(key), str) or len(provenance[key]) != 64:
            errors.append(f"provenance {key} is missing")
    if "raw_rows" in summary or "cases" in summary:
        errors.append("commit-safe summary must not contain raw sealed rows")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", nargs="?", type=Path, default=DEFAULT_SUMMARY)
    args = parser.parse_args()
    try:
        summary = json.loads(args.summary.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        print(f"A0 core-eight sealed gate: BLOCKED ({error})")
        return 1
    errors = validate(summary)
    if errors:
        print("A0 core-eight sealed gate: BLOCKED")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("A0 core-eight sealed gate: PASS (8 fields x 8,000 rows, mismatches=0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
