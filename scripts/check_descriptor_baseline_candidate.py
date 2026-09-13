#!/usr/bin/env python3
"""Validate a descriptor baseline/candidate comparison packet.

The validator does not decide which implementation is chemically better.  It
only makes an apples-to-apples comparison auditable and fail closed when the
two arms reuse an artifact, corpus, or incompatible contract.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REQUIRED_FIELDS = {
    "molecular_weight",
    "hba",
    "hbd",
    "tpsa",
    "logp",
    "molar_refractivity",
    "fsp3",
    "aromatic_ring_count",
}
ROLES = {"baseline", "candidate"}
REQUIRED_PROVENANCE = {
    "source_commit",
    "source_tree_sha256",
    "artifact_sha256",
    "rdkit_version",
    "corpus_sha256",
    "corpus_rows",
}
CONTRACT = "rdkit_descriptor_semantics_v1"


def load(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path}: top-level JSON value must be an object")
    return value


def validate_arm(arm: dict, expected_role: str) -> list[str]:
    errors: list[str] = []
    if arm.get("schema_version") != 1:
        errors.append(f"{expected_role}: schema_version must be 1")
    if arm.get("role") != expected_role:
        errors.append(f"{expected_role}: role must be {expected_role!r}")
    if arm.get("contract") != CONTRACT:
        errors.append(f"{expected_role}: contract must be {CONTRACT!r}")
    provenance_value = arm.get("provenance")
    if not isinstance(provenance_value, dict) or set(provenance_value) != REQUIRED_PROVENANCE:
        errors.append(f"{expected_role}: provenance must contain exactly the required keys")
    provenance = provenance_value
    if not isinstance(provenance, dict):
        return errors
    for key in REQUIRED_PROVENANCE - {"corpus_rows"}:
        if not isinstance(provenance.get(key), str) or not provenance[key]:
            errors.append(f"{expected_role}: provenance.{key} is missing")
    if (
        not isinstance(provenance.get("corpus_rows"), int)
        or isinstance(provenance.get("corpus_rows"), bool)
        or provenance["corpus_rows"] <= 0
    ):
        errors.append(f"{expected_role}: provenance.corpus_rows must be positive")
    fields = arm.get("fields")
    if not isinstance(fields, dict) or set(fields) != REQUIRED_FIELDS:
        errors.append(f"{expected_role}: fields must exactly match the eight-field contract")
    if not isinstance(arm.get("profile"), str) or not arm["profile"]:
        errors.append(f"{expected_role}: profile is missing")
    return errors


def validate_pair(baseline: dict, candidate: dict) -> list[str]:
    errors = validate_arm(baseline, "baseline") + validate_arm(candidate, "candidate")
    bp = baseline.get("provenance", {})
    cp = candidate.get("provenance", {})
    for key in ("rdkit_version", "corpus_sha256", "corpus_rows"):
        if bp.get(key) != cp.get(key):
            errors.append(f"baseline/candidate {key} differs")
    if bp.get("artifact_sha256") == cp.get("artifact_sha256"):
        errors.append("baseline and candidate reuse the same artifact")
    if bp.get("source_commit") == cp.get("source_commit"):
        errors.append("baseline and candidate reuse the same source commit")
    if bp.get("source_tree_sha256") == cp.get("source_tree_sha256"):
        errors.append("baseline and candidate reuse the same source tree")
    # Field-level coverage is a result to report, not a reason to reject the
    # comparison: a candidate may legitimately support rows that the baseline
    # could not parse.  The shared corpus row count above remains mandatory.
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    try:
        baseline = load(args.baseline)
        candidate = load(args.candidate)
        errors = validate_pair(baseline, candidate)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"baseline/candidate packet error: {exc}", file=sys.stderr)
        return 2
    result = {
        "schema_version": 1,
        "profile": "rdkit_descriptor_baseline_candidate",
        "baseline": str(args.baseline),
        "candidate": str(args.candidate),
        "errors": errors,
        "gate_passed": not errors,
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("baseline/candidate gate: " + ("PASS" if not errors else "BLOCKED"))
    if errors:
        print("\n".join(f"- {error}" for error in errors))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
