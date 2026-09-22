#!/usr/bin/env python3
"""Validate the committed RDKit 2026.03.6 rebaseline evidence packet."""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
PREFIX = "v1.0.19-vs-2026.03.6-2026-09-22"
ALLOWED_CLASSES = {
    "chematic_regression",
    "oracle_change",
    "contract_difference",
    "unresolved",
}


def load(name: str) -> dict:
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def first_lane(document: dict, name: str, errors: list[str]) -> dict:
    lanes = document.get("lanes")
    if not isinstance(lanes, list) or len(lanes) != 1 or not isinstance(lanes[0], dict):
        errors.append(f"{name} provenance must contain exactly one lane")
        return {}
    return lanes[0]


def validate_rows(path: Path, expected_sha256: str, errors: list[str]) -> int:
    digest = hashlib.sha256()
    count = 0
    with gzip.open(path, "rb") as handle:
        for expected_index, line in enumerate(handle):
            digest.update(line)
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                errors.append(f"rows line {expected_index + 1}: invalid JSON: {exc}")
                continue
            require(
                row.get("input_index") == expected_index,
                f"rows line {expected_index + 1}: non-contiguous input_index",
                errors,
            )
            for difference in row.get("differences", []):
                require(
                    difference.get("classification") in ALLOWED_CLASSES,
                    f"rows line {expected_index + 1}: invalid difference classification",
                    errors,
                )
            count += 1
    require(
        digest.hexdigest() == expected_sha256,
        "uncompressed rows SHA-256 changed",
        errors,
    )
    return count


def main() -> int:
    errors: list[str] = []
    contract = load(f"rdkit-rebaseline-python-binding-contract-{PREFIX}.json")
    summary = load(f"rdkit-rebaseline-python-chemistry-summary-{PREFIX}.json")
    python_provenance = load(f"rdkit-rebaseline-python-provenance-{PREFIX}.json")
    npm_provenance = load(f"rdkit-rebaseline-npm-provenance-{PREFIX}.json")

    require(contract.get("gate_passed") is True, "binding contract gate failed", errors)
    require(
        contract.get("runtime", {}).get("rdkit_version") == "2026.03.6",
        "wrong Python runtime version",
        errors,
    )
    require(
        contract.get("runtime", {}).get("backend", {}).get("value")
        == "boost_python",
        "wrong Python wrapper backend",
        errors,
    )
    require(
        contract.get("row_accounting", {})
        == {
            "input_count": 10000,
            "success_count": 10000,
            "failed_count": 0,
            "failures": [],
        },
        "binding row accounting changed",
        errors,
    )

    require(summary.get("gate_passed") is True, "chemistry evidence gate failed", errors)
    require(
        summary.get("sealed_accuracy_cohort_reused") is False,
        "sealed cohort reuse must remain false",
        errors,
    )
    require(
        summary.get("row_accounting")
        == {
            "input_count": 10000,
            "completed_count": 10000,
            "parse_failure_count": 0,
        },
        "chemistry row accounting changed",
        errors,
    )
    expected_counts = {
        "morgan_exact": 9999,
        "morgan_difference": 1,
        "cip_exact": 9880,
        "smiles_semantic_difference": 18,
        "smarts_cells": 310000,
        "smarts_differences": 14306,
    }
    for key, expected in expected_counts.items():
        require(
            summary.get("counts", {}).get(key) == expected,
            f"{key} count changed",
            errors,
        )

    rows_path = RESULTS / f"rdkit-rebaseline-python-chemistry-rows-{PREFIX}.jsonl.gz"
    row_count = validate_rows(rows_path, summary["rows"]["sha256"], errors)
    require(row_count == 10000, f"expected 10000 raw rows, found {row_count}", errors)

    python_lane = first_lane(python_provenance, "Python", errors)
    require(python_lane.get("availability") == "measured", "Python provenance lane unavailable", errors)
    require(python_lane.get("missing_dimensions") == [], "Python provenance has missing dimensions", errors)
    require(
        python_lane.get("package", {}).get("artifact_archive_sha256")
        == "e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41",
        "Python wheel SHA-256 changed",
        errors,
    )

    npm_lane = first_lane(npm_provenance, "npm", errors)
    require(npm_lane.get("availability") == "measured", "npm provenance lane unavailable", errors)
    require(npm_lane.get("missing_dimensions") == [], "npm provenance has missing dimensions", errors)
    require(
        npm_lane.get("package", {}).get("artifact_archive_sha256")
        == "3b86b72775394ae997fb96ca854c6e7c5840927d3e899c1ef117d42bca573c87",
        "npm tarball SHA-256 changed",
        errors,
    )

    for command in (
        "python-binding-contract",
        "python-chemistry",
        "python-provenance",
        "npm-provenance",
    ):
        execution = load(f"rdkit-rebaseline-{command}-execution-{PREFIX}.json")
        require(execution.get("returncode") == 0, f"{command} execution failed", errors)

    if errors:
        print("RDKit rebaseline evidence invalid:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    print(
        "RDKit rebaseline evidence OK: 10,000 complete rows; "
        "Morgan 9,999 exact; CIP 9,880 exact; SMARTS 14,306/310,000 differences"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
