#!/usr/bin/env python3
"""Validate the committed RDKit 2026.03.6 rebaseline evidence packet."""

from __future__ import annotations

import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
PREFIX = "v1.0.19-vs-2026.03.6-2026-09-22"
ISSUE632_SUMMARY = RESULTS / (
    "smiles-ez-semantic-issue632-v1.0.20-candidate-vs-rdkit-2026.03.6-2026-09-23.json"
)
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


def validate_rows(
    path: Path, expected_sha256: str, errors: list[str]
) -> tuple[int, Counter[str], Counter[tuple[str, str]], int]:
    digest = hashlib.sha256()
    classes: Counter[str] = Counter()
    operation_classes: Counter[tuple[str, str]] = Counter()
    correspondence_rows = 0
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
            correspondence = row.get("cip", {}).get("index_correspondence", {})
            if correspondence == {"atom_order": True, "bond_endpoints": True}:
                correspondence_rows += 1
            else:
                errors.append(
                    f"rows line {expected_index + 1}: CIP index correspondence is not proven"
                )
            for difference in row.get("differences", []):
                classification = difference.get("classification")
                operation = difference.get("operation")
                require(
                    classification in ALLOWED_CLASSES,
                    f"rows line {expected_index + 1}: invalid difference classification",
                    errors,
                )
                if classification in ALLOWED_CLASSES and isinstance(operation, str):
                    classes[classification] += 1
                    operation_classes[(operation, classification)] += 1
            count += 1
    require(
        digest.hexdigest() == expected_sha256,
        "uncompressed rows SHA-256 changed",
        errors,
    )
    return count, classes, operation_classes, correspondence_rows


def validate_issue632_candidate(baseline_rows_path: Path, errors: list[str]) -> None:
    evidence = json.loads(ISSUE632_SUMMARY.read_text(encoding="utf-8"))
    require(evidence.get("gate_passed") is True, "Issue #632 gate failed", errors)
    require(
        evidence.get("oracle") == {"name": "RDKit Python", "version": "2026.03.6"},
        "Issue #632 oracle changed",
        errors,
    )
    require(
        evidence.get("row_accounting")
        == {"input_count": 10000, "completed_count": 10000, "parse_failure_count": 0},
        "Issue #632 row accounting changed",
        errors,
    )

    baseline_affected: list[int] = []
    with gzip.open(baseline_rows_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if any(
                item.get("operation") == "smiles_parse_write"
                and item.get("classification") == "chematic_regression"
                for item in row.get("differences", [])
            ):
                baseline_affected.append(row["input_index"])
    require(
        evidence.get("affected_input_indices") == baseline_affected,
        "Issue #632 affected-row inventory changed",
        errors,
    )

    raw = evidence.get("raw_rows", {})
    rows_path = ROOT / str(raw.get("path", ""))
    try:
        compressed = rows_path.read_bytes()
    except OSError as exc:
        errors.append(f"Issue #632 rows unavailable: {exc}")
        return
    require(
        hashlib.sha256(compressed).hexdigest() == raw.get("compressed_sha256"),
        "Issue #632 compressed rows SHA-256 changed",
        errors,
    )

    uncompressed_digest = hashlib.sha256()
    row_count = 0
    semantic_differences = 0
    graph_differences = 0
    cip_exact = 0
    morgan_exact = 0
    smarts_differences = 0
    repaired_affected: list[int] = []
    with gzip.open(rows_path, "rb") as handle:
        for expected_index, line in enumerate(handle):
            uncompressed_digest.update(line)
            row = json.loads(line)
            require(
                row.get("input_index") == expected_index,
                f"Issue #632 row {expected_index + 1}: non-contiguous input_index",
                errors,
            )
            smiles = row.get("smiles_parse_write", {})
            if smiles.get("semantic_roundtrip") is not True:
                semantic_differences += 1
            if smiles.get("nonisomeric_roundtrip") is not True:
                graph_differences += 1
            if row.get("cip", {}).get("exact") is True:
                cip_exact += 1
            if row.get("morgan", {}).get("exact") is True:
                morgan_exact += 1
            smarts_differences += int(row.get("smarts", {}).get("difference_count", 0))
            if expected_index in baseline_affected and smiles.get("semantic_roundtrip") is True:
                repaired_affected.append(expected_index)
            row_count += 1

    require(
        uncompressed_digest.hexdigest() == raw.get("uncompressed_sha256"),
        "Issue #632 uncompressed rows SHA-256 changed",
        errors,
    )
    observed_after = {
        "smiles_semantic_difference": semantic_differences,
        "smiles_graph_difference": graph_differences,
        "cip_exact": cip_exact,
        "morgan_exact": morgan_exact,
        "smarts_differences": smarts_differences,
    }
    expected_after = evidence.get("after", {})
    for key, value in observed_after.items():
        require(
            expected_after.get(key) == value,
            f"Issue #632 {key} changed: expected {expected_after.get(key)!r}, found {value}",
            errors,
        )
    require(row_count == 10000, f"Issue #632 expected 10,000 rows, found {row_count}", errors)
    require(
        repaired_affected == baseline_affected,
        "Issue #632 did not repair every historical affected row",
        errors,
    )
    require(
        (cip_exact, morgan_exact, smarts_differences) == (9770, 9999, 14306),
        "Issue #632 candidate worsened CIP, Morgan, or SMARTS counts",
        errors,
    )


def main() -> int:
    errors: list[str] = []
    contract = load(f"rdkit-rebaseline-python-binding-contract-{PREFIX}.json")
    summary = load(f"rdkit-rebaseline-python-chemistry-summary-{PREFIX}.json")
    python_provenance = load(f"rdkit-rebaseline-python-provenance-{PREFIX}.json")
    npm_provenance = load(f"rdkit-rebaseline-npm-provenance-{PREFIX}.json")
    native_availability = load(
        "rdkit-rebaseline-native-cpp-availability-v1.0.19-vs-2026.03.6-2026-09-23.json"
    )

    require(contract.get("gate_passed") is True, "binding contract gate failed", errors)
    require(
        contract.get("runtime", {}).get("rdkit_version") == "2026.03.6",
        "wrong Python runtime version",
        errors,
    )
    require(
        contract.get("runtime", {}).get("backend", {}).get("value") == "boost_python",
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

    require(
        summary.get("gate_passed") is True, "chemistry evidence gate failed", errors
    )
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
        "cip_exact": 9770,
        "smiles_semantic_difference": 18,
        "smarts_cells": 310000,
        "smarts_differences": 14306,
        "difference_class_chematic_regression": 18,
        "difference_class_contract_difference": 9919,
        "difference_class_unresolved": 3594,
    }
    for key, expected in expected_counts.items():
        require(
            summary.get("counts", {}).get(key) == expected,
            f"{key} count changed",
            errors,
        )

    rows_path = RESULTS / f"rdkit-rebaseline-python-chemistry-rows-{PREFIX}.jsonl.gz"
    row_count, classes, operation_classes, correspondence_rows = validate_rows(
        rows_path, summary["rows"]["sha256"], errors
    )
    require(row_count == 10000, f"expected 10000 raw rows, found {row_count}", errors)
    require(
        correspondence_rows == 10000,
        f"expected 10,000 correspondence-proven CIP rows, found {correspondence_rows}",
        errors,
    )
    require(
        classes
        == {
            "chematic_regression": 18,
            "contract_difference": 9919,
            "unresolved": 3594,
        },
        f"raw-row classification counts changed: {dict(classes)}",
        errors,
    )
    require(
        operation_classes[("smiles_parse_write", "chematic_regression")] == 18,
        "expected 18 classified SMILES stereo regressions",
        errors,
    )
    require(
        operation_classes[("morgan", "contract_difference")] == 1,
        "expected one classified Morgan contract difference",
        errors,
    )
    require(
        operation_classes[("cip", "unresolved")] == 230,
        "expected 230 unresolved CIP rows",
        errors,
    )

    validate_issue632_candidate(rows_path, errors)
    require(
        operation_classes[("smarts", "unresolved")] == 3364,
        "expected 3,364 unresolved SMARTS rows",
        errors,
    )

    python_lane = first_lane(python_provenance, "Python", errors)
    require(
        python_lane.get("availability") == "measured",
        "Python provenance lane unavailable",
        errors,
    )
    require(
        python_lane.get("missing_dimensions") == [],
        "Python provenance has missing dimensions",
        errors,
    )
    require(
        python_lane.get("package", {}).get("artifact_archive_sha256")
        == "e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41",
        "Python wheel SHA-256 changed",
        errors,
    )

    npm_lane = first_lane(npm_provenance, "npm", errors)
    require(
        npm_lane.get("availability") == "measured",
        "npm provenance lane unavailable",
        errors,
    )
    require(
        npm_lane.get("missing_dimensions") == [],
        "npm provenance has missing dimensions",
        errors,
    )
    require(
        npm_lane.get("package", {}).get("artifact_archive_sha256")
        == "3b86b72775394ae997fb96ca854c6e7c5840927d3e899c1ef117d42bca573c87",
        "npm tarball SHA-256 changed",
        errors,
    )
    require(
        native_availability.get("release", {}).get("release_id") == 378318432,
        "native availability release identity changed",
        errors,
    )
    require(
        native_availability.get("query")
        == "GET /repos/rdkit/rdkit/releases/tags/Release_2026_03_6",
        "native availability query changed",
        errors,
    )
    require(
        native_availability.get("release", {}).get("tag") == "Release_2026_03_6",
        "native availability release tag changed",
        errors,
    )
    require(
        native_availability.get("release", {}).get("github_release_assets") == [],
        "native availability must preserve the observed empty GitHub asset list",
        errors,
    )
    require(
        native_availability.get("lane", {}).get("availability") == "unavailable",
        "native C++ lane must remain explicitly unavailable",
        errors,
    )
    boundaries = native_availability.get("boundaries", {})
    for key in (
        "does_not_assert_no_native_distribution_exists_elsewhere",
        "does_not_reuse_python_wheel_as_cpp_lane",
        "does_not_reuse_npm_wasm_as_cpp_lane",
        "does_not_rank_unmeasured_performance",
    ):
        require(
            boundaries.get(key) is True,
            f"native availability boundary {key} changed",
            errors,
        )
    require(
        boundaries.get("sealed_accuracy_cohort_reused") is False,
        "native availability must not reuse sealed accuracy data",
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
        "Issue #632 reduces 18 SMILES stereo regressions to zero; "
        "one typed Morgan contract difference; "
        "230 CIP and 3,364 SMARTS rows unresolved"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
