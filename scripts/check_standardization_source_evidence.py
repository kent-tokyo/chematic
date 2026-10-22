#!/usr/bin/env python3
"""Fail closed if the current-source standardization record drifts."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "validation/results/rdkit-reionize-in-memory-stereo-2026-10-22.json"
PREVIOUS_RECORD = ROOT / "validation/results/rdkit-cleanup-ring-stereo-stress-2026-10-22.json"
SUMMARY = ROOT / "benchmarks/2026-10-22-rdkit-reionize-in-memory-stereo.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_lane(record: dict, prefix: str, exact_key: str) -> None:
    compatible = record[f"{prefix}compatible_rows"]
    exact = record[exact_key]
    mismatches = record[f"{prefix}mismatches"]
    classes = record[f"{prefix}classification_counts"]
    statuses = record[f"{prefix}support_status_counts"]
    assert exact + len(mismatches) == compatible
    assert sum(classes.values()) == len(mismatches)
    assert sum(statuses.values()) == len(mismatches)
    assert all(row.get("classification") and row.get("support_status") for row in mismatches)
    assert classes.get("unclassified_semantic_difference", 0) == 0
    assert statuses.get("needs_adjudication", 0) == 0


def main() -> None:
    record = json.loads(RECORD.read_text())
    assert record["profile"] == "source_standardization_vs_rdkit_2026_09_1"
    assert record["rdkit_version"] == "2026.09.1"
    # The measured commit may be absent from a PR checkout after a squash merge.
    # Current-source hashes below still fail closed on any source drift.
    assert re.fullmatch(r"[0-9a-f]{40}", record["source_commit"])
    for relative, expected in record["source_sha256"].items():
        assert digest(ROOT / relative) == expected, relative
    assert digest(ROOT / record["corpus"]) == record["corpus_sha256"]

    # Cleanup retains its historical field names for compatibility.
    cleanup = dict(record)
    cleanup["cleanup_compatible_rows"] = record["compatible_rows"]
    cleanup["cleanup_mismatches"] = record["mismatches"]
    cleanup["cleanup_classification_counts"] = record["classification_counts"]
    cleanup["cleanup_support_status_counts"] = record["support_status_counts"]
    check_lane(record, "normalize_", "normalize_exact_rows")
    check_lane(cleanup, "cleanup_", "exact_rows")
    check_lane(record, "reionize_", "reionize_exact_rows")

    comparison = record["baseline_comparison"]
    assert not comparison["regressed_indices"]
    assert digest(PREVIOUS_RECORD) == comparison["previous_report_sha256"]
    previous = json.loads(PREVIOUS_RECORD.read_text())
    assert previous["source_commit"] == comparison["previous_source_commit"]
    old_cleanup = {row["index"] for row in previous["mismatches"]}
    new_cleanup = {row["index"] for row in record["mismatches"]}
    assert comparison["fixed_indices"] == sorted(old_cleanup - new_cleanup)
    assert comparison["regressed_indices"] == sorted(new_cleanup - old_cleanup)
    assert comparison["previous_exact"] == previous["exact_rows"]
    assert not comparison["fixed_indices"]
    assert record["exact_rows"] == previous["exact_rows"]
    old_reionize = {row["index"] for row in previous["reionize_mismatches"]}
    new_reionize = {row["index"] for row in record["reionize_mismatches"]}
    assert comparison["reionize_fixed_indices"] == sorted(old_reionize - new_reionize)
    assert comparison["reionize_regressed_indices"] == sorted(new_reionize - old_reionize)
    assert comparison["previous_reionize_exact"] == previous["reionize_exact_rows"]
    assert not new_reionize - old_reionize
    assert comparison["reionize_fixed_indices"] == [1432, 1878, 2914, 3153, 4023]
    assert record["reionize_exact_rows"] == previous["reionize_exact_rows"] + len(old_reionize)
    assert record["reionize_exact_rows"] == record["reionize_compatible_rows"]
    assert record["normalize_exact_rows"] == previous["normalize_exact_rows"]
    summary = json.loads(SUMMARY.read_text())
    assert summary["baseline_source_record_sha256"] == comparison["previous_report_sha256"]
    assert summary["baseline_dump_sha256"] == previous["source_dump_sha256"]
    assert summary["candidate_dump_sha256"] == record["source_dump_sha256"]
    assert summary["changed_indices"] == comparison["reionize_fixed_indices"]
    assert summary["baseline_reionize_exact_rows"] == previous["reionize_exact_rows"]
    assert summary["candidate_reionize_exact_rows"] == record["reionize_exact_rows"]
    assert summary["candidate_cleanup_exact_rows"] == record["exact_rows"]
    assert summary["candidate_normalize_exact_rows"] == record["normalize_exact_rows"]
    print("Current-source standardization evidence passed")


if __name__ == "__main__":
    main()
