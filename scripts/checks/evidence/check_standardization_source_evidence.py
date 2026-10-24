#!/usr/bin/env python3
"""Fail closed if current Cleanup records or historical Reionize evidence drift."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "validation/results/rdkit-cleanup-integrated-prs-stress-2026-10-24.json"
PREVIOUS_RECORD = ROOT / "validation/results/rdkit-cleanup-metal-aromaticity-stress-2026-10-22.json"
TEN_K_RECORD = ROOT / "validation/results/rdkit-cleanup-integrated-prs-10k-2026-10-24.json"
TEN_K_PREVIOUS = ROOT / "validation/results/rdkit-cleanup-metal-aromaticity-10k-2026-10-22.json"
REIONIZE_RECORD = ROOT / "validation/results/rdkit-reionize-in-memory-stereo-2026-10-22.json"
REIONIZE_PREVIOUS = ROOT / "validation/results/rdkit-cleanup-ring-stereo-stress-2026-10-22.json"
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


def check_current_record(record: dict) -> None:
    assert record["profile"] == "source_cleanup_vs_rdkit_2026_09_1"
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


def check_comparison(record: dict, previous: dict, previous_path: Path) -> None:
    comparison = record["baseline_comparison"]
    assert not comparison["regressed_indices"]
    assert digest(previous_path) == comparison["previous_report_sha256"]
    assert previous["source_commit"] == comparison["previous_source_commit"]
    all_indices = set(range(record["source_rows"]))
    previous_unscored = {
        row["index"] for field in ("excluded", "oracle_errors") for row in previous[field]
    }
    current_unscored = {
        row["index"] for field in ("excluded", "oracle_errors") for row in record[field]
    }
    previous_compatible = all_indices - previous_unscored
    current_compatible = all_indices - current_unscored
    old_cleanup = {row["index"] for row in previous["mismatches"]}
    new_cleanup = {row["index"] for row in record["mismatches"]}
    previous_exact = previous_compatible - old_cleanup
    current_exact = current_compatible - new_cleanup
    assert comparison["fixed_indices"] == sorted(old_cleanup & current_exact)
    assert comparison["regressed_indices"] == sorted(previous_exact & new_cleanup)
    assert comparison.get("newly_excluded_indices", []) == sorted(
        previous_compatible - current_compatible
    )
    assert comparison.get("newly_compatible_indices", []) == sorted(
        current_compatible - previous_compatible
    )
    assert comparison["previous_exact"] == previous["exact_rows"]
    assert record["exact_rows"] == len(current_exact)
    if previous["reionize_mismatches"] is None:
        assert previous["reionize_exact_rows"] == previous["reionize_compatible_rows"]
    old_reionize = {row["index"] for row in previous["reionize_mismatches"] or []}
    new_reionize = {row["index"] for row in record["reionize_mismatches"]}
    if "reionize_fixed_indices" in comparison:
        assert comparison["reionize_fixed_indices"] == sorted(old_reionize - new_reionize)
        assert comparison["reionize_regressed_indices"] == sorted(new_reionize - old_reionize)
        assert comparison["previous_reionize_exact"] == previous["reionize_exact_rows"]
    else:
        assert old_reionize == new_reionize == set()
        assert record["reionize_exact_rows"] == previous["reionize_exact_rows"]
    assert not new_reionize - old_reionize


def main() -> None:
    record = json.loads(RECORD.read_text())
    previous = json.loads(PREVIOUS_RECORD.read_text())
    check_current_record(record)
    check_comparison(record, previous, PREVIOUS_RECORD)
    assert record["baseline_comparison"]["fixed_indices"] == [
        403, 1471, 1497, 1630, 2478, 2602, 2677, 2787, 3154, 3238, 3501, 3737
    ]
    assert record["baseline_comparison"]["reionize_fixed_indices"] == []
    assert record["baseline_comparison"]["newly_excluded_indices"] == [
        1334, 1373, 1400, 1550, 1904, 2253, 2449, 2671, 2859, 3056, 3235, 3497,
        4019, 4066
    ]
    assert record["reionize_exact_rows"] == record["normalize_exact_rows"] == 4055
    assert previous["reionize_exact_rows"] == previous["normalize_exact_rows"] == 4069
    assert record["oracle_nonreparseable_indices"] == [1497, 1630, 2677]
    assert record["source_nonreparseable_indices"] == [1497, 1630, 2677]

    ten_k = json.loads(TEN_K_RECORD.read_text())
    ten_k_previous = json.loads(TEN_K_PREVIOUS.read_text())
    check_current_record(ten_k)
    check_comparison(ten_k, ten_k_previous, TEN_K_PREVIOUS)
    assert ten_k["exact_rows"] == 9998
    assert ten_k_previous["exact_rows"] == 9997
    assert ten_k["baseline_comparison"]["fixed_indices"] == [6991]

    # Keep the earlier Reionize evidence frozen while current source hashes
    # are checked against the new Cleanup report above.
    reionize = json.loads(REIONIZE_RECORD.read_text())
    older = json.loads(REIONIZE_PREVIOUS.read_text())
    historical = reionize["baseline_comparison"]
    assert digest(REIONIZE_PREVIOUS) == historical["previous_report_sha256"]
    assert older["source_commit"] == historical["previous_source_commit"]
    assert historical["fixed_indices"] == historical["regressed_indices"] == []
    assert historical["reionize_fixed_indices"] == [1432, 1878, 2914, 3153, 4023]
    assert historical["reionize_regressed_indices"] == []
    summary = json.loads(SUMMARY.read_text())
    assert summary["baseline_source_record_sha256"] == historical["previous_report_sha256"]
    assert summary["baseline_dump_sha256"] == older["source_dump_sha256"]
    assert summary["candidate_dump_sha256"] == reionize["source_dump_sha256"]
    assert summary["changed_indices"] == historical["reionize_fixed_indices"]
    assert summary["baseline_reionize_exact_rows"] == older["reionize_exact_rows"]
    assert summary["candidate_reionize_exact_rows"] == reionize["reionize_exact_rows"]
    assert summary["candidate_cleanup_exact_rows"] == reionize["exact_rows"]
    assert summary["candidate_normalize_exact_rows"] == reionize["normalize_exact_rows"]
    print("Current-source standardization evidence passed")


if __name__ == "__main__":
    main()
