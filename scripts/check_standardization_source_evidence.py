#!/usr/bin/env python3
"""Fail closed if the current-source standardization record drifts."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "validation/results/rdkit-standardization-current-source-2026-10-20-reionize.json"
PREVIOUS_RECORD = ROOT / "validation/results/rdkit-standardization-current-source-2026-10-20.json"


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

    assert not record["baseline_comparison"]["regressed_indices"]
    comparison = record["baseline_comparison"]
    assert digest(PREVIOUS_RECORD) == comparison["previous_report_sha256"]
    previous = json.loads(PREVIOUS_RECORD.read_text())
    assert previous["source_commit"] == comparison["previous_source_commit"]
    old_reionize = {row["index"] for row in previous["reionize_mismatches"]}
    new_reionize = {row["index"] for row in record["reionize_mismatches"]}
    assert old_reionize - new_reionize
    assert not new_reionize - old_reionize
    assert record["reionize_exact_rows"] > previous["reionize_exact_rows"]
    print("Current-source standardization evidence passed")


if __name__ == "__main__":
    main()
