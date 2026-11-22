#!/usr/bin/env python3
"""Check v1.2.2 Cleanup measurements and frozen historical Reionize evidence."""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RECORD = ROOT / "validation/results/rdkit-cleanup-metal-ligand-stable-stress-2026-11-21.json"
PREVIOUS_RECORD = ROOT / "validation/results/rdkit-cleanup-integrated-prs-stress-2026-10-24.json"
TEN_K_RECORD = ROOT / "validation/results/rdkit-cleanup-metal-ligand-stable-10k-2026-11-21.json"
TEN_K_PREVIOUS = ROOT / "validation/results/rdkit-cleanup-integrated-prs-10k-2026-10-24.json"
REIONIZE_RECORD = ROOT / "validation/results/rdkit-reionize-in-memory-stereo-2026-10-22.json"
REIONIZE_PREVIOUS = ROOT / "validation/results/rdkit-cleanup-ring-stereo-stress-2026-10-22.json"
SUMMARY = ROOT / "benchmarks/2026-10-22-rdkit-reionize-in-memory-stereo.json"
PRIORITY_BASE = ROOT / "validation/results/rdkit-cleanup-v1.2.2-source-stress-2026-11-21.json"
PRIORITY_AUDIT = ROOT / "validation/results/cleanup-47-metal-ligand-stable-audit-2026-11-21.json"
PRIORITY_BASE_AUDIT = ROOT / "validation/results/cleanup-v1.2.2-47-atom-audit-2026-11-21.json"
PRIORITY_ROOTED = ROOT / "validation/results/cleanup-rooted-smiles-31-metal-ligand-stable-2026-11-21.json"
PRIORITY_WHEEL = ROOT / "validation/results/cleanup-metal-hydrogen-metal-ligand-stable-2026-11-21.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_priority_audit(base: dict, current: dict, audit: dict, rooted: dict, wheel: dict,
                         baseline_audit: dict) -> None:
    """Cross-check the published rows, candidate source and candidate wheel."""
    def require(condition: bool, message: str) -> None:
        if not condition:
            raise ValueError(message)

    require(base["corpus_sha256"] == current["corpus_sha256"], "corpora differ")
    require(all(record["rdkit_version"] == "2026.09.1" for record in (base, current, audit, rooted, wheel)),
            "RDKit versions differ")
    require(current["source_commit"] == audit["candidate_commit"] == rooted["candidate_commit"] == wheel["candidate_commit"],
            "candidate commits differ")
    require(audit["wheel_sha256"] == rooted["wheel_sha256"] == wheel["wheel_sha256"] and
            audit["wheel"] == rooted["wheel"] == wheel["wheel"], "wheel identities differ")
    require(re.fullmatch(r"[0-9a-f]{64}", wheel["wheel_sha256"]) is not None and
            all(record["installed_package_files_verified"] > 0 for record in (audit, rooted, wheel)),
            "wheel verification is incomplete")
    require(wheel["corpus_sha256"] == base["corpus_sha256"], "wheel corpus differs")
    require(audit["source_record"] == rooted["source_record"] == wheel["source_record"] == str(PRIORITY_BASE.relative_to(ROOT)),
            "baseline source records differ")
    original = {row["index"]: row for row in base["mismatches"]}
    candidate = {row["index"]: row for row in current["mismatches"]}
    reviewed = {row["index"]: row for row in audit["rows"]}
    wheel_rows = {row["index"]: row for row in wheel["rows"]}
    require(len(wheel_rows) == len(wheel["rows"]) == 50, "wheel audit rows differ")
    require(len(original) == len(base["mismatches"]) == 47 and set(original) == set(reviewed),
            "the prioritized 47 rows differ")
    require(Counter(row["classification"] for row in original.values()) == {
        "acid_base_reionization": 15, "metal_charge_assignment": 2, "metal_disconnection_policy": 30,
    }, "original mismatch groups differ")
    require(audit["counts"] == {"chemical_difference": 32, "representation_difference_after_sanitize": 15,
                                "refusal": 0, "requires_stereo_radical_review": 0},
            "priority classifications differ")
    require(Counter(row["assessment"] for row in reviewed.values()) == {
        "chemical_difference": 32, "representation_difference_after_sanitize": 15,
    }, "priority row classifications differ")
    changed = []
    for index, row in reviewed.items():
        require(index in candidate and candidate[index]["source_error"] is None,
                f"candidate row {index} missing or refused")
        require(original[index]["input"] == candidate[index]["input"], f"input changed at {index}")
        require(index in wheel_rows and wheel_rows[index].get("output") == candidate[index]["chematic_cleanup"],
                f"wheel/source output differs at {index}")
        wheel_row = wheel_rows[index]
        require(len({wheel_row[key] for key in ("input_formula", "output_formula", "reparsed_formula", "repeated_formula")}) == 1
                and len(wheel_row["formal_charges"]) == 4 and len(set(wheel_row["formal_charges"])) == 1,
                f"wheel formula or charge differs at {index}")
        same = original[index]["chematic_cleanup"] == candidate[index]["chematic_cleanup"]
        require(row["output_matches_published_record"] == same, f"wheel/source mismatch at {index}")
        require(row["formula_preserved"] and row["chematic_second_cleanup_matches"],
                f"formula or repeated Cleanup changed at {index}")
        require(row["rdkit_second_cleanup_matches"] ==
                (row["rdkit_six_cleanup_strings"][0] == row["rdkit_six_cleanup_strings"][1]),
                f"RDKit repeat status inconsistent at {index}")
        if not same:
            changed.append(index)
    require(sorted(changed) == [2025, 2203], "unexpected published-output changes")
    baseline_rows = {row["index"]: row for row in baseline_audit["rows"]}
    require(set(baseline_rows) == set(original) and baseline_audit["rdkit_version"] == "2026.09.1",
            "published-wheel audit rows differ")
    graph_indices = {index for index, row in baseline_rows.items() if
                     not row["charge_h_connectivity_isomorphic"] or row["bond_isotope_graph_isomorphic"] is False}
    current_graph_indices = {index for index, row in reviewed.items() if
                             not row["charge_h_connectivity_isomorphic"] or row["bond_isotope_graph_isomorphic"] is False}
    require(current_graph_indices == graph_indices | {2025}, "candidate graph-difference rows differ")
    cases = {row["index"]: row for row in rooted["cases"]}
    require(len(cases) == len(rooted["cases"]) == 31 and set(cases) == graph_indices,
            "graph-difference rows differ")
    root_dependent = repeat_dependent = 0
    for index, case in cases.items():
        variants = case["variants"]
        require(not case["excluded"] and len(variants) == case["inputs"] and
                len({row["input"] for row in variants}) == case["inputs"], f"excluded or duplicate rooted input at {index}")
        chematic_outputs = {row["chematic"] for row in variants}
        rdkit_outputs = {row["rdkit"] for row in variants}
        require(chematic_outputs == {candidate[index]["chematic_cleanup"]} and
                case["chematic_distinct_outputs"] == 1 and case["chematic_refusals"] == 0,
                f"candidate rooted output differs at {index}")
        require(None not in rdkit_outputs and len(rdkit_outputs) == case["rdkit_distinct_outputs"],
                f"RDKit rooted count differs at {index}")
        if len(rdkit_outputs) > 1:
            root_dependent += 1
        elif not reviewed[index]["rdkit_second_cleanup_matches"]:
            repeat_dependent += 1
        else:
            raise ValueError(f"RDKit is stable under both checks at {index}")
    require((root_dependent, repeat_dependent) == (24, 7), "RDKit instability counts differ")


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
    assert re.fullmatch(r"[0-9a-f]{40}", record["source_commit"])
    # The original candidate remains historical evidence. A separately measured
    # v1.2.3 candidate validates current source and identical mismatch rows.
    corpus = record["corpus"]
    suffix = "stress" if corpus.endswith("stress-v1.smi") else "10k"
    current_path = ROOT / "validation" / "results" / f"rdkit-cleanup-v1.2.3-candidate-{suffix}.json"
    current = json.loads(current_path.read_text())
    assert current["corpus"] == corpus
    assert current["corpus_sha256"] == record["corpus_sha256"]
    assert current["rdkit_version"] == record["rdkit_version"]
    assert current["mismatches"] == record["mismatches"]
    for relative, expected in current["source_sha256"].items():
        assert digest(ROOT / relative) == expected, relative
    assert digest(ROOT / corpus) == record["corpus_sha256"]

    # Cleanup retains its historical field names for compatibility.
    cleanup = dict(record)
    cleanup["cleanup_compatible_rows"] = record["compatible_rows"]
    cleanup["cleanup_mismatches"] = record["mismatches"]
    cleanup["cleanup_classification_counts"] = record["classification_counts"]
    cleanup["cleanup_support_status_counts"] = record["support_status_counts"]
    check_lane(record, "normalize_", "normalize_exact_rows")
    check_lane(cleanup, "cleanup_", "exact_rows")
    check_lane(record, "reionize_", "reionize_exact_rows")
    semantic_counts = Counter(row["semantic_classification"] for row in record["mismatches"])
    assert dict(sorted(semantic_counts.items())) == record["semantic_classification_counts"]
    assert record["sanitized_equivalent_mismatches"] == sum(
        row["equivalent_after_sanitization"] is True for row in record["mismatches"]
    )
    assert record["sanitization_unavailable_mismatches"] == sum(
        row["equivalent_after_sanitization"] is None for row in record["mismatches"]
    )
    for row in record["mismatches"]:
        if row["source_error"] is not None:
            assert row["semantic_classification"] == "source_error"
        elif row["equivalent_after_sanitization"] is True:
            assert row["semantic_classification"] == "representation_after_sanitization"
        elif row["equivalent_after_sanitization"] is False:
            assert row["semantic_classification"] == "different_after_sanitization"
        else:
            assert row["semantic_classification"] == "unclassified"


def check_comparison(
    record: dict, previous: dict, previous_path: Path, expected_regressed: list[int]
) -> None:
    comparison = record["baseline_comparison"]
    # These are measured exact-SMILES regressions, not accepted as chemical
    # equivalence. Pin every row so any additional change fails closed.
    assert comparison["regressed_indices"] == expected_regressed
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
    check_priority_audit(json.loads(PRIORITY_BASE.read_text()), record,
                         json.loads(PRIORITY_AUDIT.read_text()), json.loads(PRIORITY_ROOTED.read_text()),
                         json.loads(PRIORITY_WHEEL.read_text()), json.loads(PRIORITY_BASE_AUDIT.read_text()))
    # Rows 3075 and 3402 retain explicit metal H that the earlier source and
    # RDKit dropped. They are exact-string differences with conserved formula.
    check_comparison(
        record, previous, PREVIOUS_RECORD,
        [65, 83, 101, 403, 457, 809, 911, 1294, 1331, 1366, 1367, 1368, 1387, 1413, 1432, 1437, 1471, 1497, 1556, 1577,
         1630, 1664, 1686, 1756, 1781, 1831, 1854, 1875, 1878, 1890, 2003, 2025,
         2069, 2076, 2130, 2133, 2178, 2203, 2289, 2501, 2532, 2543, 2576, 2598, 2602, 2677, 2787,
         2855, 2909, 3075, 3153, 3154, 3160, 3207, 3369, 3376, 3402, 3501, 3504, 3590, 3700,
         3720, 3755, 3849, 3910, 3944, 3963, 4023, 4027],
    )
    assert record["baseline_comparison"]["fixed_indices"] == []
    assert record["baseline_comparison"]["reionize_fixed_indices"] == []
    assert record["baseline_comparison"]["newly_compatible_indices"] == [
        1334, 1373, 1400, 1550, 1904, 2253, 2449, 2671, 2859, 3056, 3235, 3497,
        4019, 4066
    ]
    assert record["reionize_exact_rows"] == record["normalize_exact_rows"] == 4069
    assert record["exact_rows"] == 4000
    assert record["semantic_classification_counts"] == {
        "different_after_sanitization": 45,
        "representation_after_sanitization": 17,
        "source_error": 7,
    }
    mismatches = {row["index"]: row for row in record["mismatches"]}
    assert mismatches[1294]["source_error"] is None
    for index in (457, 1556, 1756, 3755, 3944):
        assert mismatches[index]["source_error"].startswith("cleanup_unstable_metal:")
    for index in (2133, 2289):
        assert mismatches[index]["source_error"].startswith("cleanup_metal_stereo:")
    for index, metal in ((2025, "[FeH+2]"), (2203, "[PtH+2]"),
                         (3075, "[CuH+]"), (3402, "[ZnH+3]")):
        assert metal in mismatches[index]["chematic_cleanup"], index
    assert record["reionize_all_rows_exact"] == 4072
    assert previous["reionize_exact_rows"] == previous["normalize_exact_rows"] == 4055
    assert record["oracle_nonreparseable_indices"] == [1497, 1630, 2677]
    assert record["source_nonreparseable_indices"] == []

    ten_k = json.loads(TEN_K_RECORD.read_text())
    ten_k_previous = json.loads(TEN_K_PREVIOUS.read_text())
    check_current_record(ten_k)
    check_comparison(ten_k, ten_k_previous, TEN_K_PREVIOUS, [6991, 7850, 7998, 9504])
    assert ten_k["exact_rows"] == 9994
    assert next(row for row in ten_k["mismatches"] if row["index"] == 7998)["source_error"].startswith("cleanup_metal_hydrogen:")
    assert next(row for row in ten_k["mismatches"] if row["index"] == 9504)["source_error"].startswith("cleanup_unstable_metal:")
    assert ten_k_previous["exact_rows"] == 9998
    assert ten_k["baseline_comparison"]["fixed_indices"] == []

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
    print("Current-source standardization evidence checked: Cleanup exact-SMILES regressions 69/4 including 7/2 refusals; Reionize unchanged")


if __name__ == "__main__":
    main()
