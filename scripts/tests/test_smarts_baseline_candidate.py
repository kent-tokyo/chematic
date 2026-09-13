from scripts.check_smarts_baseline_candidate import validate_pair


def arm(role: str, *, commit: str, tree: str, binary: str, corpus: str = "a" * 64, query: str = "b" * 64):
    return {
        "schema_version": 1,
        "role": role,
        "provenance": {
            "source_commit": commit,
            "source_tree_sha256": tree,
            "source_diff_sha256": "c" * 64,
            "binary_sha256": binary,
            "corpus_sha256": corpus,
            "query_sha256": query,
        },
        "rdkit_version": "2025.09.3",
        "rdkit_pinned_source": "fixture",
        "comparison_config": {"requires_completed_footer": True},
        "n_rows_in_dump": 2,
        "n_molecules_compared": 2,
        "n_alignment_checked": 2,
        "n_alignment_failures": 0,
        "total_cells": 4,
        "bucket_counts": {"agree_all": 4},
    }


def pair():
    return (
        arm("baseline", commit="1" * 40, tree="d" * 64, binary="e" * 64),
        arm("candidate", commit="2" * 40, tree="f" * 64, binary="0" * 64),
    )


def test_accepts_independently_built_comparable_arms():
    baseline, candidate = pair()
    assert validate_pair(baseline, candidate) == []


def test_rejects_same_binary_or_different_query_set():
    baseline, candidate = pair()
    candidate["provenance"]["binary_sha256"] = baseline["provenance"]["binary_sha256"]
    candidate["provenance"]["query_sha256"] = "9" * 64
    errors = validate_pair(baseline, candidate)
    assert "baseline and candidate reuse the same binary_sha256" in errors
    assert "baseline/candidate query_sha256 differs" in errors


def test_rejects_incomplete_row_accounting():
    baseline, candidate = pair()
    candidate["n_rows_in_dump"] = 1
    errors = validate_pair(baseline, candidate)
    assert "candidate: dump row accounting is incomplete" in errors


def test_reports_candidate_residual_regression_as_an_adoption_failure(tmp_path):
    import json
    import subprocess
    import sys
    from pathlib import Path

    baseline, candidate = pair()
    baseline["bucket_counts"] = {"agree_all": 3, "both_disagree_same_as_default": 1}
    candidate["bucket_counts"] = {"agree_all": 2, "both_disagree_same_as_default": 2}
    baseline_path = tmp_path / "baseline.json"
    candidate_path = tmp_path / "candidate.json"
    output = tmp_path / "result.json"
    baseline_path.write_text(json.dumps(baseline), encoding="utf-8")
    candidate_path.write_text(json.dumps(candidate), encoding="utf-8")
    script = Path(__file__).parents[1] / "check_smarts_baseline_candidate.py"
    result = subprocess.run(
        [sys.executable, str(script), str(baseline_path), str(candidate_path), "--output", str(output), "--require-non-regression"],
        capture_output=True,
        text=True,
    )
    packet = json.loads(output.read_text(encoding="utf-8"))
    assert result.returncode == 1
    assert packet["comparison_valid"] is True
    assert packet["adoption_allowed"] is False
    assert "residual cells regressed" in packet["adoption_errors"][0]
