"""Fail-closed output and interval checks for isolated published-wheel timing."""

import gzip
import json

import pytest

from scripts.bench_published_python_isolated_paired import ci, digest, reference_digest
from scripts.check_published_python_isolated_paired import RECORDS, ROOT, check_one, expected_outcome


def test_reference_digest_requires_contiguous_error_free_rows(tmp_path):
    path = tmp_path / "outputs.jsonl.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps({"op": "hba", "rows": [
            {"input_index": 0, "value": 2},
            {"input_index": 1, "value": 3},
        ]}) + "\n")
    assert reference_digest(path, "hba", 2) == digest(["2", "3"])
    with pytest.raises(ValueError, match="expected 3 rows"):
        reference_digest(path, "hba", 3)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(json.dumps({"op": "hba", "rows": [
            {"input_index": 0, "value": 2},
            {"input_index": 1, "error": "bad input"},
        ]}) + "\n")
    with pytest.raises(ValueError, match="invalid row 1"):
        reference_digest(path, "hba", 2)


def test_paired_interval_requires_twenty_positive_blocks():
    with pytest.raises(ValueError, match="20 positive"):
        ci([2.0] * 19)
    with pytest.raises(ValueError, match="20 positive"):
        ci([2.0] * 19 + [0.0])
    result = ci([2.0] * 20)
    assert result["median_a_over_b"] == 2.0
    assert result["bootstrap_95pct"] == [2.0, 2.0]


def test_claim_boundary_and_archived_order_tamper(tmp_path):
    assert expected_outcome(False, "parse_inclusive", False, [2.0, 3.0]) == "not_counted_output_difference"
    assert expected_outcome(True, "prepared_first_use", True, [2.0, 3.0]) == "not_counted_eager_lazy_perception"
    assert expected_outcome(True, "parse_inclusive", True, [0.9, 1.1]) == "inconclusive_interval_crosses_parity"
    name = "2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-parse-inclusive.json"
    source = ROOT / "benchmarks" / name
    record = json.loads(source.read_text(encoding="utf-8"))
    record["pairs"][1]["order"] = ["a", "b"]
    altered = tmp_path / name
    altered.write_text(json.dumps(record), encoding="utf-8")
    with pytest.raises(ValueError, match="block 1 order"):
        check_one(altered, RECORDS[name])
