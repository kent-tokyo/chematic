import json

import pytest

from scripts.check_published_rust_isolated_paired import ROOT, check_all, check_report


def record(name: str) -> dict:
    return json.loads((ROOT / "benchmarks" / name).read_text(encoding="utf-8"))


def test_all_archived_published_rust_lanes():
    assert check_all() == 0


def test_output_mismatch_cannot_earn_speed_claim():
    report = record("2026-10-03-rust-v1029-v1030-morgan-parse-inclusive-paired20.json")
    report["pairs"][0]["a"]["output_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="invalid arm a block 0"):
        check_report(report, "morgan", "parse_inclusive")


def test_order_and_interval_are_recomputed():
    report = record("2026-10-03-rust-v1029-v1030-morgan-parse-inclusive-paired20.json")
    report["pairs"][1]["order"] = ["a", "b"]
    with pytest.raises(ValueError, match="execution order"):
        check_report(report, "morgan", "parse_inclusive")
    report = record("2026-10-03-rust-v1029-v1030-morgan-parse-inclusive-paired20.json")
    report["speed_a_over_b"]["median_a_over_b"] = 99
    with pytest.raises(ValueError, match="paired statistics"):
        check_report(report, "morgan", "parse_inclusive")
