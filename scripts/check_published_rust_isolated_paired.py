#!/usr/bin/env python3
"""Recheck all archived published-Rust paired time and process-RSS blocks."""

from __future__ import annotations

import json
import math
from pathlib import Path

if __package__:
    from .bench_published_python_isolated_paired import ci, reference_digest, sha256
    from .bench_published_rust_isolated_paired import CORPUS_SHA256, MODES, ROOT, VERSIONS, validate_graph
else:
    from bench_published_python_isolated_paired import ci, reference_digest, sha256
    from bench_published_rust_isolated_paired import CORPUS_SHA256, MODES, ROOT, VERSIONS, validate_graph

PREFIX = "2026-10-03-rust-v1029-v1030"
BINARIES = {
    "1.0.29": "67ce3c57e368c3e377cc25554411078a490b85cc980cadb00ce4ce3a8a7b7b90",
    "1.0.30": "ff0133ccb25043a9ae171275ee348600f2424acc18c3766fcf2e375a5f20227e",
}


def check_report(record: dict, operation: str, mode: str) -> None:
    if (record.get("schema") != "published-rust-isolated-paired-v1"
            or record.get("operation") != operation or record.get("mode") != mode
            or record.get("corpus") != {"sha256": CORPUS_SHA256, "limit": 5000}
            or record.get("protocol", {}).get("blocks") != 20
            or record.get("protocol", {}).get("order") != "AB/BA alternating fresh processes"):
        raise ValueError(f"{operation}/{mode}: report identity or protocol changed")
    arms = record["artifacts"]
    for label, version in (("a", "1.0.29"), ("b", "1.0.30")):
        info = VERSIONS[version]
        expected = reference_digest(info["reference"], operation, 5000)
        arm = arms[label]
        if (arm["version"] != version
                or arm["lock_sha256"] != info["lock_sha256"]
                or arm["reference_archive_sha256"] != info["reference_sha256"]
                or arm["expected_output_sha256"] != expected
                or arm["binary_sha256"] != BINARIES[version]):
            raise ValueError(f"{operation}/{mode}: {label} artifact or output reference changed")
    same_output = arms["a"]["expected_output_sha256"] == arms["b"]["expected_output_sha256"]
    if record["outputs_equal_across_arms"] != same_output:
        raise ValueError(f"{operation}/{mode}: output equality claim changed")
    pairs = record["pairs"]
    if len(pairs) != 20:
        raise ValueError(f"{operation}/{mode}: 20 paired blocks required")
    for index, pair in enumerate(pairs):
        order = ["a", "b"] if index % 2 == 0 else ["b", "a"]
        if pair["block"] != index or pair["order"] != order:
            raise ValueError(f"{operation}/{mode}: execution order changed at block {index}")
        for label in order:
            row = pair[label]
            if (row["operation"] != operation or row["mode"] != mode
                    or row["input_count"] != 5000
                    or row["output_sha256"] != arms[label]["expected_output_sha256"]
                    or not all(isinstance(row[key], int) and row[key] > 0
                               for key in ("operation_ns", "peak_process_rss_bytes"))
                    or not isinstance(row["setup_ns"], int) or row["setup_ns"] < 0):
                raise ValueError(f"{operation}/{mode}: invalid arm {label} block {index}")
    speed = ci([p["a"]["operation_ns"] / p["b"]["operation_ns"] for p in pairs])
    memory = ci([p["a"]["peak_process_rss_bytes"] / p["b"]["peak_process_rss_bytes"] for p in pairs])
    if record["speed_a_over_b"] != speed or record["memory_a_over_b"] != memory:
        raise ValueError(f"{operation}/{mode}: paired statistics changed")
    expected_outcome = ("not_counted_output_difference" if not same_output else
                        "v1030_faster_on_declared_lane" if speed["bootstrap_95pct"][0] > 1 else
                        "v1029_faster_on_declared_lane" if speed["bootstrap_95pct"][1] < 1 else
                        "inconclusive_interval_crosses_parity")
    if (record["speed_outcome"] != expected_outcome
            or not all(math.isfinite(value) for value in speed["bootstrap_95pct"] + memory["bootstrap_95pct"])):
        raise ValueError(f"{operation}/{mode}: speed classification changed")


def check_all() -> int:
    corpus = ROOT / "scripts/chembl_accuracy_corpus_4999.smi"
    if sha256(corpus) != CORPUS_SHA256:
        raise ValueError("exposed corpus changed")
    for version in VERSIONS:
        validate_graph(version)
    for operation in ("hba", "morgan"):
        for mode in MODES:
            path = ROOT / "benchmarks" / f"{PREFIX}-{operation}-{mode.replace('_', '-')}-paired20.json"
            check_report(json.loads(path.read_text(encoding="utf-8")), operation, mode)
    print("six published-Rust lanes OK: exact 1.0.29/1.0.30 graphs, outputs, 20 paired blocks, separate RSS")
    return 0


if __name__ == "__main__":
    raise SystemExit(check_all())
