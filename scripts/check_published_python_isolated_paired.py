#!/usr/bin/env python3
"""Check the six archived 20-block fresh-process Python comparison lanes."""

from __future__ import annotations

import json
from pathlib import Path

if __package__:
    from .bench_published_python_isolated_paired import (
        ARTIFACT_SHA256, REFERENCE_SHA256, ci, reference_digest, sha256,
    )
else:
    from bench_published_python_isolated_paired import (
        ARTIFACT_SHA256, REFERENCE_SHA256, ci, reference_digest, sha256,
    )

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts/chembl_accuracy_corpus_4999.smi"
RESULTS = ROOT / "validation/results"
RECORDS = {
    "2026-10-03-v1.0.29-vs-v1.0.30-isolated-hba-parse-inclusive.json":
        ("hba", "parse_inclusive", ("chematic", "1.0.29"), ("chematic", "1.0.30")),
    "2026-10-03-v1.0.29-vs-v1.0.30-isolated-morgan-parse-inclusive.json":
        ("morgan", "parse_inclusive", ("chematic", "1.0.29"), ("chematic", "1.0.30")),
    "2026-10-03-v1.0.30-vs-rdkit-isolated-hba-parse-inclusive.json":
        ("hba", "parse_inclusive", ("rdkit", "2026.03.6"), ("chematic", "1.0.30")),
    "2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-parse-inclusive.json":
        ("morgan", "parse_inclusive", ("rdkit", "2026.03.6"), ("chematic", "1.0.30")),
    "2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-prepared-first-use.json":
        ("morgan", "prepared_first_use", ("rdkit", "2026.03.6"), ("chematic", "1.0.30")),
    "2026-10-03-v1.0.30-vs-rdkit-isolated-morgan-precomputed.json":
        ("morgan", "precomputed", ("rdkit", "2026.03.6"), ("chematic", "1.0.30")),
}


def expected_outcome(output_equal: bool, mode: str, different_library: bool, interval: list[float]) -> str:
    if not output_equal:
        return "not_counted_output_difference"
    if mode == "prepared_first_use" and different_library:
        return "not_counted_eager_lazy_perception"
    if interval[0] > 1:
        return "arm_b_faster_on_declared_lane"
    if interval[1] < 1:
        return "arm_a_faster_on_declared_lane"
    return "inconclusive_interval_crosses_parity"


def check_one(path: Path, spec: tuple) -> None:
    operation, mode, arm_a, arm_b = spec
    report = json.loads(path.read_text(encoding="utf-8"))
    if report["schema"] != "published-python-isolated-paired-v1" or \
            report["operation"] != operation or report["mode"] != mode or \
            report["corpus"] != {"sha256": sha256(CORPUS), "limit": 5000} or \
            report["protocol"]["blocks"] != 20 or len(report["pairs"]) != 20:
        raise ValueError(f"{path.name}: scope or denominator differs")

    expected_digests = {}
    for label, identity in (("a", arm_a), ("b", arm_b)):
        library, version = identity
        metadata = report["artifacts"][label]
        reference = RESULTS / f"v{version if library == 'chematic' else '1.0.30'}-published-python-63op-outputs.jsonl.gz"
        expected_digest = reference_digest(reference, operation, 5000)
        if (metadata["library"], metadata["version"]) != identity or \
                metadata["artifact_sha256"] != ARTIFACT_SHA256[identity] or \
                metadata["reference_archive_sha256"] != REFERENCE_SHA256[identity] or \
                metadata["expected_output_sha256"] != expected_digest or \
                sha256(reference) != REFERENCE_SHA256[identity]:
            raise ValueError(f"{path.name}: {label} artifact or reference mismatch")
        expected_digests[label] = expected_digest

    for block, pair in enumerate(report["pairs"]):
        order = ["a", "b"] if block % 2 == 0 else ["b", "a"]
        if pair["block"] != block or pair["order"] != order:
            raise ValueError(f"{path.name}: block {block} order")
        for label in ("a", "b"):
            row = pair[label]
            library, version = arm_a if label == "a" else arm_b
            if row["library"] != library or row["version"] != version or \
                    row["operation"] != operation or row["mode"] != mode or \
                    row["input_count"] != 5000 or row["output_sha256"] != expected_digests[label] or \
                    row["operation_ns"] <= 0 or row["setup_ns"] < 0 or \
                    row["peak_process_rss_bytes"] < row["rss_after_setup_highwater_bytes"] or \
                    row["rss_after_setup_highwater_bytes"] <= 0:
                raise ValueError(f"{path.name}: block {block} {label} invalid")

    speed = ci([p["a"]["operation_ns"] / p["b"]["operation_ns"] for p in report["pairs"]])
    memory = ci([p["a"]["peak_process_rss_bytes"] / p["b"]["peak_process_rss_bytes"]
                 for p in report["pairs"]])
    output_equal = expected_digests["a"] == expected_digests["b"]
    if report["speed_a_over_b"] != speed or report["memory_a_over_b"] != memory or \
            report["outputs_equal_across_arms"] != output_equal or \
            report["speed_outcome"] != expected_outcome(output_equal, mode, arm_a[0] != arm_b[0],
                                                         speed["bootstrap_95pct"]):
        raise ValueError(f"{path.name}: interval or claim classification differs")


def main() -> None:
    for name, spec in RECORDS.items():
        check_one(ROOT / "benchmarks" / name, spec)
    print("six isolated published-Python lanes OK: 20 paired blocks each, pinned outputs and separate RSS")


if __name__ == "__main__":
    main()
