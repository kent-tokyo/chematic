#!/usr/bin/env python3
"""Verify archived crates.io operation rows against published Python wheels.

The Rust runner uses one pinned Cargo graph per release. This gate checks
outputs, not timing, RDKit parity, or general reaction compatibility.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import tomllib
import argparse
from pathlib import Path

if __package__:
    from .check_published_python_version_outputs import compare
else:
    from check_published_python_version_outputs import compare

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation/results"
CORPUS = ROOT / "scripts/chembl_accuracy_corpus_4999.smi"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rows(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle]


def verify(version: str) -> dict:
    stem = f"v{version}-published-rust-63op-outputs"
    raw_path = RESULTS / f"{stem}.jsonl.gz"
    summary_path = RESULTS / f"{stem}-summary.json"
    py_path = RESULTS / f"v{version}-published-python-63op-outputs.jsonl.gz"
    gate = "published_rust_gate_v1029" if version == "1.0.29" else "published_rust_gate"
    lock_path = ROOT / "tools" / gate / "Cargo.lock"

    lock_bytes = lock_path.read_bytes()
    packages = tomllib.loads(lock_bytes.decode("utf-8"))["package"]
    chematic_packages = [p for p in packages if p["name"] == "chematic" or p["name"].startswith("chematic-")]
    if len(chematic_packages) < 14 or any(p["version"] != version or
            p.get("source") != "registry+https://github.com/rust-lang/crates.io-index" or
            len(p.get("checksum", "")) != 64 for p in chematic_packages):
        raise ValueError(f"v{version}: Cargo graph is mixed, local, or incomplete")

    compressed = raw_path.read_bytes()
    rust = rows(raw_path)
    baseline = rows(py_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary["schema"] != "published-rust-python-63op-output-slice/v1" or
            summary["crate"] != f"chematic {version} from crates.io" or
            summary["corpus_sha256"] != sha256(CORPUS.read_bytes()) or
            summary["input_count"] != 5000 or summary["adapted_operations"] != 63 or
            summary["rows_sha256"] != sha256(gzip.decompress(compressed))):
        raise ValueError(f"v{version}: summary or corpus digest mismatch")
    if len(rust) != 63 or len(baseline) != 63:
        raise ValueError(f"v{version}: expected 63 operations")

    count = 0
    for actual, expected in zip(rust, baseline, strict=True):
        if actual["op"] != expected["op"] or len(actual["rows"]) != len(expected["rows"]):
            raise ValueError(f"v{version}: operation name or denominator mismatch")
        for index, (a, b) in enumerate(zip(actual["rows"], expected["rows"], strict=True)):
            if a.get("input_index") != index or b.get("input_index") != index or \
                    (("value" in a) == ("error" in a)):
                raise ValueError(f"v{version} {actual['op']}: malformed row {index}")
            if "error" in a or "error" in b:
                raise ValueError(f"v{version} {actual['op']}: operation error at {index}")
            if a != b:
                raise ValueError(f"v{version} {actual['op']}: Rust/Python output differs at {index}")
            count += 1
    if count != 210410:
        raise ValueError(f"v{version}: expected 210410 rows, got {count}")
    return {"version": version, "operation_count": 63, "output_count": count,
            "rust_python_different_rows": 0,
            "rust_archive_sha256": sha256(compressed),
            "rust_rows_sha256": summary["rows_sha256"],
            "python_archive_sha256": sha256(py_path.read_bytes()),
            "cargo_lock_sha256": sha256(lock_bytes),
            "crate_checksums": {p["name"]: p["checksum"] for p in chematic_packages}}


def build_report() -> dict:
    releases = [verify(version) for version in ("1.0.29", "1.0.30")]
    diff = compare(RESULTS / "v1.0.29-published-rust-63op-outputs.jsonl.gz",
                   RESULTS / "v1.0.30-published-rust-63op-outputs.jsonl.gz")
    python_diff = json.loads((RESULTS / "v1.0.29-to-v1.0.30-python-63op-diff.json").read_text(encoding="utf-8"))
    if diff["changed_rows_by_operation"] != python_diff["changed_rows_by_operation"] or \
            diff["output_count"] != 210410:
        raise ValueError("published Rust and Python version differentials disagree")
    return {"schema": "published-rust-63op-cross-binding-and-version-diff/v1",
              "scope": "output identity only; no RDKit or speed claim",
              "releases": releases,
              "version_differential": diff["changed_rows_by_operation"],
              "unexpected_changed_operations": diff["unexpected_changed_operations"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true", help="regenerate the archived report")
    args = parser.parse_args()
    report = build_report()
    destination = RESULTS / "v1.0.29-to-v1.0.30-published-rust-63op-diff.json"
    if args.write_report:
        destination.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elif json.loads(destination.read_text(encoding="utf-8")) != report:
        raise ValueError("archived published Rust report differs from raw rows")
    print("checked 2 published Rust graphs, 63 operations and 210410 rows each")


if __name__ == "__main__":
    main()
