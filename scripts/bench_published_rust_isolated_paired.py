#!/usr/bin/env python3
"""Measure two exact published Rust graphs in independent paired processes.

Only HBA and compatible Morgan are supported. The emitted digest is checked
against archived published-Rust 63-operation outputs before a time ratio is
eligible for an equivalent-output claim. Peak RSS is process-wide, not an
operation-allocation estimate.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
import tomllib
from pathlib import Path

if __package__:
    from .bench_published_python_isolated_paired import ci, reference_digest, sha256
else:
    from bench_published_python_isolated_paired import ci, reference_digest, sha256

ROOT = Path(__file__).resolve().parents[1]
CORPUS_SHA256 = "1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e"
VERSIONS = {
    "1.0.29": {
        "manifest": ROOT / "tools/published_rust_gate_v1029/Cargo.toml",
        "lock": ROOT / "tools/published_rust_gate_v1029/Cargo.lock",
        "lock_sha256": "dfa39d95c4f6fca601337ed103bfff74ea048290b1b2f1bcb27d2c633f24caf4",
        "reference": ROOT / "validation/results/v1.0.29-published-rust-63op-outputs.jsonl.gz",
        "reference_sha256": "01f6563984af9889821d94a29d1765b0539703eefa38a1efe86f9391e1bb0be8",
    },
    "1.0.30": {
        "manifest": ROOT / "tools/published_rust_gate/Cargo.toml",
        "lock": ROOT / "tools/published_rust_gate/Cargo.lock",
        "lock_sha256": "0cb18ae64edaac9352c9c633c918ed3d48d3d9b6b1098457dade1a78e1dd4f88",
        "reference": ROOT / "validation/results/v1.0.30-published-rust-63op-outputs.jsonl.gz",
        "reference_sha256": "1d54bb28b09d1ba4effe17264fcf7945cd508cb2500f509806948cff5b3d1037",
    },
}
MODES = ("parse_inclusive", "prepared_first_use", "precomputed")


def validate_graph(version: str) -> None:
    info = VERSIONS[version]
    if sha256(info["lock"]) != info["lock_sha256"] or sha256(info["reference"]) != info["reference_sha256"]:
        raise ValueError(f"{version}: lockfile or published output reference changed")
    lock = tomllib.loads(info["lock"].read_text(encoding="utf-8"))
    crates = [package for package in lock["package"] if package["name"] == "chematic" or package["name"].startswith("chematic-")]
    if len(crates) < 15 or any(package["version"] != version or not package.get("checksum") for package in crates):
        raise ValueError(f"{version}: mixed or unpinned published crate graph")
    manifest = tomllib.loads(info["manifest"].read_text(encoding="utf-8"))
    if manifest["dependencies"]["chematic"]["version"] != f"={version}":
        raise ValueError(f"{version}: harness does not require its exact crate version")


def build(version: str) -> Path:
    info = VERSIONS[version]
    subprocess.run(["cargo", "build", "--release", "--locked", "--offline", "--manifest-path",
                    str(info["manifest"]), "--bin", "isolated_speed"], check=True, cwd=ROOT)
    name = "isolated_speed.exe" if sys.platform == "win32" else "isolated_speed"
    binary = info["manifest"].parent / "target/release" / name
    if not binary.is_file():
        raise ValueError(f"{version}: expected binary is missing")
    return binary


def time_command(binary: Path, corpus: Path, limit: int, operation: str, mode: str) -> list[str]:
    command = [str(binary), str(corpus), str(limit), operation, mode]
    if sys.platform == "darwin":
        return ["/usr/bin/time", "-l", *command]
    if sys.platform.startswith("linux"):
        return ["/usr/bin/time", "-v", *command]
    raise ValueError("RSS capture requires macOS or Linux /usr/bin/time")


def read_peak_rss(stderr: str) -> int:
    if sys.platform == "darwin":
        match = re.search(r"^\s*(\d+)\s+maximum resident set size\s*$", stderr, re.MULTILINE)
        scale = 1
    else:
        match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", stderr)
        scale = 1024
    if match is None or int(match.group(1)) <= 0:
        raise ValueError(f"peak RSS was not recorded: {stderr[-500:]}")
    return int(match.group(1)) * scale


def run_arm(binary: Path, corpus: Path, limit: int, operation: str, mode: str, expected: str) -> dict:
    completed = subprocess.run(time_command(binary, corpus, limit, operation, mode),
                               check=True, capture_output=True, text=True)
    value = json.loads(completed.stdout)
    if (value["output_sha256"] != expected or value["input_count"] != limit
            or value["operation"] != operation or value["mode"] != mode
            or value["operation_ns"] <= 0):
        raise ValueError("published-output or measurement mismatch")
    value["peak_process_rss_bytes"] = read_peak_rss(completed.stderr)
    return value


def measure(args: argparse.Namespace) -> dict:
    if args.blocks < 20 or args.limit < 1 or args.mode not in MODES or args.operation not in {"hba", "morgan"}:
        raise ValueError("requires at least 20 blocks and a declared operation/mode")
    if sha256(args.corpus) != CORPUS_SHA256:
        raise ValueError("exposed corpus hash changed")
    arms = {}
    for label, version in (("a", "1.0.29"), ("b", "1.0.30")):
        validate_graph(version)
        binary = build(version)
        reference = VERSIONS[version]["reference"]
        expected = reference_digest(reference, args.operation, args.limit)
        arms[label] = {"version": version, "binary": binary, "binary_sha256": sha256(binary),
                       "lock_sha256": sha256(VERSIONS[version]["lock"]),
                       "reference_archive_sha256": sha256(reference), "expected_output_sha256": expected}
        # Separate process for warmup; measured blocks never reuse its caches.
        run_arm(binary, args.corpus, args.limit, args.operation, args.mode, expected)
    pairs = []
    for block in range(args.blocks):
        order = ("a", "b") if block % 2 == 0 else ("b", "a")
        pair = {"block": block, "order": list(order)}
        for label in order:
            pair[label] = run_arm(arms[label]["binary"], args.corpus, args.limit,
                                  args.operation, args.mode, arms[label]["expected_output_sha256"])
        pairs.append(pair)
        print(f"block {block + 1}/{args.blocks}: {pair['a']['operation_ns']} / {pair['b']['operation_ns']} ns", flush=True)
    time_ratio = ci([pair["a"]["operation_ns"] / pair["b"]["operation_ns"] for pair in pairs])
    rss_ratio = ci([pair["a"]["peak_process_rss_bytes"] / pair["b"]["peak_process_rss_bytes"] for pair in pairs])
    outputs_equal = arms["a"]["expected_output_sha256"] == arms["b"]["expected_output_sha256"]
    if not outputs_equal:
        outcome = "not_counted_output_difference"
    elif time_ratio["bootstrap_95pct"][0] > 1:
        outcome = "v1030_faster_on_declared_lane"
    elif time_ratio["bootstrap_95pct"][1] < 1:
        outcome = "v1029_faster_on_declared_lane"
    else:
        outcome = "inconclusive_interval_crosses_parity"
    serializable_arms = {label: {key: value for key, value in arm.items() if key != "binary"}
                         for label, arm in arms.items()}
    return {"schema": "published-rust-isolated-paired-v1", "artifacts": serializable_arms,
            "corpus": {"sha256": sha256(args.corpus), "limit": args.limit},
            "operation": args.operation, "mode": args.mode,
            "environment": {"platform": platform.platform(), "python": sys.version.split()[0],
                            "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
                            "cargo": subprocess.check_output(["cargo", "--version"], text=True).strip()},
            "protocol": {"blocks": args.blocks, "order": "AB/BA alternating fresh processes",
                         "warmup": "one discarded process per arm, plus disjoint in-process molecule",
                         "timing": "Rust Instant around operation calls; parse only in parse_inclusive",
                         "memory": "OS whole-process peak RSS, including startup and preparation; not operation allocation"},
            "outputs_equal_across_arms": outputs_equal, "speed_a_over_b": time_ratio,
            "speed_outcome": outcome, "memory_a_over_b": rss_ratio,
            "memory_interpretation": "whole-process peak only; no per-operation allocation claim",
            "pairs": pairs}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=ROOT / "scripts/chembl_accuracy_corpus_4999.smi")
    parser.add_argument("--limit", type=int, default=5000)
    parser.add_argument("--operation", choices=("hba", "morgan"), required=True)
    parser.add_argument("--mode", choices=MODES, required=True)
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = measure(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"speed_a_over_b": result["speed_a_over_b"],
                      "memory_a_over_b": result["memory_a_over_b"],
                      "speed_outcome": result["speed_outcome"]}, sort_keys=True))


if __name__ == "__main__":
    main()
