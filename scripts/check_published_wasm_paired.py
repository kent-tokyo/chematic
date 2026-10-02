#!/usr/bin/env python3
"""Validate the archived published-WASM 20-block Node timing record offline.

This checks recorded data and arithmetic, not that the machine was idle or
that Node measurements transfer to browsers or other hosts.
"""

from __future__ import annotations

import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "benchmarks/2026-10-03-v1029-v1030-rdkitjs-node-isolated-paired20.json"
RECORD_SHA256 = "edb2b289698dc022db9791b0307f6cb001b141d961045324a524b60e0882ae37"
RUNNER = ROOT / "scripts/bench_published_wasm_paired.mjs"
CORPUS = ROOT / "scripts/descriptor_census_corpus.smi"
EXPECTED = {
    "v29": ("1.0.29", "d0af78a8d6b711a13b63985d4b36079e245441e99f1712532555dc69f11569fa",
            "4ca57d6e9b9c8cc703d854d0265550f5ef355d52c2bf56a1a7b8bfcc09f21eb4",
            "78d3a205d5b4a8914ff3fc9dab3832a5b668e2b886ecda5dc5c566a48b2ba919"),
    "v30": ("1.0.30", "fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201",
            "4ca57d6e9b9c8cc703d854d0265550f5ef355d52c2bf56a1a7b8bfcc09f21eb4",
            "e0a505c33a282724114e279201dae199b2e5dbe68c35c0ebaa474693a6eb4ae6"),
    "rdkit": ("2026.03.6", "3b86b72775394ae997fb96ca854c6e7c5840927d3e899c1ef117d42bca573c87",
              "25ea1529b8eeb106ac1f3acb7c5d338a90533248c6223614574fefd14fee772f",
              "f5981af4cf5bdcef861ffdf7410330bf9887f1195e39fbf9901b3a2f134ace2e"),
}
LANES = ["parse_only", "parse_morgan", "prepared_first_use", "prepared_reused"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def close(actual: float, expected: float, message: str) -> None:
    require(math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12), message)


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low = math.floor(position)
    return ordered[low] + (ordered[min(low + 1, len(ordered) - 1)] - ordered[low]) * (position - low)


def bootstrap(log_ratios: list[float]) -> list[float]:
    state = 0x9E3779B9
    sample = []
    for _ in range(10_000):
        total = 0.0
        for _ in log_ratios:
            state ^= (state << 13) & 0xFFFFFFFF
            state &= 0xFFFFFFFF
            state ^= state >> 17
            state &= 0xFFFFFFFF
            state ^= (state << 5) & 0xFFFFFFFF
            state &= 0xFFFFFFFF
            total += log_ratios[state % len(log_ratios)]
        sample.append(math.exp(total / len(log_ratios)))
    return [quantile(sample, 0.025), quantile(sample, 0.975)]


def check() -> dict:
    require(digest(RECORD) == RECORD_SHA256, "archived record SHA-256")
    report = json.loads(RECORD.read_text(encoding="utf-8"))
    require(report["schema"] == "published-wasm-paired-speed/v1", "schema")
    require(report["runner_sha256"] == digest(RUNNER), "runner SHA-256")
    require(report["corpus"] == {"path": "scripts/descriptor_census_corpus.smi",
                                 "sha256": digest(CORPUS), "rows": 250}, "pinned corpus")
    require(report["environment"]["platform"] == "darwin"
            and report["environment"]["arch"] == "arm64"
            and report["environment"]["node"].startswith("v24."), "recorded platform")
    for name, (version, tarball, js, wasm) in EXPECTED.items():
        artifact = report["artifacts"][name]
        require(artifact["version"] == version and artifact["tarball_sha256"] == tarball
                and artifact["js_sha256"] == js and artifact["wasm_sha256"] == wasm, f"{name} artifact")
    protocol = report["protocol"]
    require(protocol["paired_blocks"] == 20 and protocol["processes_per_block"] == 4
            and protocol["lanes"] == LANES and protocol["warmup_rows"] == 20
            and protocol["order"] == "ABBA then BAAB alternating", "protocol")
    preflight = report["preflight"]
    require(preflight["compared_rows"] == preflight["exact_atom_and_morgan_rows"] == 250, "preflight accounting")
    reference = preflight["by_artifact"]["v30"]
    require(len(reference) == 250 and all(isinstance(row[0], int) and row[0] > 0
            and isinstance(row[1], str) and len(row[1]) == 64 for row in reference), "preflight rows")
    require(all(preflight["by_artifact"][name] == reference for name in EXPECTED), "row-level output equivalence")
    require(preflight["atom_count_sum"] == sum(row[0] for row in reference), "preflight atom count")
    require(isinstance(preflight["morgan_rows_sha256"], str)
            and len(preflight["morgan_rows_sha256"]) == 64, "Morgan digest")
    comparisons = report["comparisons"]
    require([(item["left"], item["right"], item["lane"]) for item in comparisons]
            == [(left, right, lane) for lane in LANES
                for left, right in [("v29", "v30"), ("v30", "rdkit")]], "comparison matrix")
    for item in comparisons:
        left, right, lane = item["left"], item["right"], item["lane"]
        blocks = item["blocks"]
        require(len(blocks) == 20, f"{lane} block count")
        ratios = []
        samples_by_name = {left: [], right: []}
        for index, block in enumerate(blocks):
            expected_order = [left, right, right, left] if index % 2 == 0 else [right, left, left, right]
            require(block["index"] == index and block["order"] == expected_order
                    and [sample["name"] for sample in block["samples"]] == expected_order,
                    f"{lane} block {index} order")
            for sample in block["samples"]:
                require(sample["lane"] == lane and sample["parsed_rows"] == preflight["compared_rows"],
                        f"{lane} block {index} output count")
                require((sample["morgan_rows_sha256"] is None if lane == "parse_only"
                         else sample["morgan_rows_sha256"] == preflight["morgan_rows_sha256"]),
                        f"{lane} block {index} Morgan output")
                require(all(math.isfinite(sample[key]) and sample[key] > 0 for key in
                            ["elapsed_ms", "init_ms", "peak_rss_bytes", "rss_after_bytes"]),
                        f"{lane} block {index} measurement")
                samples_by_name[sample["name"]].append(sample)
            times = {name: statistics.mean(sample["elapsed_ms"] for sample in block["samples"]
                                            if sample["name"] == name) for name in (left, right)}
            ratio = times[right] / times[left]
            close(block["left_ms"], times[left], f"{lane} left timing")
            close(block["right_ms"], times[right], f"{lane} right timing")
            close(block["speed_ratio_left_over_right"], ratio, f"{lane} ratio")
            ratios.append(ratio)
        logs = [math.log(ratio) for ratio in ratios]
        summary = item["summary"]
        ci = bootstrap(logs)
        close(summary["median_speed_ratio_left_over_right"], quantile(ratios, 0.5), f"{lane} median")
        close(summary["geometric_mean_ratio"], math.exp(statistics.mean(logs)), f"{lane} geometric mean")
        for actual, expected in zip(summary["bootstrap_95_ci"], ci, strict=True):
            close(actual, expected, f"{lane} interval")
        require(summary["statistical_speed_signal"] == (ci[0] > 1)
                and summary["practical_10pct_speed_win"] == (ci[0] > 1.10)
                and summary["interval_crosses_parity"] == (ci[0] <= 1 <= ci[1]), f"{lane} claim flags")
        for name, prefix in ((left, "left"), (right, "right")):
            close(summary[f"{prefix}_peak_rss_median_bytes"],
                  quantile([sample["peak_rss_bytes"] for sample in samples_by_name[name]], 0.5),
                  f"{lane} RSS median")
            close(summary[f"{prefix}_init_median_ms"],
                  quantile([sample["init_ms"] for sample in samples_by_name[name]], 0.5),
                  f"{lane} init median")
    return {"comparisons": len(comparisons), "blocks_each": 20,
            "exact_preflight_rows": preflight["exact_atom_and_morgan_rows"]}


if __name__ == "__main__":
    print("published WASM paired gate OK:", check())
