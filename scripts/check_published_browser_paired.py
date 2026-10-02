#!/usr/bin/env python3
"""Check paired Chromium timings for pinned published CheMatic/RDKit.js artifacts.

Per-run fingerprint digests are checked before timing is interpreted. Browser
whole-process RSS is retained separately from speed and library allocations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path

if __package__:
    from .check_published_wasm_paired import bootstrap, quantile
else:
    from check_published_wasm_paired import bootstrap, quantile

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts/descriptor_census_corpus.smi"
OPERATIONS = ("parse", "parse_write", "parse_fp", "prepared_fp")
EXPECTED = {
    "schematic_wasm": "e0a505c33a282724114e279201dae199b2e5dbe68c35c0ebaa474693a6eb4ae6",
    "schematic_js": "4ca57d6e9b9c8cc703d854d0265550f5ef355d52c2bf56a1a7b8bfcc09f21eb4",
    "rdkit_wasm": "f5981af4cf5bdcef861ffdf7410330bf9887f1195e39fbf9901b3a2f134ace2e",
    "rdkit_js": "25ea1529b8eeb106ac1f3acb7c5d338a90533248c6223614574fefd14fee772f",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(source: Path) -> dict:
    document = json.loads(source.read_bytes())
    config = document["configuration"]
    artifacts = document["artifacts"]
    prepared_mode = config.get("prepared_mode", "mixed")
    if prepared_mode not in {"mixed", "split"}:
        raise ValueError("unknown prepared-molecule timing mode")
    operations = (OPERATIONS if prepared_mode == "mixed" else
                  ("parse", "parse_write", "parse_fp", "prepared_fp_first_use", "prepared_fp_reused"))
    fingerprint_operations = ("parse_fp", "prepared_fp") if prepared_mode == "mixed" else (
        "parse_fp", "prepared_fp_first_use", "prepared_fp_reused")
    if (document["schema_version"] != 2
            or document["gate"] != "wasm-vs-official-rdkit-browser-isolated"
            or config["engine"] != "chromium"
            or config["rows"] != 250 or config["warmup_rows"] != 20
            or config["repetitions"] != 20
            or config["timing_method"] != "single_batch"
            or config["process_rss"] != "summed fresh Chromium process-tree RSS sampled every 50 ms"
            or document["corpus"]["sha256"] != sha256(CORPUS)
            or document["corpus"]["rows"] != 250
            or config["fingerprint_typed_unsupported_rows"] != 0):
        raise ValueError("browser protocol or corpus differs from the frozen lane")
    if any(artifacts[name]["sha256"] != expected for name, expected in EXPECTED.items()):
        raise ValueError("published JS/WASM hash mismatch")
    package = artifacts["schematic_package"]
    if (package["name"] != "@kent-tokyo/chematic" or package["version"] != "1.0.30"
            or package["kind"] != "published"
            or package["sha256"] != "6a18eb0fa882b4a2772751c9484db8200deca21764f0a9418d8a94196d6b5b90"
            or package["tarball"]["sha256"] != "fd427a161c11b14e6cca51e78e0c35ff05bbc52faf0539a2ddc2b7090e534201"
            or artifacts["rdkit_package"]["name"] != "@rdkit/rdkit"
            or artifacts["rdkit_package"]["version"] != "2026.03.6"
            or artifacts["rdkit_package"]["sha256"] != "bb938eb6b28bf39baa9778f56aefcdea746d28060fe2c542209f877a71468fed"):
        raise ValueError("browser package identity differs from the pinned registry artifacts")
    by_arm = {}
    for arm in ("chematic", "rdkit"):
        runs = document["raw_runs"][arm]
        if len(runs) != 20 or sorted(run["repetition"] for run in runs) != list(range(20)):
            raise ValueError(f"{arm}: missing or duplicate paired repetition")
        by_arm[arm] = {run["repetition"]: run for run in runs}
        for run in runs:
            if run["arm"] != arm or run["accepted_rows"] != 250:
                raise ValueError(f"{arm}: wrong result arm or row count")
            if run["process_tree_rss"]["status"] != "measured" or run["process_tree_rss"]["samples"] < 1:
                raise ValueError(f"{arm}: missing process RSS observation")
            if not math.isfinite(run["process_tree_rss"]["peak_rss_bytes"]) or run["process_tree_rss"]["peak_rss_bytes"] <= 0:
                raise ValueError(f"{arm}: invalid process RSS")
            if not math.isfinite(run["init_ms"]) or run["init_ms"] <= 0:
                raise ValueError(f"{arm}: invalid initialization timing")
            for operation in operations:
                item = run["operations"][operation]
                if (item["input_rows"] != 250 or item["count"] != 250
                        or not math.isfinite(item["mean_ms"]) or item["mean_ms"] <= 0
                        or not math.isfinite(item["batch_elapsed_ms"]) or item["batch_elapsed_ms"] <= 0
                        or not math.isclose(item["batch_elapsed_ms"] / 250, item["mean_ms"], rel_tol=1e-10)):
                    raise ValueError(f"{arm}/{operation}: invalid timed sample")
    reference = by_arm["chematic"][0]["operations"]["parse_fp"]["output"]
    if reference != {"format": "packed-lsb-first-2048-bit", "bytes_per_row": 256,
                     "fnv1a32": "5c006dbe", "set_bits": 11307}:
        raise ValueError("250-row Morgan output digest changed")
    for repetition in range(20):
        for operation in fingerprint_operations:
            for arm in ("chematic", "rdkit"):
                if by_arm[arm][repetition]["operations"][operation]["output"] != reference:
                    raise ValueError(f"{arm}/{operation}/{repetition}: unequal fingerprint output")

    lanes = {}
    for operation in operations:
        ratios = [by_arm["rdkit"][i]["operations"][operation]["mean_ms"] /
                  by_arm["chematic"][i]["operations"][operation]["mean_ms"]
                  for i in range(20)]
        ci = bootstrap([math.log(ratio) for ratio in ratios])
        lanes[operation] = {"paired_blocks": 20,
                            "median_speed_ratio_chematic_over_rdkit": quantile(ratios, 0.5),
                            "geometric_mean_speed_ratio": math.exp(statistics.fmean(math.log(r) for r in ratios)),
                            "bootstrap_95_ci": ci,
                            "output_gate": "browser_morgan_digest_exact" if operation in fingerprint_operations
                                           else "node_row_preflight_only" if operation == "parse"
                                           else "canonical_write_not_equivalent",
                            "practical_10pct_speed_win": ci[0] > 1.10 and operation in fingerprint_operations}
    memory = {arm: {"median_peak_process_tree_rss_bytes": quantile(
                    [run["process_tree_rss"]["peak_rss_bytes"] for run in by_arm[arm].values()], 0.5),
                    "median_init_ms": quantile([run["init_ms"] for run in by_arm[arm].values()], 0.5)}
              for arm in by_arm}
    result = {"schema": "published-browser-paired20-summary/v1",
            "source_sha256": sha256(source), "corpus_sha256": sha256(CORPUS),
            "artifacts": {name: artifacts[name]["sha256"] for name in EXPECTED},
            "protocol": {"engine": "chromium", "rows": 250, "paired_blocks": 20,
                         "processes_per_block": 2, "alternating_order": "AB then BA",
                         "warmup_rows": 20, "timing_method": "single_batch"},
            "lanes": lanes, "memory": memory,
            "limits": ["Single-host Chromium, not a cross-browser or cross-host speed claim",
                       "Canonical-write outputs are not gated for equality and cannot earn a win",
                       "Parse atom-count/Morgan row identity is from the separate Node preflight on the same published artifacts",
                       "Process-tree RSS can double-count shared pages and is not library-allocated memory"]}
    if prepared_mode == "split":
        result["protocol"]["prepared_mode"] = "split: first-use on fresh prepared objects after separate warmup objects, then reused objects"
        result["limits"].append("Prepared first-use and reused lanes time the same 250 objects in sequence; neither includes parsing/preparation")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = check(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({name: lane["bootstrap_95_ci"] for name, lane in result["lanes"].items()}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
