#!/usr/bin/env python3
"""Verify the pinned three-browser Ubuntu published-artifact replication run."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

if __package__:
    from .check_published_browser_paired import check as check_speed
    from .check_v1030_published_browser_morgan_rows import check as check_rows
else:
    from check_published_browser_paired import check as check_speed
    from check_v1030_published_browser_morgan_rows import check as check_rows

ROOT = Path(__file__).resolve().parents[1]
RECORD_HASHES = {
    "published-browser-chromium.json": "a8d92068a1f22c986df95e20063e0a1c34b9856201f6fb63aef4df34fb8c92cd",
    "published-browser-chromium-summary.json": "c3cdb31481cd47ca43a81691f2dc26c02b55ca25c34e6e447f8bb24b10209815",
    "published-browser-chromium-row-direct.json": "69c670a5266d4f07e7ae3bf761da6cbef6dcbcfd55140e91632099f31f5a9124",
    "published-browser-chromium-row-prepared.json": "bc7dafe2dc0598c172687e1bc97768f4d237ab1edefc5602c1dafffafdf7e8d9",
    "published-browser-firefox.json": "9631752721d678cbecc42a7e6077de1aa0634f31e4d2b41411f74c22bcdc1115",
    "published-browser-firefox-summary.json": "0322f2c454c47d92a83a84d13a7f94f0bf8da45e4eef6eadb5b9b77217c8f3d4",
    "published-browser-firefox-row-direct.json": "3e09497571787c6681a65c8483fae130aecd8f9ec3d50b539237f51a09ed31c1",
    "published-browser-firefox-row-prepared.json": "a267ec1e47f39f59db3f54e977a7e80d8682a0ab9e3e650b3fc2d4569259623b",
    "published-browser-webkit.json": "249ec3fca431a60b8b61e183eb2459e193b7465670eb71b700cea1c70f6f03f8",
    "published-browser-webkit-summary.json": "567d6d2fbed78ca8b9a59daf0fef6aeae68b93a9b5057c9c7817e2a3741fc35b",
    "published-browser-webkit-row-direct.json": "5f12d3f810b64b17595c33d7f25730f43818f6642edef3b9307073511a19d579",
    "published-browser-webkit-row-prepared.json": "387346af953fd5983a56b12c96e5c542e3da443da39ba16b7dc6a134234446c6",
}


def check() -> dict:
    for name, digest in RECORD_HASHES.items():
        path = ROOT / "benchmarks" / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Linux browser record hash mismatch: {name}")
    results = {}
    for engine in ("chromium", "firefox", "webkit"):
        raw = ROOT / "benchmarks" / f"published-browser-{engine}.json"
        summary_path = ROOT / "benchmarks" / f"published-browser-{engine}-summary.json"
        summary = check_speed(raw)
        if (summary != json.loads(summary_path.read_bytes())
                or summary["protocol"]["engine"] != engine
                or summary["protocol"]["process_rss"] != "not_measured"
                or not summary["protocol"]["prepared_mode"].startswith("split:")):
            raise ValueError(f"Linux {engine}: timing/summary protocol mismatch")
        for operation in ("direct", "prepared"):
            path = ROOT / "benchmarks" / f"published-browser-{engine}-row-{operation}.json"
            if check_rows(path, operation, 250, engine)["exact"] != 250:
                raise ValueError(f"Linux {engine}/{operation}: row-level Morgan mismatch")
        results[engine] = {
            lane: {"median_ratio": summary["lanes"][lane]["median_speed_ratio_chematic_over_rdkit"],
                   "paired_95_ci": summary["lanes"][lane]["bootstrap_95_ci"]}
            for lane in ("parse_fp", "prepared_fp_first_use", "prepared_fp_reused")
        }
    return results


if __name__ == "__main__":
    print(json.dumps(check(), indent=2, sort_keys=True))
