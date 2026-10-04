#!/usr/bin/env python3
"""Verify the archived same-binary RDKit legacy-ring diagnostic."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

from scripts.summarize_rdkit_ring_backend import compare_triplet


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
NATIVE_PREFIX = "rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04"
RING_PREFIX = "rdkit-native-ring-backend-2026-09-1-legacy-2026-10-04"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rows(path: Path) -> tuple[bytes, list[dict]]:
    body = gzip.decompress(path.read_bytes())
    return body, [json.loads(line) for line in body.splitlines()]


def check() -> None:
    report = json.loads((RESULTS / f"{RING_PREFIX}-summary.json").read_text())
    native = json.loads((RESULTS / f"{NATIVE_PREFIX}-summary.json").read_text())
    config_path = ROOT / "validation" / "rdkit_native_rebaseline_build.json"
    config = json.loads(config_path.read_text())
    corpus_path = ROOT / config["corpus"]
    queries_path = ROOT / config["queries"]
    smiles = corpus_path.read_text().splitlines()
    queries = json.loads(queries_path.read_text())["queries"]
    require(report["schema_version"] == 1, "unexpected ring probe schema")
    require(
        report["profile"] == "rdkit_2026_09_1_legacy_ring_mode_probe_v1",
        "unexpected ring probe profile",
    )
    require(
        report["new_legacy_environment"] == {"RDK_USE_LEGACY_RING_FINDING": "1"},
        "legacy switch provenance changed",
    )
    for label in ("old", "new"):
        require(
            report[f"{label}_source_commit"] == config[f"{label}_source_commit"]
            == native["sources"][label]["commit"],
            f"{label} source commit changed",
        )
    for key, path in (
        ("build_config_sha256", config_path),
        ("corpus_sha256", corpus_path),
        ("queries_sha256", queries_path),
    ):
        require(report[key] == sha256(path.read_bytes()), f"{key} changed")
    require(len(smiles) == 10000 and len(queries) == 31, "input count changed")

    bodies = {}
    parsed = {}
    paths = {
        "old": RESULTS / f"{NATIVE_PREFIX}-old.jsonl.gz",
        "new": RESULTS / f"{NATIVE_PREFIX}-new.jsonl.gz",
        "new_with_legacy_ring_mode": RESULTS / f"{RING_PREFIX}-rows.jsonl.gz",
    }
    for label, path in paths.items():
        body, data = rows(path)
        require(sha256(body) == report["row_sha256"][label], f"{label} digest changed")
        require(len(data) == 10000, f"{label} rows incomplete")
        for index, row in enumerate(data):
            require(
                row["input_index"] == index
                and row["smiles"] == smiles[index]
                and row["status"] == "ok"
                and len(row["smarts"]) == 31,
                f"{label} row {index} incomplete or reordered",
            )
        bodies[label], parsed[label] = body, data
    require(
        bodies["old"] == bodies["new_with_legacy_ring_mode"],
        "legacy-mode rows do not reproduce old-version rows",
    )
    for label in ("old", "new"):
        require(
            report["row_sha256"][label]
            == native["sources"][label]["rows_sha256"],
            f"{label} native packet disagrees",
        )

    comparison = compare_triplet(
        parsed["old"],
        parsed["new"],
        parsed["new_with_legacy_ring_mode"],
        queries,
    )
    for key, value in comparison.items():
        require(report[key] == value, f"{key} differs from raw rows")
    require(report["old_vs_new"] == native["counts"], "old/new count changed")
    require(
        report["old_vs_new"]["changed_rows"] == 6
        and report["old_vs_new"]["smarts_changed_cells"] == 12
        and report["old_vs_new"]["smarts_cells_compared"] == 310000,
        "old/new difference scope changed",
    )
    require(
        report["old_vs_new_legacy"]["changed_rows"] == 0
        and report["old_vs_new_legacy"]["smarts_changed_cells"] == 0
        and report["legacy_reproduces_old_changed_cells"] == 12,
        "legacy mode does not explain all changed cells",
    )
    require(
        {item["input_index"] for item in report["oracle_changed_cells"]}
        == {9, 23, 28, 29, 30, 34}
        and {item["query"] for item in report["oracle_changed_cells"]}
        == {"[R2]", "[R3]"},
        "changed query cells are not the pinned ring-count cases",
    )
    print("RDKit ring backend evidence OK: 10k legacy rows reproduce old version; all 12 ring-query deltas explained")


if __name__ == "__main__":
    try:
        check()
    except (OSError, KeyError, IndexError, TypeError, ValueError) as exc:
        raise SystemExit(f"RDKit ring backend evidence invalid: {exc}") from exc
