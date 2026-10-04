"""Compare pinned RDKit old/new rows with the new binary's legacy ring mode."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from scripts.summarize_rdkit_native_rebaseline import digest, read_rows, summarize


def compare_triplet(
    old_rows: list[dict[str, Any]],
    new_rows: list[dict[str, Any]],
    legacy_rows: list[dict[str, Any]],
    queries: list[str],
) -> dict[str, Any]:
    old_new_counts, old_new_deltas = summarize(old_rows, new_rows, queries)
    new_legacy_counts, _ = summarize(new_rows, legacy_rows, queries)
    old_legacy_counts, _ = summarize(old_rows, legacy_rows, queries)
    changed_cells: list[dict[str, Any]] = []
    for delta in old_new_deltas:
        row_index = delta["input_index"]
        for difference in delta["differences"]:
            if difference["operation"] != "smarts":
                continue
            query_index = difference["query_index"]
            legacy = legacy_rows[row_index]
            legacy_set = (
                legacy["smarts"][query_index] if legacy["status"] == "ok" else None
            )
            changed_cells.append(
                {
                    "input_index": row_index,
                    "query_index": query_index,
                    "query": queries[query_index],
                    "old": difference["old"],
                    "new": difference["new"],
                    "new_with_legacy_ring_mode": legacy_set,
                    "legacy_matches_old": legacy["status"] == "ok"
                    and legacy_set == difference["old"],
                }
            )
    return {
        "old_vs_new": old_new_counts,
        "new_vs_new_legacy": new_legacy_counts,
        "old_vs_new_legacy": old_legacy_counts,
        "oracle_changed_cells": changed_cells,
        "legacy_reproduces_old_changed_cells": sum(
            cell["legacy_matches_old"] for cell in changed_cells
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("build-config", "old-rows", "new-rows", "legacy-rows", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.build_config.read_text(encoding="utf-8"))
    root = args.build_config.resolve().parents[1]
    corpus = root / config["corpus"]
    queries_path = root / config["queries"]
    smiles = corpus.read_text(encoding="utf-8").splitlines()
    queries = json.loads(queries_path.read_text(encoding="utf-8"))["queries"]
    if (
        len(smiles) != config["expected_rows"]
        or len(queries) != config["expected_queries"]
    ):
        raise ValueError("pinned corpus or queries changed")
    old = read_rows(args.old_rows, smiles, len(queries))
    new = read_rows(args.new_rows, smiles, len(queries))
    legacy = read_rows(args.legacy_rows, smiles, len(queries))
    comparison = compare_triplet(old, new, legacy, queries)
    changed = comparison["oracle_changed_cells"]
    if (
        len(changed) != 12
        or {cell["input_index"] for cell in changed} != {9, 23, 28, 29, 30, 34}
        or {cell["query"] for cell in changed} != {"[R2]", "[R3]"}
    ):
        raise ValueError("old/new ring-query delta changed from pinned 12 cells")
    report = {
        "schema_version": 1,
        "profile": "rdkit_2026_09_1_legacy_ring_mode_probe_v1",
        "old_source_commit": config["old_source_commit"],
        "new_source_commit": config["new_source_commit"],
        "new_legacy_environment": {"RDK_USE_LEGACY_RING_FINDING": "1"},
        "build_config_sha256": digest(args.build_config),
        "corpus_sha256": digest(corpus),
        "queries_sha256": digest(queries_path),
        "row_sha256": {
            "old": digest(args.old_rows),
            "new": digest(args.new_rows),
            "new_with_legacy_ring_mode": digest(args.legacy_rows),
        },
        **comparison,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        "Legacy ring probe:",
        report["legacy_reproduces_old_changed_cells"],
        "/ 12 changed cells match the old version",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
