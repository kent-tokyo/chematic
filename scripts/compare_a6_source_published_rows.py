#!/usr/bin/env python3
"""Compare fixed A6 source and published row outcomes on the same host."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def load_lane(rows_path: Path, metadata_path: Path, kind: str) -> tuple[list[dict], dict]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    configuration = metadata["configuration"]
    if (
        metadata["package"]["kind"] != kind
        or metadata["result"]["termination"] != "completed"
        or metadata["result"]["emitted_rows"] != 265
        or configuration["tiers"] != "AB"
        or configuration["arms"] != ["chematic_pipeline_v2_mmff94_strict_stereo_safe"]
        or configuration["force_field_max_iterations"] != 200
        or configuration["random_seed"] != 20260801
    ):
        raise ValueError(f"{kind} lane is not the fixed 265-row A6 run")
    digest = hashlib.sha256(rows_path.read_bytes()).hexdigest()
    if digest != metadata["result"]["output_sha256"]:
        raise ValueError(f"{kind} row SHA-256 differs from runner metadata")
    rows = [json.loads(line) for line in rows_path.read_text(encoding="utf-8").splitlines()]
    if len(rows) != 265 or [row.get("row_index") for row in rows] != list(range(265)):
        raise ValueError(f"{kind} lane has missing, duplicate or reordered rows")
    return rows, metadata


def compare(source: list[dict], published: list[dict]) -> dict:
    if len(source) != len(published):
        raise ValueError("source and published lanes have different row counts")
    lost = []
    gained = []
    for source_row, published_row in zip(source, published):
        if (
            source_row["row_index"] != published_row["row_index"]
            or source_row["name"] != published_row["name"]
            or source_row["smiles"] != published_row["smiles"]
            or source_row["arm"] != published_row["arm"]
        ):
            raise ValueError(f"corpus mismatch at row {source_row['row_index']}")
        source_ok = source_row["status"] == "success"
        published_ok = published_row["status"] == "success"
        if published_ok and not source_ok:
            lost.append(
                {
                    "row_index": source_row["row_index"],
                    "source_status": source_row["status"],
                    "source_cause": source_row.get("failure_cause"),
                }
            )
        elif source_ok and not published_ok:
            gained.append(published_row["row_index"])
    return {
        "published_statuses": dict(sorted(Counter(r["status"] for r in published).items())),
        "source_statuses": dict(sorted(Counter(r["status"] for r in source).items())),
        "lost_published_successes": lost,
        "gained_successes": gained,
        "source_preserves_published_successes": not lost,
        "full_a6_quality_gate": "not_run",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-rows", required=True, type=Path)
    parser.add_argument("--source-metadata", required=True, type=Path)
    parser.add_argument("--published-rows", required=True, type=Path)
    parser.add_argument("--published-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    source, source_meta = load_lane(args.source_rows, args.source_metadata, "source_candidate")
    published, published_meta = load_lane(
        args.published_rows, args.published_metadata, "published"
    )
    report = {
        "schema": "a6-source-vs-published-status/v1",
        "source_wheel_sha256": source_meta["package"]["wheel"]["sha256"],
        "published_wheel_sha256": published_meta["package"]["wheel"]["sha256"],
        "source_rows_sha256": source_meta["result"]["output_sha256"],
        "published_rows_sha256": published_meta["result"]["output_sha256"],
        **compare(source, published),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                key: report[key]
                for key in (
                    "published_statuses",
                    "source_statuses",
                    "lost_published_successes",
                    "source_preserves_published_successes",
                )
            },
            sort_keys=True,
        )
    )
    return 0  # Diagnoses status retention, never substitutes for full A6 scoring.


if __name__ == "__main__":
    raise SystemExit(main())
