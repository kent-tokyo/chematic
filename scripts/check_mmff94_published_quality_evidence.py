#!/usr/bin/env python3
"""Check the two v1.0.26 published-wheel MMFF94 quality repetitions (#637)."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

from summarize_public_package_3d import paired_speedup, quality_summary, timing_summary


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation/results/mmff94-public-v1026-20260926"
CHEMATIC_SHA256 = "6ebe5337b1f52d529b3d956e0822e695c1983369fb97c694cdd4c32cd73a629f"
RDKIT_SHA256 = "e16c467cb254a223e59a0cf81358c6b39da15a99d2909170d693e95778fddb41"
ARMS = {
    "chematic": "chematic_pipeline_v2_mmff94_strict_stereo_safe",
    "rdkit": "rdkit_etkdgv3_mmff94",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def expected_inputs() -> list[tuple[str, str, str]]:
    result: list[tuple[str, str, str]] = []
    for tier, filename in (
        ("A", "pipeline_v2_vs_rdkit_etkdgv3_tier_a.json"),
        ("B", "pipeline_v2_vs_rdkit_etkdgv3_tier_b.json"),
    ):
        manifest = load_json(ROOT / "validation/manifests" / filename)
        result.extend(
            (tier, row["name"], row["smiles"]) for row in manifest["molecules"]
        )
    require(len(result) == 265, "manifest no longer contains 265 declared inputs")
    return result


def validate_run(
    suffix: str, inputs: list[tuple[str, str, str]]
) -> tuple[datetime, datetime]:
    files = {
        engine: {
            "rows": RESULTS / f"{engine}-3d{suffix}.jsonl",
            "meta": RESULTS / f"{engine}-3d{suffix}.meta.json",
        }
        for engine in ARMS
    }
    rows_by_engine: dict[str, list[dict]] = {}
    times: dict[str, datetime] = {}
    for engine, paths in files.items():
        rows, meta = load_jsonl(paths["rows"]), load_json(paths["meta"])
        expected_version = "1.0.26" if engine == "chematic" else "2026.03.6"
        expected_sha = CHEMATIC_SHA256 if engine == "chematic" else RDKIT_SHA256
        require(
            meta["package"]["runtime_version"] == expected_version,
            f"{engine}: runtime version",
        )
        require(
            meta["package"]["wheel"]["sha256"] == expected_sha, f"{engine}: wheel hash"
        )
        require(
            meta["result"]["output_sha256"] == digest(paths["rows"]),
            f"{engine}: row hash",
        )
        require(
            meta["result"]["selected_molecules"] == 265
            and meta["result"]["emitted_rows"] == 265
            and meta["result"]["termination"] == "completed",
            f"{engine}: incomplete runner accounting",
        )
        require(meta["configuration"]["random_seed"] == 20260801, f"{engine}: seed")
        require(meta["configuration"]["arms"] == [ARMS[engine]], f"{engine}: arm")
        require(len(rows) == 265, f"{engine}: expected 265 rows")
        require(
            [(row["tier"], row["name"], row["smiles"]) for row in rows] == inputs
            and [row["row_index"] for row in rows] == list(range(265)),
            f"{engine}: input row identity/order changed",
        )
        statuses = Counter(row["status"] for row in rows)
        expected = (
            {"success": 265}
            if engine == "chematic"
            else {"success": 264, "internal_error": 1}
        )
        require(dict(statuses) == expected, f"{engine}: terminal statuses changed")
        rows_by_engine[engine] = rows
        times[engine] = datetime.fromisoformat(meta["collected_at_utc"])

    scored_path = RESULTS / f"common-scored{suffix}.jsonl"
    scored = load_jsonl(scored_path)
    scored_keys = [
        (row["engine"], row["tier"], row["name"])
        for row in scored
        if row["status"] == "scored"
    ]
    successful_keys = [
        (engine, row["tier"], row["name"])
        for engine, rows in rows_by_engine.items()
        for row in rows
        if row["status"] == "success"
    ]
    require(
        len(scored_keys) == 529 and len(set(scored_keys)) == 529,
        "scorer: duplicate/missing rows",
    )
    require(
        set(scored_keys) == set(successful_keys),
        "scorer: successful inputs not accounted for",
    )
    require(
        all(
            row["independently_sound"] is True
            and row["gross_clash_count"] == 0
            and row["stereo"]["violated"] == 0
            for row in scored
            if row["status"] == "scored"
        ),
        "scorer: unsound, stereo-violating, or clashing success",
    )
    require(
        Counter(row["status"] for row in scored) == {"scored": 529, "paired_rmsd": 264},
        "scorer: unexpected status or paired-row count",
    )

    summary = load_json(RESULTS / f"summary{suffix}.json")
    for engine, paths in files.items():
        require(
            summary["inputs"][engine]["sha256"] == digest(paths["rows"]),
            f"summary: {engine} hash",
        )
    require(
        summary["inputs"]["scored"]["sha256"] == digest(scored_path),
        "summary: scorer hash",
    )
    pair = summary["arm_pairs"]["mmff94_stereo_safe"]
    for engine, rows in rows_by_engine.items():
        require(
            pair["timing"][engine] == timing_summary(rows),
            f"summary: {engine} timing does not match raw rows",
        )
        require(
            pair["quality"][engine]
            == quality_summary(scored, engine, ARMS[engine], 265),
            f"summary: {engine} quality does not match scored rows",
        )
    require(
        pair["timing"]["paired"]
        == paired_speedup(rows_by_engine["chematic"], rows_by_engine["rdkit"]),
        "summary: paired timing does not match raw rows",
    )
    require(
        pair["quality"]["chematic"]["usable"] == 265, "summary: CheMatic quality count"
    )
    require(pair["quality"]["rdkit"]["usable"] == 264, "summary: RDKit quality count")
    require(
        pair["timing"]["paired"]["common_successes"] == 264, "summary: pairing count"
    )
    return times["chematic"], times["rdkit"]


def main() -> int:
    inputs = expected_inputs()
    first_c, first_r = validate_run("", inputs)
    second_c, second_r = validate_run("-r2", inputs)
    require(
        first_c < first_r < second_r < second_c,
        "two runs were not recorded in reversed order",
    )
    print(
        "MMFF94 published-quality evidence OK: 2 x 265 complete rows, 265/265 vs 264/265 usable"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
