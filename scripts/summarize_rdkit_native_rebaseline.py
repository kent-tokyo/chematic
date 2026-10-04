#!/usr/bin/env python3
"""Check complete RDKit C++ source-build rows and summarize old/new deltas."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path, smiles: list[str], query_count: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            row = json.loads(line)
            if index >= len(smiles) or row.get("input_index") != index:
                raise ValueError(f"{path}: missing, excess, or reordered row {index}")
            if row.get("smiles") != smiles[index]:
                raise ValueError(f"{path}: SMILES mismatch at row {index}")
            if row.get("status") != "ok":
                raise ValueError(
                    f"{path}: failed native row {index}: {row.get('status')}"
                )
            if not isinstance(row.get("canonical"), str):
                raise ValueError(f"{path}: missing canonical at row {index}")
            if not isinstance(row.get("cip_atoms"), dict) or not isinstance(
                row.get("cip_bonds"), dict
            ):
                raise ValueError(f"{path}: missing CIP maps at row {index}")
            bits = row.get("morgan_on_bits")
            if (
                not isinstance(bits, list)
                or any(
                    not isinstance(bit, int) or bit < 0 or bit >= 2048 for bit in bits
                )
                or bits != sorted(set(bits))
            ):
                raise ValueError(f"{path}: invalid Morgan bit set at row {index}")
            matches = row.get("smarts")
            if not isinstance(matches, list) or len(matches) != query_count:
                raise ValueError(f"{path}: missing SMARTS cells at row {index}")
            rows.append(row)
    if len(rows) != len(smiles):
        raise ValueError(f"{path}: expected {len(smiles)} rows, got {len(rows)}")
    return rows


def morgan_binary_sha256(on_bits: list[int]) -> str:
    """Hash RDKit's 2048-bit little-endian BitVectToBinaryText representation."""
    data = bytearray(2048 // 8)
    for bit in on_bits:
        data[bit // 8] |= 1 << (bit % 8)
    return hashlib.sha256(data).hexdigest()


def verify_python_baseline(path: Path, rows: list[dict[str, Any]]) -> int:
    """Cross-check the old C++ oracle against the pinned old Python wheel oracle."""
    count = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index >= len(rows):
                raise ValueError(f"{path}: excess Python baseline row {index}")
            reference = json.loads(line)
            native = rows[index]
            checks = {
                "input_index": (reference.get("input_index"), index),
                "smiles": (reference.get("smiles"), native["smiles"]),
                "status": (reference.get("status"), "completed"),
                "canonical": (
                    reference.get("smiles_parse_write", {}).get("rdkit_canonical"),
                    native["canonical"],
                ),
                "cip_atoms": (
                    reference.get("cip", {}).get("rdkit_atoms"),
                    native["cip_atoms"],
                ),
                "cip_bonds": (
                    reference.get("cip", {}).get("rdkit_bonds"),
                    native["cip_bonds"],
                ),
                "morgan_sha256": (
                    reference.get("morgan", {}).get("rdkit_sha256"),
                    morgan_binary_sha256(native["morgan_on_bits"]),
                ),
            }
            for field, (expected, actual) in checks.items():
                if expected != actual:
                    raise ValueError(
                        f"{path}: row {index} {field} differs from old C++ oracle"
                    )
            count += 1
    if count != len(rows):
        raise ValueError(f"{path}: expected {len(rows)} Python rows, got {count}")
    return count


def summarize(
    old_rows: list[dict[str, Any]], new_rows: list[dict[str, Any]], queries: list[str]
) -> tuple[dict[str, int], list[dict[str, Any]]]:
    if len(old_rows) != len(new_rows):
        raise ValueError("old/new row counts differ")
    counts: Counter[str] = Counter()
    deltas: list[dict[str, Any]] = []
    fields = ("canonical", "cip_atoms", "cip_bonds", "morgan_on_bits")
    for old, new in zip(old_rows, new_rows, strict=True):
        if old["input_index"] != new["input_index"] or old["smiles"] != new["smiles"]:
            raise ValueError("old/new input order differs")
        differences: list[dict[str, Any]] = []
        for field in fields:
            if old[field] != new[field]:
                counts[f"{field}_changed_rows"] += 1
                differences.append(
                    {"operation": field, "old": old[field], "new": new[field]}
                )
        for query_index, query in enumerate(queries):
            counts["smarts_cells_compared"] += 1
            if old["smarts"][query_index] != new["smarts"][query_index]:
                counts["smarts_changed_cells"] += 1
                differences.append(
                    {
                        "operation": "smarts",
                        "query_index": query_index,
                        "query": query,
                        "old": old["smarts"][query_index],
                        "new": new["smarts"][query_index],
                    }
                )
        counts["rows"] += 1
        if differences:
            counts["changed_rows"] += 1
            deltas.append(
                {
                    "input_index": old["input_index"],
                    "smiles": old["smiles"],
                    "differences": differences,
                }
            )
    for key in (
        "changed_rows",
        "smarts_changed_cells",
        *(f"{f}_changed_rows" for f in fields),
    ):
        counts.setdefault(key, 0)
    return dict(counts), deltas


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "build-config",
        "old-source",
        "new-source",
        "old-binary",
        "new-binary",
        "old-rows",
        "new-rows",
        "python-baseline-rows",
        "summary",
        "delta",
    ):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.build_config.read_text(encoding="utf-8"))
    if config.get("schema_version") != 1:
        raise ValueError("unsupported build config")
    root = args.build_config.resolve().parents[1]
    corpus = root / config["corpus"]
    queries_path = root / config["queries"]
    smiles = corpus.read_text(encoding="utf-8").splitlines()
    queries = json.loads(queries_path.read_text(encoding="utf-8"))["queries"]
    if (
        len(smiles) != config["expected_rows"]
        or len(queries) != config["expected_queries"]
    ):
        raise ValueError("corpus or query count changed")
    sources = {}
    for label in ("old", "new"):
        source = getattr(args, f"{label}_source")
        commit = subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
        ).strip()
        if commit != config[f"{label}_source_commit"]:
            raise ValueError(f"{label} source commit changed: {commit}")
        sources[label] = {
            "commit": commit,
            "binary_sha256": digest(getattr(args, f"{label}_binary")),
            "rows_sha256": digest(getattr(args, f"{label}_rows")),
        }
    old_rows = read_rows(args.old_rows, smiles, len(queries))
    new_rows = read_rows(args.new_rows, smiles, len(queries))
    python_baseline_verified_rows = verify_python_baseline(
        args.python_baseline_rows, old_rows
    )
    counts, deltas = summarize(old_rows, new_rows, queries)
    args.delta.parent.mkdir(parents=True, exist_ok=True)
    with args.delta.open("w", encoding="utf-8") as handle:
        for delta in deltas:
            handle.write(
                json.dumps(delta, sort_keys=True, separators=(",", ":")) + "\n"
            )
    result = {
        "schema_version": 1,
        "profile": config["profile"],
        "interpretation": "version-pinned independent C++ source builds; not distributed binaries, Python wrappers, or a speed comparison",
        "build_config_sha256": digest(args.build_config),
        "runner_source_sha256": digest(root / "tools/rdkit_native_rebaseline/main.cpp"),
        "cmake_flags": config["cmake_flags"],
        "corpus_sha256": digest(corpus),
        "queries_sha256": digest(queries_path),
        "python_baseline_rows_sha256": digest(args.python_baseline_rows),
        "python_baseline_verified_rows": python_baseline_verified_rows,
        "sources": sources,
        "toolchain": {
            "platform": platform.platform(),
            "compiler": subprocess.check_output(
                ["c++", "--version"], text=True
            ).splitlines()[0],
            "cmake": subprocess.check_output(
                ["cmake", "--version"], text=True
            ).splitlines()[0],
        },
        "counts": counts,
        "delta_sha256": digest(args.delta),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(counts, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
