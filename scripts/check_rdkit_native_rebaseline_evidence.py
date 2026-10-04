#!/usr/bin/env python3
"""Independently verify the pinned RDKit C++ old/new source-build packet."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

from summarize_rdkit_native_rebaseline import summarize, verify_python_baseline


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"
PREFIX = "rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04"
NPM_DELTA = RESULTS / (
    "rdkit-rebaseline-npm-oracle-delta-2026.03.6-to-2026.09.1-2026-10-04.jsonl.gz"
)
PYTHON_BASELINE = RESULTS / (
    "rdkit-rebaseline-python-chemistry-v1.0.33-vs-2026.03.6-2026-10-04.jsonl.gz"
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_gzip_jsonl(path: Path) -> tuple[list[dict], str]:
    body = gzip.decompress(path.read_bytes())
    return [json.loads(line) for line in body.splitlines()], digest(body)


def normalize_atom_sets(value: list) -> list[tuple[int, ...]]:
    return sorted(
        tuple(json.loads(item) if isinstance(item, str) else item) for item in value
    )


def check() -> None:
    summary = json.loads((RESULTS / f"{PREFIX}-summary.json").read_text())
    build = ROOT / "validation" / "rdkit_native_rebaseline_build.json"
    corpus = ROOT / "validation" / "benchmark_corpora" / "rdkit-js-browser-10k-v1.smi"
    queries_path = ROOT / "validation" / "rdkit_rebaseline_smarts_queries.json"
    queries = json.loads(queries_path.read_text())["queries"]
    smiles = corpus.read_text().splitlines()
    require(summary["schema_version"] == 1, "native summary schema changed")
    require(
        summary["profile"] == "rdkit_native_2026_03_6_to_2026_09_1_exposed_10k_v1",
        "native profile changed",
    )
    require(
        summary["build_config_sha256"] == digest(build.read_bytes()),
        "build packet changed",
    )
    require(
        summary["runner_source_sha256"]
        == digest(
            (ROOT / "tools" / "rdkit_native_rebaseline" / "main.cpp").read_bytes()
        ),
        "native runner source changed",
    )
    require(summary["corpus_sha256"] == digest(corpus.read_bytes()), "corpus changed")
    require(
        summary["queries_sha256"] == digest(queries_path.read_bytes()),
        "queries changed",
    )
    require(len(smiles) == 10000 and len(queries) == 31, "comparison inputs incomplete")
    require(
        summary["python_baseline_rows_sha256"] == digest(PYTHON_BASELINE.read_bytes())
        and summary["python_baseline_verified_rows"] == 10000,
        "old Python cross-check provenance changed",
    )
    config = json.loads(build.read_text())
    require(summary["cmake_flags"] == config["cmake_flags"], "build flags changed")
    for label in ("old", "new"):
        require(
            summary["sources"][label]["commit"] == config[f"{label}_source_commit"],
            f"{label} source commit changed",
        )
        require(
            len(summary["sources"][label]["binary_sha256"]) == 64,
            f"{label} compiled probe hash missing",
        )
    for name in ("compiler", "cmake", "platform"):
        require(bool(summary["toolchain"].get(name)), f"{name} provenance missing")

    rows = {}
    for label in ("old", "new"):
        path = RESULTS / f"{PREFIX}-{label}.jsonl.gz"
        data, uncompressed_hash = load_gzip_jsonl(path)
        require(
            uncompressed_hash == summary["sources"][label]["rows_sha256"],
            f"{label} raw rows hash changed",
        )
        require(len(data) == 10000, f"{label} rows incomplete")
        for index, row in enumerate(data):
            require(
                row["input_index"] == index and row["smiles"] == smiles[index],
                f"{label} row {index} input correspondence changed",
            )
            require(row["status"] == "ok", f"{label} row {index} not comparable")
            require(len(row["smarts"]) == 31, f"{label} row {index} SMARTS incomplete")
        rows[label] = data
    require(
        verify_python_baseline(PYTHON_BASELINE, rows["old"]) == 10000,
        "old Python oracle differs",
    )

    counts, derived_deltas = summarize(rows["old"], rows["new"], queries)
    require(counts == summary["counts"], "native summary counts disagree with rows")
    require(
        (
            counts["rows"],
            counts["smarts_cells_compared"],
            counts["changed_rows"],
            counts["smarts_changed_cells"],
        )
        == (10000, 310000, 6, 12),
        "native old/new difference scope changed",
    )
    for field in (
        "canonical_changed_rows",
        "cip_atoms_changed_rows",
        "cip_bonds_changed_rows",
        "morgan_on_bits_changed_rows",
        "noncomparable_rows",
        "parse_status_changed_rows",
        "smarts_old_parse_failure_cells",
        "smarts_new_parse_failure_cells",
    ):
        require(counts[field] == 0, f"native {field} is nonzero")

    native_deltas, delta_hash = load_gzip_jsonl(RESULTS / f"{PREFIX}-delta.jsonl.gz")
    require(delta_hash == summary["delta_sha256"], "native delta hash changed")
    require(native_deltas == derived_deltas, "native delta file differs from raw rows")
    npm_rows, _ = load_gzip_jsonl(NPM_DELTA)
    npm_deltas = {row["input_index"]: row for row in npm_rows if row["differences"]}
    require(len(npm_rows) == 10000, "npm oracle row packet incomplete")
    require(
        {row["input_index"] for row in native_deltas} == npm_deltas.keys(),
        "native and npm oracle changed-row sets differ",
    )
    compared = 0
    for native in native_deltas:
        npm = npm_deltas[native["input_index"]]
        require(native["smiles"] == npm["smiles"], "native/npm input differs")
        by_key = {
            (item["operation"], item.get("query_index")): item
            for item in npm["differences"]
        }
        require(
            len(native["differences"]) == len(by_key), "native/npm delta count differs"
        )
        for item in native["differences"]:
            key = (item["operation"], item.get("query_index"))
            require(
                key in by_key,
                f"native/npm query differs at row {native['input_index']}",
            )
            other = by_key[key]
            require(item["query"] == other["query"], "native/npm SMARTS text differs")
            for version in ("old", "new"):
                require(
                    normalize_atom_sets(item[version])
                    == normalize_atom_sets(other[version]),
                    f"native/npm {version} SMARTS atom sets differ",
                )
            compared += 1
    require(compared == 12, "native/npm changed-cell comparison incomplete")
    print(
        "RDKit native source rebaseline OK: 10k old/new rows, old Python baseline, 6 rows/12 SMARTS cells identical to npm delta"
    )


if __name__ == "__main__":
    try:
        check()
    except (OSError, KeyError, IndexError, TypeError, ValueError) as exc:
        raise SystemExit(f"RDKit native rebaseline evidence invalid: {exc}") from exc
