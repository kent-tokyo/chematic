#!/usr/bin/env python3
"""Compare RDKit-compatible top-k search across Rust, Python, and Node/WASM."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path
from validation_provenance import collect

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "descriptor_census_corpus.smi"
REPORT = ROOT / "validation" / "results" / "rdkit-search-cross-binding-parity-v1.0.13.json"
CASE_COUNT = 500
LIBRARY_SIZE = 4_500
CHUNK_SIZE = 1_000
K_VALUES = (1, 10, 100)


def load_jsonl(output: str) -> list[dict]:
    return [json.loads(line) for line in output.splitlines() if line.strip()]


def run(command: list[str], cases: list[dict]) -> list[dict]:
    result = subprocess.run(
        command,
        cwd=ROOT,
        input="".join(json.dumps(case) + "\n" for case in cases),
        text=True,
        capture_output=True,
        check=True,
    )
    return load_jsonl(result.stdout)


def python_rows(cases: list[dict]) -> list[dict]:
    import chematic

    rows = []
    indexes = {}
    for case in cases:
        key = tuple(case["db"])
        index = indexes.get(key)
        if index is None:
            index = chematic.PreparedFingerprintIndex.from_smiles(case["db"], fp="rdkit_ecfp4")
            indexes[key] = index
        hits = index.search(case["query"], k=case["k"])
        rows.append({"results": [{"index": i, "tanimoto": score} for i, score in hits]})
    return rows


def compare(left: list[dict], right: list[dict]) -> dict:
    mismatches = [
        {"case": i, "left": a, "right": b}
        for i, (a, b) in enumerate(zip(left, right, strict=True))
        if a != b
    ]
    return {
        "cases": len(left),
        "matches": len(left) - len(mismatches),
        "mismatches": len(mismatches),
        "examples": mismatches[:3],
    }


def merge_chunks(rows: list[dict], chunk_count: int, case_count: int, k: int) -> list[dict]:
    """Merge local chunk hits into global top-k results."""
    merged = []
    for query_index in range(case_count):
        hits = []
        for chunk_index in range(chunk_count):
            row = rows[chunk_index * CASE_COUNT + query_index]
            hits.extend(
                (item["index"] + chunk_index * CHUNK_SIZE, item["tanimoto"])
                for item in row["results"]
            )
        hits.sort(key=lambda item: (-item[1], item[0]))
        merged.append({"results": [{"index": i, "tanimoto": score} for i, score in hits[:k]]})
    return merged


def main() -> int:
    smiles = [line.strip() for line in CORPUS.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(smiles) < LIBRARY_SIZE + CASE_COUNT:
        raise SystemExit("corpus is smaller than the full library/query split")
    chunks = [
        smiles[start : min(start + CHUNK_SIZE, LIBRARY_SIZE)]
        for start in range(0, LIBRARY_SIZE, CHUNK_SIZE)
    ]
    pairwise = {}
    for k in K_VALUES:
        cases = [
            {"query": smiles[LIBRARY_SIZE + query], "db": chunk, "k": k}
            for chunk in chunks
            for query in range(CASE_COUNT)
        ]
        rows = {
            "rust": run(["cargo", "run", "-p", "chematic-fp", "--release", "--offline", "--example", "rdkit_search_binding_dump"], cases),
            "python": python_rows(cases),
            "node_wasm": run(["node", "scripts/binding_dump.mjs", "rdkit-search"], cases),
        }
        merged = {
            name: merge_chunks(value, len(chunks), CASE_COUNT, k) for name, value in rows.items()
        }
        pairwise[str(k)] = {
            "rust_vs_python": compare(merged["rust"], merged["python"]),
            "rust_vs_node_wasm": compare(merged["rust"], merged["node_wasm"]),
            "python_vs_node_wasm": compare(merged["python"], merged["node_wasm"]),
        }
    report = {
        "schema_version": 2,
        "target_version": "1.0.13",
        "contract": "same RDKit-compatible Morgan profile, database indices, full-precision scores, and tie ordering",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "configuration": {
            "queries": CASE_COUNT,
            "library_rows": LIBRARY_SIZE,
            "chunk_size": CHUNK_SIZE,
            "chunks": len(chunks),
            "k_values": list(K_VALUES),
            "score_comparison": "JSON f64 values; no display rounding",
        },
        "provenance": collect(ROOT),
        "pairwise": pairwise,
        "gate_passed": all(
            comparison["mismatches"] == 0
            for by_binding in pairwise.values()
            for comparison in by_binding.values()
        ),
        "boundary": [
            "This is binding parity for the RDKit-compatible chematic profile, not an independent RDKit accuracy claim.",
            "The independent RDKit oracle gate remains the separate benchmark evidence.",
        ],
    }
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=REPORT.parent, prefix=f".{REPORT.name}.", delete=False
    ) as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
        temporary_path = handle.name
    os.replace(temporary_path, REPORT)
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
