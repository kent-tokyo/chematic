#!/usr/bin/env python3
"""Measure native MACCS 166-bit byte parity across Rust, Python, and WASM."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "chembl_accuracy_corpus_4999.smi"
REPORT = ROOT / "validation" / "results" / "maccs-cross-binding-parity-5000-v1.0.9.json"


def parse_jsonl(value: str) -> list[dict]:
    return [json.loads(line) for line in value.splitlines() if line.strip()]


def run(command: list[str], corpus: str) -> list[dict]:
    return parse_jsonl(subprocess.run(command, cwd=ROOT, input=corpus, text=True, capture_output=True, check=True).stdout)


def python_records(corpus: str) -> list[dict]:
    import chematic

    rows = []
    for index, raw in enumerate(corpus.splitlines()):
        smiles = raw.strip()
        if not smiles:
            continue
        try:
            mol = chematic.from_smiles(smiles)
            rows.append({"index": index, "smiles": smiles, "status": "ok", "maccs_hex": mol.maccs().hex()})
        except Exception as error:
            rows.append({"index": index, "smiles": smiles, "status": "error", "error": str(error)})
    return rows


def compare(left: list[dict], right: list[dict]) -> dict:
    differences = [i for i, (a, b) in enumerate(zip(left, right, strict=True)) if a != b]
    return {
        "rows": len(left),
        "matches": len(left) - len(differences),
        "mismatches": len(differences),
        "examples": [{"index": i, "left": left[i], "right": right[i]} for i in differences[:10]],
    }


def main() -> int:
    corpus = CORPUS.read_text(encoding="utf-8")
    expected_rows = sum(bool(line.strip()) for line in corpus.splitlines())
    # MACCS evaluates a relatively large SMARTS key set. Run the independent
    # bindings concurrently so this lane measures parity without serializing
    # three identical O(N * key-set) traversals.
    jobs = {
        "rust": (["cargo", "run", "-p", "chematic-chem", "--release", "--offline", "--example", "maccs_binding_dump"],),
        "node_wasm": (["node", "scripts/maccs_binding_dump.mjs"],),
    }
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            name: executor.submit(run, command[0], corpus) for name, command in jobs.items()
        }
        futures["python"] = executor.submit(python_records, corpus)
        bindings = {name: future.result() for name, future in futures.items()}
    if any(len(rows) != expected_rows for rows in bindings.values()):
        raise RuntimeError({name: len(rows) for name, rows in bindings.items()})
    pairwise = {
        "rust_vs_python": compare(bindings["rust"], bindings["python"]),
        "rust_vs_node_wasm": compare(bindings["rust"], bindings["node_wasm"]),
        "python_vs_node_wasm": compare(bindings["python"], bindings["node_wasm"]),
    }
    report = {
        "schema_version": 1,
        "target_version": "1.0.9",
        "contract": "same chematic source SMILES, native MACCS 166-bit 21-byte LSB-first representation",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "rows": expected_rows, "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "binding_status_counts": {name: {status: sum(row["status"] == status for row in rows) for status in ("ok", "error")} for name, rows in bindings.items()},
        "pairwise": pairwise,
        "gate_passed": all(pair["mismatches"] == 0 for pair in pairwise.values()),
        "boundary": ["This is native chematic binding parity, not an RDKit accuracy claim.", "RDKit MACCS key mapping remains a separate external comparison."],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "rows": expected_rows, "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
