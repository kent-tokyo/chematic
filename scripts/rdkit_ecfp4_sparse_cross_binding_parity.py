#!/usr/bin/env python3
"""Measure RDKit-compatible ECFP4 sparse-count parity across bindings."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "chembl_accuracy_corpus_4999.smi"
REPORT = ROOT / "validation" / "results" / "rdkit-ecfp4-sparse-cross-binding-parity-5000-v1.0.9.json"


def jsonl(value: str) -> list[dict]:
    return [json.loads(line) for line in value.splitlines() if line.strip()]


def run(command: list[str], corpus: str) -> list[dict]:
    result = subprocess.run(command, cwd=ROOT, input=corpus, text=True, capture_output=True, check=True)
    return jsonl(result.stdout)


def python_records(corpus: str) -> list[dict]:
    import chematic

    rows = []
    for index, raw in enumerate(corpus.splitlines()):
        smiles = raw.strip()
        if not smiles: continue
        try:
            mol = chematic.from_smiles(smiles)
            _, sparse_counts, _, _ = mol.rdkit_ecfp4_detail()
            pairs = [[int(identifier), int(count)] for identifier, count in sorted(sparse_counts.items())]
            rows.append({"index": index, "smiles": smiles, "status": "ok", "sparse_counts": pairs})
        except Exception as error:
            rows.append({"index": index, "smiles": smiles, "status": "error", "error": str(error)})
    return rows


def compare(left: list[dict], right: list[dict]) -> dict:
    differences = []
    matches = 0
    for index, (a, b) in enumerate(zip(left, right, strict=True)):
        same = a["status"] == b["status"]
        if same and a["status"] == "ok": same = a["sparse_counts"] == b["sparse_counts"]
        if same: matches += 1
        elif len(differences) < 10: differences.append({"index": index, "left": a, "right": b})
    return {"rows": len(left), "matches": matches, "mismatches": len(left) - matches, "examples": differences}


def main() -> int:
    corpus = CORPUS.read_text(encoding="utf-8")
    expected_rows = sum(bool(line.strip()) for line in corpus.splitlines())
    bindings = {
        "rust": run(["cargo", "run", "-p", "chematic-chem", "--release", "--offline", "--example", "rdkit_ecfp4_sparse_binding_dump"], corpus),
        "python": python_records(corpus),
        "node_wasm": run(["node", "scripts/rdkit_ecfp4_sparse_binding_dump.mjs"], corpus),
    }
    if any(len(rows) != expected_rows for rows in bindings.values()): raise RuntimeError({name: len(rows) for name, rows in bindings.items()})
    pairwise = {"rust_vs_python": compare(bindings["rust"], bindings["python"]), "rust_vs_node_wasm": compare(bindings["rust"], bindings["node_wasm"]), "python_vs_node_wasm": compare(bindings["python"], bindings["node_wasm"])}
    report = {
        "schema_version": 1,
        "target_version": "1.0.9",
        "contract": "same chematic source SMILES, RDKit-compatible ECFP4 raw sparse identifier/count pairs",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "rows": expected_rows, "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "binding_status_counts": {name: {status: sum(row["status"] == status for row in rows) for status in ("ok", "error")} for name, rows in bindings.items()},
        "pairwise": pairwise,
        "gate_passed": all(pair["mismatches"] == 0 for pair in pairwise.values()),
        "boundary": ["This checks sparse-count binding parity for chematic's RDKit-compatible implementation, not an independent RDKit oracle run.", "Raw/folded bitInfo provenance remains a separate gate."],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "rows": expected_rows, "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__": sys.exit(main())
