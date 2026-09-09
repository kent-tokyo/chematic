#!/usr/bin/env python3
"""Measure native ECFP4 byte parity across Rust, Python, and Node/WASM."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "chembl_accuracy_corpus_4999.smi"
REPORT = ROOT / "validation" / "results" / "ecfp4-cross-binding-parity-5000-v1.0.9.json"


def records(output: str) -> list[dict]:
    return [json.loads(line) for line in output.splitlines() if line.strip()]


def run(command: list[str], corpus: str) -> list[dict]:
    result = subprocess.run(command, cwd=ROOT, input=corpus, text=True, capture_output=True, check=True)
    return records(result.stdout)


def python_records(corpus: str) -> list[dict]:
    import chematic

    result = []
    for index, raw in enumerate(corpus.splitlines()):
        smiles = raw.strip()
        if not smiles:
            continue
        try:
            mol = chematic.from_smiles(smiles)
            result.append({
                "index": index,
                "smiles": smiles,
                "status": "ok",
                "ecfp4_hex": mol.ecfp4().hex(),
            })
        except Exception as error:
            result.append({"index": index, "smiles": smiles, "status": "error", "error": str(error)})
    return result


def compare(left: list[dict], right: list[dict]) -> dict:
    mismatches = []
    for lrow, rrow in zip(left, right, strict=True):
        if lrow == rrow:
            continue
        if len(mismatches) < 10:
            mismatches.append({"index": lrow.get("index"), "left": lrow, "right": rrow})
    return {"rows": len(left), "matches": len(left) - len(mismatches) if not mismatches else sum(a == b for a, b in zip(left, right)), "mismatches": sum(a != b for a, b in zip(left, right)), "examples": mismatches}


def main() -> int:
    corpus = CORPUS.read_text(encoding="utf-8")
    expected_rows = sum(bool(line.strip()) for line in corpus.splitlines())
    bindings = {
        "rust": run(["cargo", "run", "-p", "chematic-chem", "--release", "--offline", "--example", "ecfp4_binding_dump"], corpus),
        "python": python_records(corpus),
        "node_wasm": run(["node", "scripts/ecfp4_binding_dump.mjs"], corpus),
    }
    for name, rows in bindings.items():
        if len(rows) != expected_rows:
            raise RuntimeError(f"{name} emitted {len(rows)} rows, expected {expected_rows}")
    pairwise = {
        "rust_vs_python": compare(bindings["rust"], bindings["python"]),
        "rust_vs_node_wasm": compare(bindings["rust"], bindings["node_wasm"]),
        "python_vs_node_wasm": compare(bindings["python"], bindings["node_wasm"]),
    }
    report = {
        "schema_version": 1,
        "target_version": "1.0.9",
        "contract": "same chematic source SMILES, native ECFP4 2048-bit LSB-first bytes",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "rows": expected_rows, "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "binding_status_counts": {name: {status: sum(row["status"] == status for row in rows) for status in ("ok", "error")} for name, rows in bindings.items()},
        "pairwise": pairwise,
        "gate_passed": all(pair["mismatches"] == 0 for pair in pairwise.values()),
        "boundary": [
            "This is native chematic binding parity, not an RDKit fingerprint-accuracy claim.",
            "RDKit-compatible ECFP4 and sparse/explanation parity remain separate gates.",
            "The existing held-out parity ledger remains unchanged.",
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "rows": expected_rows, "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
