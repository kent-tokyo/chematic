#!/usr/bin/env python3
"""Measure RDKit-compatible RDKFingerprint parity across Rust, Python, and WASM."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "chembl_accuracy_corpus_4999.smi"
REPORT = ROOT / "validation" / "results" / "rdkit-rdk-cross-binding-parity-5000-v1.0.9.json"


def jsonl(value: str) -> list[dict]:
    return [json.loads(line) for line in value.splitlines() if line.strip()]


def run(command: list[str], corpus: str) -> list[dict]:
    result = subprocess.run(command, cwd=ROOT, input=corpus, text=True, capture_output=True, check=True)
    return jsonl(result.stdout)


def rust_records(corpus: str) -> list[dict]:
    result = subprocess.run(
        [
            "cargo", "run", "-p", "chematic-cli", "--release", "--offline", "--",
            "batch-fingerprints", "--input", str(CORPUS), "--algorithm", "rdkit_rdk",
        ], cwd=ROOT, text=True, capture_output=True, check=True,
    )
    batch = json.loads(result.stdout)
    rows = []
    for index, record in enumerate(batch["records"]):
        if record.get("error") is not None:
            rows.append({"index": index, "smiles": record["input_smiles"], "status": "error"})
            continue
        bits = record["fingerprint"]["set_bits"]
        payload = bytearray(256)
        for bit in bits:
            payload[bit // 8] |= 1 << (bit % 8)
        rows.append({"index": index, "smiles": record["input_smiles"], "status": "ok", "rdk_hex": payload.hex()})
    return rows


def python_records(corpus: str) -> list[dict]:
    import chematic

    rows = []
    for index, raw in enumerate(corpus.splitlines()):
        smiles = raw.strip()
        if not smiles:
            continue
        try:
            mol = chematic.from_smiles(smiles)
            rows.append({"index": index, "smiles": smiles, "status": "ok", "rdk_hex": mol.rdkit_rdk_fp().hex()})
        except Exception as error:
            rows.append({"index": index, "smiles": smiles, "status": "error", "error": str(error)})
    return rows


def compare(left: list[dict], right: list[dict]) -> dict:
    differences = []
    matches = 0
    for a, b in zip(left, right, strict=True):
        same = a["status"] == b["status"] and (a["status"] != "ok" or a["rdk_hex"] == b["rdk_hex"])
        if same:
            matches += 1
        elif len(differences) < 10:
            differences.append({"index": a["index"], "left": a, "right": b})
    return {"rows": len(left), "matches": matches, "mismatches": len(left) - matches, "examples": differences}


def main() -> int:
    corpus = CORPUS.read_text(encoding="utf-8")
    expected_rows = sum(bool(line.strip()) for line in corpus.splitlines())
    bindings = {
        "rust": rust_records(corpus),
        "python": python_records(corpus),
        "node_wasm": run(["node", "scripts/rdkit_rdk_binding_dump.mjs"], corpus),
    }
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
        "contract": "same chematic source SMILES, RDKit-compatible RDKFingerprint 2048-bit LSB-first bytes",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "rows": expected_rows, "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "binding_status_counts": {name: {status: sum(row["status"] == status for row in rows) for status in ("ok", "error")} for name, rows in bindings.items()},
        "pairwise": pairwise,
        "gate_passed": all(pair["mismatches"] == 0 for pair in pairwise.values()),
        "boundary": ["This is cross-binding parity for chematic's RDKit-compatible RDKFingerprint implementation, not an independent RDKit oracle run."],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "rows": expected_rows, "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
