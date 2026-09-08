#!/usr/bin/env python3
"""Measure the 5,000-row Rust/Python/Node-WASM descriptor binding contract.

This is a same-source binding agreement report, not an RDKit accuracy claim.
The corpus and descriptor semantics are pinned in the report so the result is
reproducible without changing the existing external held-out ledger.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "chembl_accuracy_corpus_4999.smi"
REPORT = ROOT / "validation" / "results" / "descriptor-cross-binding-parity-5000-v1.0.9.json"
FIELDS = ("mw", "tpsa", "hbd", "hba", "heavy_atoms")
FLOAT_FIELDS = {"mw", "tpsa"}
FLOAT_TOLERANCE = 1e-9


def read_records(lines: str) -> list[dict]:
    return [json.loads(line) for line in lines.splitlines() if line.strip()]


def run_jsonl(command: list[str], corpus: str) -> list[dict]:
    result = subprocess.run(command, cwd=ROOT, input=corpus, text=True, capture_output=True, check=True)
    return read_records(result.stdout)


def python_records(corpus: str) -> list[dict]:
    import chematic

    records = []
    for index, raw in enumerate(corpus.splitlines()):
        smiles = raw.strip()
        if not smiles:
            continue
        try:
            mol = chematic.from_smiles(smiles)
            records.append({
                "index": index,
                "smiles": smiles,
                "status": "ok",
                "descriptors": {
                    "mw": mol.mw,
                    "tpsa": mol.tpsa,
                    "hbd": mol.hbd,
                    "hba": mol.hba,
                    "heavy_atoms": mol.heavy_atoms,
                },
            })
        except Exception as error:  # Keep invalid-row accounting binding-stable.
            records.append({"index": index, "smiles": smiles, "status": "error", "error": str(error)})
    return records


def compare(left: list[dict], right: list[dict]) -> dict:
    mismatches = []
    matches = 0
    for lrow, rrow in zip(left, right, strict=True):
        ok = lrow.get("index") == rrow.get("index") and lrow.get("smiles") == rrow.get("smiles")
        ok = ok and lrow.get("status") == rrow.get("status")
        if ok and lrow["status"] == "ok":
            for field in FIELDS:
                lv, rv = lrow["descriptors"][field], rrow["descriptors"][field]
                if field in FLOAT_FIELDS:
                    ok = ok and abs(lv - rv) <= FLOAT_TOLERANCE
                else:
                    ok = ok and lv == rv
        if ok:
            matches += 1
        elif len(mismatches) < 10:
            mismatches.append({"index": lrow.get("index"), "left": lrow, "right": rrow})
    return {"rows": len(left), "matches": matches, "mismatches": len(left) - matches, "examples": mismatches}


def main() -> int:
    corpus = CORPUS.read_text(encoding="utf-8")
    corpus_lines = [line for line in corpus.splitlines() if line.strip()]
    corpus_sha256 = hashlib.sha256(CORPUS.read_bytes()).hexdigest()
    rust = run_jsonl(["cargo", "run", "-p", "chematic-chem", "--release", "--offline", "--example", "descriptor_binding_dump"], corpus)
    node = run_jsonl(["node", "scripts/descriptor_binding_dump.mjs"], corpus)
    python = python_records(corpus)
    bindings = {"rust": rust, "python": python, "node_wasm": node}
    expected_rows = len(corpus_lines)
    for name, records in bindings.items():
        if len(records) != expected_rows:
            raise RuntimeError(f"{name} emitted {len(records)} rows, expected {expected_rows}")

    report = {
        "schema_version": 1,
        "target_version": "1.0.9",
        "contract": "same chematic source SMILES, Rust/PyO3/Node-WASM binding agreement",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "rows": expected_rows, "sha256": corpus_sha256},
        "descriptor_fields": list(FIELDS),
        "float_tolerance": FLOAT_TOLERANCE,
        "binding_status_counts": {
            name: {status: sum(row["status"] == status for row in records) for status in ("ok", "error")}
            for name, records in bindings.items()
        },
        "pairwise": {
            "rust_vs_python": compare(rust, python),
            "rust_vs_node_wasm": compare(rust, node),
            "python_vs_node_wasm": compare(python, node),
        },
        "gate_passed": all(
            pair["mismatches"] == 0
            for pair in (
                compare(rust, python),
                compare(rust, node),
                compare(python, node),
            )
        ),
        "boundary": [
            "This is binding parity, not an RDKit accuracy or chemical-validity oracle.",
            "The existing v1.0.8 held-out parity ledger remains unchanged.",
            "Full held-out parity for fingerprints, topology, torsion, and standardization remains open.",
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "gate_passed": report["gate_passed"], "rows": expected_rows}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
