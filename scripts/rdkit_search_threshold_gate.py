#!/usr/bin/env python3
"""Verify inclusive RDKit-compatible search thresholds across all bindings."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from validation_provenance import collect

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "scripts" / "descriptor_census_corpus.smi"
REPORT = ROOT / "validation" / "results" / "rdkit-search-threshold-gate-v1.0.13.json"
LIBRARY_SIZE = 32
QUERY_COUNT = 8
K = 128
THRESHOLD_EPSILON = 1e-12


def load_rows(output: str) -> list[dict]:
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
    return load_rows(result.stdout)


def rdkit_reference(query: str, db: list[str], threshold: float) -> list[dict]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    query_fp = generator.GetFingerprint(Chem.MolFromSmiles(query))
    scored = []
    for index, smiles in enumerate(db):
        score = DataStructs.TanimotoSimilarity(
            query_fp, generator.GetFingerprint(Chem.MolFromSmiles(smiles))
        )
        if score >= threshold:
            scored.append({"index": index, "tanimoto": score})
    scored.sort(key=lambda row: (-row["tanimoto"], row["index"]))
    return scored[:K]


def normalize(row: dict) -> list[dict]:
    return row["results"]


def main() -> int:
    smiles = [line.strip() for line in CORPUS.read_text(encoding="utf-8").splitlines() if line.strip()]
    db = smiles[:LIBRARY_SIZE]
    queries = smiles[LIBRARY_SIZE : LIBRARY_SIZE + QUERY_COUNT]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    cases = []
    expected = []
    labels = []
    for query in queries:
        query_fp = generator.GetFingerprint(Chem.MolFromSmiles(query))
        scores = sorted(
            {
                DataStructs.TanimotoSimilarity(
                    query_fp, generator.GetFingerprint(Chem.MolFromSmiles(smiles))
                )
                for smiles in db
            },
            reverse=True,
        )
        thresholds = [0.0, 1.0, 0.5]
        for score in scores[:3]:
            thresholds.extend([max(0.0, score - THRESHOLD_EPSILON), score, min(1.0, score + THRESHOLD_EPSILON)])
        for threshold in dict.fromkeys(thresholds):
            case = {"query": query, "db": db, "k": K, "threshold": threshold}
            cases.append(case)
            expected.append(rdkit_reference(query, db, threshold))
            labels.append({"threshold": threshold, "query": query})

    bindings = {
        "rust": run(
            ["cargo", "run", "-p", "chematic-fp", "--release", "--offline", "--example", "rdkit_search_binding_dump"],
            cases,
        ),
        "node_wasm": run(["node", "scripts/rdkit_search_binding_dump.mjs"], cases),
    }
    import chematic

    index = chematic.PreparedFingerprintIndex.from_smiles(db, fp="rdkit_ecfp4")
    bindings["python"] = [
        {"results": [{"index": i, "tanimoto": score} for i, score in index.search_threshold(case["query"], case["threshold"], case["k"])]}
        for case in cases
    ]

    comparisons = {}
    for name, rows in bindings.items():
        mismatches = []
        for index, (actual, want) in enumerate(zip(rows, expected, strict=True)):
            if normalize(actual) != want:
                mismatches.append({"case": index, "label": labels[index], "actual": normalize(actual), "expected": want})
        comparisons[name] = {"cases": len(cases), "matches": len(cases) - len(mismatches), "mismatches": len(mismatches), "examples": mismatches[:3]}

    report = {
        "schema_version": 1,
        "target_version": "1.0.13",
        "contract": "inclusive score >= threshold, descending score then original index, full-precision scores",
        "corpus": {"path": str(CORPUS.relative_to(ROOT)), "sha256": hashlib.sha256(CORPUS.read_bytes()).hexdigest()},
        "configuration": {"library_rows": LIBRARY_SIZE, "query_rows": QUERY_COUNT, "cases": len(cases), "k": K, "epsilon": THRESHOLD_EPSILON},
        "provenance": collect(ROOT),
        "comparisons": comparisons,
        "gate_passed": all(item["mismatches"] == 0 for item in comparisons.values()),
        "boundary": ["RDKit 2025.09.3 oracle in this environment", "zero-score candidates are included at threshold 0.0", "invalid thresholds are covered by binding/unit tests"],
    }
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=REPORT.parent, prefix=f".{REPORT.name}.", delete=False
    ) as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
        temporary_path = handle.name
    os.replace(temporary_path, REPORT)
    print(json.dumps({"report": str(REPORT.relative_to(ROOT)), "cases": len(cases), "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
