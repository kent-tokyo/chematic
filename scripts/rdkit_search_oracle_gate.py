#!/usr/bin/env python3
"""Compare the search result against an independent RDKit Morgan oracle."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from validation_provenance import collect

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPUS = ROOT / "scripts" / "descriptor_census_corpus.smi"
DEFAULT_OUTPUT = ROOT / "validation" / "results" / "rdkit-search-oracle-v1.0.13.json"
LIBRARY_SIZE = 4_500
QUERY_COUNT = 500
CHUNK_SIZE = 1_000
K_VALUES = (1, 10, 100)
SCORE_TOLERANCE = 1e-12


def run_rust(cases: list[dict]) -> list[dict]:
    result = subprocess.run(
        [
            "cargo", "run", "-p", "chematic-fp", "--release", "--offline",
            "--example", "rdkit_search_binding_dump",
        ],
        cwd=ROOT,
        input="".join(json.dumps(case) + "\n" for case in cases),
        text=True,
        capture_output=True,
        check=True,
    )
    return [json.loads(line) for line in result.stdout.splitlines() if line.strip()]


def rdkit_fingerprints(smiles: list[str]):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
        includeChirality=False,
        useBondTypes=True,
        includeRedundantEnvironments=False,
    )
    return [generator.GetFingerprint(Chem.MolFromSmiles(smi)) for smi in smiles]


def oracle_rows(cases: list[dict], fps: dict[str, object]) -> list[dict]:
    rows = []
    for case in cases:
        query_fp = fps[case["query"]]
        scored = [
            (DataStructs.TanimotoSimilarity(query_fp, fps[smi]), index)
            for index, smi in enumerate(case["db"])
        ]
        scored.sort(key=lambda item: (-item[0], item[1]))
        rows.append({
            "results": [
                {"index": index, "tanimoto": score}
                for score, index in scored[: case["k"]]
            ]
        })
    return rows


def compare(actual: list[dict], expected: list[dict]) -> dict:
    mismatches = []
    max_score_delta = 0.0
    for case_index, (got, want) in enumerate(zip(actual, expected, strict=True)):
        got_results = got["results"]
        want_results = want["results"]
        if len(got_results) != len(want_results):
            mismatches.append({"case": case_index, "reason": "length", "actual": got, "oracle": want})
            continue
        for rank, (got_hit, want_hit) in enumerate(zip(got_results, want_results, strict=True)):
            delta = abs(float(got_hit["tanimoto"]) - float(want_hit["tanimoto"]))
            max_score_delta = max(max_score_delta, delta)
            if got_hit["index"] != want_hit["index"] or delta > SCORE_TOLERANCE:
                mismatches.append({
                    "case": case_index,
                    "rank": rank,
                    "actual": got_hit,
                    "oracle": want_hit,
                    "score_delta": delta,
                })
                break
    return {
        "cases": len(actual),
        "matches": len(actual) - len(mismatches),
        "mismatches": len(mismatches),
        "max_score_delta": max_score_delta,
        "score_tolerance": SCORE_TOLERANCE,
        "examples": mismatches[:5],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.corpus = args.corpus.resolve()
    args.output = args.output.resolve()

    smiles = [line.strip() for line in args.corpus.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(smiles) < LIBRARY_SIZE + QUERY_COUNT:
        raise SystemExit("corpus is smaller than the declared library/query split")
    library = smiles[:LIBRARY_SIZE]
    queries = smiles[LIBRARY_SIZE : LIBRARY_SIZE + QUERY_COUNT]
    unique_smiles = sorted(set(library + queries))
    molecules = {smi: Chem.MolFromSmiles(smi) for smi in unique_smiles}
    if any(mol is None for mol in molecules.values()):
        raise SystemExit("RDKit failed to parse a corpus molecule")
    fps = dict(zip(unique_smiles, rdkit_fingerprints(unique_smiles), strict=True))
    chunks = [library[start : start + CHUNK_SIZE] for start in range(0, len(library), CHUNK_SIZE)]

    by_k = {}
    for k in K_VALUES:
        cases = [
            {"query": query, "db": chunk, "k": k}
            for chunk in chunks
            for query in queries
        ]
        actual = run_rust(cases)
        expected = oracle_rows(cases, fps)
        by_k[str(k)] = compare(actual, expected)

    report = {
        "schema_version": 1,
        "profile": "rdkit_morgan_ecfp4_search_oracle",
        "chematic_version": "1.0.13",
        "rdkit_version": rdBase.rdkitVersion,
        "corpus": {
            "path": str(args.corpus.relative_to(ROOT)),
            "sha256": hashlib.sha256(args.corpus.read_bytes()).hexdigest(),
            "library_rows": len(library),
            "query_rows": len(queries),
            "chunk_size": CHUNK_SIZE,
        },
        "configuration": {
            "radius": 2,
            "fp_size": 2048,
            "include_chirality": False,
            "use_bond_types": True,
            "k_values": list(K_VALUES),
        },
        "provenance": collect(ROOT),
        "comparison": by_k,
        "gate_passed": all(item["mismatches"] == 0 for item in by_k.values()),
        "boundary": [
            "RDKit computes the independent fingerprint and exhaustive similarity ordering.",
            "This does not establish accuracy on an unused corpus; the corpus is a declared regression split.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=args.output.parent, prefix=f".{args.output.name}.", delete=False
    ) as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
        temporary_path = handle.name
    os.replace(temporary_path, args.output)
    try:
        report_path = args.output.relative_to(ROOT)
    except ValueError:
        report_path = args.output
    print(json.dumps({"report": str(report_path), "gate_passed": report["gate_passed"]}))
    return 0 if report["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
