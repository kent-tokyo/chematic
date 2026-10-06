#!/usr/bin/env python3
"""Gate ``find_matches_rdkit_parity(profile="2026.09.1")`` on RDKit 2026.09.1.

The oracle is the version-pinned native C++ RDKit 2026.09.1 build's SMARTS
cells for the exposed 10k corpus and the pinned 31-query matrix
(``rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-new.jsonl.gz``;
RDKit 2026.09.1 has no Python wheel). It differs from 2026.03.6 in 12
``[R2]``/``[R3]`` cells of six bis-quinolinium macrocycles, which the
2026.03.6 profile refuses as ``ring_model_ambiguous``. The 2026.09.1 profile
counts relevant cycles and must give every one of the 310,000 cells exactly,
with no refusal.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import chematic

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
QUERIES = ROOT / "validation/rdkit_rebaseline_smarts_queries.json"
RESULTS = ROOT / "validation/results"
ORACLE = RESULTS / "rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-new.jsonl.gz"
SUMMARY = RESULTS / "rdkit-native-rebaseline-2026-03-6-to-2026-09-1-2026-10-04-summary.json"
CHANGED_ROWS = (9, 23, 28, 29, 30, 34)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def match_sets(value: object) -> frozenset[tuple[int, ...]]:
    return frozenset(tuple(sorted(match)) for match in value)  # type: ignore[union-attr]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scope", default="source-built Python wheel; not a published registry artifact")
    args = parser.parse_args()
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    if (summary["sources"]["new"]["rows_sha256"] != hashlib.sha256(gzip.open(ORACLE).read()).hexdigest()
            or summary["queries_sha256"] != sha256_text(QUERIES)
            or summary["corpus_sha256"] != sha256_text(CORPUS)):
        raise ValueError("RDKit 2026.09.1 oracle provenance changed")
    queries = json.loads(QUERIES.read_text(encoding="utf-8"))["queries"]
    corpus = CORPUS.read_text(encoding="utf-8").splitlines()
    counts: Counter[str] = Counter()
    unexpected = []
    with gzip.open(ORACLE, "rt", encoding="utf-8") as stream:
        for index, smiles in enumerate(corpus):
            row = json.loads(next(stream))
            if row["input_index"] != index or row["smiles"] != smiles or len(row["smarts"]) != len(queries):
                raise ValueError(f"oracle row {index} is not aligned")
            mol = chematic.from_smiles(smiles)
            for position, query in enumerate(queries):
                expected = match_sets(row["smarts"][position])
                for profile in ("2026.09.1", "2026.03.6"):
                    result = mol.find_matches_rdkit_parity(query, profile=profile)
                    if result["status"] == "ok":
                        ok = match_sets(result["matches"]) == expected
                        counts[f"{profile}:exact" if ok else f"{profile}:wrong"] += 1
                        if not ok and profile == "2026.09.1":
                            unexpected.append([index, query, "wrong_confident"])
                    else:
                        counts[f"{profile}:{result['status']}:{result['reason']}"] += 1
                        if profile == "2026.09.1":
                            unexpected.append([index, query, result["status"], result["reason"]])
        if next(stream, None) is not None:
            raise ValueError("oracle has extra rows")
    report = {
        "schema": "python-opt-in-smarts-parity-310k-rdkit-2026.09.1/v1",
        "scope": args.scope,
        "chematic_version": chematic.__version__,
        "oracle": {"path": str(ORACLE.relative_to(ROOT)), "sha256": sha256(ORACLE),
                   "rdkit": "2026.09.1 native C++ build", "commit": summary["sources"]["new"]["commit"]},
        "input_cells": len(corpus) * len(queries),
        "counts": dict(sorted(counts.items())),
        "unexpected": unexpected[:50],
        "unexpected_total": len(unexpected),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"counts": report["counts"], "unexpected": len(unexpected)}))
    if unexpected or counts["2026.09.1:exact"] != len(corpus) * len(queries):
        raise ValueError("2026.09.1 profile did not give every cell exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
