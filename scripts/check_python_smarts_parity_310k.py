#!/usr/bin/env python3
"""Gate the opt-in Python SMARTS profile on the pinned 310,000-cell oracle."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import chematic
from rdkit import rdBase

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
QUERIES = ROOT / "validation/rdkit_rebaseline_smarts_queries.json"
ORACLE = ROOT / "validation/results/v1.0.30-rdkit-smarts-all-cells.jsonl.gz"
CLASSIFICATION = ROOT / "validation/results/rdkit-rebaseline-residual-classification-v1.0.27-issue634-vs-2026.03.6-2026-09-28.json"
PUBLISHED = ROOT / "validation/results/v1.0.30-published-python-chemistry-rows.jsonl.gz"
REFUSED_ROWS = (9, 23, 28, 29, 30, 34)
REFUSED_QUERIES = ("[R1]", "[R2]", "[R3]")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def match_sets(value: object) -> frozenset[tuple[int, ...]]:
    if not isinstance(value, list) or any(
        not isinstance(match, list) or any(not isinstance(atom, int) for atom in match)
        for match in value
    ):
        raise ValueError("match sets must be lists of atom-index lists")
    return frozenset(tuple(sorted(match)) for match in value)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--module-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != "2026.03.6":
        parser.error(f"expected RDKit 2026.03.6, got {rdBase.rdkitVersion}")
    module_path = Path(chematic.__file__).resolve()
    if not module_path.is_relative_to(args.module_root.resolve()):
        parser.error(f"imported {module_path}, not the isolated wheel")
    query_doc = json.loads(QUERIES.read_text(encoding="utf-8"))
    queries = query_doc["queries"]
    if len(queries) != 31 or len(set(queries)) != 31:
        raise ValueError("query list is not the pinned 31-query matrix")
    corpus = CORPUS.read_text(encoding="utf-8").splitlines()
    if len(corpus) != 10_000:
        raise ValueError("corpus row count changed")
    classification = json.loads(CLASSIFICATION.read_text(encoding="utf-8"))
    originally_wrong = {
        (row["input_index"], cell["query"])
        for row in classification["smarts"]["rows"] for cell in row["cells"]
    }
    if len(originally_wrong) != 200 or classification["smarts"]["differing_cells"] != 200:
        raise ValueError("published baseline no longer pins 200 residual cells")
    published_boolean = set()
    published_count = 0
    with gzip.open(PUBLISHED, "rt", encoding="utf-8") as stream:
        for index, line in enumerate(stream):
            published_count += 1
            row = json.loads(line)
            if row["input_index"] != index:
                raise ValueError("published row order changed")
            for diff in row["smarts"]["differences"]:
                if bool(diff["rdkit"]) != bool(diff["chematic"]):
                    published_boolean.add((index, diff["query"]))
    if (published_count != 10_000 or len(published_boolean) != 43
            or not published_boolean <= originally_wrong):
        raise ValueError("published Boolean residual accounting changed")
    expected_refused = {(i, q) for i in REFUSED_ROWS for q in REFUSED_QUERIES}
    counters: Counter[str] = Counter()
    refused = []
    unexpected = []
    with gzip.open(ORACLE, "rt", encoding="utf-8") as stream:
        manifest = json.loads(next(stream))
        if (manifest.get("rdkit_version") != rdBase.rdkitVersion
                or manifest.get("corpus_sha256") != sha256(CORPUS)
                or manifest.get("queries_sha256") != sha256(QUERIES)
                or manifest.get("input_count") != len(corpus)
                or manifest.get("query_count") != len(queries)):
            raise ValueError("RDKit oracle provenance changed")
        for index, smiles in enumerate(corpus):
            oracle = json.loads(next(stream))
            if (oracle["input_index"] != index or oracle["smiles"] != smiles
                    or len(oracle["matches"]) != len(queries)):
                raise ValueError(f"oracle row {index} is not aligned")
            mol = chematic.from_smiles(smiles)
            for position, query in enumerate(queries):
                cell = (index, query)
                expected = match_sets(oracle["matches"][position])
                result = mol.find_matches_rdkit_parity(query)
                status = result["status"]
                if status == "ok":
                    if result["reason"] is not None:
                        raise ValueError(f"successful cell has a reason: {cell}")
                    actual = match_sets(result["matches"])
                    if actual != expected:
                        unexpected.append([index, query, "wrong_confident"])
                    else:
                        counters["exact"] += 1
                        counters["fixed_original_residual"] += cell in originally_wrong
                        counters["fixed_original_boolean"] += cell in published_boolean
                elif status in {"typed_unsupported", "typed_refusal"}:
                    if result["matches"] is not None or not isinstance(result["reason"], str):
                        raise ValueError(f"refusal is not typed: {cell}")
                    refused.append([index, query, status, result["reason"]])
                    counters["refused_original_residual"] += cell in originally_wrong
                    counters["refused_original_correct"] += cell not in originally_wrong
                else:
                    raise ValueError(f"unknown result status {status!r} at {cell}")
        if next(stream, None) is not None:
            raise ValueError("RDKit oracle has extra rows")
    observed_refused = {(i, q) for i, q, _, _ in refused}
    report = {
        "schema": "python-opt-in-smarts-parity-310k/v1",
        "scope": "source-built Python wheel; not a published registry artifact",
        "rdkit_version": rdBase.rdkitVersion,
        "wheel_sha256": sha256(args.wheel),
        "module_path_within_wheel": str(module_path.relative_to(args.module_root.resolve())),
        "inputs": {name: sha256(path) for name, path in {
            "corpus": CORPUS, "queries": QUERIES, "oracle": ORACLE,
            "classification": CLASSIFICATION, "published_baseline": PUBLISHED,
        }.items()},
        "counts": dict(sorted(counters.items())),
        "input_cells": len(corpus) * len(queries),
        "refused": refused,
        "unexpected": unexpected[:50],
        "unexpected_total": len(unexpected),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if (unexpected or observed_refused != expected_refused
            or any(status != "typed_unsupported" or reason != "ring_model_ambiguous"
                   for _, _, status, reason in refused)
            or counters["exact"] != 309_982
            or counters["fixed_original_residual"] != 183
            or counters["fixed_original_boolean"] != 43
            or counters["refused_original_residual"] != 17
            or counters["refused_original_correct"] != 1):
        raise ValueError(f"source-wheel SMARTS gate did not meet the predeclared profile: {report['counts']}; "
                         f"refused={len(refused)}, unexpected={len(unexpected)}")
    print(json.dumps({"input_cells": report["input_cells"], "counts": report["counts"],
                      "refused": len(refused), "unexpected": len(unexpected)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
