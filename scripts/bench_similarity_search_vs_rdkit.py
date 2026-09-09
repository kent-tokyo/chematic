#!/usr/bin/env python3
"""Benchmark exact top-k similarity search against RDKit on one corpus.

The two lanes use the same SMILES split, query order, radius, bit width, k,
warm-up policy, and timing loop. Fingerprint definitions are still named
separately: chematic's native ECFP4 and RDKit Morgan are not bit-identical.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CORPUS = ROOT / "scripts" / "descriptor_census_corpus.smi"


def percentile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * p
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    fraction = position - low
    return ordered[low] + (ordered[high] - ordered[low]) * fraction


def summary(samples_us: list[float]) -> dict[str, float]:
    return {
        "p50_us": round(percentile(samples_us, 0.50), 3),
        "p95_us": round(percentile(samples_us, 0.95), 3),
        "mean_us": round(statistics.mean(samples_us), 3),
        "min_us": round(min(samples_us), 3),
        "max_us": round(max(samples_us), 3),
    }


def time_queries(search, queries, warmup: int) -> tuple[dict[str, float], list[list[int]]]:
    for query in queries[:warmup]:
        search(query)
    samples = []
    result_ids = []
    for query in queries:
        started = time.perf_counter_ns()
        result = search(query)
        samples.append((time.perf_counter_ns() - started) / 1_000)
        result_ids.append([int(item[0]) for item in result])
    return summary(samples), result_ids


def overlap(left: list[list[int]], right: list[list[int]]) -> dict[str, float]:
    recalls = []
    jaccards = []
    for a, b in zip(left, right):
        aset, bset = set(a), set(b)
        recalls.append(len(aset & bset) / len(bset) if bset else 1.0)
        union = aset | bset
        jaccards.append(len(aset & bset) / len(union) if union else 1.0)
    return {
        "mean_recall_of_rdkit_top_k": round(statistics.mean(recalls), 6),
        "p05_recall_of_rdkit_top_k": round(percentile(recalls, 0.05), 6),
        "mean_top_k_jaccard": round(statistics.mean(jaccards), 6),
        "queries": len(recalls),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--library", type=int, default=4500)
    parser.add_argument("--queries", type=int, default=500)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--warmup", type=int, default=10)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    corpus = args.corpus if args.corpus.is_absolute() else ROOT / args.corpus
    output = args.output if args.output.is_absolute() else ROOT / args.output
    smiles = [line.strip() for line in corpus.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(smiles) < args.library + args.queries:
        raise SystemExit("corpus is smaller than library + query split")
    library = smiles[: args.library]
    queries = smiles[args.library : args.library + args.queries]

    import chematic
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    started = time.perf_counter_ns()
    chematic_index = chematic.PreparedFingerprintIndex.from_smiles(library, fp="ecfp4")
    chematic_build_ms = (time.perf_counter_ns() - started) / 1_000_000
    started = time.perf_counter_ns()
    rdkit_molecules = [Chem.MolFromSmiles(value) for value in library]
    if any(molecule is None for molecule in rdkit_molecules):
        raise SystemExit("RDKit rejected an input in the fixed library split")
    rdkit_fps = [generator.GetFingerprint(molecule) for molecule in rdkit_molecules]
    rdkit_build_ms = (time.perf_counter_ns() - started) / 1_000_000

    def chematic_search(query: str):
        return chematic_index.search(query, k=args.k)

    def rdkit_search(query_smiles: str):
        query_molecule = Chem.MolFromSmiles(query_smiles)
        if query_molecule is None:
            raise ValueError("RDKit rejected a query in the fixed corpus split")
        similarities = DataStructs.BulkTanimotoSimilarity(
            generator.GetFingerprint(query_molecule), rdkit_fps
        )
        return sorted(enumerate(similarities), key=lambda item: (-item[1], item[0]))[: args.k]

    chematic_timing, chematic_ids = time_queries(chematic_search, queries, args.warmup)
    rdkit_timing, rdkit_ids = time_queries(rdkit_search, queries, args.warmup)

    result = {
        "schema_version": 1,
        "target_version": getattr(chematic, "__version__", "unknown"),
        "corpus": {
            "path": str(corpus.relative_to(ROOT)) if corpus.is_relative_to(ROOT) else str(corpus),
            "sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(),
            "total_rows": len(smiles),
            "library_rows": args.library,
            "query_rows": args.queries,
            "split": "library=first rows; queries=following rows",
        },
        "configuration": {
            "fingerprint": "ECFP4/Morgan radius=2, nBits=2048",
            "k": args.k,
            "warmup_queries": args.warmup,
            "ordering": "similarity descending, original index ascending on ties",
        },
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "chematic": getattr(chematic, "__version__", "unknown"),
            "rdkit": rdBase.rdkitVersion,
        },
        "lanes": {
            "chematic_prepared_ecfp4": {"build_ms": round(chematic_build_ms, 3), "query": chematic_timing},
            "rdkit_morgan": {"build_ms": round(rdkit_build_ms, 3), "query": rdkit_timing},
        },
        "ranking_overlap": overlap(chematic_ids, rdkit_ids),
        "interpretation": "Speed lanes are implementation-specific because native fingerprint definitions differ; overlap is a retrieval-agreement metric, not ground truth.",
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
