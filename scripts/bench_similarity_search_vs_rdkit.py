#!/usr/bin/env python3
"""Benchmark exact top-k search in separate fingerprint-compatibility lanes."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CORPUS = ROOT / "scripts" / "descriptor_census_corpus.smi"


def percentile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * p
    low = int(position)
    high = min(low + 1, len(ordered) - 1)
    return ordered[low] + (ordered[high] - ordered[low]) * (position - low)


def summary(samples_us: list[float]) -> dict[str, float]:
    return {"p50_us": round(percentile(samples_us, .5), 3),
            "p95_us": round(percentile(samples_us, .95), 3),
            "mean_us": round(statistics.mean(samples_us), 3),
            "min_us": round(min(samples_us), 3), "max_us": round(max(samples_us), 3)}


def time_queries(search, queries, warmup: int):
    for query in queries[:warmup]:
        search(query)
    samples, result_ids = [], []
    for query in queries:
        started = time.perf_counter_ns()
        result = search(query)
        samples.append((time.perf_counter_ns() - started) / 1_000)
        result_ids.append([int(item[0]) for item in result])
    return summary(samples), result_ids


def overlap(left, right, label):
    recalls, jaccards = [], []
    for a, b in zip(left, right):
        aset, bset = set(a), set(b)
        common = len(aset & bset)
        recalls.append(common / len(bset) if bset else 1.0)
        union = aset | bset
        jaccards.append(common / len(union) if union else 1.0)
    return {"reference": label, "mean_top10_recall": round(statistics.mean(recalls), 6),
            "p05_top10_recall": round(percentile(recalls, .05), 6),
            "mean_top10_jaccard": round(statistics.mean(jaccards), 6),
            "queries": len(recalls)}


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

    import chematic
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import rdFingerprintGenerator

    library_input = smiles[:args.library]
    query_input = smiles[args.library:args.library + args.queries]
    native_index = chematic.PreparedFingerprintIndex.from_smiles(library_input, fp="ecfp4")
    compatible_index = chematic.PreparedFingerprintIndex.from_smiles(library_input, fp="rdkit_ecfp4")
    native_failed = set(native_index.failed_indices())
    compatible_failed = set(compatible_index.failed_indices())
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    rdkit_failed = {i for i, value in enumerate(library_input) if Chem.MolFromSmiles(value) is None}
    failed_library = native_failed | compatible_failed | rdkit_failed
    valid_indices = [i for i in range(len(library_input)) if i not in failed_library]
    valid_smiles = [library_input[i] for i in valid_indices]
    native_molecules = [chematic.from_smiles(value) for value in valid_smiles]
    rdkit_molecules = [Chem.MolFromSmiles(value) for value in valid_smiles]
    native_fps = [molecule.ecfp4() for molecule in native_molecules]
    rdkit_fps = [generator.GetFingerprint(molecule) for molecule in rdkit_molecules]

    valid_queries, failed_queries = [], []
    for i, value in enumerate(query_input):
        if Chem.MolFromSmiles(value) is None:
            failed_queries.append(i)
            continue
        try:
            chematic.from_smiles(value)
            compatible_index.search(value, k=1)
        except (ValueError, RuntimeError):
            failed_queries.append(i)
        else:
            valid_queries.append(value)
    if not valid_queries:
        raise SystemExit("no queries survived the shared validation scope")

    started = time.perf_counter_ns()
    native_index = chematic.PreparedFingerprintIndex.from_smiles(library_input, fp="ecfp4")
    native_build_ms = (time.perf_counter_ns() - started) / 1_000_000
    started = time.perf_counter_ns()
    compatible_index = chematic.PreparedFingerprintIndex.from_smiles(library_input, fp="rdkit_ecfp4")
    compatible_build_ms = (time.perf_counter_ns() - started) / 1_000_000
    started = time.perf_counter_ns()
    rdkit_fps = [generator.GetFingerprint(molecule) for molecule in rdkit_molecules]
    rdkit_build_ms = (time.perf_counter_ns() - started) / 1_000_000

    def native_prepared(query):
        return native_index.search(query, k=args.k)

    def native_reference(query):
        query_fp = chematic.from_smiles(query).ecfp4()
        scores = [(i, chematic.tanimoto(query_fp, fp)) for i, fp in zip(valid_indices, native_fps)]
        return sorted(((i, s) for i, s in scores if s > 0.0), key=lambda x: (-x[1], x[0]))[:args.k]

    def compatible_prepared(query):
        return compatible_index.search(query, k=args.k)

    def rdkit_search(query):
        query_fp = generator.GetFingerprint(Chem.MolFromSmiles(query))
        similarities = DataStructs.BulkTanimotoSimilarity(query_fp, rdkit_fps)
        scores = [(i, s) for i, s in zip(valid_indices, similarities) if s > 0.0]
        return sorted(scores, key=lambda x: (-x[1], x[0]))[:args.k]

    native_timing, native_ids = time_queries(native_prepared, valid_queries, args.warmup)
    native_ref_timing, native_ref_ids = time_queries(native_reference, valid_queries, args.warmup)
    compatible_timing, compatible_ids = time_queries(compatible_prepared, valid_queries, args.warmup)
    rdkit_timing, rdkit_ids = time_queries(rdkit_search, valid_queries, args.warmup)
    native_lane = overlap(native_ids, native_ref_ids, "native byte-fingerprint reference")
    compatible_lane = overlap(compatible_ids, rdkit_ids, "RDKit Morgan")
    cross_lane = overlap(native_ids, rdkit_ids, "RDKit Morgan diagnostic")
    result = {
        "schema_version": 2, "target_version": getattr(chematic, "__version__", "unknown"),
        "corpus": {"path": str(corpus.relative_to(ROOT)) if corpus.is_relative_to(ROOT) else str(corpus),
                   "sha256": hashlib.sha256(corpus.read_bytes()).hexdigest(), "total_rows": len(smiles),
                   "library_rows": args.library, "query_rows": args.queries,
                   "valid_library_rows": len(valid_indices), "valid_query_rows": len(valid_queries),
                   "split": "library=first rows; queries=following rows"},
        "configuration": {"fingerprint_native": "chematic ECFP4, 2048 bits",
                           "fingerprint_compatible": "RDKit Morgan radius=2, nBits=2048", "k": args.k,
                           "warmup_queries": args.warmup,
                           "ordering": "similarity descending, original input index ascending on ties"},
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "machine": platform.machine(), "chematic": getattr(chematic, "__version__", "unknown"),
                        "rdkit": rdBase.rdkitVersion},
        "failures": {"library": {"chematic_native": sorted(native_failed),
                                   "chematic_compatible": sorted(compatible_failed), "rdkit": sorted(rdkit_failed),
                                   "shared_excluded_count": len(failed_library)},
                     "queries": {"shared_excluded_indices": failed_queries, "count": len(failed_queries)}},
        "lanes": {
            "native_native": {"implementation": "PreparedFingerprintIndex(ecfp4) vs byte-fingerprint brute-force",
                              "prepared_build_ms": round(native_build_ms, 3), "prepared_query": native_timing,
                              "reference_query": native_ref_timing, "top10_recall": native_lane},
            "rdkit_compatible": {"implementation": "PreparedFingerprintIndex(rdkit_ecfp4) vs RDKit Morgan + BulkTanimotoSimilarity",
                                 "prepared_build_ms": round(compatible_build_ms, 3), "prepared_query": compatible_timing,
                                 "reference_build_ms": round(rdkit_build_ms, 3), "reference_query": rdkit_timing,
                                 "top10_recall": compatible_lane},
            "cross_profile": {"implementation": "chematic native ECFP4 vs RDKit-compatible Morgan",
                               "top10_recall": cross_lane,
                               "interpretation": "diagnostic cross-profile retrieval agreement, not a compatibility score"}},
        "gate": {"name": "rdkit_compatible_top10_recall", "target": 0.99,
                 "measured": compatible_lane["mean_top10_recall"],
                 "passed": compatible_lane["mean_top10_recall"] >= 0.99,
                 "scope": "shared valid library/query inputs"}}
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
