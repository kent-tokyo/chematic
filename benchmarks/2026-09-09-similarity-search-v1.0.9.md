# Similarity search comparison — 2026-09-09

This is a same-host local benchmark of exact top-k similarity search using a
4,500-molecule library and 500 fixed queries from the checked-in 5,000-row
ChEMBL-derived corpus. Both lanes include library SMILES parsing and
fingerprint construction in the build measurement, and query SMILES parsing,
fingerprint construction, full-library comparison, and top-k selection in the
query measurement. Ten queries warm each lane before the 500 timed queries.

The fingerprint definitions are named rather than treated as identical:
chematic uses native ECFP4 and RDKit uses Morgan radius 2. Both use 2,048 bits.

## Result

Measured with chematic v1.0.9 and RDKit 2025.09.3 on Python 3.13.6, macOS
26.5.2 arm64. `k=10`; timings are per query.

| Lane | Library build | Query p50 | Query p95 |
|---|---:|---:|---:|
| chematic `PreparedFingerprintIndex` + native ECFP4 | 888.616 ms | **265.688 µs** | **475.139 µs** |
| RDKit Morgan + `BulkTanimotoSimilarity` | 654.179 ms | 1,550.167 µs | 1,609.174 µs |

Under this exact protocol, the prepared chematic search query p50 is 5.83×
lower than the RDKit lane. This is an operation-specific result: it includes
different fingerprint implementations and must not be generalized to all
RDKit or all chematic workflows.

## Ranking agreement

Across the same 500 queries, the mean overlap with RDKit's top-10 candidate set
was 72.6% when measured as recall of the RDKit set; the 5th-percentile recall
was 30.0%, and mean top-10 Jaccard similarity was 60.8%. This is retrieval
agreement, not chemical ground truth or a claim that either ranking is correct.

The different bit assignments and aromaticity/invariant semantics are the
expected primary confounder. A fair bit-identical follow-up should use a
verified RDKit-compatible fingerprint on both sides, then repeat this protocol
without changing the split or tie-breaking rule.

## Reproduce

```bash
/tmp/chematic-sim-bench-venv/bin/python \
  scripts/bench_similarity_search_vs_rdkit.py \
  --output benchmarks/2026-09-09-similarity-search-v1.0.9.json
```

The machine-readable result is
[`2026-09-09-similarity-search-v1.0.9.json`](2026-09-09-similarity-search-v1.0.9.json).
The corpus SHA-256, configuration, engine versions, raw timing summary, and
ranking-overlap summary are recorded there.
