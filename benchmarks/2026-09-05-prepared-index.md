# Prepared fingerprint index

Measured 2026-09-05 on Apple arm64, macOS 26.5.2, from the local v1.0.6
working tree. This is source-level evidence, not a published-artifact claim.

`PreparedFingerprintIndex` computes the selected fingerprint for the database
once and reuses it for subsequent exact top-k queries. The one-shot comparison
rebuilds the database fingerprints for every query. Both paths use ECFP4, the
same ten-molecule fixture, and `k=5`.

| Lane | Median per ten queries |
| --- | ---: |
| Prepared index | 460.48 us |
| Rebuild database per query | 3.362 ms |
| Reuse speedup | **7.30x** |

The criterion ranges were 418.36–515.06 us and 3.2658–3.4635 ms. The prepared
index is an exact reuse optimization: the Rust tests compare its results with
the existing one-shot API, and the Python contract test checks self-hit order,
invalid-SMILES original-index mapping, and unknown fingerprint rejection.

## Reproduce

```bash
CARGO_TARGET_DIR=/private/tmp/chematic-binding-target \
  cargo bench -p chematic-fp --bench ecfp_bench --offline -- \
  prepared_search_10mol --noplot

CARGO_TARGET_DIR=/private/tmp/chematic-binding-target \
  cargo bench -p chematic-fp --bench ecfp_bench --offline -- \
  rebuild_search_10mol --noplot

/private/tmp/chematic-speed-venv-2/bin/python -m pytest -p no:asyncio \
  crates/chematic-py/tests/test_similarity_index.py -q
```

The current scope intentionally does not claim parallel/tiled similarity
matrix acceleration; that remains a separate roadmap item.
