# Similarity search comparison — 2026-09-11

This is a same-host local benchmark of exact top-10 search using the fixed
4,500-molecule library and 500-query split from the checked-in 5,000-row
corpus. Two query rows and 20 library rows were excluded from the shared valid
scope because the new RDKit-compatible preprocessing rejected them. Failures
are recorded separately in the JSON result.

## Result

Measured with chematic v1.0.12 and RDKit 2026.03.6 on Python 3.13.6, macOS
26.5.2 arm64. `k=10`; ties use score descending, then original input index
ascending.

| Lane | Search/reference | Query p50 | Top-10 mean recall |
|---|---|---:|---:|
| native/native | chematic prepared native ECFP4 / native byte-FP brute force | 219.666 µs / 1,562.062 µs | 99.5984% |
| RDKit-compatible | chematic prepared `rdkit_ecfp4` / RDKit Morgan | 406.438 µs / 1,585.522 µs | **100.0000%** |
| cross-profile | chematic native ECFP4 / RDKit Morgan | — | 72.6506% (diagnostic) |

The RDKit-compatible gate passes the roadmap target of at least 99% mean
top-10 recall on the 4480-library / 498-query valid-input scope. The previous
72.6% result was a cross-profile comparison and is not a compatibility score.
The native/native lane is an independent implementation check; its small
residual is retained as evidence rather than rounded to 100%.

## Failure accounting

The compatible profile rejected library indices
`1206, 1213, 1214, 1231, 1238, 1287, 1370, 1371, 1829, 1830, 1843, 1975,
2007, 2018, 2079, 2080, 2103, 4385, 4471, 4498` during preprocessing. Two
query rows were excluded from the shared scope. No RDKit library parse failures
were observed. See the machine-readable record for the complete lists.

## Reproduce

```bash
.venv-v1gate/bin/pip install rdkit
.venv-v1gate/bin/maturin build --release -m crates/chematic-py/Cargo.toml --offline -o target/chematic-wheels
.venv-v1gate/bin/pip install --force-reinstall target/chematic-wheels/chematic-1.0.12-*.whl
.venv-v1gate/bin/python scripts/bench_similarity_search_vs_rdkit.py \
  --output benchmarks/2026-09-11-similarity-search-v1.0.12.json
```

The machine-readable result is
[`2026-09-11-similarity-search-v1.0.12.json`](2026-09-11-similarity-search-v1.0.12.json).
