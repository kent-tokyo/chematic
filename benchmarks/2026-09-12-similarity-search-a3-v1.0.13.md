# Similarity search A3 rerun — v1.0.13

This is the post-A3 rerun of the fixed 5,000-row corpus. It uses 4,500
library molecules, 500 queries, ECFP4/Morgan radius 2, 2,048 bits, `k=10`,
and deterministic input-index tie breaking. The machine-readable result is
[`2026-09-12-similarity-search-a3-v1.0.13.json`](2026-09-12-similarity-search-a3-v1.0.13.json).

Environment: chematic 1.0.13, RDKit 2026.03.6, CPython 3.13.6, macOS
26.5.2 arm64. Corpus SHA-256:
`d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`.

The previous run excluded 20 library rows and 2 query rows because the
RDKit-parity aromaticity path rejected valid explicit-aromatic fused inputs.
The new validated preservation path recovers all 4,500 library rows and all
500 queries. No rows are excluded in this rerun.

| Lane | Result |
|---|---:|
| Native ECFP4 vs native reference mean top-10 recall | 1.000 |
| RDKit-compatible Morgan vs RDKit mean top-10 recall | 1.000 |
| RDKit-compatible Morgan p05 top-10 recall | 1.000 |
| RDKit-compatible prepared-query p50 | 415.458 µs |
| RDKit-compatible prepared-query p95 | 860.562 µs |
| RDKit-compatible prepared-query mean | 492.520 µs |
| Cross-profile native vs Morgan diagnostic recall | 0.727 |

The A3 compatible retrieval gate passes (`target >= 0.99`, measured `1.0`).
The 0.727 cross-profile figure is retained as a diagnostic only: native
ECFP4 and RDKit Morgan are different profiles, so it is not an RDKit
compatibility score. Cross-binding parity in source-built Python, Node/WASM,
and Rust remains open.
