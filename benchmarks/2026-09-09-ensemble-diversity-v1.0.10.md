# v1.0.10 deterministic ensemble diversity

Local evidence from `distance_geometry_v2_gap_check` on macOS arm64.

Command:

```text
cargo run --release -p chematic-3d --example distance_geometry_v2_gap_check
```

The gate reproduced all 58 corpus molecules bit-identically for the same
seed, and all 58 produced different coordinates for adjacent seeds. Pairwise
Kabsch-aligned RMSD across eight seeds was positive for each flexible-corpus
probe:

| Molecule | Minimum Å | Mean Å | Maximum Å |
|---|---:|---:|---:|
| decane | 0.7099 | 1.4885 | 2.0473 |
| diphenhydramine | 0.9670 | 1.8131 | 2.4630 |
| hexadecane | 1.0724 | 2.1435 | 2.8972 |

This closes the local deterministic-diversity measurement slice only. It does
not close MMFF94/UFF soundness, external-oracle parity, or class-level failure
rate gates.

Machine-readable evidence: [`2026-09-09-ensemble-diversity-v1.0.10.json`](2026-09-09-ensemble-diversity-v1.0.10.json).
