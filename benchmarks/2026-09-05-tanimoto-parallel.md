# Parallel Tanimoto matrix

Measured 2026-09-05 on Apple arm64, macOS 26.5.2, from the local v1.0.6
working tree. This is source-level evidence, not a published-artifact claim.

The new `tanimoto_matrix_parallel` implementation parallelizes independent
query rows and retains the existing dense row-major `Vec<f32>` contract.
For a 256×256 ECFP4 matrix from the checked-in ten-molecule fixture repeated
to 256 rows and columns:

| Lane | Median |
| --- | ---: |
| Existing serial matrix | 1.5248 ms |
| Parallel row-wise matrix | 1.2638 ms |
| Speedup | **1.21x** |

The Rust regression test compares every output value and position against the
serial implementation, including empty-input behavior. The Python bulk
matrix functions now use the parallel core path; the public dense shape and
invalid-SMILES filtering contracts are unchanged.

## Reproduce

```bash
CARGO_TARGET_DIR=/private/tmp/chematic-binding-target \
  cargo bench -p chematic-fp --bench ecfp_bench --offline -- \
  tanimoto_matrix_serial_256x256 --noplot

CARGO_TARGET_DIR=/private/tmp/chematic-binding-target \
  cargo bench -p chematic-fp --bench ecfp_bench --offline -- \
  tanimoto_matrix_parallel_256x256 --noplot
```

This does not claim threshold-aware pruning or tiled storage; those remain
separate follow-up work for larger matrices.
