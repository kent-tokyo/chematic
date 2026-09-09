# v1.0.10 3D class-level failure rates

The existing `distance_geometry_v2_gap_check` was re-run on the frozen 58-
molecule corpus. A row is successful when its engine status is `ok`.

Command:

```text
cargo run --release -p chematic-3d --example distance_geometry_v2_gap_check
```

| Category | Checked | New embedder failures | Legacy DFS failures |
|---|---:|---:|---:|
| rigid_ring | 8 | 0 | 2 |
| fused_aromatic | 7 | 0 | 0 |
| flexible_chain | 6 | 0 | 4 |
| macrocycle | 3 | 0 | 0 |
| stereocenter_implicit_h | 10 | 0 | 0 |
| stereocenter_quaternary | 6 | 0 | 0 |
| alkene_ez | 8 | 0 | 0 |
| druglike | 7 | 0 | 1 |
| druglike_rigid | 1 | 0 | 0 |
| druglike_stress | 2 | 0 | 1 |
| **total** | **58** | **0** | **8** |

This closes the local status-class measurement slice only. It does not close
MMFF94/UFF soundness, external-oracle parity, symmetry-aware TFD, or statistical
generalization.

Machine-readable evidence: [`2026-09-09-3d-class-failure-rates-v1.0.10.json`](2026-09-09-3d-class-failure-rates-v1.0.10.json).
