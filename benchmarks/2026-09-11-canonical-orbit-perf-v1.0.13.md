# Canonical SMILES orbit-pruning benchmark — 2026-09-11

This is a current-source remeasurement for issue #372. It compares the
orbit-pruned canonical search with the retained exhaustive implementation on
the checked-in Tier A/Tier B fixtures. Both paths receive the same parsed
molecules; this is not a cross-engine chemistry-compatibility measurement.

Environment: chematic workspace version 1.0.13, release build, Rust offline
build, arm64 macOS. Command:

```bash
CANONICAL_ORBIT_PERF_RUN_LEGACY=1 \
cargo run -p chematic-smiles --example canonical_orbit_perf --release --offline \
  --features canonical-search-instrumentation
```

All 20 fixtures parsed, produced non-empty output, and had zero output
mismatches between the exhaustive and orbit-pruned implementations.

The Tier A set now includes an explicitly named minimized RENKIN-shaped
`CC(C)(C)OC(=O)NCC` Boc/tBu witness. In the verbose run it reduced exhaustive
leaves from 6 to 1; the aggregate Tier A result was 7.07x geometric-mean
speedup with 99.8% fewer leaves. This is a minimized local witness, not the
held-out `uspto50k_test.smi` target from RENKIN, so downstream end-to-end
acceptance remains separate.

| Tier | Fixtures | Geomean speedup (exhaustive/orbit-pruned) | Leaf reduction |
|---|---:|---:|---:|
| High symmetry | 14 | 7.07x | 99.8% |
| Low symmetry control | 6 | 1.42x | 98.0% |

The instrumented search visited 62 nodes and wrote 13 leaves for Tier A, and
12 nodes and 3 leaves for Tier B. No search budget was exhausted.

This supports the existing exact local-twin/orbit-pruning optimization. It
does not close issue #372 by itself: the exact RENKIN witness is not present
in this repository, so the downstream 2x acceptance target remains
unmeasured. The low-symmetry control also shows that this is not a universal
speedup claim.
