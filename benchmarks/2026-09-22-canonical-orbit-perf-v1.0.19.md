# Canonical SMILES orbit pruning — exact RENKIN witness

This record closes the missing measurement in issue #372. It compares the
current exact local-twin/orbit-pruned canonical search with the retained
pre-pruning exhaustive implementation on the exact RENKIN issue #128 target 2,
not only on the earlier minimized Boc/tBu proxy.

## Scope

- Source: local v1.0.19 candidate based on `main` at `ab1b7176`; production
  canonicalization code is unchanged by this measurement.
- Witness: line 5 of RENKIN `data/uspto50k_test.smi`, recorded in
  `canonical_orbit_perf.rs` as `RENKIN exact target 2`.
- Operation: parse once, then generate canonical SMILES with the exhaustive
  reference and current orbit-pruned engines.
- Environment: arm64 macOS, Rust release build.
- This is a canonicalization microbenchmark. It is not a rerun of RENKIN's
  complete reaction search.

## Command

```bash
CANONICAL_ORBIT_PERF_VERBOSE=1 \
cargo run --release -p chematic-smiles \
  --features canonical-search-instrumentation \
  --example canonical_orbit_perf
```

The command was run five times after compilation. Every run passed the full
Tier A/B old-vs-new output differential.

## Exact witness result

| Metric | Exhaustive reference | Current orbit-pruned |
|---|---:|---:|
| Search leaves | 24 | 1 |
| Search nodes | — | 5 |
| Automorphism tests | — | 1 |
| Pruned children | — | 5 |
| Median of five elapsed samples | 102 µs | 21 µs |

The paired median ratio is **4.86x**. Individual old/new samples in
microseconds were `102/21`, `101/23`, `100/23`, `108/18`, and `108/21`.
Timing at this scale is environment-sensitive, so the exact structural result
(24 leaves to 1) is the stronger mechanism evidence; both exceed issue #372's
preferred 2x target.

Across the complete 15-fixture high-symmetry tier, the measured geometric-mean
speedup was 6.77x to 6.95x across the recorded runs, with 6,216 exhaustive
leaves reduced to 15 and no budget exhaustion. The six-fixture low-symmetry
control was not slower in these runs (1.08x to 1.25x), but that small timing
difference is not presented as a general speed claim.

## Correctness boundary

- Exact witness: exhaustive output equals current output; output is non-empty.
- Tier A/B: 0 parse failures, 0 empty outputs, 0 old-vs-new mismatches.
- Existing canonical, stereo, isotope, charge, map, dative-bond, and
  disconnected-graph regression suites remain the release gate.
- The result establishes that the residual downstream witness cost was the
  already-addressed local symmetry pattern. It does not claim that every part
  of RENKIN's reaction-template pipeline is 4.86x faster.

Machine-readable evidence is in
[`canonical-orbit-perf-v1.0.19.json`](../validation/results/canonical-orbit-perf-v1.0.19.json).
