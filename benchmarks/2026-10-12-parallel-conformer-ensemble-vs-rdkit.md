# Parallel conformer ensemble candidate versus RDKit 2026.09.1

This record covers a narrow development milestone: parallelizing independent
CheMatic conformer attempts without changing their deterministic result. It is
not a broad claim that CheMatic's 3D stack is more accurate or universally
faster than RDKit.

## Provenance

- CheMatic source: `2f6f9516bd77e06d6d7d17453e26edf381cbc034`
- RDKit: `2026.09.1`, isolated CPython 3.13.6 environment
- Host: Mac16,12, arm64, 10 logical CPUs, macOS 27.0.1 (26A434)
- Rust: 1.97.0
- Molecules: aspirin, caffeine, and ibuprofen (the exact SMILES are in the raw
  JSON)
- Work per molecule: 10 seeded attempts, seed 42, MMFF maximum 200 iterations
- Timing: 21 alternating blocks within each engine; lower is faster

Raw samples are in
[`2026-10-12-parallel-conformer-ensemble-chematic.json`](2026-10-12-parallel-conformer-ensemble-chematic.json)
and
[`2026-10-12-parallel-conformer-ensemble-rdkit.json`](2026-10-12-parallel-conformer-ensemble-rdkit.json).

## Result

| Lane | Median for 3 molecules |
|---|---:|
| CheMatic sequential | 62.417 ms |
| CheMatic parallel, 10 Rayon threads | 17.649 ms |
| RDKit single thread | 216.854 ms |
| RDKit all threads | 68.300 ms |

The CheMatic parallel path is 3.54x faster than its sequential path. Against
the recorded RDKit lane, the median ratio is 3.87x for both engines' native
parallel settings and 3.47x for their sequential settings.

Before every CheMatic timing run, the benchmark asserts exact equality between
the sequential and parallel results: attempt indices, seeds, success/failure,
energies, force-field provenance, duplicate dispositions, RMSD bit patterns,
and every output coordinate are identical.

## What this does and does not prove

- It proves that CheMatic can parallelize the independent generation phase
  without changing its established deterministic result on these inputs.
- It is a strong candidate speed result, but it is not yet a formal
  equivalent-output win against RDKit. The engines do not produce the same
  conformer ensemble, and the cross-engine runs were not interleaved in one
  harness.
- The set contains only three ordinary drug-like molecules. It does not cover
  macrocycles, metals, highly flexible molecules, cages, failed embeddings, or
  peak memory.
- This is source-candidate evidence. It is not evidence from a published wheel
  or crate, and it has not been repeated on a second host.
- Older A6 evidence keeps energy, gradient, convergence, stereo, clash, and
  best-of-N quality separate. Those quality gates must be rerun against RDKit
  2026.09.1 before a broader superiority claim.

## Reproduce

```bash
SCHEMATIC_BENCH_REPEATS=21 cargo run --release -p chematic-3d \
  --example conformer_ensemble_parallel_bench

python scripts/benchmarks/bench_rdkit_conformer_ensemble.py --repeats 21
```

The next acceptance packet will add the frozen A6 corpus plus an independent
holdout, paired quality measures, peak memory, a second host, and exact public
artifacts.
