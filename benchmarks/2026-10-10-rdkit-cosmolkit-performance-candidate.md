# RDKit/COSMolKit performance candidate — 2026-10-10

This record evaluates commit `e3904b5acf0571c1f188ae159056f3c1d039a622`
against RDKit 2026.03.1 and COSMolKit 0.5.0rc15. It is source-candidate
evidence, not a published-package or universal speed claim.

## Contract

- Host: macOS 27.0.1, Apple silicon, 10 logical CPUs, Python 3.13.6.
- Input: the first 1,000 non-empty rows of
  `scripts/chembl_accuracy_corpus_4999.smi` (SHA-256
  `1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`).
- Artifact: CPython 3.13 arm64 wheel built from the clean commit above,
  SHA-256 `3526dfb908b96a078512b059a52c2381800e8d468917d3af500c6ff1816e91cb`.
- Timing: 21 process-isolated blocks with rotating engine order. The table
  reports medians. Confidence intervals are deterministic paired bootstrap
  intervals over the 21 per-block speed ratios.
- Operation-only creates fresh parsed molecules before each measured call but
  excludes that parsing time. Pipeline includes parsing. For MMFF it includes
  parsing, hydrogen/coordinate preparation, energy and gradient; the common
  RDKit seed-42 coordinates are prepared once outside all timed lanes.
- Output equivalence is a prerequisite. A faster result with different output
  is not counted as a win.

Raw timings: [JSON](2026-10-10-rdkit-cosmolkit-pipeline-candidate.json).

## Result

All seven measured pipelines, including the six operations that were slower
than COSMolKit in the v1.0.41 record, were faster than both comparators in all
21 blocks.

| Operation | chematic pipeline ms | COSMolKit ms | RDKit ms | speedup vs COSMolKit, 95% CI | speedup vs RDKit, 95% CI |
|---|---:|---:|---:|---:|---:|
| parse | 3.169 | 124.923 | 66.821 | 39.59x [38.44, 40.20] | 21.22x [20.57, 21.27] |
| ring count | 6.142 | 125.040 | 66.891 | 20.43x [19.97, 20.64] | 10.90x [10.78, 11.02] |
| TPSA | 7.004 | 125.968 | 69.404 | 18.07x [17.57, 18.13] | 9.93x [9.72, 10.00] |
| Labute ASA | 6.906 | 125.588 | 68.413 | 18.08x [17.84, 18.31] | 9.92x [9.74, 9.99] |
| chiral Morgan radius 2 | 26.322 | 138.802 | 88.067 | 5.25x [5.15, 5.32] | 3.34x [3.26, 3.37] |
| amide/amine reaction | 29.128 | 138.014 | 92.826 | 4.73x [4.69, 4.75] | 3.19x [3.16, 3.21] |
| MMFF energy + gradient | 357.351 | 753.857 | 663.543 | 2.11x [2.10, 2.12] | 1.86x [1.85, 1.86] |

MMFF also wins as a prepared operation: 1.53x [1.52, 1.54] versus COSMolKit
and 1.61x [1.61, 1.61] versus RDKit, with 21/21 winning blocks.

## Accuracy and regression gates

The candidate wheel was compared directly with RDKit on ChEMBL and NCI:

| Gate | ChEMBL | NCI |
|---|---:|---:|
| ring count | 5,000 / 5,000 | 4,991 / 4,991 |
| non-chiral Morgan bits | 5,000 / 5,000 | 4,991 / 4,991 |
| chiral Morgan bits | 5,000 / 5,000 | 4,991 / 4,991 |

The complete row-level result is
[`rdkit-fastpath-parity-e3904b5a.json`](../validation/results/rdkit-fastpath-parity-e3904b5a.json).
The bit-only Morgan implementation also equals the detailed bit-info path on
4,999/4,999 ChEMBL rows and 4,998/4,998 parsed NCI rows. Changed-crate Rust
tests, Python/WASM clippy, formatting, and the benchmark aggregation tests pass.
The broader 161-operation accuracy boundary remains the separately recorded
v1.0.41 source comparison; it is not widened by this performance result.

## What changed

- Folded Morgan callers no longer build sparse counts and bit provenance they
  do not return. Acyclic and non-stereo inputs avoid unnecessary allocations.
- One-cycle SMILES and acyclic parser provenance avoid redundant ring
  perception while mutation still invalidates the shortcut.
- RDKit-order MMFF setup no longer constructs a complete general MMFF model
  and then rebuilds the same force-field tables. Fixed-key maps use a faster
  deterministic hasher.
- The benchmark now preserves raw block samples, reports operation and
  pipeline lanes separately, and records paired confidence intervals.

## Remaining performance work

This does **not** establish a complete prepared-operation win. COSMolKit and
RDKit precompute ring/perception state while parsing, so operation-only ring
count, TPSA, Labute ASA, chiral Morgan and the reaction remain slower. Their
pipeline wins must not be presented as hot-call wins. The next gate is to make
each equivalent prepared operation's paired 95% lower bound exceed 1.0 without
moving expensive work into parsing solely to improve the benchmark. Public
PyPI/npm/crates.io and additional-host reruns also remain open.

