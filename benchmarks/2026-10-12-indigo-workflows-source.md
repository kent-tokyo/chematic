# CheMatic / Indigo workflow comparison (source candidate)

Date: 2026-10-12

This record compares three bounded workflow surfaces. It is not a general
claim that CheMatic is faster or more accurate than Indigo across all
reactions, standardization policies, or R-group decomposition modes.

## Provenance

- CheMatic: `1.1.0`, source commit
  `9bfc4ef9ea22cbce5db7aa23a77056d214538656`
- Indigo: `1.46.0` (`epam.indigo`)
- RDKit observer: `2026.03.6`
- Host: macOS 27.0.1, arm64, CPython 3.13.6
- Raw record: [`2026-10-12-indigo-workflows-source.json`](2026-10-12-indigo-workflows-source.json)

The runner records SHA-256 values for every corpus. Timings use 21
alternating-order blocks. A speed ratio is emitted only for rows where both
engines exactly satisfy the declared observer.

## Results

| Workflow | Accuracy: CheMatic | Accuracy: Indigo | Timed exact rows | CheMatic median | Indigo median | Indigo / CheMatic median | Paired 95% interval | Blocks won |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Single-reactant, single-product transform | 15/15 | 13/15 | 13 | 0.305 ms | 3.180 ms | **10.66x** | 10.17–11.13x | 21/21 |
| Largest-fragment parent selection | 39/39 | 38/39 | 38 | 1.867 ms | 2.347 ms | **1.27x** | 1.24–1.28x | 21/21 |
| Bounded labelled R-group decomposition | 16/16 | 16/16 | 16 | 1.528 ms | 3.989 ms | **2.60x** | 2.59–2.63x | 21/21 |

CheMatic's two reaction advantages are stereochemical rows: Indigo inverted
one identity transform and returned no match for the other. Indigo's one
fragment-selection miss chose explicit hydrogen over the expected propane
parent. These rows remain in the accuracy denominator and are excluded only
from the equivalent-output timing subset.

## Boundaries

- The reaction lane is restricted to one reactant and one single-molecule
  product because Indigo's `transform` API is not a product-set enumeration
  surface. Multi-product and multi-reactant rows remain covered by CheMatic's
  separate 83-row and BioTransformer-derived gates.
- The standardization lane measures the checked-in Phase-1 parent-selection
  policy, not complete equivalence between CheMatic's full standardization
  pipeline and Indigo's configurable standardizer.
- The R-group lane covers terminal substituents on four declared single cores.
  It does not cover multiple cores, MCS alignment, tautomer/core enumeration,
  non-terminal groups, or GA scoring.
- This is a source-built result on one host. Published packages, a second host,
  memory use, larger frozen corpora, and prepared-object APIs remain open.

## Reproduce

```bash
.venv/bin/python validation/indigo_workflows/benchmark.py \
  --blocks 21 \
  --output benchmarks/2026-10-12-indigo-workflows-source.json
```

The runner refuses comparator-version drift and reports full-corpus accuracy
separately from the common-correct timing subset.
