# Open-work ledger

The [CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md)
tracks completed changes. The [validation report](validation.md) records
published v1.2.2 and links the dated evidence.

## Release and compatibility

- Check the v1.2.3 source and bindings against the v1.2.2 tag. Verify [Rust Cleanup compatibility](cleanup-migration.md) against v1.2.2.
- Qualify the exact release artifacts across Rust, Python and WASM. Preserve
  comparator versions and typed refusals when reporting results.

## Chemistry correctness

- **MCS:** v1.2.1 fixes the v1.2.0 ring/chain input-order regression,
  equal-atom bond-count ties, complete-ring valid-best selection, and typed
  resource limits. Rust MCS queries now retain observed atom and bond
  alternatives under relaxed comparators. Python and WASM now expose these
  alternatives as SMARTS for configurations without charge, isotope, or
  chirality matching. Encoding those constraints in the result, full FMCS
  parity, the explicit-H policy difference, and a newly selected independent
  holdout remain open.
- **Standardization:** adjudicate the remaining Cleanup graph differences
  and the hydrogen-total change in the 575-row metal audit. The
  [Cleanup audit](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-21-cleanup-v122-regression-audit.md)
  records the published wheel and the subsequently tested candidate.
- **CTfile:** one chematic readback failure and ten stereo-only rows remain;
  16 RDKit oracle failures stay outside the chematic failure count.
- **Reaction/SMARTS:** retain the 80 exact reaction rows plus three jointly
  invalid inputs, the 310,000-cell SMARTS gate, autoconf probes,
  BioTransformer-derived cases, and minimized fuzz regressions.
- **R-group:** terminal wildcards and selected labelled/symmetric rows are
  covered. Multi-core alignment, tautomer/core enumeration, and general
  non-terminal R-groups remain outside the contract.
- **File interchange:** expand V3000, MOL2, CML, CIF, SDF, and CDXML semantic
  tests with stratified real files before broad competitor claims.

Closed issues remain represented by regression fixtures; closure does not
broaden the supported domain.

## Performance and memory

- Extend the v1.2.2 single-host [descriptor and MCS comparison](benchmark.md)
  to a second host. Measure ASA, QED, reactions and 3D on exact packages where
  comparable. Record parsing, first use, prepared reuse, startup and peak RSS.
- Do not rank rows with unequal outputs, comparator failures, skips, or
  unsupported work. Preparation cost cannot be hidden at another boundary.

## Independent and 3D evaluation

- Frozen holdout v1 was scored once against immutable v1.1.8. SDF graph,
  title, property, and stereo columns are 64/64; MCS atom/bond counts are
  31/64; 3D embed status and stereo graph are 96/96, conformer count 95/96,
  and strict pair-distance coordinates 4/96.
- Those inspected rows are historical evaluation data. Use newly selected
  development rows for MCS and multi-conformer work; do not rewrite the frozen
  score.
- Evaluate typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time, and memory separately on unused data.

## Recurring evidence refreshes

- Re-run the official RDKit.js comparison when its pinned package changes.
- Re-run public-package chemistry gates for every release candidate.
- Update competitor-watch claims only from pinned primary sources.

## Dependency classes

- **`local-open`:** implementation, tests, documentation, and source gates in
  this repository.
- **`toolchain-open`:** cross-platform wheels, npm/WASM builds, and browser
  runners that require hosted toolchains.
- **`external-open`:** comparator packages, upstream review, and independent
  adjudication.
- **`data-sealed`:** one-time evaluation data; inspected rows cannot become a
  sealed set again.
- **`historical`:** retained provenance that does not establish current
  behavior.

## Accuracy packages

| Package | State |
|---|---|
| A0 Evaluation contract | Complete |
| A1 Perception and descriptors | Active |
| A2 Stereo and identity | Active |
| A3 Fingerprints and retrieval | Maintenance |
| A4 Workflows and interchange | Active |
| A5 Independent adjudication | External |
| A6 3D and force fields | Experimental |

## Completion rule

Every gate records its support domain, versions, corpus, options, failure
policy, and counts for success, failure, refusal, unsupported, and skipped
inputs. A refusal is not a match, source evidence is not package evidence, and
registry availability is not chemistry validation.
