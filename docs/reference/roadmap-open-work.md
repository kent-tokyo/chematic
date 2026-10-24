# Open-work ledger

Updated 2026-10-23 for **v1.2.1**. Completed
work belongs in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md);
dated evidence belongs in the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Release and compatibility

- Re-run changed chemistry lanes from exact PyPI, npm/WASM, and crates.io
  artifacts. Source, package, browser, and operating-system evidence stay
  separate.
- Compare every public Rust crate with the immediately previous release and
  retain consumer fixtures where static semver tools cannot observe behavior.
- Preserve RDKit 2026.03.6/2026.09.1 ring differences and typed CIP
  abstentions instead of treating one oracle version as universal truth.

## Chemistry correctness

- **MCS:** v1.2.1 fixes the v1.2.0 ring/chain input-order regression,
  equal-atom bond-count ties, complete-ring valid-best selection, and typed
  resource limits. Rust MCS queries now retain observed atom and bond
  alternatives under relaxed comparators. Python and WASM now expose these
  alternatives as SMARTS for configurations without charge, isotope, or
  chirality matching. Encoding those constraints in the result, full FMCS
  parity, the explicit-H policy difference, and a newly selected independent
  holdout remain open.
- **Standardization:** current source matches RDKit 2026.09.1 on Normalize
  4,069/4,069, Cleanup 4,057/4,069, and Reionize 4,069/4,069 comparable rows.
  Classify the remaining 12 Cleanup rows before adding more transforms.
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

- Re-run ASA, TPSA, chiral Morgan, QED, reaction, and 3D lanes from exact
  public artifacts. Split parsing, first use, prepared reuse, startup, and peak
  RSS.
- Add a second host to current macOS source measurements.
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
