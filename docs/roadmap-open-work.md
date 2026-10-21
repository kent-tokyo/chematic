# Open-work ledger

Updated 2026-10-21 for **v1.2.0**. Completed work
belongs in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md);
dated evidence belongs in the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Release and compatibility

- Re-run changed chemistry lanes from exact PyPI, npm/WASM and crates.io
  artifacts. Keep source, package, browser and operating-system evidence
  separate.
- Previous-release gates now use v1.1.9 for changed public Rust crates, Python
  declarations/dependencies and generated WASM declarations. Add consumer
  fixtures where static semver tools cannot observe behavior.
- Preserve RDKit 2026.03.6/2026.09.1 ring differences and typed CIP
  abstentions instead of treating one oracle version as universal truth.

## Chemistry correctness

- **MCS:** the source candidate fixes the v1.2.0 ring/chain input-order
  regression and exposes native search limits distinctly from timeouts across
  Python and WASM. Full FMCS feature/output parity, generalized SMARTS output,
  the existing explicit-H policy difference and a newly selected independent
  holdout remain open.
- **Standardization:** current source matches RDKit 2026.09.1 on Normalize
  4,065/4,069, Cleanup 3,933/4,069 and Reionize 3,150/4,069 rows. Prioritize
  the remaining bounded acid/base and metal disconnection/charge/bond/stereo
  rows; keep policy and representation differences explicit.
- **CTfile:** one chematic readback failure and ten stereo-only rows remain;
  16 RDKit oracle failures stay outside the chematic failure count.
- **Reaction/SMARTS:** retain the 80 exact reaction rows plus three jointly
  invalid inputs, the 310,000-cell SMARTS gate, autoconf probes,
  BioTransformer-derived cases and minimized fuzz regressions.
- **R-group:** the terminal-wildcard subset is covered across bindings.
  Multi-core alignment, tautomer/core enumeration and non-terminal R-groups
  remain outside the current contract.
- **File interchange:** expand V3000, MOL2, CML, CIF, SDF and CDXML semantic
  tests to stratified real files before broad competitor claims.

Issues #734, #754, #739, #769, #779 and #785 are closed. Their fixtures remain
regression tests; closure does not broaden the supported domain.

## Performance and memory

- The lazy Python source candidate shortens parse and parse-inclusive work by
  deferring descriptor, aromatic-view and ring preparation. First-use-only
  costs rise and prepared reuse is approximately flat. Verify a clean hosted
  wheel and second host before release claims.
- Re-run ASA, TPSA, chiral Morgan, QED and reaction lanes from exact public
  artifacts. Split parsing, first use, prepared reuse, startup and peak RSS.
- Do not rank rows with unequal outputs, comparator failures, skips or
  unsupported work, and do not hide preparation cost at another boundary.

## Independent and 3D evaluation

- Frozen holdout v1 was scored once against immutable v1.1.8. SDF graph,
  title, property and stereo columns are 64/64; MCS atom/bond counts are
  31/64; 3D embed status and stereo graph are 96/96, conformer count 95/96,
  and strict pair-distance coordinates 4/96.
- Treat those inspected rows as historical evaluation data. Use newly selected
  development rows for MCS and multi-conformer improvements; never rewrite the
  frozen score. The current MCS candidate was diagnosed against these inspected
  rows only after implementation, so its 34/64 exact rows and five typed
  resource limits are not a new sealed score.
- Retain typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time and memory. Add alternating blocks, a second host, exact
  public artifacts and unused quality data before general 3D claims.

## Recurring evidence refreshes

- Re-run the official RDKit.js comparison when its pinned package changes.
- Re-run public-package chemistry gates for every release candidate.
- Update competitor-watch claims only from pinned primary sources.

## Dependency classes

- **`local-open`:** implementation, regression tests, documentation and source
  gates that can be completed in this repository.
- **`toolchain-open`:** cross-platform wheels, npm/WASM builds and browser
  runners that require hosted toolchains.
- **`external-open`:** comparator packages, upstream review and independent
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
policy and counts for success, failure, refusal, unsupported and skipped
inputs. A refusal is not a match, source evidence is not package evidence and
registry availability is not chemistry validation.
