# Open-work ledger

Updated 2026-10-19 for **v1.1.8**. This file contains unresolved work and
current denominators only. Completed changes belong in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md);
dated measurements belong in the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Release and compatibility

- The public v1.1.8 macOS arm64 wheel completed the 10,000-row Labute ASA
  lane on one host. A broad cross-platform v1.1.8 artifact packet is still
  open; v1.1.7 remains the latest complete hosted chemistry workflow.
- Keep PyPI, npm/WASM, crates.io, browser, source, and operating-system results
  separate. Extend previous-release consumer tests where semver tooling cannot
  observe behavioral contracts.
- Preserve the versioned RDKit 2026.03.6/2026.09.1 ring disagreements and
  typed CIP abstentions instead of selecting one oracle as universal truth.

## Chemistry correctness

- **Standardization:** current source matches RDKit 2026.09.1 on Normalize
  4,065/4,069, Cleanup 3,923/4,069, and Reionize 3,150/4,069 rows. Implement
  high-value missing assignments and functional groups; keep independent
  ionization and metal-stereo policy differences explicit.
- **CTfile:** one chematic readback failure and ten stereo-only rows remain;
  16 RDKit oracle failures stay outside the chematic failure count.
- **Reaction/SMARTS:** retain the 80 exact reaction rows plus three jointly
  invalid inputs, the 310,000-cell SMARTS gate, autoconf probes,
  BioTransformer-derived cases, and minimized fuzz regressions.
- **R-group:** the terminal-wildcard subset is covered across bindings.
  Multi-core, MCS alignment, tautomer/core enumeration, and non-terminal
  R-groups remain outside the current contract.
- **File interchange:** extend V3000, MOL2, CML, CIF, SDF, CDXML, and related
  semantic tests to stratified real files before broad competitor claims.

Issues #734, #754, #739, #769, #779, and #785 are closed. Their fixtures
remain regression tests; closure does not broaden the supported domain.

## Performance and memory

- The public v1.1.8 ASA lane establishes exact output agreement and a bounded
  one-host result. COSMolKit first-use ASA is a tie; repeated and
  parse-inclusive lanes favor chematic. Six unequal COSMolKit ring-count rows
  remain unranked.
- Re-run ASA, TPSA, chiral Morgan, QED, and reaction findings from exact public
  artifacts on a second host. Split parsing, first use, prepared reuse,
  startup, throughput, and peak RSS.
- Do not hide preparation at another boundary or rank rows with unequal
  outputs, comparator errors, skips, or unsupported work.

## Independent and 3D evaluation

- Candidate `f5d6cd6f` predates the hash-pinned NCI/PubChem sources. MCS
  comparison (64), SDF attributes (64), and stereo/3D (96) are frozen and
  unscored. Evaluate once without changing the candidate; never tune on the
  observed outcomes.
- Retain typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time, and memory. Development sets and inspected holdouts cannot be
  promoted to independent evidence.
- Add 21 alternating blocks, a second host, peak memory, exact public
  artifacts, and unused quality data before a general 3D speed-and-quality
  claim.

## Recurring evidence refreshes

- Re-run relevant chemistry gates from exact registry artifacts.
- Pin comparator versions and hashes; report every failure class.
- Keep strict/reporting paths for lossy MOL, V3000, CDXML, Markush/polymer,
  Standard InChI, and canonical identity.
- Expand independent corpora for multi-reactant reactions, multi-core
  R-groups, unseen file formats, and OGB workflows.

## Dependency classes

- **`local-open`:** implementation, regression tests, documentation, source gates.
- **`toolchain-open`:** cross-platform wheels, npm/WASM builds, browser runners.
- **`external-open`:** comparator packages, upstream review, adjudication.
- **`data-sealed`:** one-time evaluation; inspected rows cannot be sealed again.
- **`historical`:** retained provenance that does not establish current behavior.

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

Every gate records support domain, versions, corpus, options, failure policy,
and counts for success, failure, refusal, unsupported, and skipped inputs. A
refusal is not a match, a source pass is not a package result, and registry
availability is not chemistry validation.
