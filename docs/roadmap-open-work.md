# Open-work ledger

Updated 2026-10-18 for **v1.1.7**.
This file lists unresolved work only. Completed changes belong in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md); measurements belong in the
[validation report](validation.md) and [benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Active follow-up

- **v1.1.7 artifact packet:** run the hosted Linux/macOS/Windows chemistry
  workflow after publication and check in artifact/corpus hashes and failure
  counts. The successful v1.1.5 hosted run remains the latest complete package
  lane; npm/WASM and crates.io remain separate evidence channels.
- **Descriptor performance:** repeated QED now has a mutation-safe cache.
  Complete the hosted Linux 20-block run, compare it with macOS, and keep
  first-use chiral Morgan and Labute ASA deficits visible. Do not turn a hot
  cache win into a first-use or toolkit-wide claim.
- **Standardization implementation:** every current Normalize, Reionize and
  Cleanup residual is classified from the same source revision. Prioritize 36
  missing Reionize assignments, aromatic N-oxide kekulization and common
  functional groups. Independent ionization and metal-stereo policy
  differences remain declared contracts.
- **CTfile:** current source has 4,045 exact readbacks, ten stereo-only rows,
  one chematic readback failure and 16 RDKit oracle failures. Fix the chematic
  rows without absorbing oracle failures into the success count.
- **Compatibility:** changed public crates, Python declarations/dependencies
  and WASM declarations now have previous-release gates. Extend feature-matrix
  consumer fixtures where `cargo-semver-checks` cannot observe behavior.
- **Equivalent-work performance:** close prepared-operation deficits without
  moving work into parsing. Separate parse-inclusive, first-use, prepared,
  startup, and peak-memory results; add a second host before broad claims.
- **RDKit version contract:** retain the versioned 2026.03.6/2026.09.1 ring
  disagreements and typed CIP abstentions. Do not present one oracle version
  as a universal chemical truth.
- **File interchange:** expand V3000, MOL2, CML, CIF, SDF, and CDXML semantic
  tests to stratified real-file corpora before broad Open Babel claims.
- **Reaction and standardization:** keep the 83-row reaction and fuzz gates,
  then expand multi-reactant, explicit-H, radical, charge, stereo, and
  rejection strata. Broader Cleanup parity remains open.
- **R-group:** maintain the terminal-wildcard subset across Rust, Python, and
  WASM. Multi-core, MCS alignment, tautomer/core enumeration, and non-terminal
  R-groups remain outside the current contract.
- **CDK/Indigo comparison:** expand beyond the bounded 15 reaction, 39 parent,
  and 16 labelled R-group rows; add prepared-object, memory, second-host, and
  exact-package lanes.
- **Frozen trust holdouts:** candidate `f5d6cd6f` predates the acquired,
  hash-pinned NCI/PubChem sources. MCS comparison (64), SDF attributes (64),
  and stereo/3D (96) are unscored. Evaluate them once without changing the
  candidate, and do not use their outcomes for tuning.
- **A6 quality:** topology bounds are aligned on the development corpus, while
  minimization and basin-diversity tails remain. Retain quality tails,
  convergence, time, and memory when the new stereo/3D holdout is evaluated.
- **A6 publication:** add 21 alternating blocks, a second host, peak memory,
  exact public artifacts, and an unused quality set before claiming a general
  speed-and-quality win.
- **Independent workflow corpora:** add multi-reactant reactions, multi-core
  R-groups, and unseen SDF/V3000/MOL2/CJSON/CIF files while retaining refusals
  and semantic losses in the denominator. Extend the OGB lane beyond the known
  EGFR sample before making a broad throughput claim.

Issues #734 and #754 are closed, as are #739, #769, #779, and #785. Their
fixtures remain regression tests; closed status does not broaden the supported
domain.

## Recurring evidence refreshes

- Rerun chemistry, reaction, SMARTS, CIP, fingerprint, format, and 3D gates
  from exact registry artifacts on their declared platforms.
- Compare equal work and equal acceptable outputs in alternating repeated
  blocks; split parsing, perception, prepared reuse, and memory.
- Pin new comparator artifacts and hashes. Keep source, package, and browser
  lanes separate.
- Retain strict/reporting paths for lossy MOL, V3000, CDXML, Markush/polymer,
  Standard InChI, and canonical identity.

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
