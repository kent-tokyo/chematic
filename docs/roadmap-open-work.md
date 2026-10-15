# Open-work ledger

Updated 2026-10-16 for **v1.1.5**.
This file lists unresolved work only. Completed changes belong in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md); measurements belong in the
[validation report](validation.md) and [benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Active follow-up

- **v1.1.5 artifact packet:** after publication, a complete checked-in
  Linux/macOS/Windows Python, npm/WASM, and crates.io chemistry packet is still
  required before replacing the current v1.1.3 hosted chemistry packet and
  v1.0.42 broad performance packet.
- **Rust compatibility:** keep the v1.1.3 downstream compile fixture and
  `cargo-semver-checks` release lane. Public enum variants and struct fields
  remain semver-sensitive even when chemistry behavior is additive.
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
- **A6 quality:** topology bounds are aligned on the frozen corpus, while
  minimization and basin-diversity tails remain. The 265-row set and holdouts
  v1–v4 are development data. A fresh candidate must be frozen before new
  outcome-blind data are selected.
- **A6 publication:** add 21 alternating blocks, a second host, peak memory,
  exact public artifacts, and an unused quality set before claiming a general
  speed-and-quality win.
- **Independent workflow corpora:** add multi-reactant reactions, multi-core
  R-groups, and real MOL2/CJSON/CIF files while retaining refusals and semantic
  losses in the denominator.

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
