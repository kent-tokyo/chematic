# Open-work ledger

Updated 2026-10-16 for **v1.1.5**.
This file lists unresolved work only. Completed changes belong in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md); measurements belong in the
[validation report](validation.md) and [benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Active follow-up

- **v1.1.5 artifact packet:** the hosted published-wheel workflow completed all
  11 Linux/macOS/Windows, RDKit-comparison, A6, and autoconf jobs. Check in a
  durable packet with artifact/corpus hashes and failure counts; npm/WASM and
  crates.io remain separate evidence channels. The v1.0.42 broad performance
  packet is not superseded by this correctness run.
- **QED and LogP performance:** reproduce the externally reported v1.1.4 to
  v1.1.5 slowdown with 21 alternating blocks before tuning. Preserve exact
  outputs and split parse, first-use, prepared reuse, and host conversion.
- **Standardization:** reproduce the reported 254-row lane and classify
  Normalize, Reionize, and Cleanup residuals as bug, semantic equivalent,
  policy difference, or typed unsupported. Only the bug class is parity work.
- **Rust compatibility:** extend the v1.1.3 downstream fixture and
  `cargo-semver-checks` from selected surfaces to all 20 public crates and
  relevant feature sets. Cover enums, error types, structs, and traits.
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
  v1–v4 are development data. Freeze a new outcome-blind holdout before any
  further 3D tuning and retain quality tails, convergence, time, and memory.
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
