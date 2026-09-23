# Roadmap open-work ledger

Updated 2026-09-24. This is a compact dependency and evidence ledger for work
that remains open after v1.0.22. It does not repeat completed implementation
history; use CHANGELOG, dated validation artifacts, and Git history for that.

## Immediate queue

| Order | Work | Current evidence | Exit |
|---:|---|---|---|
| 1 | #632 SMILES E/Z semantics | Released source `9808f54f`: RDKit 2026.03.6 fixed 10k lane moves from 18 semantic differences to 0; graph differences remain 0 | Complete the 28 x 1,024 relabel audit without widening the stable-key contract; do not count the interrupted audit as passed |
| 2 | #634 CIP residuals | 9,770/10,000 correspondence-correct exact labels; 230 classified rows remain | Each row is fixed, explicitly unsupported, or justified as an oracle/representation boundary; order invariance remains intact |
| 3 | #635 SMARTS residuals | 14,306 differing cells across 310,000; 3,364 rows affected | Differences are grouped by semantic family with truth-table regressions and explicit unsupported boundaries |
| 4 | A6 MMFF94 | Public v1.0.20 passes 265/265 stereo/clash quality but measures 0.944x RDKit speed; current source has a separate speed candidate | Energy/term residuals, timeout/convergence accounting, broader conformer quality, and a published-package rerun pass separately |
| 5 | Next RDKit rebaseline | Historical 2026.03.6 Python/npm packet is reproducible | Run only after the next official stable artifact is pinned; retain old/new results side by side |

## Accuracy packages

| Package | State | Remaining work |
|---|---|---|
| A0 Evaluation contract | Complete | Preserve frozen/exposed cohort rules and failure accounting for every later package |
| A1 Perception and descriptors | Open | Additional descriptor families, aromaticity residuals, potential-center coverage, independent holdouts |
| A2 Stereo and identity | Active | #632, #634, spelling/permutation/file round trips, stable-key scope, enhanced stereo |
| A3 Fingerprints and retrieval | Open | Option coverage, top-k invariants, named compatible-Morgan parity, cross-binding checks |
| A4 Workflows and interchange | Open | #635, reactions, V3000/query semantics, attachment metadata, batch accounting |
| A5 Independent adjudication | External/open | Independently reviewed gold data and non-maintainer assessment without exposed-cohort reuse |
| A6 3D and force fields | Active | Energy/terms, timeout/convergence, conformer quality, publication-level speed and reproducibility |

An open package can include completed sub-gates. Only A0 currently satisfies
all declared exits.

## Cross-cutting work

| Area | Remaining exit |
|---|---|
| Parser security | Broaden malformed parser-state coverage while retaining process/time/memory bounds and typed outcomes |
| BatchResult | Preserve original index and terminal outcome across Rust/Python/Node/WASM, including cancellation and unprocessed ranges |
| Browser runtime | Worker cancellation, structured limits/errors, package startup/memory evidence, and multi-browser reproducibility |
| Interchange | Separate opaque preservation from editable semantics for V3000 coordination/haptic/polymer and rich CDXML objects |
| Reactions/SMARTS | Broader semantic truth tables and curated precision/recall/timeout reports |
| Release operations | Version/docs/metadata synchronization, package builds, registry publication, GitHub release, Pages and channel verification |
| External review | Non-maintainer review, security audit, and later advisory/backport rehearsal |

## Dependency classification

Every open item uses one of these labels:

| Class | Meaning | Current examples |
|---|---|---|
| `local-open` | Can be implemented and checked in this repository with available source and fixtures | #632 focused checks, #634/#635 classification, A6 numerical gates, documentation synchronization |
| `toolchain-open` | Repository work is defined, but a specific browser, compiler target, package, or runtime must be available | rebuilt Python/Node/WASM artifacts, browser-memory lanes, cross-platform execution |
| `data-sealed` | Evaluation data must remain unused until a frozen candidate and overlap audit are recorded | future independent descriptor, stereo, retrieval, or conformer holdouts |
| `external-open` | Requires publication, credentials, an official future artifact, or a non-maintainer decision | next RDKit stable rebaseline, registry release, independent review/security audit |
| `historical` | Retained for traceability but not a current completion gate | abandoned additional 1.10x SMILES parse stretch and rejected sealed candidates |

`local-open` and available `toolchain-open` work may proceed autonomously.
`data-sealed` inputs cannot be inspected for tuning. `external-open` work is not
complete merely because a local packet exists.

## Evidence boundaries

- **Published release:** v1.0.22 channel verification is prepared in
  `validation/results/release-channel-verification-v1.0.22.json`; v1.0.21 is
  the latest fully verified historical channel record until each v1.0.22 channel
  is observed.
- **Source operation matrix:**
  `benchmarks/2026-09-24-python-op-matrix-vs-rdkit-perf-branch.md` records
  strict output-comparable outcomes; it is not a published-package result.
- **Browser performance:** the current published-package record is
  `benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md`.
- **A0 acceptance:** the commit-safe summary is
  `validation/results/a0-core-eight-sealed-acceptance-20260922.json`.
- **RDKit rebaseline:** the pinned 2026.03.6 artifacts and complete exposed
  rows are under `validation/results/rdkit-rebaseline-*`.
- **#632 release-source diagnostic:**
  `validation/results/smiles-ez-semantic-issue632-v1.0.20-candidate-vs-rdkit-2026.03.6-2026-09-23.json`
  is a bounded source diagnostic; the interrupted long audit remains open.
- **A6 source candidate:**
  `benchmarks/2026-09-23-mmff94-stereo-safe-performance.md` is not a
  registry-package result.

## Completion rules

An item is complete only when:

1. its support domain, comparator, versions, options, corpus, and failure
   policy are declared;
2. success, failure, refusal, and skipped inputs account for the whole input;
3. deterministic regressions cover the corrected behavior;
4. affected bindings are rebuilt and checked where the API is shared;
5. the evidence is committed and linked from the maintained documentation;
6. local, source-candidate, published-package, and public-channel states are
   named accurately.

No row is silently dropped, no exposed cohort becomes sealed again, and no
safe refusal is counted as an exact match.
