# chematic roadmap

> Updated 2026-10-23. Release line: **v1.2.0**. Changes after the tag remain
> **Unreleased** until the next release.

chematic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap favors stable contracts
and reproducible evidence over feature-count parity.

## Current position

- Stable selected paths cover SMILES/SMARTS, descriptors, fingerprints,
  search, and common MOL/SDF workflows through one Rust core.
- v1.2.0 public wheels passed the hosted three-platform runtime smoke and the
  11-job chemistry gate. Later source fixes are not public-package evidence.
- Current source matches the pinned RDKit standardization corpus on Normalize
  4,069/4,069, Cleanup 4,057/4,069, and Reionize 4,069/4,069 rows.
- MCS input-order and resource-limit fixes are in Unreleased source. They still
  need a new independent holdout and released-artifact verification.
- The R17-R25 refactoring program is complete. CI now holds the production
  debt ratchet at eight file and 50 function allowances.
- 3D/MMFF94, pKa/ADMET screening, IUPAC naming, broad Markush/polymer support,
  and rich CDXML editing remain experimental or bounded.

Exact denominators belong in [validation](docs/validation.md), dated evidence
in the [benchmark index](benchmarks/README.md), and unfinished work in the
[open-work ledger](docs/roadmap-open-work.md).

## Priority order

### 1. Release current source with complete evidence

- Run changed Rust, Python, WASM/Node, compatibility, security, and package
  gates from the exact release commit.
- Re-run the changed chemistry lanes from exact PyPI, npm, and crates.io
  artifacts. Keep source, package, browser, and operating-system evidence
  separate.
- Freeze a new MCS holdout only after the implementation is fixed; previously
  inspected rows remain historical evidence.

**Exit:** versions, hashes, artifacts, corpora, options, and terminal outcomes
are recorded, and all release checks pass on the tagged commit.

### 2. Close correctness residuals

- Classify the remaining 12 Cleanup rows as implementation defects,
  representation differences, unsupported chemistry, or comparator policy.
- Fix the one chematic CTfile readback failure and reduce ten stereo-only rows;
  keep 16 RDKit oracle failures outside the chematic failure count.
- Extend reaction, SMARTS, R-group, and file-interchange tests with stratified
  real-world inputs while preserving typed refusals and loss reports.

**Exit:** every residual is an exact match, documented policy difference,
typed unsupported case, or reproducible bug with an owner.

### 3. Prove equivalent-work performance

- Require exact or toleranced output agreement before timing.
- Separate parse-inclusive, first-use, prepared/reuse, startup, RSS, and batch
  lanes. Do not rank operations with unequal output or skipped work.
- Re-run current ASA, TPSA, chiral Morgan, QED, reaction, and 3D findings on a
  second host and exact public artifacts.

**Exit:** each performance claim uses balanced repeated blocks, paired
confidence bounds, equivalent work, exact artifacts, and explicit host scope.

### 4. Expand workflows and 3D quality

- Compare graph, charge, isotope, bond order, stereo, coordinates, properties,
  and metadata before file-I/O speed.
- Measure 3D typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time, and memory separately on unused data.
- Keep browser workflows local-first, cancellable, resource-bounded, and
  explicit about failure.

**Exit:** supported rows preserve declared semantics or return a typed
loss/refusal; no known wrong-confident result remains in the declared domain.

## Product phases

| Phase | Scope | State |
|---|---|---|
| P0 | Release integrity, typed failures, evidence provenance | Active |
| P1 | Parsing, reactions, SMARTS, descriptors, interchange | Maintenance |
| P2 | Fingerprints, search, stereo, identity | Active |
| P3 | 2D/3D, force fields, conformer quality | Experimental |
| P4 | Python, WASM/Node, CLI, MCP, browser products | Active |
| P5 | Operational hardening and cross-platform publication | Active |
| P6 | Advanced chemistry breadth | Selective |

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

## Completion rules

- Record versions, hashes, corpus, options, host, and failure policy.
- Count success, failure, refusal, unsupported, and skipped inputs separately.
- A typed refusal is safer than a wrong answer, but it is not a match.
- A source pass is not a published-package result.
- Performance requires equivalent work and acceptable output.
- Coverage, availability, and feature presence are not chemistry validation.
- Sealed data is evaluated once; inspected rows cannot become sealed again.

## Detailed references

- [Open-work ledger](docs/roadmap-open-work.md)
- [Validation report](docs/validation.md)
- [Benchmark method](docs/benchmark.md)
- [Compatibility scope](docs/compatibility-scope.md)
- [Format capability matrix](docs/format-capabilities.md)
- [Refactoring closeout](docs/refactoring-plan.md)
- [Historical roadmap snapshots](docs/archive/)
