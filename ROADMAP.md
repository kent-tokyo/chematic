# chematic roadmap

> Updated 2026-10-23. Release line: **v1.2.1**. Hosted package evidence keeps
> its original release label until the tagged workflows finish.

chematic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap favors reliable contracts
and reproducible evidence over feature-count parity.

## Current position

- Stable selected paths cover SMILES/SMARTS, descriptors, fingerprints,
  search, and common MOL/SDF workflows through one Rust core.
- v1.2.1 is the release line. v1.2.0 public wheels passed the hosted
  three-platform runtime smoke and 11-job chemistry workflow; those results do
  not become v1.2.1 package evidence.
- Current source matches the pinned RDKit standardization corpus on Normalize
  4,069/4,069, Cleanup 4,057/4,069, and Reionize 4,069/4,069 comparable rows.
- v1.2.1 fixes bounded MCS order dependence and resource-limit reporting. Full
  FMCS breadth and released-artifact holdout verification remain open.
- 3D/MMFF94, pKa/ADMET screening, IUPAC naming, broad Markush/polymer support,
  and rich CDXML editing remain experimental or bounded.

Exact counts belong in [validation](docs/reference/validation.md), dated
evidence in the [benchmark index](benchmarks/README.md), and unresolved details
in the [open-work ledger](docs/reference/roadmap-open-work.md).

## Priority order

### P0 — Release integrity

- Run changed Rust, Python, WASM/Node, compatibility, security, and package
  checks from the exact release commit.
- Re-run affected chemistry lanes from exact PyPI, npm, and crates.io
  artifacts. Keep source, package, browser, and OS evidence separate.
- Compare every public Rust crate with the preceding release.

**Exit:** versions, hashes, artifacts, corpora, options, and terminal outcomes
are recorded, and release checks pass on the tagged commit.

### P1 — Correctness residuals

- Classify the remaining Cleanup rows as defects, representation differences,
  unsupported chemistry, or comparator policy.
- Resolve the remaining CTfile readback and stereo-only differences.
- Extend reaction, SMARTS, R-group, and interchange tests with stratified
  real-world inputs while preserving typed refusals and loss reports.

**Exit:** every residual is an exact match, a documented policy difference, a
typed unsupported case, or a reproducible bug with an owner.

### P2 — Equivalent-work performance

- Require exact or toleranced output agreement before timing.
- Separate parse-inclusive, first-use, prepared/reuse, startup, RSS, and batch
  lanes.
- Re-run ASA, TPSA, chiral Morgan, QED, reaction, and 3D comparisons on a
  second host and exact public artifacts.

**Exit:** promoted claims use balanced repeated blocks, paired confidence
bounds, equivalent work, and explicit artifact/host scope.

### P3 — Workflow and 3D quality

- Check graph, charge, isotope, bond order, stereo, coordinates, properties,
  and metadata before ranking file-I/O speed.
- Measure 3D typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time, and memory independently on unused data.
- Keep browser workflows local-first, cancellable, resource-bounded, and
  explicit about failure.

**Exit:** supported rows preserve declared semantics or return a typed
loss/refusal; no known wrong-confident result remains in the declared domain.

## Product phases

| Phase | Focus | State |
|---|---|---|
| P0 | Release integrity and evidence | Active |
| P1 | Parsing, reactions, SMARTS, descriptors, interchange | Maintenance |
| P2 | Fingerprints, search, stereo, identity | Active |
| P3 | 2D/3D, force fields, conformer quality | Experimental |
| P4 | Bindings, CLI, MCP, browser products | Active |
| P5 | Cross-platform publication and hardening | Active |
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

## Evidence rules

- Record versions, hashes, corpus, options, host, and failure policy.
- Count success, failure, refusal, unsupported, and skipped inputs separately.
- A typed refusal is safer than a wrong answer, but it is not a match.
- A source pass is not a published-package result.
- Performance requires equivalent work and acceptable output.
- Coverage and feature presence are not chemistry validation.
- Sealed data is evaluated once; inspected rows cannot become sealed again.

## Detailed references

- [Open-work ledger](docs/reference/roadmap-open-work.md)
- [Validation report](docs/reference/validation.md)
- [Benchmark method](docs/reference/benchmark.md)
- [Compatibility scope](docs/reference/compatibility-scope.md)
- [Format capability matrix](docs/reference/format-capabilities.md)
- [Refactoring guardrails](docs/plans/refactoring-plan.md)
- [Historical snapshots](docs/archive/)
