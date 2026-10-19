# chematic roadmap

> Updated 2026-10-19. Release line: **v1.1.8**. Work merged after the
> release remains **Unreleased** until the next tag.

chematic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap prioritizes stable
contracts and reproducible evidence over feature-count parity.

## Current position

- Stable selected paths cover SMILES/SMARTS, descriptors, fingerprints,
  search, and common MOL/SDF workflows through one Rust core.
- The public v1.1.8 macOS arm64 wheel reproduces the 10,000-row Labute ASA
  result on one host. This is one-operation evidence, not a cross-platform or
  toolkit-wide superiority claim.
- Broad hosted package chemistry is established for v1.1.7; the durable
  checked-in cross-platform packet remains v1.1.3; the broadest published
  RDKit package benchmark remains v1.0.42.
- Reaction, SMARTS, fragment, and labelled R-group results are bounded by
  declared corpora. Standardization and CTfile residuals are classified but
  not fully closed.
- Hash-pinned MCS, SDF-attribute, and stereo/3D holdouts are frozen and
  unscored. They must not be used for tuning before their one-time evaluation.
- 3D/MMFF94, pKa/ADMET screening, IUPAC naming, broad Markush/polymer
  expansion, and rich CDXML editing remain experimental or bounded.

Exact denominators belong in [validation](docs/validation.md), dated results
in the [benchmark index](benchmarks/README.md), and unresolved work in the
[open-work ledger](docs/roadmap-open-work.md).

## Priority order

### 1. Finish release-grade evidence

- Evaluate the frozen holdouts once from a fixed candidate.
- Re-run relevant chemistry gates from exact PyPI, npm, and crates.io
  artifacts; keep source, package, browser, and platform lanes separate.
- Record versions, hashes, corpora, options, host, and all failure classes.
- Extend semver, downstream-consumer, Python stub, and WASM declaration gates
  to every changed public surface.

**Exit:** every public claim points to a reproducible artifact record, and the
frozen holdout candidate is unchanged after results are inspected.

### 2. Close correctness residuals

- Implement high-value standardization gaps without silently adopting
  comparator-specific policy differences.
- Fix chematic CTfile readback failures and reduce stereo-only differences;
  keep comparator oracle failures separate.
- Maintain reaction, SMARTS, autoconf, fuzz, CIP-abstention, and identity
  regression gates across Rust, Python, and WASM.
- Preserve typed refusals and loss reports instead of guessed successes.

**Exit:** each residual is an exact match, a documented policy difference, a
typed unsupported case, or a reproducible bug with an owner.

### 3. Prove equivalent-work performance

- Require exact or toleranced output agreement before timing.
- Separate parse-inclusive, first-use, prepared/hot, startup, RSS, and batch
  lanes; reject skipped or asymmetric work.
- Re-run current ASA, TPSA, chiral Morgan, QED, and reaction findings on a
  second host and exact public artifacts.
- Do not rank operations whose outputs differ.

**Exit:** a performance claim uses balanced repeated blocks, paired confidence
bounds, equivalent outputs, exact artifacts, and an explicit host scope.

### 4. Expand workflows and file interchange

- Extend V3000, MOL2, CML, CIF, SDF, CDXML, standardization, reaction, and
  R-group tests to stratified real-world corpora.
- Compare graph, charge, isotope, bond order, stereo, coordinates, properties,
  and metadata before speed.
- Preserve batch accounting:
  `input = success + refusal + unsupported + invalid + internal_error`.

**Exit:** every supported row preserves its declared semantics or returns a
typed loss/refusal. See the [Open Babel plan](docs/openbabel-interchange-plan.md).

### 5. Improve 3D and product surfaces

- Measure typing, energy, gradient, convergence, stereo, clash, conformer
  quality, time, and memory separately on unused data.
- Keep browser workflows local-first, cancellable, resource-bounded, and
  explicit about failure.
- Add features only with a support boundary, fixtures, and maintenance owner.

**Exit:** no known wrong-confident result remains in the declared domain; any
speed claim also meets the published quality gate.

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
- [Open Babel interchange plan](docs/openbabel-interchange-plan.md)
- [Refactoring plan](docs/refactoring-plan.md)
- [Historical roadmap snapshots](docs/archive/)
