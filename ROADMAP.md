# chematic roadmap

> Updated 2026-10-18. Release line: **v1.1.7**. Work merged after the
> release remains **Unreleased** until the next tag.

chematic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap prioritizes stable
contracts and reproducible evidence over feature-count parity.

## Current position

- Stable selected paths cover SMILES/SMARTS, descriptors, fingerprints,
  search, and common MOL/SDF workflows through one Rust core.
- v1.1.7 adds checked scientific-file conversion and export, custom
  fingerprint and force-field inputs, source-preserving SDF/PDB/MOL2
  workflows, and typed alignment, information-theory, and MIF helpers.
- Issues #734 and #754 are closed, but their reaction, SMARTS, autoconf, and
  fuzz corpora remain regression gates. Closure is not universal parity.
- The bounded workflow gates cover 80 exact reaction rows plus three jointly
  invalid inputs, 39/39 largest-fragment decisions, and 16/16 labelled
  R-group rows. Broader standardization and R-group chemistry remain open.
- The v1.1.5 published-wheel chemistry workflow completed all 11 Linux,
  macOS, Windows, RDKit-comparison, A6, and autoconf jobs. Its hosted run is
  package evidence; a durable checked-in v1.1.5 packet is still required.
- The broadest published-package RDKit performance packet remains v1.0.42.
  Source or exploratory timings are not relabelled as package evidence.
- The processed Linux Codecov lane reports 96.93% Rust line coverage.
  Coverage proves execution, not chemical correctness.
- Current-source RDKit 2026.09.1 classification covers every Normalize,
  Reionize, Cleanup and CTfile stress-row residual. The remaining differences
  are bounded support gaps or declared policy/representation differences, not
  untracked matches.
- Post-candidate NCI/PubChem MCS, SDF-attribute and stereo/3D holdouts are
  frozen under a hash-checked provenance manifest. They are unscored
  evaluation data and must not be tuned on.
- 3D, MMFF94, pKa/ADMET screening, IUPAC naming, broad Markush/polymer
  expansion, and rich CDXML editing remain experimental or bounded.

Exact denominators belong in [validation](docs/validation.md), dated results
in the [benchmark index](benchmarks/README.md), and unresolved work in the
[open-work ledger](docs/roadmap-open-work.md).

## Priority order

### Next candidate — immediate sequence

1. **Evaluate the frozen holdouts once.** Run MCS comparison, SDF attributes
   and stereo/3D quality from a frozen release candidate; report every row and
   do not change that candidate after inspecting results.
2. **Implement the classified standardization gaps by value.** Prioritize the
   36 missing Reionize assignments, aromatic N-oxide handling and common
   functional groups. Keep 883 independent-ionization and 26 metal-stereo
   policy differences explicit instead of copying them silently.
3. **Close the remaining CTfile rows.** Fix the one current chematic readback
   failure and reduce ten stereo-only rows without counting 16 RDKit oracle
   failures against the writer.
4. **Finish two-host performance evidence.** Repeated QED caching is complete;
   retain first-use output identity and use the hosted Linux 20-block lane to
   adjudicate chiral Morgan, Labute ASA and QED alongside macOS.

**Exit:** the release candidate passes all changed-crate/binding compatibility
gates, evaluates the frozen rows exactly once, and has separate two-host
first-use, reuse, parse-inclusive and RSS evidence.

### P0 — Maintain evidence from exact v1.1.x artifacts

1. Preserve the successful v1.1.5 published-wheel chemistry run and check in
   its versions, hashes, counts, and failure accounting as a durable packet.
2. Run reaction, SMARTS, fingerprint, format, browser, and A6 gates from exact
   PyPI, npm, and crates.io artifacts when a claim depends on that channel.
3. Record artifact and corpus hashes, comparator builds, options, platform,
   and failure accounting.
4. Keep Linux, macOS, Windows, browser, and native results separate.

**Exit:** every public claim points to a reproducible artifact record.

### P0 — Protect API and chemistry compatibility

1. Run `cargo-semver-checks` against v1.1.7 for every changed publishable crate
   and retain focused downstream compile fixtures where Rust semver tooling
   cannot express the contract.
2. Treat public enum variants and struct fields as semver-significant.
3. Keep RDKit 2026.03.6 and 2026.09.1 ring-count differences versioned rather
   than selecting one oracle as universally correct.
4. Preserve typed refusals and explicit loss reports across bindings.
5. Reject Python removals/new required arguments/dependency narrowing and WASM
   TypeScript declaration removals before merge.

**Exit:** supported downstream fixtures compile, chemistry contracts are
versioned, and unsupported inputs never become guessed successes.

### P0 — Close equivalent-work performance deficits

1. Require exact or toleranced output agreement before timing.
2. Separate parse-inclusive, first-use, prepared/hot, startup, and memory
   lanes.
3. Prioritize reproducible deficits in QED/LogP first-use, ring perception,
   Labute ASA/TPSA, chiral Morgan, and reaction application.
4. Reject rows with comparator errors, skipped work, or asymmetric APIs.

**Exit:** a declared win uses at least 21 alternating blocks, wins every block,
and has a paired 95% lower bound above 1.0 on two recorded hosts.

### P0/A4 — Expand semantic file interchange

1. Extend the V3000, MOL2, CML, CIF, SDF, and CDXML comparison to stratified
   corpora.
2. Compare graph, charge, isotope, bond order, stereo, coordinates,
   properties, and format metadata before speed.
3. Preserve batch accounting:
   `input = success + refusal + unsupported + invalid + internal_error`.
4. Keep cold start, parse, write, round trip, peak RSS, and large-file
   throughput separate.

**Exit:** every supported cell round-trips under its declared observer or
returns a typed loss/refusal. See the [Open Babel plan](docs/openbabel-interchange-plan.md).

### P1 — Maintain reaction, SMARTS, standardization, and R-group contracts

1. Keep the 83-row reaction gate, 310,000-cell SMARTS gate, autoconf probes,
   BioTransformer-derived cases, and property-fuzz minimizers in CI.
2. Test product graph, atom origin, maps, explicit H, stereo, and rejection
   diagnostics independently.
3. Expand standardization and labelled R-group corpora with atom-order
   invariant fixtures and typed refusals.
4. Add second-host and exact-package reruns before promoting source timing
   results against RDKit or Indigo.

**Exit:** supported rows are deterministic across Rust, Python, and WASM;
unsupported chemistry remains observable.

### P1/A2/A6 — Maintain stereo, identity, and 3D quality

1. Retain permutation, spelling, round-trip, CIP-abstention, MOL stereo, and
   stable-key gates.
2. Measure typing, energy, gradient, convergence, stereo, clash, conformer
   quality, time, and memory separately.
3. Treat the 265-row A6 set and inspected holdouts v1–v4 as development data.
4. Diagnose remaining minimization and basin-diversity tails without tuning
   on inspected holdouts; freeze a candidate before selecting new data.
5. Promote prepared force fields and parallel conformer generation only after
   two-host, 21-block, peak-memory, exact-artifact, and no-regression gates.

**Exit:** no known wrong-confident result remains in the declared domain. A
speed claim also requires equal acceptable quality and published artifacts.

### P2 — Product and ecosystem work

- Keep browser workflows local-first, cancellable, resource-bounded, and
  explicit about failure.
- Maintain shared semantics across Python, Rust, Node/WASM, CLI, and MCP
  without forcing identical host-language API shapes.
- Add requested workflow features only with a documented support boundary,
  fixtures, and maintenance owner.

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
- Performance requires equivalent work and equivalent acceptable output.
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
