# chematic roadmap

> Updated 2026-10-12. Release line: **v1.1.1**. Source, published packages,
> and historical measurements are kept separate.

CheMatic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap favors reliable contracts
and reproducible evidence over feature-count parity.

## Current position

- Stable selected paths cover SMILES/SMARTS, descriptors, fingerprints,
  search, and common MOL/SDF workflows across one Rust core.
- Issues #734 and #754 are closed. Their reaction, SMARTS, autoconf, and fuzz
  corpora remain regression gates; closure is not a universal SMIRKS-parity
  claim.
- v1.1.1 exposes typed reaction rejection diagnostics across Rust, Python, and
  WASM, and input-aligned `Core`/`R1`/`R2` columns across Rust, Python, and
  WASM/Node for the documented terminal-wildcard subset.
- The broadest published-package RDKit performance packet is still v1.0.42.
  Source results and older package results are not relabeled as v1.1.1
  evidence.
- Rust line coverage on the processed Linux Codecov lane is 96.93%. Coverage
  demonstrates executed code, not chemical correctness or compatibility.
- 3D, MMFF94, pKa/ADMET screening, IUPAC naming, broad Markush/polymer
  expansion, and rich CDXML editing remain experimental or bounded.
- The source candidate has deterministic native parallel conformer generation:
  on the recorded three-molecule lane it is exactly equal to the sequential
  result and 3.54x faster. Its median is 3.87x faster than RDKit 2026.09.1's
  all-thread lane, but latest-oracle quality, independent holdout, memory,
  second-host, and public-artifact gates remain open.

Exact denominators belong in [validation](docs/validation.md) and dated
[benchmark records](benchmarks/README.md). Unresolved work is tracked in the
[open-work ledger](docs/roadmap-open-work.md).

## Priority order

### P0 — Publish v1.1.x evidence from exact artifacts

1. Rebuild the chemistry, reaction, SMARTS, fingerprint, format, and A6 gates
   from exact PyPI, npm, and crates.io artifacts.
2. Record artifact hashes, comparator builds, corpus hashes, options, platform,
   and failure accounting.
3. Keep Linux, macOS, Windows, browser, and native lanes separate. A source
   pass is not a package pass.

**Exit:** every public claim points to a reproducible artifact record; old
measurements remain explicitly historical.

### P0 — Preserve correctness while closing prepared-operation deficits

1. Require exact or toleranced output agreement before timing an operation.
2. Measure parse-inclusive, first-use, prepared/hot, memory, and startup lanes
   separately.
3. Prioritize remaining equivalent-work deficits in ring perception,
   Labute ASA/TPSA, chiral Morgan, and reaction application.
4. Reject benchmark rows with comparator errors, skipped work, or asymmetric
   APIs instead of reporting a speed ratio.

**Exit:** each declared win uses at least 21 alternating blocks, wins every
block, and has a paired 95% lower bound above 1.0 on two recorded hosts.

### P0/A4 — Expand semantic file interchange

1. Extend the bounded Open Babel comparison beyond the initial V3000, MOL2,
   CML, and CDXML fixtures to stratified corpora.
2. Compare graph, charge, isotope, bond order, stereo, coordinates, properties,
   and format-specific metadata before measuring speed.
3. Preserve batch accounting:
   `input = success + refusal + unsupported + invalid + internal_error`.
4. Keep cold start, parse, write, round trip, peak RSS, and large-file
   throughput as separate results.

**Exit:** every supported cell either round-trips under the declared observer
or reports a typed loss/refusal. See the detailed
[Open Babel plan](docs/openbabel-interchange-plan.md).

### P1 — Maintain reaction, SMARTS, and R-group contracts

1. Keep the 83-row reaction gate, 310,000-cell SMARTS gate, autoconf probes,
   BioTransformer-derived cases, and property-fuzz minimizers in CI.
2. Preserve product graph, atom origin, template-map, explicit-H, stereo, and
   rejection-diagnostic contracts independently.
3. Expand R-group coverage only with atom-order-invariant fixtures and typed
   refusals. Multi-core/MCS alignment, tautomer/core enumeration,
   non-terminal R-groups, and GA scoring stay out of the bounded API until
   independently specified.
4. Maintain the pinned Indigo 1.46 workflow gate. The current source packet
   covers 15 single-result reactions, 39 largest-fragment decisions, and 16
   labelled R-group rows. Expand each into a stratified frozen corpus, add a
   second host, split prepared from parse-inclusive work, and rerun exact
   registry artifacts before making a broad Indigo comparison claim.

**Exit:** supported rows are deterministic across Rust, Python, and WASM;
unsupported chemistry never becomes a guessed successful result.

### P1/A2/A6 — Maintain stereo, identity, and 3D quality

1. Retain permutation, spelling, round-trip, CIP-abstention, MOL stereo, and
   canonical stable-key gates.
2. Keep force-field typing, per-term energy, gradient, convergence, stereo,
   clash, and conformer quality as separate measurements.
3. Add new 3D breadth only after failure reasons and quality evidence are
   available on the existing supported domain.
4. Rebaseline the frozen 265-row A6 set against RDKit 2026.09.1, then evaluate
   a previously unused holdout. Report typing, energy, gradient, convergence,
   stereo, clash, and best-of-N quality independently.
5. Promote parallel conformer generation only after the same result is shown
   from exact public artifacts on two hosts, with at least 21 alternating
   blocks, peak memory, and no quality regression. Add macrocycles, metals,
   fullerenes, cages, highly flexible molecules, and failed embeddings to the
   stress lane.

**Exit:** no known wrong-confident result remains in the declared domain;
unresolved cases are observable and typed. A speed claim additionally requires
equal acceptable quality, favorable paired uncertainty, two hosts, and exact
published artifacts.

### P2 — Product and ecosystem work

- Keep browser workflows local-first, bounded, cancellable, and explicit
  about memory/time limits.
- Maintain Python, Rust, Node/WASM, CLI, and MCP contracts from the shared
  core without promising identical API shapes where the host language differs.
- Add requested workflow features only when their support boundary, fixtures,
  and maintenance owner are clear.

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
- A local/source pass is not a published-package result.
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
