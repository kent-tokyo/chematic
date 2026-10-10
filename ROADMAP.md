# chematic roadmap

> Updated 2026-10-10. Release line: **v1.0.42**. Source, published-package, and
> public-channel evidence are tracked separately.

CheMatic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap favors explicit support
boundaries and reproducible evidence over feature-count parity.

## Current position

- **Core:** selected SMILES/SMARTS, descriptor, fingerprint, search, and
  MOL/SDF paths are stable and share one Rust implementation across bindings.
- **Compatibility:** the v1.0.41 source comparison covers 161 operations on
  ChEMBL 5k and RDKit.js 10k. It is not a claim of universal or
  published-package parity.
- **Performance:** a clean source candidate beats RDKit 2026.03.1 and
  COSMolKit 0.5.0rc15 for the seven recorded parse-inclusive pipelines on one
  Apple-silicon host. MMFF also wins prepared; five other prepared-operation
  lanes and public-artifact replication remain open.
- **Reactions:** the published 83-row gate is 80 exact graph/origin/map rows
  and three inputs invalid in both engines. Broader SMIRKS behavior remains
  intentionally bounded.
- **3D:** A6 now passes its 265-row geometry, stereo, and clash gate on Linux,
  macOS, and Windows CI. Force-field and conformer quality remain experimental;
  quality evidence takes precedence over speed claims.
- **Release channels:** v1.0.42 is verified through the recorded GitHub, PyPI,
  crates.io, npm, docs.rs, and Pages channels. Availability alone is not
  chemistry validation.

See [validation](docs/validation.md) for denominators and
[benchmark records](benchmarks/README.md) for exact versions, hashes, and
commands.

## Priority order

### 1. P0/A1 — Close the v1.0.42 correctness and API regressions

- The source candidate fixes the four false RDKit-compatible ring counts caused
  by treating every single SMILES closure as a graph cycle. Keep coordination,
  disconnected-component closure, and mutation/cache regressions in Rust,
  Python, and WASM.
- Classify the v1.0.42 explicit-hydrogen reaction residuals by deletion,
  bond-breaking, radical, and re-sanitization causes. The first confirmed
  radical-producing defect is fixed in source; rerun the complete pinned
  corpus before changing the published baseline.
- Preserve the Rust reaction report's existing struct-literal shape. Add
  diagnostics through separate opt-in report types rather than public fields
  on an existing struct.
- Prefer checked descriptor and legacy-MMFF APIs where an unsupported result
  could otherwise be mistaken for `0.0`; keep legacy calls only for source
  compatibility.

**Exit:** the four ring regressions pass in every exposed binding, no confirmed
explicit-H implementation defect remains in the pinned corpus, old Rust report
struct literals compile, and unsupported calculations have a checked path.

### 2. P0 — Refresh public-artifact evidence

- The v1.0.42 macOS arm64 CPython 3.13 wheel is now hash-pinned and measured
  for the chemistry packet and 265-row A6 coordinates. Its explicit-H reaction
  counts are retained as historical v1.0.42 behavior, not silently replaced by
  the source candidate.
- Run the same expected-value packet on the published Linux, macOS, and Windows
  artifacts. Keep platform output, comparator execution, and source-candidate
  results separate.
- Rerun chemistry and performance gates from exact PyPI/npm/crates.io artifacts
  rather than carrying source results forward.

**Exit:** every published claim names its artifact hash, comparator, corpus,
operation, options, and failure policy.

### 3. P0/A3 — Complete the equivalent-work performance gate

- Preserve the existing RDKit agreement denominators before accepting a speed
  change; a typed refusal or changed result is not a performance win.
- Keep parsing, first use, prepared/hot calls, 3D preparation, and memory as
  separate lanes. Do not move costly perception into parsing only to improve a
  downstream benchmark.
- Close the five remaining prepared-operation deficits: ring count, TPSA,
  Labute ASA, chiral Morgan, and amide/amine reaction.
- The first follow-up removed redundant TPSA/Labute allocations and made the
  ordinary SSSR count O(number of SMILES ring closures), but did not close the
  RDKit symmetrized-ring or isolated descriptor gates. Continue from measured
  profiles rather than shifting their cost into parsing.
- Rerun from exact PyPI/npm/crates.io artifacts and on a second host.

**Exit:** for every declared equivalent-output lane, 21 or more alternating
blocks, all blocks faster, and the paired 95% speedup lower bound above 1.0
against both RDKit and COSMolKit. Published-package claims require exact
artifact hashes.

### 4. P1/A4 — Close the remaining SMARTS/SMIRKS contract

Issues [#734](https://github.com/kent-tokyo/chematic/issues/734) and
[#754](https://github.com/kent-tokyo/chematic/issues/754) remain open.

- Keep the 83-row reaction gate and the original 57 fixtures unchanged.
- Classify every fuzz difference as exact, typed refusal, unsupported dialect,
  invalid in both engines, or implementation defect.
- Submit the minimized upstream probes and rerun the published macOS/Windows
  package lanes.
- Do not count refusals or jointly invalid inputs as matches.

**Exit:** no wrong-confident supported result in the pinned corpora, all rows
accounted for, and package evidence recorded separately from source evidence.

### 5. P2/A6 — Maintain stereo, identity, and 3D quality

- Preserve atom-order, spelling, and round-trip invariance gates for CIP/E/Z.
- Keep lossy MOL output observable through report/strict APIs.
- Measure MMFF/UFF typing, energy, convergence, stereo, clash, and conformer
  quality independently.
- Keep canonical SMILES distinct from the fail-closed stable identity key.

**Exit:** no silent information loss and no wrong-confident stereo result in
the declared domain.

### 6. External — Rebaseline new RDKit distributions

Pin the exact distributed Python, npm/WASM, and native artifacts before a
rebaseline. Preserve old results as historical records; do not substitute a
source build for an unavailable package.

## Product phases

| Phase | Area | Status |
|---|---|---|
| P0 | Core chemistry and reproducible comparison | Active: prepared-operation and artifact performance gates |
| P1 | Parsing, files, reactions, SMARTS | Active: #734, #754 |
| P2 | Stereo, identity, canonicalization | Active maintenance |
| P3 | Browser, Node, Python, agents | RDKit WASM parity implemented; maintain cross-binding gates |
| P4 | Medicinal chemistry workflows | Bounded implementation |
| P5 | 3D, force fields, materials | Experimental; quality-gated |
| P6 | Release operations and ecosystem trust | Automated with recurring review |

Phase numbers identify areas, not a strict delivery sequence. Security or
silent-corruption regressions always take priority.

## Accuracy packages

| Package | State |
|---|---|
| A0 Evaluation contract | Complete |
| A1 Perception and descriptors | Open |
| A2 Stereo and identity | Active |
| A3 Fingerprints and retrieval | Active: output-preserving performance |
| A4 Workflows and interchange | Active |
| A5 Independent adjudication | External |
| A6 3D and force fields | Active |

## Completion rules

Every completed item must include:

1. a declared support domain and typed failure policy;
2. regression tests for the corrected behavior;
3. affected binding and package checks;
4. committed evidence with versions, hashes, and denominators;
5. synchronized README, CHANGELOG, validation, and release metadata.

## References

- [Open-work ledger](docs/roadmap-open-work.md)
- [Compatibility scope](docs/compatibility-scope.md)
- [Validation report](docs/validation.md)
- [Benchmark methodology](docs/benchmark.md)
- [Accuracy plan](docs/rdkit-accuracy-plan.md)
- [Trust Release rules](docs/trust-release-plan.md)
- [Release history](CHANGELOG.md)
