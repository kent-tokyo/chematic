# chematic roadmap

> Updated 2026-10-10. Release line: **v1.0.42**. Source, published-package, and
> public-channel evidence are tracked separately.

CheMatic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap favors explicit support
boundaries and reproducible evidence over feature-count parity.

## Current position

- **Core:** selected SMILES/SMARTS, descriptor, fingerprint, search, and
  MOL/SDF paths are stable and share one Rust implementation across bindings.
- **Compatibility:** the source comparison covers 161 operations on ChEMBL 5k
  and RDKit.js 10k. It is not a claim of universal or published-package
  parity. Current work is rebaselined to RDKit 2026.09.1 and COSMolKit
  0.5.0rc22; older comparator results remain historical.
- **Performance:** the current macOS arm64 candidate wins every measured
  parse-inclusive pipeline, but does not yet win every prepared operation.
  RDKit-compatible TPSA is statistically tied with RDKit in the pinned
  21-block run; ring count, Labute ASA, chiral Morgan, and amide/amine reaction
  still lose to at least one comparator. Public-artifact replication remains
  open.
- **CDK/Indigo:** a source-only smoke gate now pins CDK 2.13 and Indigo
  1.46.0. CheMatic is exact on all 15 measured operations of the 24-row common
  matrix. Every rankable speed median and paired 95% interval favours
  CheMatic; one high-repetition `[#8]` lane wins 20/21 blocks, so the strict
  every-block exit remains open. This is not broad or package evidence.
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

The target is a bounded, reproducible win over the pinned distributed
artifacts, not a universal claim over every chemistry operation or workload.

- Freeze the primary comparator set at RDKit 2026.09.1 and COSMolKit
  0.5.0rc22 until the acceptance packet closes. A comparator upgrade starts a
  new packet; it does not rewrite historical results.
- Accuracy comes first: on every declared supported row, chematic must equal
  the pinned oracle or return the documented typed refusal. A faster wrong
  result, comparator exception, or skipped row is not a win.
- Maintain two performance contracts: parse-inclusive pipelines and prepared
  operation-only calls. Parsing cannot absorb expensive perception solely to
  make the prepared lane look faster.

- Preserve the existing RDKit agreement denominators before accepting a speed
  change; a typed refusal or changed result is not a performance win.
- Keep parsing, first use, prepared/hot calls, 3D preparation, and memory as
  separate lanes. Do not move costly perception into parsing only to improve a
  downstream benchmark.
- Close the remaining prepared-operation deficits in this order: ring count,
  Labute ASA, chiral Morgan, and amide/amine reaction against COSMolKit. Keep
  TPSA above the RDKit lower-confidence-bound gate while closing its COSMolKit
  deficit.
- The current follow-up caches the shared RDKit-model safety predicate,
  removes a redundant cleanup scan, and avoids Labute bond-order allocation.
  On the pinned 1,000-row/21-block run, all six parse-inclusive pipelines win
  every block. Prepared TPSA is statistically tied with RDKit, and the
  amide/amine reaction wins against RDKit; neither yet beats COSMolKit.
- Fix the benchmark before trusting it: API drift or any per-row comparator
  error invalidates that operation's speed ratio and confidence interval.
- Rerun from exact PyPI/npm/crates.io artifacts and on a second host.

#### CDK 2.13 and Indigo 1.46.0 lane

- Keep the checked-in common-schema adapters independent and public-API-only;
  do not copy LGPL CDK or Apache Indigo implementation code.
- Promote the 24-row smoke corpus to a frozen broad corpus plus adversarial
  metal, isotope, stereo, aromatic, disconnected, and hypervalent strata.
- Separate semantic-common operations from comparator-specific definitions.
  In particular, report the RDKit-model/SSSR ring-count difference without
  treating cubane's basis choice as an accuracy failure.
- Expand exact-output speed gates from parse, formula, canonical stability,
  and six SMARTS queries to MOL/SDF round trip, descriptors, fingerprints,
  reactions, and batch file I/O only where all engines perform equivalent
  work.
- Add cold start and peak RSS, repeat on a second host, then rerun from exact
  published CheMatic artifacts before making a public claim.

**Exit:** CheMatic has no wrong-confident row on the frozen common matrix and
is at least as accurate as both comparators. Every declared speed win passes
21 rotating blocks, wins every block, and has a paired 95% lower bound above
1.0 in both pipeline and prepared lanes; package, host, and memory evidence
are recorded separately.

**Exit:** on the frozen common-operation matrix, chematic is at least as
accurate as each comparator row by row, with no wrong-confident result. For
every equivalent-output performance lane, 21 or more alternating blocks, all
blocks faster, and the paired 95% speedup lower bound above 1.0 against both
RDKit and COSMolKit. Published-package claims require exact artifact hashes,
a second host, and separate peak-memory results.

### 4. P0/A4 — Beat Open Babel on bounded production file I/O

- Use the versioned contract in
  [`validation/openbabel_file_io_contract_v1.json`](validation/openbabel_file_io_contract_v1.json)
  for SDF, MOL V2000/V3000, MOL2, CML, CDXML, PDB, and mmCIF. Do not use raw
  format count as the success metric.
- First close semantic round trips in both directions: graph, charge, isotope,
  bond order, stereo, coordinates, record properties, and format-specific
  metadata. Any loss must be reported or refused; silent loss is a defect.
- Current source progress: V3000, MOL2, CML, and CDXML have a common-observer
  single-fixture semantic gate. The new loss-aware MOL2 record API preserves
  partial charges, atom types, residue fields, status bits, formal charges,
  and opaque sections. Expand this from one fixture to stratified corpora before
  making a broad accuracy claim.
- Preserve `input = success + typed refusal + unsupported + invalid + internal
  error` for batch work. Require zero internal errors, panics, and crashes on
  malformed and resource-limit corpora.
- Compare cold start, parse, write, round trip, and peak RSS as separate lanes.
  Require equivalent outputs, at least 21 alternating blocks, every block
  faster, and a paired 95% speedup lower bound above 1.0 before declaring a
  format win.
- Current source progress: the fresh-process CLI round-trip lane passes the
  strict 21-block rule for all four semantic-gated fixtures. The same-process
  hot-loop gate also passes all 12 parse, write, and round-trip lanes for
  V3000, MOL2, CML, and CDXML. Large-file throughput, peak RSS, broad semantic
  corpora, PDB/mmCIF, and published artifacts remain open.
- Publish source and package evidence separately. A win on the Tier-A profile
  is not a claim to match Open Babel's total format breadth.

**Exit:** every Tier-A format has bidirectional semantic evidence, zero silent
loss, complete typed accounting, bounded malformed-input behavior, and a
reproducible equivalent-work comparison against pinned Open Babel 3.2.1.

See [the comparison contract](docs/openbabel-file-io.md).

### 5. P1/A4 — Maintain the closed SMARTS/SMIRKS contract

Issues [#734](https://github.com/kent-tokyo/chematic/issues/734) and
[#754](https://github.com/kent-tokyo/chematic/issues/754) are closed. Their
corpora remain active regression gates rather than release-completion claims.

- Keep the 83-row reaction gate and the original 57 fixtures unchanged.
- Classify every fuzz difference as exact, typed refusal, unsupported dialect,
  invalid in both engines, or implementation defect.
- Submit the minimized upstream probes and rerun the published macOS/Windows
  package lanes.
- Do not count refusals or jointly invalid inputs as matches.

**Maintenance rule:** no wrong-confident supported result in the pinned
corpora, all rows accounted for, and package evidence recorded separately from
source evidence.

### 6. P2/A6 — Maintain stereo, identity, and 3D quality

- Preserve atom-order, spelling, and round-trip invariance gates for CIP/E/Z.
- Keep lossy MOL output observable through report/strict APIs.
- Measure MMFF/UFF typing, energy, convergence, stereo, clash, and conformer
  quality independently.
- Keep canonical SMILES distinct from the fail-closed stable identity key.

**Exit:** no silent information loss and no wrong-confident stereo result in
the declared domain.

### 7. External — Rebaseline new RDKit distributions

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
