# chematic roadmap

> Updated 2026-10-21. Release line: **v1.2.0**. Work merged after the
> release remains **Unreleased** until the next tag.

chematic is a safe, typed, local-first chemistry kernel for Rust, Python,
Node, WebAssembly, and agent workflows. The roadmap prioritizes stable
contracts and reproducible evidence over feature-count parity.

## Current position

- Stable selected paths cover SMILES/SMARTS, descriptors, fingerprints,
  search, and common MOL/SDF workflows through one Rust core.
- The v1.2.0 candidate adds checked format and geometry surfaces, stronger
  fail-closed chemistry paths, and output-preserving module refactors. It must
  earn its own published-wheel evidence after tagging; older package, browser,
  holdout and source records remain separate evidence channels.
- Reaction, SMARTS, fragment and labelled R-group results remain bounded by
  declared corpora. Cleanup now matches 3,933/4,069 pinned stress rows; other
  standardization and CTfile residuals remain classified but not fully closed.
- The frozen post-candidate NCI/PubChem holdouts have been scored once against
  immutable v1.1.8. SDF graph/title/property/stereo columns are 64/64; MCS
  atom/bond counts are 31/64; 3D embed status/stereo graph are 96/96,
  conformer count is 95/96 and strict pair-distance coordinates are 4/96.
  These inspected rows are historical evaluation evidence, not tuning data.
- 3D, MMFF94, pKa/ADMET screening, IUPAC naming, broad Markush/polymer
  expansion, and rich CDXML editing remain experimental or bounded.

Exact denominators belong in [validation](docs/validation.md), dated results
in the [benchmark index](benchmarks/README.md), and unresolved work in the
[open-work ledger](docs/roadmap-open-work.md).

## Priority order

### 1. Finish release-grade evidence

1. **Freeze a new independent MCS evaluation after the correctness fix.** The
   source candidate now returns three atoms/two bonds for `C1CC1` versus `CCC`
   in either input order and reports search-node limits explicitly. Keep the
   inspected holdout-v1 result historical, develop on newly selected rows, and
   freeze a disjoint atom/bond/status corpus only after the candidate is fixed.
2. **Implement classified standardization gaps by value.** Cleanup is now
   3,933/4,069 after ten further fixes with zero regressions. Prioritize the
   remaining 26 bounded acid/base rows and 101 metal disconnection/charge/
   bond/stereo rows. Keep independent-ionization and representation/policy
   differences explicit instead of copying them silently.
3. **Close the remaining CTfile rows.** Fix the one current chematic readback
   failure and reduce ten stereo-only rows without counting 16 RDKit oracle
   failures against the writer.
4. **Finish two-host performance evidence.** Repeated QED caching is complete;
   retain first-use output identity and use the hosted Linux 20-block lane to
   adjudicate chiral Morgan, Labute ASA and QED alongside macOS.
   The [terminal-TPSA source follow-up](benchmarks/2026-10-18-terminal-tpsa.md)
   improves both fixed first-use windows. The
   [Python ASA/ring follow-up](benchmarks/2026-10-18-python-first-use-ring-asa.md)
   closes the measured RDKit first-use ASA and ring-count deficits on one
   macOS source lane. The later
   [10k ASA follow-up](benchmarks/2026-10-19-python-labute-asa-cosmolkit.md)
   also clears the fixed COSMolKit ASA lane in 20/20 blocks while preserving
   all values. First-use TPSA plus one chiral-Morgan window remain open.
   COSMolKit's six unequal ring counts stay unranked instead of being treated
   as a speed deficit. Keep source evidence separate from public artifacts and
   state the single-macOS-host limit explicitly.

**Exit:** the release candidate passes all changed-crate/binding compatibility
gates, preserves the immutable frozen score, and has separate two-host
first-use, reuse, parse-inclusive and RSS evidence.

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
