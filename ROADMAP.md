# chematic roadmap

> Updated 2026-10-03. Release line: **v1.0.31**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

| Gate | Verified on published v1.0.30 artifacts | Remaining exit |
|---|---|---|
| P0.1 outputs | Python, npm and pinned Rust artifacts reran 10,000 chemistry inputs, 310,000 SMARTS cells and 57 legacy reactions. Python/Rust covered 63 operation outputs; only HBA and its bundle changed from v1.0.29. Published HBA: 5,000/5,000 vs pinned RDKit 2026.03.6. | **Audit complete, parity not complete:** 200 SMARTS failures, five typed CIP abstentions and four npm-unexposed operations remain. See the [acceptance policy](benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md). |
| P0.2 speed | The published Python 63-operation matrix has 20 paired blocks; 20 operations pass exact-output and confidence-interval gates on one host. Published npm/WASM Morgan also has Node and Chromium lanes; Ubuntu 24.04 Chromium, Firefox and WebKit each have 250/250 bit-exact direct/prepared rows and favorable 20-block speed intervals (smallest lower bound 1.63×). | Equal-work perception and library-memory boundaries, other operations/bindings and more corpora. Process RSS is not library allocation; changed HBA output cannot earn an equivalent-output speed win. |
| P1 reactions | Original 57 fixtures match. Published v1.0.30 had 73/83 graph+origin+template-map matches. The checked v1.0.31-source profile has 76 exact on all three axes, three typed unsupported, one diagnosed refusal and three jointly invalid rows; see the [source audit](benchmarks/2026-10-03-source-reaction-83-checked-graph-origin-map.md). | Publish/rebuild the source fix and rerun applicable Rust/Python/npm provenance gates. Typed refusal and invalid input are not exact matches; no general SMIRKS parity claim. |
| P1 SMARTS | Published v1.0.30 has 200 differing **match-set** cells / 310,000, but only 43 hit/no-hit differences. The opt-in source profile corrects 154 match sets / 37 Booleans, retains 29 / 6 wrong-confident cells and refuses 18 cells (including formerly correct ones); see the [source audit](benchmarks/2026-10-03-source-smarts-310k-profiles.md). | Resolve remaining Boolean and `[k6]`/ring-set residuals, adjudicate six charged-polycycle refusal rows, then expose and measure a deliberate opt-in published profile. Native SSSR remains unchanged. |
| P2 A6 quality | The [published v1.0.31 MMFF94 stereo-safe rerun](benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md) has 265/265 independently sound, stereo-clean and gross-clash-free outputs, but only 100/265 force-field convergence. On the separate same-coordinate lane, 262/262 comparable total energies are within 1 kcal/mol; the historical two >5 kcal/mol residuals no longer appear. A separate v1.0.25 census still has 2,908 heavy-atom type differences. | Keep typing, convergence, per-term energy and independent conformer-quality gates separate; account for two embedding failures and one declared unsupported row. Current geometry/stereo and bounded total-energy results do not complete A6. |

The [v1.0.30 published-artifact packet](benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md),
[benchmark index](benchmarks/README.md), and [validation report](docs/validation.md)
hold versions, hashes, denominators and residuals. v1.0.31 is published on
all six tracked channels, but its WASM formula and E/Z JSON fixes still need
package-output checks. Exposed data are
not sealed evidence; A0's frozen descriptor result and Experimental 3D/MMFF94
have separate gates.

## Priority order and acceptance gates

1. **P1 / A4 — Close the reaction 83-row residuals by outcome, not by forcing
   83 exact matches.** First reconcile the ten named rows against product graph,
   distinct-product multiplicity, atom origins and template-map labels; retain
   raw embedding counts as a separate diagnostic. Adjudicate the one map-only
   row and the symmetric-origin row before changing code. Convert genuine
   out-of-contract semantics to typed unsupported; fix supported-domain graph
   or provenance defects with row-specific regression tests. Then rerun the
   checked-source profile on rebuilt published Rust/Python/npm artifacts,
   including the npm E/Z JSON path and three new-product-map supplement.
   **Exit:** all 83 rows accounted for on each exposed axis, zero wrong
   confident supported-domain results, typed reasons for refusals, and no
   regression in the original 57. The three jointly invalid rows remain
   classified as such; reaction yield/selectivity is outside this gate.
2. **P1 / A4 — Reduce the 200 SMARTS match-set residuals by frequency.** Start
   with `[R2]` (64 cells), `[R3]` (62), `[R1]` (59), then `[k6]` (10),
   `[k5]` (1) and the four other organometallic queries (one each).
   The first three cover 185/200 cells, including two organometallic cells.
   For each family, record the independent graph/ring interpretation, RDKit
   behavior, supported domain and a specification-difference / unsupported /
   implementation-defect verdict before editing. Keep native SSSR unchanged;
   use an explicit opt-in RDKit-style profile only where its semantics are
   justified. **Exit:** rerun the full pinned 310,000-cell match-set matrix on
   affected published bindings, report both match-set and Boolean deltas by
   query and row, add focused
   regression tests, and permit no new wrong confident matches or silent
   acceptance of unsupported organometallic semantics. This is not a general
   RDKit SMARTS parity claim.
3. **P2 / A6 — Gate 3D quality before speed.** In order, adjudicate the 2,908
   heavy-atom MMFF94 typing differences and parameter/charge coverage; compare
   per-term energies and analytic gradients at identical explicit-H coordinates;
   account for convergence, timeout and cancellation; check stereo retention,
   bond geometry and gross clashes; then assess conformer quality against an
   independent reference on matched seeds and a held-out cohort. Unsupported
   chemistry must return a typed reason, not plausible coordinates. The
   published v1.0.31 geometry/stereo/clash rerun is complete but 165/265 force
   fields remain un-converged. **Exit:** predeclared tolerances and complete
   outcome accounting for each remaining dimension on current published
   artifacts; make no 3D speed or RDKit-equivalence claim before quality exits.
4. **P0.2 — Complete paired speed evidence after quality gates.** Extend
   equivalent-output, alternating-order, ≥20-block measurements beyond the
   finished lanes. Separate parse, perception, prepared/reused calls and memory
   allocations. **Exit:** pinned scripts/raw observations, paired intervals
   and an explicit no-claim outcome whenever outputs or work differ.
5. **P1 / A2 — Preserve CIP safety.** Keep five typed abstentions stable under
   atom-order, spelling and file round trips; adjudicate phosphorus independently
   and define a lone-pair convention before assigning R/S. **Exit:** complete
   outcome accounting and no wrong confident label.
6. **P2 — Test representation limits separately.** Keep CDXML opaque
   preservation distinct from semantic editing, Markush/polymer expansion
   bounded, Standard InChI separate from native identifiers, and canonical
   SMILES separate from the fail-closed stable key. **Exit:** cross-binding
   fixtures with no silent information loss.
7. **External — Rebaseline RDKit only after pinning a new official release.**
   Retain 2026.03.6 evidence and do not tune on sealed data.

Silent corruption or a security regression takes precedence. Source,
published-package and sealed results remain separate evidence classes.

## Accuracy packages

| Package | State | Exit focus |
|---|---|---|
| A0 Evaluation contract | Complete | Preserve frozen/exposed-cohort and failure-accounting rules |
| A1 Perception and descriptors | Open | Descriptor families, aromaticity boundaries, potential centers, independent holdouts |
| A2 Stereo and identity | Active | #632 audit, CIP adjudication, permutation/spelling/file round trips, stable-key scope |
| A3 Fingerprints and retrieval | Open | Compatible-Morgan options, top-k invariants, cross-binding parity |
| A4 Workflows and interchange | Open | SMARTS ring counts, reactions, V3000/query semantics, attachment metadata, batch accounting |
| A5 Independent adjudication | External | Non-maintainer gold data and review without exposed-cohort reuse |
| A6 3D and force fields | Active | Heavy-atom typing, convergence, conformer quality, publication-level speed |

See [the accuracy plan](docs/rdkit-accuracy-plan.md) for acceptance criteria
and [the open-work ledger](docs/roadmap-open-work.md) for the dependency state.

The [Trust Release plan](docs/trust-release-plan.md) tracks the cross-cutting
version, runtime, security and review workstreams without repeating this queue.

## Product phases

| Phase | Area | Status |
|---|---|---|
| P0 | Core chemistry and reproducible comparison | Core stable; v1.0.30 artifact accounting complete, paired-speed work open |
| P1 | Safe parsing and files | Stable bounded paths; malformed-state coverage continues |
| P2 | Stereo, identity, canonicalization | Active: #632 evidence and CIP adjudication |
| P3 | Browser, Node, Python, agents | Selected stable surface; runtime controls remain |
| P4 | Reactions, SMARTS, medicinal chemistry | Bounded implementation; ring-count semantics remain |
| P5 | 3D, force fields, materials | Experimental; A6 gates remain |
| P6 | Release operations and ecosystem trust | Partly automated; review/rebaseline recur |

Phase numbers describe areas, not a strict delivery order.

## Reference documents

- [Open-work ledger](docs/roadmap-open-work.md)
- [Trust Release plan](docs/trust-release-plan.md)
- [Compatibility scope](docs/compatibility-scope.md)
- [Validation report](docs/validation.md)
- [RDKit accuracy plan](docs/rdkit-accuracy-plan.md)
- [Benchmark methodology](docs/benchmark.md)
- [Release history](CHANGELOG.md)
