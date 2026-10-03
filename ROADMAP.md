# chematic roadmap

> Updated 2026-10-03. Release line: **v1.0.31**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

| Gate | Evidence by artifact and version | Remaining exit |
|---|---|---|
| P0.1 outputs | Python, npm and pinned Rust artifacts reran 10,000 chemistry inputs, 310,000 SMARTS cells and 57 legacy reactions. Python/Rust covered 63 operation outputs; only HBA and its bundle changed from v1.0.29. Published HBA: 5,000/5,000 vs pinned RDKit 2026.03.6. | **Audit complete, parity not complete:** 200 SMARTS failures, five typed CIP abstentions and four npm-unexposed operations remain. See the [acceptance policy](benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md). |
| P0.2 speed | The published Python 63-operation matrix has 20 paired blocks; 20 operations pass exact-output and confidence-interval gates on one host. Published npm/WASM Morgan also has Node and Chromium lanes; Ubuntu 24.04 Chromium, Firefox and WebKit each have 250/250 bit-exact direct/prepared rows and favorable 20-block speed intervals (smallest lower bound 1.63×). | Equal-work perception and library-memory boundaries, other operations/bindings and more corpora. Process RSS is not library allocation; changed HBA output cannot earn an equivalent-output speed win. |
| P1 reactions | Original 57 fixtures match. Published v1.0.30 had 73/83 graph+origin+template-map matches. The checked v1.0.31-source profile has 76 exact on all three axes, three typed unsupported, one diagnosed refusal and three jointly invalid rows. [Post-merge source and binding gates](benchmarks/2026-10-03-reaction-83-postmerge-binding-gate.md) confirm the same graph/status counts in Linux/macOS CI source wheels and 3/3 new-product-map supplement matches in Rust source. | Rebuild/publish and rerun Rust/Python/npm artifact gates. Python/WASM still do not expose product atom origins or template maps. Typed refusal and invalid input are not exact matches; no general SMIRKS parity claim. |
| P1 SMARTS | Published v1.0.30 has 200 differing **match-set** cells / 310,000, but only 43 hit/no-hit differences. The [opt-in source candidate](benchmarks/2026-10-03-source-smarts-organometallic-ring-view.md) corrects 183 match sets and all 43 Booleans, with zero wrong-confident cells and 18 typed refusals; the native matcher is unchanged. The [ring-family diagnosis](benchmarks/2026-10-03-smarts-ring-family-diagnostic.md) finds three source families matching the original RDKit ring sets and three differing; RDKit ring sets themselves varied under atom renumbering on five of six structures. | Keep all 18 typed refusals until an order-independent supported contract or runtime safety criterion is proved. Then expose and measure a deliberate opt-in published profile. Source improvements and typed refusals are not exact parity. |
| P2 A6 quality | Published v1.0.31: [265/265 geometry/stereo/clash](benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md), 100/265 converged; 262/262 comparable total energies within 1 kcal/mol. The [typing census](benchmarks/2026-10-03-a6-published-v1031-typing-census.md) has 2,908 heavy-atom differences; its largest bucket spans 1,322 atoms on 317 rows, only two with different symmetrized ring sets. A 400-iteration diagnostic converges on 164/265 without changing the default. | Compare MMFF-specific aromatic-ring acceptance/order and bond flags; classify the 101 non-converged rows, remaining types, per-term energy and unsupported cases; measure independent conformer quality. No 3D speed claim yet. |

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
   including the npm E/Z JSON path. The three new-product-map supplement
   cases pass in latest Rust source, but remain to be rerun on a published
   checked-profile crate; Linux/macOS source-wheel graph/status gates pass.
   **Exit:** all 83 rows accounted for on each exposed axis, zero wrong
   confident supported-domain results, typed reasons for refusals, and no
   regression in the original 57. The three jointly invalid rows remain
   classified as such; reaction yield/selectivity is outside this gate.
2. **P1 / A4 — Reduce SMARTS residuals without broadening silent support.**
   The published baseline remains 200 differing match sets / 43 Booleans.
   Current opt-in source profiling corrects 183 match sets and all 43 Boolean
   errors without a wrong-confident cell. The Fe-containing residual was a
   narrowly verified hypervalent carbanion-to-Fe single bond that RDKit
   sanitizes to dative before ring perception. The source fix does not provide
   general organometallic sanitization. Independently adjudicate the 18
   charged-polycycle refusals by ring-set identity (not just count). The
   six-row diagnosis shows RDKit atom-order dependence and three source
   ring-family differences; define a defensible order-independent contract
   before relaxing the guard. Document the supported boundary; preserve
   native SSSR. **Exit:** rerun the full pinned 310,000-cell matrix on affected
   published bindings with row-level match sets, Booleans and refusal counts;
   zero new wrong-confident cells. The 18 source-only refusals do not count
   as exact parity.
3. **P2 / A6 — Gate 3D quality before speed.** The published v1.0.31
   geometry/stereo/clash and same-coordinate total-energy checks are done,
   but full MMFF equivalence is open. Investigate the largest typing bucket
   first (1,322 RDKit type-37/CheMatic type-2 carbons across 317 molecules)
   through MMFF-specific aromatic-ring acceptance and bond flags. The
   compact-ring source candidate reduces that bucket to 16, but published
   artifacts and complex-ring boundaries are not yet cleared. Then adjudicate
   other heavy types and 22 CheMatic refusal/error rows. Classify the 101
   rows still not converged under an exploratory 400-iteration limit (97 at
   cap, four early stops); do not silently raise the shipped limit without independent
   quality and runtime checks. Bound per-term energy/gradient differences,
   timeout/cancellation, parameter/charge coverage and independent conformer
   quality on matched seeds and a held-out cohort. **Exit:** predeclared
   tolerances and complete outcome accounting on published artifacts, with
   typed unsupported for unhandled chemistry. No 3D speed claim before exit.
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
