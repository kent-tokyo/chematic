# chematic roadmap

> Updated 2026-10-03. Release line: **v1.0.30**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

| Gate | Verified on exposed data | Still open |
|---|---|---|
| P0.1 published outputs | The v1.0.30 Python wheel, npm tarball, and pinned crates.io graph each rerun 10,000 chemistry inputs, 310,000 SMARTS cells, and the legacy 57 reactions. Python and Rust each cover all 63 operation outputs (210,410 rows/version); only HBA and its bundle change from v1.0.29, on the same 1,359 inputs. Published HBA agrees with RDKit on 5,000/5,000 rows. The 69 initially unresolved chemistry classifications are now mapped to 200 SMARTS cells and five typed CIP abstentions by full published-row identity. npm covers 59/63 operation APIs. | Decide the compatibility policy for the 200 SMARTS residual cells and representation differences. Retain four npm-unexposed APIs as explicit gaps. Do not treat the source-only formula fix as a v1.0.30 result. |
| P0.2 paired speed | Published Python has 20 alternating blocks for the 63-operation matrix; 20 operations clear its exact-output and interval gate on one host. Six isolated follow-up lanes split parse-inclusive, first-use and precomputed calls and record whole-process peak RSS separately. | Equal-work perception boundaries, other operations and bindings, matched memory accounting, and independent-host replication. The old “21 faster” count remains historical. |
| P1 reactions | The original 57 fixtures still match. The 83-case exposed extension covers maps, H, charges, aromaticity, products, stereo and invalid input. An unpublished checked-source profile has 76 semantic matches, three typed unsupported, one typed refusal and three jointly invalid rows, with no wrong-confident output inside that declared profile. | The published v1.0.30 artifacts still have confident differences. Recheck the source fixes in a future artifact; finish product multiplicity, map labels and atom origins in each applicable binding. Do not claim general SMIRKS parity. |

The [published-artifact packet](benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md)
and [isolated time/memory record](benchmarks/2026-10-03-v1.0.30-isolated-python-time-memory.md)
hold the hashes, denominators, raw outputs, intervals and exact residuals.
An exposed corpus is not sealed evidence. The 9,995/10,000 CIP agreement with
five typed abstentions and Experimental 3D/MMFF94 status remain separate
boundaries; A0's frozen descriptor gate is a separate result too.

Exact versions, corpora, and limits are in [validation](docs/validation.md) and
the [benchmark index](benchmarks/README.md). Completed detail belongs in dated
records, the [CHANGELOG](CHANGELOG.md), or Git history rather than this plan.

## Priority order and acceptance gates

1. **P0.1 — Re-run all applicable gates on published v1.0.30 artifacts.** Pin
   each PyPI wheel, npm tarball, and crates.io crate graph by version and
   hash, plus RDKit 2026.03.6, corpus hashes, options, runtime, and hardware.
   Re-run the exposed 10,000-row chemistry, 310,000-cell SMARTS, 57-reaction,
   and 5,000-row/63-operation suites in every binding that exposes the
   relevant API. Record unsupported or unexposed cells explicitly; do not
   pretend all artifacts expose identical operations. Compare each overlapping
   v1.0.29→v1.0.30 output row, treating the named HBA correction as the only
   *expected* change, not assuming everything else stayed unchanged.
   **Exit:** artifact-specific raw outputs and hashes, full input/outcome
   accounting, a machine-checked cross-version diff, and independent
   adjudication of every non-HBA delta. Keep the published 5,000/5,000 HBA
   result as a completed sub-gate, not a substitute for this packet. Exposed
   inputs do not become sealed data.
   The Rust 63-operation output differential is now a completed sub-gate;
   finish the remaining cross-binding outcome adjudication and full packet
   before closing P0.1.
2. **P0.2 — Make speed comparisons symmetric and repeatable.** Measure both
   v1.0.29 versus v1.0.30 and RDKit versus chematic on pinned, published
   artifacts and equivalent work. Balance execution order (for example ABBA
   and BAAB), use at least 20 paired blocks after warm-up, and retain every
   observation, environment detail, and a paired confidence interval. Separate
   parse-inclusive from pre-parsed calls, prepared from unprepared queries,
   perception/precomputation from the timed operation, and peak memory from
   throughput. Gate each speed claim on its stated output-equivalence scope;
   a mismatching HBA profile cannot earn an equivalent-work win.
   **Exit:** reproducible scripts and raw results, effect sizes and uncertainty
   by operation/binding, and an explicit no-claim outcome where the interval
   crosses parity or the work differs. The historical “21 operations faster”
   remains a reference diagnostic until replaced.
3. **P1 — Expand reaction compatibility beyond 57 fixtures.** Freeze and
   publish a stratified exposed corpus covering atom maps, aromaticity,
   implicit/explicit H, charges, multiple reactants/products, stereo, and
   malformed or unsupported inputs. Compare against pinned RDKit by product
   graph, multiplicity, maps, and provenance rather than SMILES spelling.
   Classify each row as semantic match, typed refusal, or typed unsupported;
   count confident mismatches, duplicate/missing products, crashes, and
   timeouts as failures rather than hiding them in those categories.
   **Exit:** 100% row accounting, zero wrong confident output in the declared
   supported domain, and no regression in the original 57 fixtures. This is
   transformation compatibility, not reaction yield or selectivity accuracy.
   The checked source profile reaches zero wrong-confident rows on the exposed
   83-case corpus by explicitly declining mismatched chiral-template semantics
   and reporting filtered valence products. Current source retains both
   symmetric ether-cleavage origin assignments in a regression test; the
   published v1.0.30 crate does not. Before closing P1, publish and rerun
   that provenance lane, separately gate template map labels and carried/new
   atom origins on supported cases, verify the npm E/Z serialization fix in
   the next published artifact, and extend the supported-domain corpus. Do
   not count an opt-in refusal profile as general RDKit reaction parity.
   Recount #635's 200/310,000 SMARTS cells separately; keep native SSSR as
   default and any RDKit-style ring-count profile opt-in.
4. **P1 — Preserve and then narrow the CIP abstention boundary.** The v1.0.28-source
   9,995/10,000 exact result and five typed abstentions are historical: four
   `oracle_unstable` phosphorus cases and one `lone_pair_center`. First gate
   zero wrong confident labels and stable abstention reasons under atom-order,
   SMILES-spelling, and V2000/V3000 round trips. Independently adjudicate the
   phosphorus cases; define an explicit lone-pair stereo convention before
   converting any abstention to R/S. **Exit:** complete row accounting and
   independent gold for each newly resolved case. Five safe abstentions are a
   valid bounded contract, not complete CIP compatibility.
5. **P2 — Complete A6 before changing the Experimental 3D label.** Address
   MMFF94 heavy-atom typing residuals, missing parameters, analytic gradients,
   convergence/timeouts, stereo retention, and conformer quality on matched
   seeds and coordinates. Keep the v1.0.26 262/262 same-coordinate energy and
   265-row wheel quality packets historical. Predeclare quality-equivalent
   scope and non-inferiority thresholds before any speed comparison.
   **Exit:** no silent success for untyped or nonconverged inputs, full outcome
   accounting, independent geometry/stereo scoring, and published-artifact
   replication. Do not infer ETKDG or broad MMFF94 parity from finite energy.
6. **P2 — Keep representation limits explicit and test each separately.** For
   CDXML, distinguish exact-source/opaque preservation from semantic editing;
   for Markush/polymer, require explicit bounded expansion or typed refusal;
   for Standard InChI, distinguish the optional standard implementation from
   native identifiers; for canonical SMILES, test idempotency, permutation and
   file round trips, and stable-key no-false-merge/abstention separately from
   string equality. **Exit:** a cross-binding capability table and round-trip
   fixtures with no silent information loss; expand the supported domain only
   after its independent gate passes.
7. **External dependency — Rebaseline a later RDKit release only after an
   official artifact is pinned.** Keep 2026.03.6 and the new oracle results
   side by side; do not rewrite older evidence or tune on sealed data.

A confirmed silent-corruption or security regression takes precedence. A0
sealed results, historical source timings, and published-package measurements
remain separate evidence classes throughout this plan.

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

## Trust Release workstreams

| ID | Outcome |
|---|---|
| T0 | Keep source, package, public, and sealed evidence distinct; account for every input |
| T1 | Publish operation-specific compatibility profiles against pinned oracles |
| T2 | Keep versions, metadata, README, validation, and release channels synchronized |
| T3 | Provide bounded browser/agent runtime controls and deterministic batch outcomes |
| T4 | Maintain malformed-input corpus, resource bounds, and typed refusals |
| T5 | Extend stereo permutation, spelling, round-trip, and upstream-regression gates |
| T6 | Rebaseline dependencies and publish independent review findings faithfully |

## Product phases

| Phase | Area | Status |
|---|---|---|
| P0 | Core chemistry and reproducible comparison | Core stable; v1.0.30 artifact and paired-speed gates open |
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
