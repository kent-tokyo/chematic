# chematic roadmap

> Updated 2026-10-02. Release line: **v1.0.29**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

- **v1.0.29 corrects #679/#680.** SMIRKS product-side atomic-number atoms
  produce literal atoms; repeated SMARTS charge signs represent one charge.
  The published Python wheel has now been measured on the exposed 10,000-row
  chemistry/310,000-cell SMARTS lane, 57 reaction fixtures, and a 5,000-row
  operation matrix. This is a *Python-only partial P0 packet*, not broad
  reaction/SMARTS parity or a paired, quality-equivalent speed claim. The
  named RDKit HBA profile is corrected after the wheel measurement: that
  source candidate must not be reported as a v1.0.29 package result.
- **v1.0.28 closes the classified #634 CIP residuals at release source.**
  The pinned, exposed 10,000-row RDKit 2026.03.6 lane has 9,995 exact rows;
  the remaining five are typed abstentions. Lone-pair stereocentres are not
  yet assigned R/S labels. This is not complete CIP parity.
- **v1.0.27 added Rust atom tags and release-contract checks.** Caller labels
  survive selected in-memory transformations; SMILES does not carry them.
  The v1.0.26 MMFF94 published-wheel and source timing records stay pinned to
  their measured version, not promoted to v1.0.29 results.
- **A0 core-eight descriptors is complete.** A frozen candidate passed the
  declared 2,000-row development and one-time 8,000-row sealed evaluations.
- **Browser Parse + compatible Morgan has a scoped public-package win.** The
  v1.0.20 npm record is bit-exact on 9,999 supported rows and faster on its
  declared browser lanes; it is not a general performance claim.
- **v1.0.26 advances #632/#634/#635 at release source.** Its #632 long
  relabel audit has zero divergent components. Against pinned RDKit 2026.03.6,
  CIP reaches 9,994/10,000 exact labels and SMARTS has 200/310,000 differing
  cells; remaining differences are classified. Python/WASM SMARTS now use the
  perceived aromatic view, with an explicit behavior change.
- **v1.0.26 advances #637 at release source.** On shared coordinates,
  262/262 comparable MMFF94 rows are within 1 kcal/mol (maximum 0.32). This is
  same-coordinate energy/term evidence, not convergence, conformer-quality,
  or published-package evidence.
- **#637 evidence-regeneration scope complete:** two v1.0.26/RDKit 2026.3.6 Mac arm64
  runs retain all 265 rows; independent geometry/stereo scoring finds
  265/265 versus 264/265 usable. This is a named-cohort quality result,
  separate from the source-built energy packet. Broader A6 conformer-quality,
  timeout/convergence, and heavy-atom typing gates remain open.
- **v1.0.26 improves selected hot paths at release source.** The dated
  output-differential record preserves the measured output boundary while
  reducing SMARTS, ring, and fingerprint work. Shared-VM source timing is not
  a universal or published-package performance claim.

Exact versions, corpora, and limits are in [validation](docs/validation.md) and
the [benchmark index](benchmarks/README.md). Completed detail belongs in dated
records, the [CHANGELOG](CHANGELOG.md), or Git history rather than this plan.

## Priority order and acceptance gates

1. **P0 — Measure the published v1.0.29 artifacts.** Pin the PyPI wheel, npm
   tarball, Rust crate, their hashes, RDKit 2026.03.6, preprocessing, options,
   runtime, and hardware. Re-run the applicable *exposed* versioned 10k
   chemistry, 310k-cell SMARTS, reaction, and 5k operation-matrix lanes; mark
   non-applicable operations explicitly. Report exact or graph-semantic
   agreement, every unsupported/refused/failed/timed-out row, and a paired
   v1.0.28→v1.0.29 delta for the changed operations. Time only equivalent
   operations with the same work included (parse-inclusive and prepared lanes
   separately), repeated samples, raw observations, and uncertainty intervals.
   **Exit:** reproducible raw data, hashes, complete denominators, binding-
   specific results, and a machine-checked report; no unmeasured speed or
   whole-library parity claim. Do not reuse an exposed cohort as sealed data.
   **Progress (2026-10-02):** the published arm64 Python wheel and RDKit
   2026.03.6 have a machine-checked exposed 10k/310k packet, paired
   v1.0.28→v1.0.29 row comparison, 57/57 bounded reaction cases, and a
   5k/63-operation diagnostic matrix. CIP is 9,995 exact + five typed
   abstentions; Morgan is 9,999 exact + one typed refusal; SMARTS retains
   200 classified cells. The changed #679 products pass the small fixture,
   but the 10k lane does not test reaction equivalence. The v1.0.29 named
   HBA profile has a measured 3,641/5,000 agreement regression against this
   RDKit version; the native HBA profile is 5,000/5,000, and the unreleased
   source correction now shares that rule. An importable dev-profile source
   wheel independently reaches 5,000/5,000 on the named API; release-profile
   wheel CI and published-artifact replication are pending. npm/Rust artifact
   lanes, full reaction scope, equivalence-controlled paired timing, and cross-host
   uncertainty remain open. See the [dated record](benchmarks/2026-10-02-v1.0.29-python-accuracy.md).
2. **P1 — Expand #679/#680 into a SMARTS/SMIRKS semantic gate.** Add repeated
   charges, mapped/unmapped product atoms, explicit/implicit H, aromatic and
   aliphatic reactants, atom order, stereo, and multi-reactant cases to a
   versioned query/reaction corpus. Compare supported matches and products to
   pinned RDKit by graph, maps, and provenance rather than SMILES spelling;
   count duplicate or silently missing products. Recount #635's 200/310,000
   SMARTS residual cells (194 symmetrized-ring `[Rn]`/`[kn]`, six ferrocene)
   separately from new #679/#680 cases. Native SSSR remains the default; any
   RDKit-style ring-count mode must be opt-in with its own contract.
   **Exit:** no wrong result on the declared supported fixtures, no regression
   in previously passing cells, and every remaining difference classified as
   fixed, intentional model difference, or typed unsupported outcome. Do not
   call this general SMIRKS/SMARTS parity.
3. **P1 — Preserve and then narrow the CIP abstention boundary.** Re-run the
   published v1.0.29 lane before quoting a v1.0.29 result. The v1.0.28-source
   9,995/10,000 exact result and five typed abstentions are historical: four
   `oracle_unstable` phosphorus cases and one `lone_pair_center`. First gate
   zero wrong confident labels and stable abstention reasons under atom-order,
   SMILES-spelling, and V2000/V3000 round trips. Independently adjudicate the
   phosphorus cases; define an explicit lone-pair stereo convention before
   converting any abstention to R/S. **Exit:** complete row accounting and
   independent gold for each newly resolved case. Five safe abstentions are a
   valid bounded contract, not complete CIP compatibility.
4. **P2 — Complete A6 before changing the Experimental 3D label.** Address
   MMFF94 heavy-atom typing residuals, missing parameters, analytic gradients,
   convergence/timeouts, stereo retention, and conformer quality on matched
   seeds and coordinates. Keep the v1.0.26 262/262 same-coordinate energy and
   265-row wheel quality packets historical. Predeclare quality-equivalent
   scope and non-inferiority thresholds before any speed comparison.
   **Exit:** no silent success for untyped or nonconverged inputs, full outcome
   accounting, independent geometry/stereo scoring, and published-artifact
   replication. Do not infer ETKDG or broad MMFF94 parity from finite energy.
5. **P2 — Keep representation limits explicit and test each separately.** For
   CDXML, distinguish exact-source/opaque preservation from semantic editing;
   for Markush/polymer, require explicit bounded expansion or typed refusal;
   for Standard InChI, distinguish the optional standard implementation from
   native identifiers; for canonical SMILES, test idempotency, permutation and
   file round trips, and stable-key no-false-merge/abstention separately from
   string equality. **Exit:** a cross-binding capability table and round-trip
   fixtures with no silent information loss; expand the supported domain only
   after its independent gate passes.
6. **External dependency — Rebaseline a later RDKit release only after an
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
| P0 | Core chemistry and reproducible comparison | Stable; rebaseline is recurring |
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
