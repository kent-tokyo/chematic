# chematic roadmap

> Updated 2026-09-23. Published release: **v1.0.21**. Development target:
> **v1.0.22**. Published-package results and current-source candidates are
> reported separately.

CheMatic's priority is a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. The roadmap favors compatible,
measurable behavior over feature-count competition with RDKit, Indigo, Open
Babel, CDK, or other Rust chemistry libraries.

## Current position

- **A0 core-eight descriptors is complete.** Candidate `5e9211a6` passed the
  frozen 2,000-row development set and the one-time 8,000-row sealed holdout
  for all eight declared fields. Those rows are now exposed and cannot be
  reused as a future holdout.
- **Published v1.0.20 Parse + compatible Morgan is ahead of official
  RDKit.js on the measured browser lanes.** The registry-installed package
  recorded 1.398x parse-inclusive and 3.511x prepared speedups, with
  9,999/9,999 supported rows bit-exact. This is not a universal performance
  claim.
- **Published v1.0.20 MMFF94 is quality-complete for the fixed 265-row lane,
  but not faster than RDKit.** It produced 265/265 sound, stereo-clean,
  clash-free outputs at 0.944x RDKit speed. Current source `e9f178fa` closes
  that measured speed gap on the fixed set and a post-freeze 100-row holdout;
  publication and the remaining A6 numerical gates are separate.
- **Issue #632 is released with a bounded source diagnostic.** Commit `9808f54f`
  reduces the pinned
  RDKit 2026.03.6 SMILES semantic differences from 18/10,000 to 0/10,000 while
  preserving zero graph differences. CIP, Morgan, and SMARTS counts are
  unchanged. The long 28-component x 1,024-relabel audit was interrupted, so
  this is not a completed long-audit or published-package remeasurement claim.
- **The next RDKit residuals are classified, not hidden.** CIP has 230/10,000
  correspondence-correct differences and SMARTS has 3,364 affected rows
  (14,306/310,000 cells). They are tracked by #634 and #635. The one Morgan
  difference is a typed unsupported Fe(II) coordination case.

Detailed evidence lives in [validation](docs/validation.md), the
[benchmark index](benchmarks/README.md), and the compact
[open-work ledger](docs/roadmap-open-work.md). Historical plans remain in Git
history and `docs/archive/`.

## Priority order

1. **Finish #632 verification.** Rerun the long permutation audit when
   practical and preserve fail-closed behavior outside the measured coupled-E/Z
   domain. Its release-source 10k diagnostic is not a substitute for that gate.
2. **Resolve #634 CIP residuals.** Classify representation, oracle,
   unsupported, and implementation differences before changing labels. Never
   turn an unresolved center into a guessed result.
3. **Resolve #635 SMARTS residuals.** Work by semantic family, add small truth
   tables, and keep parse support separate from match semantics.
4. **Finish A6 MMFF94 correctness.** Close same-coordinate energy/term
   residuals, timeout accounting, and broader conformer-quality gates before
   adding new 3D breadth.
5. **Prepare the next RDKit rebaseline.** Keep 2026.03.6 historical and run a
   new packet only against an official, version-pinned artifact.
6. **Complete release evidence for v1.0.22.** Cross-binding checks,
   documentation consistency, package builds, security gates, and public
   package reruns remain distinct acceptance steps.

A confirmed silent-corruption or security regression takes precedence over
this order.

## Accuracy packages

- [x] **A0 — Evaluation contract:** frozen candidates, overlap audit, exposed
  development data, one-time sealed evaluation, complete failure accounting,
  and commit-safe summaries. Core-eight adoption is complete.
- [ ] **A1 — Perception and descriptors:** additional descriptor families,
  aromaticity boundaries, potential stereocenters, and independent holdouts.
- [ ] **A2 — Stereo and identity:** finish #632 verification, then permutation,
  spelling, round-trip, stable-key, enhanced-stereo, and CIP residual gates.
- [ ] **A3 — Fingerprints and retrieval:** preserve named compatible-Morgan
  exactness, top-k contracts, option coverage, and cross-binding parity.
- [ ] **A4 — Workflows and interchange:** reaction/query semantics, V3000
  round trips, attachment metadata, parser accounting, and bounded batch
  outcomes across bindings.
- [ ] **A5 — Independent adjudication:** obtain non-maintainer or independently
  sourced gold evaluation without reusing exposed cohorts.
- [ ] **A6 — 3D and force fields:** energy/term parity, convergence and timeout
  accounting, stereo/clash safety, conformer quality, and published-package
  speed. Keep quality and speed as separate gates.

The acceptance definitions are maintained in
[`docs/rdkit-accuracy-plan.md`](docs/rdkit-accuracy-plan.md). An unchecked
package may contain completed slices; it remains open until every declared
exit passes.

## Trust Release workstreams

| ID | Workstream | Remaining release-relevant exit |
|---|---|---|
| T0 | Evidence discipline | Keep source, package, public, and sealed claims distinct; retain every failure and unsupported row |
| T1 | Compatibility contract | Finish #632's long audit, close #634/#635, and publish operation-specific profiles against pinned oracles |
| T2 | Release and documentation | Synchronize versions, CHANGELOG, READMEs, validation, package metadata, migration notes, and public channels |
| T3 | Browser and agent runtime | Worker/cancellation limits, deterministic batch accounting, typed errors, and reproducible browser comparisons |
| T4 | Parser security | Fixed malformed-input corpus, bounded time/memory, process isolation, and explicit refusal reasons |
| T5 | Stereo torture suite | Atom-order, spelling, file round-trip, selected-center, macrocycle, and upstream-regression cases |
| T6 | Maintenance and independent review | Rebaseline new official dependencies, conduct external review, and publish findings without rewriting history |

See [`docs/trust-release-plan.md`](docs/trust-release-plan.md) for sequencing
and release gates, not for a second backlog.

## Product phases

| Phase | Product outcome | Status |
|---|---|---|
| P0 | Core chemistry and reproducible comparison | Stable core; version-pinned rebaseline remains recurring work |
| P1 | Safe parsing and files | Stable bounded paths; exhaustive malformed-state coverage remains open |
| P2 | Stereo, identity, and canonicalization | Active; #632 long-audit verification and #634 follow |
| P3 | Browser, Node, Python, and agent delivery | Stable selected surface; runtime controls and broader parity remain open |
| P4 | Reactions, SMARTS, and medicinal chemistry | Bounded implementation; #635 and quality reports remain open |
| P5 | 3D, force fields, and materials | Experimental; A6 correctness gates remain open |
| P6 | Release operations and ecosystem trust | Automated in part; external review and future rebaseline are recurring |

Phase numbers describe product areas, not a promise to finish strictly in
numeric order. Current priority is driven by correctness risk and evidence
readiness.

## Release gate for v1.0.22

A v1.0.22 candidate may be proposed when:

- #632's long-running regression is complete or is explicitly deferred and not
  claimed;
- every changed chemistry path has a deterministic regression and binding
  boundary where exposed;
- workspace tests, clippy, documentation/evidence checks, parser-security
  checks, and dependency/license checks pass in their declared environments;
- release metadata, READMEs, CHANGELOG, validation, benchmark indexes, and
  package versions agree;
- source-candidate measurements are not described as registry-package or
  public-channel results.

Release publication, tag creation, registry upload, and public-channel
verification are separate steps. Passing local gates alone does not prove that
a release is available to users.

## Contracts and reference documents

- [Compatibility scope](docs/compatibility-scope.md)
- [Validation report](docs/validation.md)
- [RDKit accuracy plan](docs/rdkit-accuracy-plan.md)
- [Benchmark methodology](docs/benchmark.md)
- [Open-work ledger](docs/roadmap-open-work.md)
- [Trust Release plan](docs/trust-release-plan.md)
- [Release history](CHANGELOG.md)

Completed implementation details belong in CHANGELOG entries, dated evidence,
or Git history. They should not be copied back into this roadmap unless they
change current priority or an acceptance boundary.
