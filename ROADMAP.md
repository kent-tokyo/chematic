# chematic roadmap

> Updated 2026-09-25. Current release: **v1.0.25**. Next development target:
> **v1.0.26** (only after its declared gates pass).

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

- **A0 core-eight descriptors is complete.** A frozen candidate passed the
  declared 2,000-row development and one-time 8,000-row sealed evaluations.
- **Browser Parse + compatible Morgan has a scoped public-package win.** The
  v1.0.20 npm record is bit-exact on 9,999 supported rows and faster on its
  declared browser lanes; it is not a general performance claim.
- **Current RDKit residuals are visible.** #632 has a released 10k diagnostic
  with zero remaining SMILES semantic differences, while its long relabel audit
  remains open. #634 CIP and #635 SMARTS retain classified residuals.
- **v1.0.25 improves provenance and HBA compatibility.** It adds atom-output
  and source-atom APIs, corrects the named RDKit-compatible Python HBA profile,
  and fixes plain SMILES ring-closure spelling.

Exact versions, corpora, and limits are in [validation](docs/validation.md) and
the [benchmark index](benchmarks/README.md). Completed detail belongs in dated
records, the [CHANGELOG](CHANGELOG.md), or Git history rather than this plan.

## Priority order

1. **Finish #632 verification.** Complete the long permutation audit without
   widening the stable-key contract.
2. **Resolve #634 CIP residuals.** Classify each row as fixed, unsupported,
   representation-boundary, or oracle-boundary; never guess a label.
3. **Resolve #635 SMARTS residuals.** Work by semantic family, with small truth
   tables and explicit parse-versus-match boundaries.
4. **Finish A6 MMFF94 correctness.** Close energy/term residuals, timeout and
   convergence accounting, and broader conformer-quality gates before new 3D
   breadth. Quality and speed remain separate gates.
5. **Prepare a pinned RDKit rebaseline.** Keep 2026.03.6 historical and only
   compare a newly pinned official artifact side by side with it.

A confirmed silent-corruption or security regression takes precedence.

## Accuracy packages

| Package | State | Exit focus |
|---|---|---|
| A0 Evaluation contract | Complete | Preserve frozen/exposed-cohort and failure-accounting rules |
| A1 Perception and descriptors | Open | Descriptor families, aromaticity boundaries, potential centers, independent holdouts |
| A2 Stereo and identity | Active | #632, #634, permutation/spelling/file round trips, stable-key scope |
| A3 Fingerprints and retrieval | Open | Compatible-Morgan options, top-k invariants, cross-binding parity |
| A4 Workflows and interchange | Open | #635, reactions, V3000/query semantics, attachment metadata, batch accounting |
| A5 Independent adjudication | External | Non-maintainer gold data and review without exposed-cohort reuse |
| A6 3D and force fields | Active | Energy/terms, convergence, conformer quality, publication-level speed |

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
| P2 | Stereo, identity, canonicalization | Active: #632 and #634 |
| P3 | Browser, Node, Python, agents | Selected stable surface; runtime controls remain |
| P4 | Reactions, SMARTS, medicinal chemistry | Bounded implementation; #635 remains |
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
