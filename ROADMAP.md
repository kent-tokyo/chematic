# chematic roadmap

> Updated 2026-09-25. Current release: **v1.0.26**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

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
- **v1.0.26 improves selected hot paths at release source.** The dated
  output-differential record preserves the measured output boundary while
  reducing SMARTS, ring, and fingerprint work. Shared-VM source timing is not
  a universal or published-package performance claim.

Exact versions, corpora, and limits are in [validation](docs/validation.md) and
the [benchmark index](benchmarks/README.md). Completed detail belongs in dated
records, the [CHANGELOG](CHANGELOG.md), or Git history rather than this plan.

## Priority order

1. **Verify v1.0.26 delivery across release channels.** Record GitHub Release,
   npm, PyPI, crates.io, docs.rs, and Pages independently; do not infer one
   channel from another.
2. **Adjudicate the remaining CIP center.** Resolve row 4480 atom 3
   independently before adopting either engine's label; never guess a label.
3. **Decide optional RDKit-style SMARTS ring counts.** #635's 310,000-cell
   classification is complete; `[Rn]`/`[kn]` still have a documented native
   SSSR-versus-symmetrized-ring boundary. Any parity option must be explicit.
4. **Finish remaining A6 MMFF94 gates.** Heavy-atom typing residuals,
   timeout/convergence accounting, and broader conformer quality remain ahead
   of new 3D breadth.
5. **Prepare a pinned RDKit rebaseline.** Retain 2026.03.6 as historical and
   compare any new official artifact side by side with it.

A confirmed silent-corruption or security regression takes precedence.

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
