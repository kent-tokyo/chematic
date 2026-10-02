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
| P1 reactions | Original 57 fixtures match; the exposed 83-case extension classifies products and raw matches across Python/Rust/npm. Published Rust graph+origin+template-map agrees on 73/83 rows, plus one map-only difference. An unpublished checked-source profile has 76 semantic matches, three typed unsupported, one refusal and three jointly invalid rows. | Published artifacts still have confident differences. Recheck source fixes after publication; broaden newly created mapped atoms/multiple products and gate origins in Python/npm. No general SMIRKS parity claim. |

The [v1.0.30 published-artifact packet](benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md),
[benchmark index](benchmarks/README.md), and [validation report](docs/validation.md)
hold versions, hashes, denominators and residuals. The WASM formula and E/Z
JSON fixes are included in v1.0.31 source but need published-package checks
after release. Exposed data are
not sealed evidence; A0's frozen descriptor result and Experimental 3D/MMFF94
have separate gates.

## Priority order and acceptance gates

1. **P0.2 — Complete paired speed evidence.** Extend equivalent-output,
   alternating-order, ≥20-block measurements beyond the finished lanes.
   Separate parse, perception, prepared/reused calls and memory allocations.
   **Exit:** pinned scripts/raw observations, paired intervals and an explicit
   no-claim outcome whenever outputs or work differ.
2. **P1 — Close reaction compatibility gaps.** Extend the 83-case corpus across
   maps, H, charge, aromaticity, stereo, multiple products and invalid input.
   Compare graph, distinct-product multiplicity, maps and atom origins, not
   SMILES spelling. Publish and rerun the checked-source fixes, expand the
   three-case new-product-map supplement, and expose/gate applicable Python/npm
   provenance; verify the npm E/Z JSON fix in a future package. **Exit:** every
   row classified, no wrong confident result in the declared supported domain,
   and no regression in the original 57. This says nothing about yield or
   selectivity. Count #635 SMARTS residuals separately.
3. **P1 — Preserve CIP safety and reduce SMARTS residuals.** Keep the five
   typed CIP abstentions stable under atom-order, spelling and file round trips;
   adjudicate phosphorus independently and define a lone-pair convention
   before assigning R/S. Address 200/310,000 SMARTS failures without silently
   changing native SSSR; any RDKit-style ring-count profile is opt-in.
   **Exit:** complete outcome accounting and no wrong confident label/match.
4. **P2 — Complete A6 before promoting 3D.** Address MMFF94 typing,
   parameters, gradients, convergence, stereo and conformer quality on matched
   seeds/coordinates. **Exit:** typed failure for unsupported inputs and
   independent published-artifact quality evidence before speed claims.
5. **P2 — Test representation limits separately.** Keep CDXML opaque
   preservation distinct from semantic editing, Markush/polymer expansion
   bounded, Standard InChI separate from native identifiers, and canonical
   SMILES separate from the fail-closed stable key. **Exit:** cross-binding
   fixtures with no silent information loss.
6. **External — Rebaseline RDKit only after pinning a new official release.**
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
