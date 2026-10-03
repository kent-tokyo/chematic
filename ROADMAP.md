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
| P1 reactions | Original 57 fixtures match. Published v1.0.30 had 73/83 graph+origin+template-map matches. The checked v1.0.31-source profile has 76 exact on all three axes, three typed unsupported, one diagnosed refusal and three jointly invalid rows. [Post-merge source and binding gates](benchmarks/2026-10-03-reaction-83-postmerge-binding-gate.md) confirm graph/status counts in Linux/macOS CI source wheels and 3/3 new-product-map supplement matches in Rust source. A [local Python source extension](benchmarks/2026-10-03-reaction-83-python-provenance-source.md) now exposes origins/maps and reproduces the 76/3/1/3 accounting on all three axes. | Verify the new Python provenance gate on Linux/macOS release-profile source wheels and WASM metadata in Node CI. Then rebuild/publish and rerun Rust/Python/npm artifact gates. Typed refusal and invalid input are not exact matches; no general SMIRKS parity claim. |
| P1 SMARTS | Published v1.0.30 has 200 differing **match-set** cells / 310,000, but only 43 hit/no-hit differences. The opt-in Python source profile passed [local and Linux/macOS source-wheel CI](benchmarks/2026-10-03-python-source-smarts-optin-310k.md): 309,982 exact match sets, 18 typed unsupported and zero wrong-confident cells; 183 original residuals and all 43 Boolean errors are corrected. The native matcher is unchanged. The [ring-family diagnosis](benchmarks/2026-10-03-smarts-ring-family-diagnostic.md) finds three source families matching the original RDKit ring sets and three differing; RDKit ring sets themselves varied under atom renumbering on five of six structures. | Rerun a deliberately versioned **published** profile and independent corpus. Keep all 18 typed refusals until an order-independent supported contract or runtime safety criterion is proved. Source improvements and typed refusals are not exact parity. |
| P2 A6 quality | Published v1.0.31: [265/265 geometry/stereo/clash on macOS](benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md), 100/265 converged; 262/262 comparable total energies within 1 kcal/mol. The published [typing census](benchmarks/2026-10-03-a6-published-v1031-typing-census.md) has 2,908 heavy-atom differences. Local and Linux/macOS CI [source-wheel type gates](benchmarks/2026-10-03-a6-source-wheel-mmff94-typing.md) have 451 heavy and eight H differences. The [published 400-iteration diagnosis](benchmarks/2026-10-03-a6-v1031-nonconvergence-classification.md) separates 97 cap hits from four early stops. [Same-host termination diagnostics](benchmarks/2026-10-03-a6-source-mmff94-termination.md): macOS source/published 265/265 success; Linux/Python 3.9 source/published 263/265, with identical typed stereo failures at rows 53 and 246. Published macOS/Python 3.9 succeeds on both rows. | Linux failures are pre-existing on that host, not an introduced source regression. Both rows cross from satisfied stereo before minimization to one violated center after it; diagnose the post-minimization repair and platform/toolchain boundary under [#739](https://github.com/kent-tokyo/chematic/issues/739). Then rerun independent geometry/stereo/clash, per-term energy, missing torsions and conformer quality on published and candidate artifacts. No 3D speed claim yet. |

The [v1.0.30 published-artifact packet](benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md),
[benchmark index](benchmarks/README.md), and [validation report](docs/validation.md)
hold versions, hashes, denominators and residuals. v1.0.31 is published on
all six tracked channels, but its WASM formula and E/Z JSON fixes still need
package-output checks. Exposed data are
not sealed evidence; A0's frozen descriptor result and Experimental 3D/MMFF94
have separate gates.

## Priority order and acceptance gates

1. **P1 / A4 — Close the reaction 83-row gate by classified outcome.**
   Current checked source has 76 exact graph/origin/map rows, three typed
   unsupported, one diagnosed refusal and three jointly invalid. The local
   Python extension reproduces this; Linux/macOS release-wheel and WASM
   metadata CI are pending. Next, rerun the checked profile on published
   Rust/Python/npm artifacts, including npm E/Z JSON and the three
   new-product-map cases. **Exit:** every row classified on all three axes,
   zero wrong-confident supported results and no regression in the original
   57. Refusals and invalid rows are not counted as matches.
2. **P1 / A4 — Publish a bounded SMARTS compatibility profile.** The
   published v1.0.30 baseline has 200 match-set and 43 Boolean differences.
   Source-wheel opt-in checks have 309,982/310,000 exact sets and 18 typed
   unsupported cells, with no wrong-confident result. Adjudicate the 18
   charged-polycycle cases under an order-independent ring contract; preserve
   native SSSR. **Exit:** a full pinned published-binding rerun and an
   independent corpus with row-level exact/refused/failed accounting.
3. **P2 / A6 — Gate 3D quality before speed.** Published v1.0.31 macOS
   passes 265/265 geometry/stereo/clash but converges on only 100/265. The
   source typing improvement is not published, and same-host Linux wheels
   have two pre-existing typed stereo failures (rows 53/246; [#739](https://github.com/kent-tokyo/chematic/issues/739)).
   **Exit:** resolve the platform boundary and separately bound atom typing,
   per-term energy, convergence, timeout/cancellation, geometry/stereo/clash
   and independent conformer quality on published and candidate artifacts.
   No 3D speed claim before this exit.
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
