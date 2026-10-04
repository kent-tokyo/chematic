# chematic roadmap

> Updated 2026-10-04. Release line: **v1.0.33**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

| Gate | Current evidence | Next exit |
|---|---|---|
| P0.1 outputs | Published v1.0.30 artifact audit covers 10k chemistry, 310k SMARTS, 57 reactions and 63 Python/Rust operations. HBA is 5,000/5,000 on its declared lane. | Retain 200 SMARTS differences, five CIP abstentions and four npm gaps as open; see [acceptance policy](benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md). |
| P0.2 speed | Published Python: 20 exact-output operations pass paired-interval gates on one host; three-browser Morgan has a separate scoped result. | Equal-work and memory lanes on more artifacts/corpora; see [validation](docs/validation.md). |
| P1 reactions | Published v1.0.30: 73/83 exact. Checked source: 76 exact, 3 typed unsupported, 1 refusal, 3 jointly invalid. Python Linux/macOS release-profile wheels and WASM Node CI passed the [source gate](benchmarks/2026-10-03-reaction-83-python-provenance-source.md). | Rerun Rust/Python/npm **published** artifacts and npm E/Z JSON; refusals are not matches. |
| P1 SMARTS | Published v1.0.30: 200/310k match-set differences. Opt-in Python source wheels: 309,982 exact, 18 typed unsupported, zero wrong-confident. | Independent corpus and published-profile rerun; retain native SSSR. |
| P2 A6 quality | Published v1.0.31 macOS: 265/265 geometry/stereo/clash, 100/265 converged. Linux/Python 3.9 has two typed stereo failures on published and source wheels. | Resolve [#739](https://github.com/kent-tokyo/chematic/issues/739); rerun per-term energy, convergence and independent conformer quality before speed claims. |

Detailed versions, hashes and denominators are in [validation](docs/validation.md)
and the [benchmark index](benchmarks/README.md). v1.0.32 publication and
published-package reruns are tracked separately. The earlier v1.0.31 WASM
formula and E/Z JSON fixes still lack a registry-package output rerun.
Exposed data are not sealed evidence.

## Priority order and acceptance gates

1. **P1 / A4 — Close the reaction 83-row gate by classified outcome.**
   Checked source has 76 exact graph/origin/map rows, three typed unsupported,
   one diagnosed refusal and three jointly invalid. Linux/macOS release-profile
   source wheels and WASM Node CI passed. Next, rerun on published
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
7. **External — Complete RDKit 2026.09.1 rebaseline.** The official npm/WASM
   old/new lanes and published Python 2026.03.6 baseline are pinned and measured
   on exposed 10k inputs; the **2026.09.1** Python/nanobind and independent
   native lanes remain unmeasured until versioned artifacts/build packets exist.
   Preserve the old baseline and never tune on sealed data.

Silent corruption or a security regression takes precedence. Source,
published-package and sealed results remain separate evidence classes.

The [accuracy plan](docs/rdkit-accuracy-plan.md) defines A0–A6 acceptance
criteria; the [open-work ledger](docs/roadmap-open-work.md) tracks dependencies.

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
