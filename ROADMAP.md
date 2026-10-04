# chematic roadmap

> Updated 2026-10-04. Release line: **v1.0.34**. Release-source and
> published-package results are kept separate.

CheMatic prioritizes a safe, typed, local-first chemistry kernel for Rust,
Python, Node, WebAssembly, and agent workflows. Compatibility, explicit
boundaries, and reproducible evidence take priority over feature-count races.

## Current position

| Gate | Current evidence | Next exit |
|---|---|---|
| P0.1 outputs | Published v1.0.30 artifact audit covers 10k chemistry, 310k SMARTS, 57 reactions and 63 Python/Rust operations. HBA is 5,000/5,000 on its declared lane. | Retain 200 SMARTS differences, five CIP abstentions and four npm gaps as open; see [acceptance policy](benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md). |
| P0.2 speed | Published Python: 20 exact-output operations pass paired-interval gates on one host; three-browser Morgan has a separate scoped result. | Equal-work and memory lanes on more artifacts/corpora; see [validation](docs/validation.md). |
| P1 reactions | Published v1.0.34 PyPI Linux wheel, npm package and crates.io crate: 80 exact graph/origin/map, 3 jointly invalid, no typed unsupported (v1.0.30: 73/83). | Rerun the macOS and Windows wheels; refusals are not matches. |
| P1 SMARTS/SMIRKS dialect ([#734](https://github.com/kent-tokyo/chematic/issues/734), [#754](https://github.com/kent-tokyo/chematic/issues/754)) | Contributor's [xsmarts-autoconf](https://github.com/swamidasslab/xsmarts-autoconf): v1.0.34 65/83 behaviour flags give RDKit 2026.03.6's value (v1.0.30: 42/83); unreleased source 75/83, with no wrong-confident flag left and a decided outcome for each of the 8 others. Sweep findings fixed or decided; tweak kinds 16 → 4. [Record](benchmarks/2026-10-04-xsmarts-autoconf-v1034.md). | Release the source fixes (isotope default and `random_smiles`/`write` stereo are user-visible); pinned autoconf lane runs in CI. |
| P1 SMARTS | Published v1.0.30: 200/310k match-set differences. Opt-in profile: 309,982/310k exact with 18 typed refusals, which are exactly the cells whose RDKit answer changes with atom order; independent ChEMBL 4,625 rows: 143,375/143,375 exact on PyPI v1.0.34 and source. [Record](benchmarks/2026-10-04-stereo-integrity-smarts-a6-followups.md). | Published-binding 310k rerun of the opt-in profile; retain native SSSR. |
| P2 A6 quality | Published v1.0.31 macOS: 265/265 geometry/stereo/clash, 100/265 converged. Published Linux wheels: 263/265 ([#739](https://github.com/kent-tokyo/chematic/issues/739)); source with `libm` math: 265/265 on Linux. | Rerun macOS and the external geometry/stereo/clash scorer on the `libm` build; per-term energy, convergence and independent conformer quality before speed claims. |

Detailed versions, hashes and denominators are in [validation](docs/validation.md)
and the [benchmark index](benchmarks/README.md). v1.0.32 publication and
published-package reruns are tracked separately. The v1.0.31 WASM formula
and E/Z JSON fixes are confirmed on the published v1.0.34 npm package.
Exposed data are not sealed evidence.

## Priority order and acceptance gates

1. **P1 / A4 — Release the SMARTS/SMIRKS dialect fixes (#734, #754).**
   On published v1.0.34 the autoconf catalog and sweep found confident
   answers RDKit does not give (isotope primitives, aromatic `v`,
   directional bonds, Kekulé reactants, explicit-H radicals, dropped
   `[c±]` products, `i`/`^6`/`^7`). The source aligns 75/83 flags, settles
   the other 8 (sanitized views agree, implicit-H equivalence from #734
   item 5, an adapter difference, native SSSR), and pins the tool in CI.
   The same work fixed `random_smiles`/`write` stereo (39% of stereo rows
   wrong) and `remove_hydrogens` on aromatic NH. **Exit:** released and the
   published package rerun at 75/83 with the same 8 expected differences.
2. **P1 / A4 — Close the reaction 83-row gate on published artifacts.**
   Published v1.0.34 PyPI Linux wheel (CPython 3.9), npm package and
   crates.io crate: 80 exact graph/origin/map rows and three jointly invalid,
   as in source; npm `formula()` and E/Z JSON confirmed. Remaining: the
   macOS and Windows wheels (not runnable on the Linux evidence host).
   **Exit:** every row classified on all three axes on each published
   artifact, zero wrong-confident supported results and no regression in the
   original 57. Refusals and invalid rows are not counted as matches.
3. **P1 / A4 — Broaden reaction evidence.** Run the contributor's
   approximately 1,200 BioTransformer rules when the rule file (it needs the
   BioTransformer 3.0 jar) is available, and keep the remaining sweep
   decisions (`~` product bonds on charged N, `[nH+2]`, native stereo for
   created centres) under review. **Exit:** a published rule denominator
   with exact / differing / refused counts.
4. **P1 / A4 — Publish a bounded SMARTS compatibility profile.** The
   published v1.0.30 baseline has 200 match-set and 43 Boolean differences,
   mostly the `R<n>` ring-count basis (SSSR versus RDKit's symmetrized
   rings). The opt-in profile is exact on 309,982/310,000 pinned cells and on
   all 143,375 cells of an independent ChEMBL corpus (PyPI v1.0.34 and
   source); its 18 typed refusals are exactly the cells whose RDKit answer
   changes with input atom order, so they stay refusals. `^n` follows RDKit's
   hybridization model (99.94% per-atom agreement). **Exit:** a pinned
   published-binding 310k rerun with row-level exact/refused/failed
   accounting; native SSSR is preserved.
5. **P2 / A6 — Gate 3D quality before speed.** Published v1.0.31 macOS
   passes 265/265 geometry/stereo/clash but converges on only 100/265.
   Published Linux wheels failed stereo on rows 53/246
   ([#739](https://github.com/kent-tokyo/chematic/issues/739)); the cause was
   host-libm last-bit differences amplified by the minimiser, and the source
   now uses the `libm` crate (as WASM already did): Linux 265/265.
   **Exit:** rerun macOS and the external scorer on the `libm` build, then
   bound atom typing, per-term energy, convergence, timeout/cancellation,
   geometry/stereo/clash and independent conformer quality on published and
   candidate artifacts. No 3D speed claim before this exit.
6. **P0.2 — Complete paired speed evidence after quality gates.** Extend
   equivalent-output, alternating-order, ≥20-block measurements beyond the
   finished lanes. Separate parse, perception, prepared/reused calls and memory
   allocations. **Exit:** pinned scripts/raw observations, paired intervals
   and an explicit no-claim outcome whenever outputs or work differ.
7. **P1 / A2 — Preserve CIP safety.** On the 1,687 stereo rows of the
   exposed 10k lane, Accurate-mode labels, the five typed abstentions and
   bond-keyed E/Z are identical across random atom orders, an H round trip
   and a canonical reparse. Remaining: adjudicate phosphorus independently
   and define a lone-pair convention before assigning R/S. **Exit:**
   complete outcome accounting and no wrong confident label.
8. **P2 — Test representation limits separately.** Keep CDXML opaque
   preservation distinct from semantic editing, Markush/polymer expansion
   bounded, Standard InChI separate from native identifiers, and canonical
   SMILES separate from the fail-closed stable key. MOL writing without
   coordinates drops tetrahedral stereo and writes E/Z as wedge codes only
   chematic reads back; a coordinate-and-wedge writer is open. **Exit:**
   cross-binding fixtures with no silent information loss.
9. **External — Complete RDKit 2026.09.1 rebaseline.** The official npm/WASM
   old/new lanes, published Python 2026.03.6 baseline, and independent C++
   source-build old/new lanes are pinned and measured on exposed 10k inputs.
   The source-build delta reproduces the npm delta: six rows and 12 SMARTS
   cells, all `[R2]`/`[R3]`; old C++ agrees with the old Python wheel on
   canonical SMILES, CIP and Morgan for all 10k rows. The **2026.09.1**
   distributed Python/nanobind lane and distributed native binary remain
   unavailable/unmeasured. Preserve the old baseline, keep source builds
   separate from distributed packages, and never tune on sealed data.

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
