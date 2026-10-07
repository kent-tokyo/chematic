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
| P1 SMARTS/SMIRKS dialect ([#734](https://github.com/kent-tokyo/chematic/issues/734), [#754](https://github.com/kent-tokyo/chematic/issues/754)) | Contributor's [xsmarts-autoconf](https://github.com/swamidasslab/xsmarts-autoconf): published v1.0.35 75/83 behaviour flags give RDKit 2026.03.6's value (v1.0.30: 42/83), with a decided outcome for each of the 8 others. BioTransformer public rules: the unreleased source gives RDKit's product sets for 1,051/1,051 implicit-H and 6,604/6,652 explicit-H pairs where RDKit has a product (v1.0.35: 775 and 3,807). [Record](benchmarks/2026-10-05-biotransformer-corpus-followups-2.md). | Release the corpus fixes (canonical E/Z writer fixes, RDKit-profile ring model, dummy atoms, E/Z rules, kekulization, MOL wedge reading, template cache); rerun the corpus on the published package. |
| P1 SMARTS | Published v1.0.30: 200/310k match-set differences. Opt-in profile: 309,982/310k exact with 18 typed refusals, which are exactly the cells whose RDKit answer changes with atom order; independent ChEMBL 4,625 rows: 143,375/143,375 exact on PyPI v1.0.34 and source; the published PyPI v1.0.34 Linux wheel passes the 310k opt-in gate. [Record](benchmarks/2026-10-04-stereo-integrity-smarts-a6-followups.md), [310k](benchmarks/2026-10-05-mol-stereo-cip-mmff-followups.md). | macOS/Windows wheel reruns; retain native SSSR. |
| P2 A6 quality | Published v1.0.31 macOS: 265/265 geometry/stereo/clash, 100/265 converged. Published Linux wheels: 263/265 ([#739](https://github.com/kent-tokyo/chematic/issues/739)); source with `libm` math: 265/265 on Linux from a clean commit, 265/265 sound/stereo-clean/clash-free on the external scorer ([record](benchmarks/2026-10-05-mol-stereo-loss-and-clean-a6.md)). | Rerun macOS and Windows (CI jobs) on the `libm` build; per-term energy, convergence and independent conformer quality before speed claims. |

Detailed versions, hashes and denominators are in [validation](docs/validation.md)
and the [benchmark index](benchmarks/README.md). v1.0.32 publication and
published-package reruns are tracked separately. The v1.0.31 WASM formula
and E/Z JSON fixes are confirmed on the published v1.0.34 npm package.
Exposed data are not sealed evidence.

## Priority order and acceptance gates

1. **P1 / A4 — BioTransformer-corpus SMIRKS fixes (#734): released.**
   v1.0.36 ships them. The published PyPI Linux wheel gives the source's
   rows on BioTransformer's public rules: RDKit's product sets for all 1,051
   implicit-H and 6,604/6,652 explicit-H pairs where RDKit has a product,
   none refused, the rest typed (implicit-H equivalence, RDKit truncation,
   products RDKit cannot re-sanitize)
   ([record](benchmarks/2026-10-06-v1036-published-reruns-and-followups.md)).
   RDKit.js 2026.03.6 and 2026.09.1 give the same product sets on 8,627 of
   the corpus's 8,706 rows with products, every change a stereo difference
   in the MOL-block read-back, so the 2026.03.6 profile stands. **Next:**
   rerun against the Python 2026.09.1 wheel once it is published.
2. **P1 / A4 — Close the reaction 83-row gate on published artifacts.**
   Published v1.0.36 PyPI Linux wheel (CPython 3.9), npm package and
   crates.io crate: 80 graph/origin/map rows and three jointly invalid, as
   in source. Remaining: the macOS and Windows wheels, not runnable on the
   Linux evidence host; `published-wheel-chemistry-gates.yml` runs this gate
   and those of items 3–5 on them and checks every count against the Linux
   numbers ([record](benchmarks/2026-10-06-734-754-followups-batch11.md)); on
   the published Linux aarch64 wheel it gives every x86-64 count
   ([record](benchmarks/2026-10-06-734-754-followups-batch12.md)).
   **Exit:** every row classified on all three axes on
   each published artifact, zero wrong-confident supported results and no
   regression in the original 57. Refusals and invalid rows are not counted
   as matches.
3. **P1 / A4 — Broaden reaction evidence.** The BioTransformer public
   rule denominator is measured and its stereo gaps are closed (E/Z through
   `remove_hydrogens`, created and copied double bonds; canonical writer
   fix for atom2-anchored direction stashes; a macrocycle alkene's marker on
   a ring-closure bond, so RDKit reads all 572 explicit-H E/Z molecules
   back; products independent of hash-map order). Closed or typed: the
   naphthalene rows RDKit cannot re-sanitize give RDKit's products after one
   sanitize (BTMR1031, typed in the harness), cyclophosphazene P gets
   RDKit's (spelling-dependent) CIP label, and the `[nH+2]` and native
   created-centre stereo cases are decided
   ([record](benchmarks/2026-10-04-xsmarts-autoconf-v1034.md)); released in
   v1.0.36 and rerun on the published wheel with the same rows. **Exit:**
   the macOS and Windows wheels.
4. **P1 / A4 — Publish a bounded SMARTS compatibility profile.** The
   published v1.0.30 baseline has 200 match-set and 43 Boolean differences,
   mostly the `R<n>` ring-count basis (SSSR versus RDKit's symmetrized
   rings). The opt-in profile is exact on 309,982/310,000 pinned cells and on
   all 143,375 cells of an independent ChEMBL corpus (PyPI v1.0.34 and
   source); its 18 typed refusals are exactly the cells whose RDKit answer
   changes with input atom order, so they stay refusals. `^n` follows RDKit's
   hybridization model (every atom of the exposed 10k, ChEMBL 5k and
   Kekulé-written 15k corpora agrees). The published PyPI
   v1.0.36 Linux wheel passes the pinned 310k gate (309,982 exact, 18 typed
   refusals, none unexpected). **Exit:** the same rerun on the macOS and
   Windows wheels; native SSSR is preserved.
5. **P2 / A6 — Gate 3D quality before speed.** Published v1.0.31 macOS
   passes 265/265 geometry/stereo/clash but converges on only 100/265.
   Published Linux wheels failed stereo on rows 53/246
   ([#739](https://github.com/kent-tokyo/chematic/issues/739)); the cause was
   host-libm last-bit differences amplified by the minimiser, and the source
   now uses the `libm` crate (as WASM already did): Linux 265/265, and
   265/265 sound, stereo-clean and clash-free on the external scorer.
   Source MMFF94 typing matches RDKit on every heavy atom of the exposed 10k,
   ChEMBL 5k and Kekulé-written 15k corpora (v1.0.36: 81, 5 and 39 differ),
   reading RDKit's canonical Kekulé structure and ring list
   ([record](benchmarks/2026-10-06-734-754-followups-batch11.md)), and every
   per-term energy is within 1e-12 kcal/mol of RDKit's on the 262
   same-coordinate rows ([record](benchmarks/2026-10-06-734-754-followups-batch12.md))
   and on all 14,684 single-fragment rows of both corpora; the 3D pipeline
   minimizes with a port of RDKit's BFGS (A6: 135/265 converged, was 100;
   [record](benchmarks/2026-10-07-734-754-followups-batch13.md)). With the
   1,000-iteration pipeline default 264/265 converge (a minimization stopped
   by a stereo constraint re-embeds from another seed;
   [record](benchmarks/2026-10-07-734-754-followups-batch15.md); row 0036 moves to a
   lower basin after iteration 700 and converges at 1,077, as RDKit's own run
   needs about 1,400: decided,
   [record](benchmarks/2026-10-07-734-754-followups-batch16.md)); relaxed by RDKit's
   MMFF94, 155/265 source conformers are within 1 kcal/mol of RDKit's best
   of ten ETKDGv3 conformers (RDKit's first conformer: 125)
   ([record](benchmarks/2026-10-07-734-754-followups-batch14.md)).
   **Exit:** rerun macOS on the `libm` build, then bound convergence,
   timeout/cancellation, geometry/stereo/clash and independent conformer
   quality on published artifacts. No 3D speed claim before this exit.
6. **P0.2 — Complete paired speed evidence after quality gates.** Extend
   equivalent-output, alternating-order, ≥20-block measurements beyond the
   finished lanes. Separate parse, perception, prepared/reused calls and memory
   allocations. **Exit:** pinned scripts/raw observations, paired intervals
   and an explicit no-claim outcome whenever outputs or work differ.
7. **P1 / A2 — Preserve CIP safety.** On the 1,687 stereo rows of the
   exposed 10k lane, Accurate-mode labels, the typed abstentions and
   bond-keyed E/Z are identical across random atom orders, an H round trip
   and a canonical reparse. Source: acyclic phosphorus and acyclic lone-pair
   centres are labelled as RDKit's CIPLabeler does, and saturated ring
   lone-pair centres (bridgehead amines) get RDKit's label for the parsed
   spelling, marked spelling-dependent; none abstain on exposed 10k. P on
   unsaturated rings gets RDKit's label, marked Kekulé-dependent.
   Pseudo-asymmetric lone-pair centres get RDKit's `r`/`s`
   ([record](benchmarks/2026-10-07-734-754-followups-batch13.md)), and a
   carbon whose rule 5 runs through a lone-pair centre gets RDKit's label
   ([record](benchmarks/2026-10-07-734-754-followups-batch14.md)). **Exit:**
   complete outcome accounting and no wrong confident label.
8. **P2 — Test representation limits separately.** Keep CDXML opaque
   preservation distinct from semantic editing, Markush/polymer expansion
   bounded, Standard InChI separate from native identifiers, and canonical
   SMILES separate from the fail-closed stable key. The source MOL writer
   draws stereo (one checked wedge per centre, E/Z by geometry, macrocycle
   E/Z included); RDKit reads 1,686 of 1,687 exposed-10k and all 1,670
   ChEMBL-5k stereo rows back unchanged and none inverted. What a block
   cannot carry (the bridgehead amine chematic reads as OpenSMILES does) is reported by
   `write_mol_with_stereo_report`, Python warns (`StereoLossWarning`), and
   Python `strict=True` / WASM `to_mol_block_strict` refuse it; the default
   writer keeps returning the block (decided). SVG and depiction data draw
   the same stereo depiction (RDKit reads the drawings as the input on
   1,686/1,687 and 1,670/1,670 stereo rows; PyPI v1.0.36: 286 and 317;
   [record](benchmarks/2026-10-07-734-754-followups-batch15.md)); every Python MOL block
   carries 2D coordinates, reaction SVGs and similarity maps draw stereo, and
   branch points are drawn 120° apart (15,000 rows as drawn: 14,757 clean,
   25 with a clash, 242 with a crossing; RDKit 14,504, 276 and 402;
   [record](benchmarks/2026-10-07-734-754-followups-batch16.md)). **Exit:** cross-binding
   fixtures with no silent information loss.
9. **External — Complete RDKit 2026.09.1 rebaseline.** The official npm/WASM
   old/new lanes, published Python 2026.03.6 baseline, and independent C++
   source-build old/new lanes are pinned and measured on exposed 10k inputs.
   The source-build delta reproduces the npm delta: six rows and 12 SMARTS
   cells, all `[R2]`/`[R3]`; old C++ agrees with the old Python wheel on
   canonical SMILES, CIP and Morgan for all 10k rows. The **2026.09.1**
   distributed Python/nanobind lane and distributed native binary remain
   unavailable/unmeasured (no PyPI 2026.9.1 wheel at the 2026-10-06 check).
   The six imine E/Z residuals are default (legacy) CIP labels; accurate mode
   gives rdCIPLabeler's Z, and agrees on all 780 exposed-10k E/Z bonds. The
   12 `[R2]`/`[R3]` cells are relevant-cycle ring counts in 2026.09.1
   ([record](benchmarks/2026-10-06-734-754-followups-batch12.md)); the opt-in
   `profile="2026.09.1"` counts them so (310,000/310,000 native cells,
   [record](benchmarks/2026-10-07-734-754-followups-batch13.md)). Preserve the old baseline, keep source builds
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
