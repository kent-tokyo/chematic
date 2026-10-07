# Open-work ledger

Updated 2026-10-08 for the **v1.0.37 release line**. The published v1.0.37
Linux CPython 3.9 wheel, npm package and crate give the gate counts below
(macOS and Windows wheels: dispatch the repaired
`published-wheel-chemistry-gates.yml`).
Published chemistry/benchmark evidence is versioned:
the P0.1 comparison remains v1.0.30, while the A6 quality packet includes
v1.0.31. Unreleased source results are labeled separately. The
[roadmap](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md) sets priority and exit criteria; this ledger names
dependencies. Completed work and raw evidence stay in the
[CHANGELOG](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md), [validation report](validation.md), and
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md).

## Queue

| Order | Open work | Next verifiable exit |
|---:|---|---|
| 1 · P1/A4 | [#734](https://github.com/kent-tokyo/chematic/issues/734)/[#754](https://github.com/kent-tokyo/chematic/issues/754) autoconf: published v1.0.36 sdist 75/83 flags give RDKit 2026.03.6's value, the other 8 decided (sanitized views agree, implicit-H equivalence, adapter difference, SSSR); tool pinned in CI ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-05-biotransformer-rule-corpus.md)). | Run `published-wheel-chemistry-gates.yml` with `version=1.0.37` on the macOS and Windows wheels. Its first runs failed before any chemistry ran: RDKit 2026.03.6 has no macOS x86-64 wheel, and the rule-table step looked for files the BioTransformer repository stores under other paths. The workflow now records chematic's outputs on each platform without RDKit and compares them on Linux with one hash-pinned RDKit wheel ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-08-734-754-followups-batch20.md)); published v1.0.37 Linux x86-64 gives every count in `validation/published-wheel-chemistry-gates-expected-v1.0.37.json`. Keep the 8 decided differences under review. |
| 2 · P1/A4 | Published v1.0.36 PyPI Linux wheel (CPython 3.9), npm package and crates.io crate: 80 graph/origin/map rows and three jointly invalid on the 83 reaction fixtures, as in source. | Rerun the macOS and Windows wheels. Preserve the original 57; invalid rows are not matches. |
| 3 · P1/A4 | BioTransformer public rules (974 × 400 molecules, RDKit profile): where RDKit has a product, published v1.0.36 (PyPI Linux wheel) and source give the same product sets for 1,051/1,051 implicit-H and 6,604/6,652 explicit-H pairs and refuses none (v1.0.35: 775 and 3,807); the explicit-H rest is implicit-H equivalence (35 rows), RDKit truncation (12) and one naphthalene row whose products RDKit gives after a single sanitize but cannot re-sanitize ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch10.md)). | RDKit.js 2026.03.6 vs 2026.09.1: 8,627/8,706 rows the same, every change a MOL read-back stereo difference; rerun on the Python 2026.09.1 wheel when published ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-v1036-published-reruns-and-followups.md)). |
| 4 · P1/A4 | Published v1.0.30 SMARTS: 200/310,000 match-set differences and 43 Boolean differences. Opt-in profile: 309,982/310k exact; its 18 typed refusals are exactly the cells whose RDKit answer depends on atom order; independent ChEMBL 4,625 rows: 143,375/143,375 exact on PyPI v1.0.34 and source ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-stereo-integrity-smarts-a6-followups.md)); published PyPI v1.0.34 (CPython 3.9) passes the 310k opt-in gate against the archived oracle ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-05-mol-stereo-cip-mmff-followups.md)). | Rerun on the macOS and Windows wheels; native SSSR stays unchanged. |
| 5 · P2/A6 | [Published v1.0.31 macOS](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md): 265/265 geometry/stereo/clash, 100/265 converged. Published Linux wheels fail stereo on rows 53/246 ([#739](https://github.com/kent-tokyo/chematic/issues/739)); source with `libm` math and the new MMFF typing: Linux 265/265, and 265/265 sound, stereo-clean and clash-free on the external scorer. Source MMFF typing (RDKit's canonical Kekulé structure and ring list): 0 heavy atoms differ from RDKit on exposed 10k, ChEMBL 5k and the Kekulé-written 15k (v1.0.36: 81, 5, 39; [record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch11.md)). Source MMFF energies: every term within 1e-12 kcal/mol of RDKit on the 262 same-coordinate rows (max 0.32 before; out-of-plane terms, angle constants, the 10 Å vdW cutoff, torsion lookup and empirical rule), A6 265/265 sound, stereo-clean, clash-free, 100 converged ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch12.md)). Every term equals RDKit's on all 14,684 single-fragment rows of both corpora; RDKit's interfragment rule is opt-in (102/102 salts); RDKit's BFGS minimizer is ported and the pipeline uses it: A6 135 converged, 265/265 sound, stereo-clean, clash-free ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch13.md)). The pipeline's default budget is 1,000 iterations (A6 264/265 converged; a minimization stopped by a stereo constraint re-embeds, [record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch15.md)); relaxed by RDKit's MMFF94, 155/265 source conformers are within 1 kcal/mol of RDKit's best of ten ETKDGv3 conformers (RDKit's first: 125; [record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch14.md)). | Rerun macOS on the `libm` build; conformer quality on published artifacts before 3D speed claims. A6 row 0036 crosses to a lower basin after iteration 700 and converges at 1,077 (RDKit's own run needs about 1,400): default budget unchanged, decided ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch16.md)). |
| 6 · P0.2 | Extend paired performance beyond the completed Python 63-operation and Node/Chromium/Firefox/WebKit Morgan lanes. | Equal-output, equal-work ≥20-block intervals for more operations/corpora; distinguish parsing, perception, prepared reuse and library memory from process RSS. The historical “21 faster” count is not a current package gate. |
| 7 · P1/A2 | Source CIP: acyclic phosphorus and acyclic lone-pair centres are labelled as RDKit's CIPLabeler does (multiple bonds at the centre not expanded; lone pair as lowest phantom); exposed 10k 4,395 agree, none abstain or differ; saturated ring lone-pair centres (bridgehead amines, cyclic phosphines) get RDKit's label for the parsed spelling, marked spelling-dependent (520/520 RDKit spellings; [record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch12.md)). | Pseudo-asymmetric lone-pair centres get RDKit's `r`/`s` ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch13.md)); a carbon whose rule 5 runs through such a centre gets RDKit's label (112 spellings: 3,136/3,136 pseudo-asymmetric centres; [record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch14.md)). P on unsaturated rings gets RDKit's label, marked Kekulé-dependent (ChEMBL 5k: 4,186 agree, none abstain). |
| 8 · P2 | CDXML, Markush/polymer, Standard InChI and canonical identity are bounded. Source MOL writer draws stereo (wedges + E/Z geometry, macrocycles included): RDKit reads 1,686/1,687 exposed-10k and 1,670/1,670 ChEMBL-5k stereo rows back unchanged, never inverted ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch11.md)); the remaining loss (a bridgehead amine) is reported by `write_mol_with_stereo_report`, warned about in Python (`StereoLossWarning`; decided: default unchanged) and refused by Python `strict=True` / WASM `to_mol_block_strict` ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch12.md)). Depictions (SVG, depict data, grids) draw the same stereo: RDKit reads the drawings as the input on 1,686/1,687 and 1,670/1,670 stereo rows (PyPI v1.0.36: 286, 317); bridged layout rows with a clash 212 → 5 of 531 (RDKit 102) ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch15.md)). Branch points drawn 120° apart, Python MOL blocks always carry 2D coordinates, reaction SVGs and similarity maps draw stereo; as drawn, 26 of 15,000 rows have a clash (RDKit 276; batch 17) ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch16.md)); source (batch 20): 14,767 of 15,000 rows clean, 23 with a clash and 232 with a crossing (RDKit 14,504, 276, 402), adamantane and porphyrin templates ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-08-734-754-followups-batch20.md)). | Per-format preservation, semantic round-trip, typed-refusal and no-false-merge gates across bindings. |
| 9 · external | Official RDKit.js 2026.09.1 is pinned on exposed 10k: Morgan unchanged on 9,999 supported CheMatic rows; graph-checked CIP retains six pre-existing imine E/Z mismatches, five abstentions and one unproven Fe row. Published CheMatic WASM SMARTS has 194 old-oracle versus 195 new-oracle mismatches in 310k cells, with 31 index-unproven cells per lane; old/new RDKit SMARTS differs in 12 `[R2]`/`[R3]` cells. [npm record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-npm-rebaseline.md). A [same-binary C++ legacy-ring probe](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-native-source-rebaseline.md) now attributes all 12 old/new cells on this corpus to RDKit's selectable ring backend; it does not adjudicate chemical correctness. The corrected [published Python 2026.03.6 baseline](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-python-baseline.md) also retains six imine E/Z mismatches and five abstentions. PyPI 2026.9.1 was unavailable at the recorded check. | The six E/Z rows are default (legacy) CIP labels: accurate mode gives rdCIPLabeler's Z (780/780 exposed-10k E/Z bonds). The 12 new-default `[R2]`/`[R3]` cells are relevant-cycle ring counts (all 12 reproduced; [record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-06-734-754-followups-batch12.md)); the 2026.03.6 profile keeps refusing them, and `profile="2026.09.1"` gives all 310,000 cells of the native 2026.09.1 grid ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-07-734-754-followups-batch13.md)). Recheck wheel/conda availability, hash the exact 2026.09.1 Python artifact, identify its wrapper backend, rerun the complete 10k Python chemistry/typed-error/overhead packet and classify residuals against the pinned old baseline. Keep native C++ separate. |

The [independent native C++ source-build packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-native-source-rebaseline.md)
now reproduces the same six changed rows and 12 `[R2]`/`[R3]` cells. Its old
oracle matches the distributed Python 2026.03.6 baseline on canonical SMILES,
CIP and Morgan across all 10k rows. Compiled binaries/libraries were not
retained, and this source-build result does not satisfy the still-open
2026.09.1 distributed Python/nanobind or native-binary gates. The
[separate Python 2026.09.1 command catalog](https://github.com/kent-tokyo/chematic/blob/main/validation/rdkit_rebaseline_python_2026_09_execution.json)
is ready for an exact wheel when published; no new-version Python result is
inferred from it.

P0.1's **published-artifact accounting is complete**, not strict RDKit
parity: [the acceptance policy](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-p0-acceptance-policy.md)
retains the 200 SMARTS failures, five CIP abstentions and four npm API gaps.
The [published-artifact packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-02-v1.0.30-published-artifact-gates.md)
contains the versioned denominators and raw-row identities.

## Accuracy package disposition

| Package | State | Queue link |
|---|---|---|
| A0 Evaluation contract | Complete | Preserve frozen/exposed-cohort accounting |
| A1 Perception and descriptors | Open | Representation and descriptor boundaries |
| A2 Stereo and identity | Active | Queue 7: CIP abstentions and adjudication |
| A3 Fingerprints and retrieval | Open | Separate compatible-Morgan/holdout gates |
| A4 Workflows and interchange | Active | Queues 1–4: dialect parity, reactions and SMARTS |
| A5 Independent adjudication | External | Gold data and non-maintainer review |
| A6 3D and force fields | Active | Queue 5: quality before speed |

## Cross-cutting dependencies

- **`local-open`:** parser security and batch-result accounting, browser cancellation
  and typed resource limits, stereo regression fixtures, release-document
  synchronization. These can proceed without a new oracle.
- **`toolchain-open`:** rebuilt Python/npm/WASM artifacts, browser lanes and
  cross-platform builds need their named runtimes; a source-only fix is not a
  published-package result.
- **`data-sealed`:** freeze a candidate and audit overlap before one-time
  evaluation. Exposed or inspected rows never become sealed again.
- **`external-open`:** the RDKit 2026.09.1 Python/native distribution artifacts,
  registry publication for future releases,
  non-maintainer review and independent security assessment require their
  respective outside artifacts or people.
- **`historical`:** source-only A/B timings, rejected candidates and older
  release-channel records remain available for provenance, not current claims.

## Completion rule

For every gate, declare the supported domain, comparator, versions, corpus,
options and failure policy. Account for success, failure, refusal and skipped
inputs; test the corrected behavior; rebuild affected bindings; and link
committed evidence. Name source, published-package and public-channel state
separately. A refusal is not an exact match, and a passing local check is not
a verified release.
