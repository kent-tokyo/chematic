# Open-work ledger

Updated 2026-10-05 for the **v1.0.35 release line**. v1.0.35 publication is
verified; published-package chemistry reruns remain open.
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
| 1 · P1/A4 | [#734](https://github.com/kent-tokyo/chematic/issues/734)/[#754](https://github.com/kent-tokyo/chematic/issues/754) autoconf: v1.0.34 65/83 flags give RDKit 2026.03.6's value; unreleased source 75/83, the other 8 decided (sanitized views agree, implicit-H equivalence, adapter difference, SSSR), tool pinned in CI; `random_smiles`/`write` stereo and `remove_hydrogens` aromatic NH fixed ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-xsmarts-autoconf-v1034.md)). | Release (isotope default is a behaviour change); rerun the published package at 75/83 with the same 8 expected differences. |
| 2 · P1/A4 | Published v1.0.34 PyPI Linux wheel (CPython 3.9), npm package and crates.io crate: 80 exact and three jointly invalid, as in source; npm `formula()` 5,000/5,000 and E/Z JSON confirmed. | Rerun the macOS and Windows wheels. Preserve the original 57; invalid rows are not matches. |
| 3 · P1/A4 | Broader reaction evidence: the contributor's approximately 1,200 BioTransformer rules are untested (the rule file needs the BioTransformer 3.0 jar). | Rule-file denominator with exact/differing/refused counts; keep the `~`, `[nH+2]` and native-stereo sweep decisions under review. |
| 4 · P1/A4 | Published v1.0.30 SMARTS: 200/310,000 match-set differences and 43 Boolean differences. Opt-in profile: 309,982/310k exact; its 18 typed refusals are exactly the cells whose RDKit answer depends on atom order; independent ChEMBL 4,625 rows: 143,375/143,375 exact on PyPI v1.0.34 and source ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-stereo-integrity-smarts-a6-followups.md)); published PyPI v1.0.34 (CPython 3.9) passes the 310k opt-in gate against the archived oracle ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-05-mol-stereo-cip-mmff-followups.md)). | Rerun on the macOS and Windows wheels; native SSSR stays unchanged. |
| 5 · P2/A6 | [Published v1.0.31 macOS](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-a6-published-v1031-mmff94-quality.md): 265/265 geometry/stereo/clash, 100/265 converged. Published Linux wheels fail stereo on rows 53/246 ([#739](https://github.com/kent-tokyo/chematic/issues/739)); source with `libm` math and the new MMFF typing: Linux 265/265, and 265/265 sound, stereo-clean and clash-free on the external scorer. Source MMFF typing: 89 heavy atoms differ from RDKit on exposed 10k (451 before this batch; 78 in one fullerene cage) ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-05-mol-stereo-cip-mmff-followups.md)). | Rerun macOS on the `libm` build; separately gate atom types, termination, per-term energy and independent conformer quality on published and candidate artifacts before 3D speed claims. |
| 6 · P0.2 | Extend paired performance beyond the completed Python 63-operation and Node/Chromium/Firefox/WebKit Morgan lanes. | Equal-output, equal-work ≥20-block intervals for more operations/corpora; distinguish parsing, perception, prepared reuse and library memory from process RSS. The historical “21 faster” count is not a current package gate. |
| 7 · P1/A2 | Source CIP: acyclic phosphorus and acyclic lone-pair centres are labelled as RDKit's CIPLabeler does (multiple bonds at the centre not expanded; lone pair as lowest phantom); one abstention left on exposed 10k (a bridgehead amine), 4,394 labels agree, none differ; stable across orders and round trips ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-05-mol-stereo-cip-mmff-followups.md)). | Ring lone-pair centres: RDKit's reading depends on how a ring-closure digit on the centre was written; decide whether to keep the typed abstention. P on unsaturated rings stays `OracleUnstable`. |
| 8 · P2 | CDXML, Markush/polymer, Standard InChI and canonical identity are bounded. Source MOL writer draws stereo (wedges + E/Z geometry, macrocycles included): RDKit reads 1,681/1,687 exposed-10k and 1,663/1,670 ChEMBL-5k stereo rows back unchanged, never inverted; the remaining losses (cage centres, one ring E/Z) are reported by `write_mol_with_stereo_report` and refused by Python `strict=True` / WASM `to_mol_block_strict` ([record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-05-mol-stereo-loss-and-clean-a6.md)). | Decide whether the default writer refuses lossy stereo; cage layouts; per-format preservation, semantic round-trip, typed-refusal and no-false-merge gates across bindings. |
| 9 · external | Official RDKit.js 2026.09.1 is pinned on exposed 10k: Morgan unchanged on 9,999 supported CheMatic rows; graph-checked CIP retains six pre-existing imine E/Z mismatches, five abstentions and one unproven Fe row. Published CheMatic WASM SMARTS has 194 old-oracle versus 195 new-oracle mismatches in 310k cells, with 31 index-unproven cells per lane; old/new RDKit SMARTS differs in 12 `[R2]`/`[R3]` cells. [npm record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-npm-rebaseline.md). A [same-binary C++ legacy-ring probe](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-native-source-rebaseline.md) now attributes all 12 old/new cells on this corpus to RDKit's selectable ring backend; it does not adjudicate chemical correctness. The corrected [published Python 2026.03.6 baseline](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-04-rdkit-2026-09-1-python-baseline.md) also retains six imine E/Z mismatches and five abstentions. PyPI 2026.9.1 was unavailable at the recorded check. | Adjudicate the six E/Z and new-default ring-model residuals separately; recheck wheel/conda availability, hash the exact 2026.09.1 Python artifact, identify its wrapper backend, rerun the complete 10k Python chemistry/typed-error/overhead packet and classify residuals against the pinned old baseline. Keep native C++ separate. |

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
