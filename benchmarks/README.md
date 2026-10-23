# chematic benchmark records

This directory stores dated evidence. Every result is scoped to the named
source or package, comparator, corpus, host, runtime, options, and timed
boundary. Methodology is documented in [`docs/benchmark.md`](../docs/benchmark.md);
the [complete inventory](archive-index.md) links every retained record.

## Start here

Read the boundary column before using a number. `Published` means the named
registry artifact was measured. `Source` means a pinned revision was built
locally. `Historical` means the result remains reproducible evidence but is
not a statement about v1.2.0 or later source.

| Question | Current record | Boundary |
|---|---|---|
| Do XML reference spellings change CML chemistry? | [CML character-reference repair](../docs/benchmarks/2026-10-22-cml-character-references.md) | Original 14-input source packet: seven valid field-preserving outputs and seven explicit rejections. Valid coverage matches Open Babel 3.2.1; rejections are not successful coverage or broad superiority. |
| Are standard CJSON names and conflicting declarations handled explicitly? | [CJSON ⇄ MOL2 name preservation](../docs/benchmarks/2026-10-21-cjson-name-preservation.md) | Authored source development packet: seven retained outputs and four explicit refusals. Legacy names count as retained meaning; no independent holdout or broad superiority claim. |
| How are symmetric unlabelled R-group sites aligned? | [source comparison](2026-10-22-rdkit-rgroup-symmetric-columns.md) | A first-row tie correction brings one exposed 11-batch graph-and-label panel to 11/11 against RDKit 2026.09.1. The separate 20-row symmetric-ring residual and broader R-group accuracy remain open. |
| Is MCS input-order independent and are incomplete searches explicit? | [tetrahedral neighbor parity](2026-10-22-mcs-stereo-neighbor-parity.md) · [stereo count record](2026-10-22-mcs-stereo-neighbor-parity.json) · [complete-ring valid-best correction](2026-10-22-mcs-complete-rings-valid-best.md) · [equal-atom bond-count correction](2026-10-22-mcs-equal-atom-bond-count.md) · [MCS non-induced/status follow-up](2026-10-21-mcs-non-induced-order-status.md) | Source fixes recover 14/14 selected two-input stereo counts across both molecule orders, a valid common chain under complete-ring constraints, an equal-atom higher-bond tie, and the earlier ring/chain case. A mixed-stereo three-input count still differs from RDKit; full FMCS parity, a new frozen holdout and broad performance evidence remain open. |
| Where does source CIP still trail RDKit? | [CIP lazy-preparation audit](2026-10-19-cip-lazy-preparation.md) | 24 balanced blocks; source CIP first-use 1.74× and parse-inclusive 1.62× faster. Standalone CIP/SMILES/SMARTS gaps remain; stress mismatches stay unranked. |
| What is the broadest cross-platform public-wheel chemistry packet? | [v1.1.3 published artifact packet](2026-10-14-v113-published-artifact-packet.md) | Four exact PyPI wheels; npm/WASM, crates.io, speed, and newer releases are separate. |
| What is the broadest published RDKit performance comparison? | [v1.0.42 / RDKit 2026.09.1](2026-10-11-v1042-published-rdkit-packet.md) | Exact PyPI/npm artifacts on hosted Linux x64 and macOS arm64. |
| What are the latest Python descriptor and Morgan source results? | [terminal TPSA](2026-10-18-terminal-tpsa.md) · [Morgan view/ring scratch](2026-10-18-morgan-borrowed-view-ring-scratch.md) | Source wheels on one macOS arm64 host; exact retention gates, not package or toolkit-wide superiority. |
| What do the current parse/first-use/reuse/RSS lanes show? | [published PyPI 1.1.8 Labute ASA](2026-10-19-v118-published-labute-asa-cosmolkit.md) · [lazy-preparation boundary](2026-10-19-python-lazy-preparation.md) · [10k source follow-up](2026-10-19-python-labute-asa-cosmolkit.md) | Published v1.1.8 matches all 10,000 ASA values on one host. The later source candidate shortens parse and parse-inclusive work while moving preparation cost into first use; prepared reuse is approximately flat. Clean-wheel and second-host evidence remain open. |
| How do reaction, standardization, and labelled R-group workflows compare with RDKit? | [current automatic-column source check](2026-10-22-rdkit-rgroup-auto-columns.md) · [16-row machine record](2026-10-22-rdkit-rgroup-auto-columns.json) · [previous auto-H check](2026-10-18-rdkit-rgroup-auto-hydrogen-followup.md) · [previous R-group hydrogen rows](2026-10-18-rdkit-rgroup-hydrogen-followup.json) · [previous source check](2026-10-18-rdkit-workflows-followup.md) · [83 reaction rows](2026-10-18-rdkit-reaction-83-followup.json) · [stereo safety](2026-10-18-rdkit-reaction-stereo-safety-followup.json) · [expanded workflow packet](2026-10-12-rdkit-workflows-expanded-source.md) | Source-built Python extension: 16/16 selected automatic R-group rows and 16/16 mapped-hydrogen rows match RDKit 2026.09.1 by graph and label. The reaction gate has 80/83 exact rows; three inputs are invalid in both engines. A separate 20-row symmetric-ring R-group batch remains 0/20 exact labelled output, so broad superiority remains open. |
| How did earlier Cleanup work progress? | [2026-10-19 source check](2026-10-19-rdkit-cleanup-acid-metal-followup.md) · [4,069-row source result](2026-10-19-rdkit-cleanup-acid-metal-followup.json) · [previous source check](2026-10-18-rdkit-cleanup-alkali-oxygen-followup.md) · [25-row coordination-tag fix](2026-10-18-rdkit-cleanup-nontetra-followup.md) · [earlier source sequence](2026-10-18-rdkit-cleanup-source.md) | Historical source candidate: 3,933/4,069 exact. Use the current standardization row below for the latest Cleanup count. |
| How closely do Normalize, Reionize, Cleanup and CTfile readback match RDKit? | [current Cleanup source](2026-10-22-rdkit-cleanup-metal-aromaticity.md) · [previous Cleanup source](2026-10-22-rdkit-cleanup-square-planar.md) · [current Reionize source](2026-10-22-rdkit-reionize-in-memory-stereo.md) · [machine record](2026-10-22-rdkit-reionize-in-memory-stereo.json) · [previous standardization source](2026-10-21-rdkit-reionize-catalog-source.md) · [older source summary](2026-10-18-trust-gates-and-residuals.md) · [CTfile classification](../validation/results/rdkit-molblock-readback-current-2026-10-18.json) | RDKit 2026.09.1, 4,069 comparable source rows: Normalize 4,069, Cleanup 4,057, Reionize 4,069 exact. The Cleanup metal-aromaticity fix has no regressions on a separate 10,000-row corpus. CTfile is a separate lane. |
| Does Cleanup preserve total formal charge? | [historical v1.1.8 charge audit](2026-10-18-rdkit-cleanup-charge-invariant-audit.json) · [current RDKit-contract change](2026-10-19-rdkit-cleanup-acid-metal-followup.md) | The historical audit is superseded for three rows: current source follows pinned RDKit metal representations even where total formal charge changes. This is a versioned compatibility choice, not a general charge model. |
| How do bounded workflows compare with Indigo? | [Indigo workflow packet](2026-10-12-indigo-workflows-source.md) | One-host source lanes with 15 reactions, 39 parent decisions, and 16 labelled R-group rows. |
| What is the current 3D/MMFF boundary? | [third ETKDG seed batch](2026-10-22-etkdg-third-seed-batch.md) · [lowest-unfinished ETKDG retry](2026-10-22-etkdg-lowest-unfinished-retry.md) · [symmetric conformer clustering](2026-10-22-conformer-symmetric-clustering.md) · [previous ETKDG retry](2026-10-22-etkdg-unfinished-minimum-retry.md) · [earlier convergence retry](2026-10-22-etkdg-quality-convergence-retry.md) · [pooled ETKDG return set](2026-10-22-etkdg-quality-pooled.md) · [distinct-seed ETKDG development](2026-10-22-etkdg-quality-seed-budget.md) · [sixfold ETKDG gradient development](2026-10-21-etkdg-sixfold-gradient-ensemble.md) · [frozen holdout score](2026-10-19-trust-holdout-v1-score.md) · [RDKit 2026.09 bounds](2026-10-12-a6-rdkit-2026-09-bounds-followup.md) · [gradient rescue](2026-10-13-a6-mmff-gradient-rescue-development.md) | Source clustering merges symmetry-equivalent conformers in focused RDKit-checked fixtures. A bounded MMFF retry for the lowest unfinished candidate raises converged returned candidates from 190 to 191 of 240 on an already used 12-molecule development slice, with two strict converged-energy wins versus RDKit. The third seed batch closes the previous 0.308 kcal/mol development gap and matches RDKit 60 on 12/12 within 1e-6 kcal/mol, with extra MMFF retry work. Equal-compute quality, timing, memory, and diversity remain open. |
| How does file interchange compare with Open Babel? | [comparison program](../docs/openbabel-file-io.md) · [SDF follow-up](2026-10-14-openbabel-sdf-read-followup.md) · [V3000 packet](2026-10-14-openbabel-v3000-speed.md) · [stereo preservation](2026-10-15-openbabel-stereo-paired.md) · [historical v1.0.42 same-process baseline](2026-10-10-openbabel-file-io-same-process-v1.0.42.md) | Open Babel 3.2.1 is the file-I/O comparator. Semantic equality is required before speed ranking; diverse residual rows stay unranked. Newer records are one-host source evidence, while the v1.0.42 baseline is historical and must not be presented as current-package performance. |
| How does COSMolKit compare? | [0.5 comparison](2026-10-10-cosmolkit-050.md) · [seeded ETKDG](2026-10-09-cosmolkit-parity-2.md) | Source-built operation and unusual-SMILES evidence; not a published-package or universal claim. |
| How does browser/WASM Morgan compare? | [three-browser packet](2026-10-03-v1030-published-linux-three-browser-paired20.md) | Published npm artifacts, 250 exact rows, and 20 blocks in Chromium, Firefox, and WebKit. |
| What is the broad historical Python operation matrix? | [v1.0.30 63-operation matrix](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | Historical published-wheel result against RDKit 2026.03.6. |

## Recent focused follow-ups

- [SDF named field preservation](../docs/benchmarks/2026-10-22-sdf-data-item-framing.md):
  12 original field cases in 48 format/line-ending variants; exact field strings
  retained on 48/48 versus 20/48 before and 36/48 with Open Babel 3.2.1 defaults.
  Named-field projection only; property order/header metadata remain outside scope.
- [Checked CML source retention](../docs/benchmarks/2026-10-22-cml-source-conversion.md):
  six original documents, 6/6 declared source-field preservation versus 3/6
  for Open Babel 3.2.1 defaults and 3/6 with property/metadata output enabled.
  Exact source retention is opaque; interpreted chemistry is not verified.

- [CML supplied 2D coordinate precision](../docs/benchmarks/2026-10-22-cml-coordinate-precision.md):
  eight original drawing layouts; source reader/writer preserves 8/8 within
  1e-10 Å absolute, versus 2/8 before and 3/8 with Open Babel 3.2.1 defaults.
  Numeric preservation evidence, not geometry accuracy or toolkit-wide parity.

- [Symmetry-aware conformer clustering](2026-10-22-conformer-symmetric-clustering.md):
  merges atom-permutation-only duplicates in focused CF3, acetate and nitro
  fixtures while retaining a distinct geometry; RDKit 2026.09.1 agrees on the
  swapped-pair RMSD. Force-field and search quality are separate.

- [ETKDG unfinished-minimum retry](2026-10-22-etkdg-unfinished-minimum-retry.md):
  retries the lowest-energy candidate when it did not converge, even when
  another candidate did. Two of 12 development molecules gain converged
  minima below RDKit 40 with extra optimization work; the 0.308065 kcal/mol
  gap on a third remains.

- [ETKDG quality convergence retry](2026-10-22-etkdg-quality-convergence-retry.md):
  adds a bounded MMFF retry only when the returned pool has no converged
  candidate. One of 12 development molecules gains a converged candidate;
  the extra optimization work is explicit, and another molecule's 0.308065
  kcal/mol converged-energy deficit remains.

- [Pooled ETKDG quality return set](2026-10-22-etkdg-quality-pooled.md):
  retains converged conformers across both candidate batches; 7/12 development
  rows improve returned convergence, none regress, and all 240 returned
  conformers pass the pinned stereochemistry and same-coordinate energy checks.
  One converged-energy deficit of 0.308065 kcal/mol remains versus RDKit 40.

- [Distinct-seed ETKDG quality profile](2026-10-22-etkdg-quality-seed-budget.md):
  uses the second candidate batch for new seeds when corrected and compatibility
  gradients coincide; 4/12 lower-energy results versus the prior profile, 0
  regressions on that metric, and 3/12 deficits versus RDKit 40. This is a
  development slice with no converged strict win over RDKit.

- [Sixfold ETKDG gradient ensemble](2026-10-21-etkdg-sixfold-gradient-ensemble.md):
  source quality profile improves best MMFF94 energy on 2/6 selected
  sulfonamides at 40 total candidates versus RDKit 2026.09.1's 40 candidates;
  four agree within 1e-6 kcal/mol. RDKit master has since fixed the derivative.

- [MCS non-induced and status follow-up](2026-10-21-mcs-non-induced-order-status.md):
  fixes the `C1CC1`/`CCC` input-order regression, separates timeout from native
  resource limits across Python/WASM, and records the focused correctness and
  timing boundary. The old holdout is treated as inspected historical data.

- [Reionize ordered acid/base catalog](2026-10-21-rdkit-reionize-catalog-source.md):
  25 additional exact rows and no regressions on the 4,069-row RDKit 2026.09.1
  source corpus; five in-memory stereo serialization differences remain.

- [Conformer symmetry RMSD correction](2026-10-21-conformer-symmetry-rmsd.md):
  terminal-group and isotope matching agree with RDKit on eight exposed cases;
  3D geometry and MMFF94 energy quality are separate gates.

- [Python lazy preparation boundary](2026-10-19-python-lazy-preparation.md):
  [single-operation JSON](2026-10-19-python-lazy-preparation-frontier.json),
  [two-descriptor JSON](2026-10-19-python-lazy-preparation-bundle.json),
  [prepared-reuse JSON](2026-10-19-python-lazy-preparation-reuse.json), and
  [source manifest](2026-10-19-python-lazy-preparation-source-manifest.json).
  One-host source evidence separates parse-only, first-use, parse-inclusive,
  prepared reuse, RSS, and exact output retention.

- [COSMolKit rotor and Layered speed audit](2026-10-19-rotor-layered.md): current-source screening, exact old/new outputs, all 24 formal orders and all 72 comparisons; remaining rotor and MQN boundaries are explicit.

- [Morgan, TPSA and Python ring-count caches](2026-10-19-python-ring-cache.md): fifth combined source candidate passes its three primary improvements and all 21 control floors; all four rejected prior candidates and all 63 comparisons are preserved.

- [Compact core cache rejected pilot](2026-10-19-morgan-tpsa-compact-cache-pilot.md): substantial scalar first-use gains remain diagnostic; ring repeated/parse-inclusive conditions prevent adoption. All eight blocks and 21 lanes retained.

- [Numeric cache placement rejected pilot](2026-10-19-morgan-tpsa-hot-scalars-pilot.md): all eight diagnostic blocks and 21 lanes retained; TPSA/ring repeated-use regressions prevented formal adoption.

- [Morgan hashing and TPSA environment typing](2026-10-19-morgan-tpsa-environments.md): combined source gates fail; both first-use primaries, all 21 source floors, all 63 comparisons and the rejected standalone epoch remain auditable.

- [Direct Morgan connectivity hashing](2026-10-19-morgan-direct-connectivity.md): source adoption fails on the full 10,000-row corpus; all 63 comparisons, source floors, external non-wins and raw numeric accounting retained.

- [Published PyPI 1.1.8 Labute ASA comparison](2026-10-19-v118-published-labute-asa-cosmolkit.md)
  ([raw JSON](2026-10-19-python-labute-asa-cosmolkit-public-wheel.json)):
  exact 10,000-row agreement, public-wheel and installed-extension hashes,
  strict 20/20 wins over RDKit in first-use, repeated-use, and parse-inclusive
  lanes. Repeated and parse-inclusive use also beat COSMolKit 20/20; first use
  is correctly reported as a tie. One macOS arm64 wheel and host only.

- [10k Labute ASA versus RDKit and COSMolKit](2026-10-19-python-labute-asa-cosmolkit.md)
  ([raw JSON](2026-10-19-python-labute-asa-cosmolkit.json)):
  exact 10,000-row agreement and strict 20/20 source wins for first use,
  repeated use, and parse-inclusive use against both comparators. Ring count
  remains a separate unranked lane because COSMolKit differs from the pinned
  RDKit result on six rows. The record now starts with a user-facing summary,
  API example, availability boundary, and guidance on what it does not prove.
  One-host source evidence only.

- [Python first-use Labute ASA and ring count](2026-10-18-python-first-use-ring-asa.md):
  historical 1,000-row output agreement against RDKit 2026.09.1. First-use ASA is
  1.908x and ring count 1.419x faster with 20/20 block wins. The bounded Python
  parser prepares equivalent ring/aromatic state; Rust remains lazy. The later
  10k record supersedes its current-candidate ASA conclusion. Retained inputs:
  [focused raw JSON](2026-10-18-python-first-use-ring-asa.json) and
  [earlier trust-lane JSON](2026-10-18-trust-performance-macos-arm64.json).

- [RDKit Cleanup source sequence](2026-10-18-rdkit-cleanup-source.md):
  [initial full-corpus source run](2026-10-18-rdkit-cleanup-full-corpus-source.json),
  [intermediate source result](2026-10-18-rdkit-cleanup-source.json),
  [pre-metal-scope baseline](2026-10-18-rdkit-cleanup-before-metal-scope.json),
  [integrated source gate](2026-10-18-rdkit-cleanup-source-gate.json),
  [N-oxide result](2026-10-18-rdkit-cleanup-n-oxide-followup.json),
  [coordination-tag report](2026-10-18-rdkit-cleanup-nontetra-followup.md)
  and [result](2026-10-18-rdkit-cleanup-nontetra-followup.json),
  [pyridone report](2026-10-18-rdkit-cleanup-pyridone-followup.md) and
  [result](2026-10-18-rdkit-cleanup-pyridone-followup.json),
  [isolated-metal report](2026-10-18-rdkit-cleanup-isolated-metals-followup.md)
  and [result](2026-10-18-rdkit-cleanup-isolated-metals-followup.json), and
  [alkali/oxygen result](2026-10-18-rdkit-cleanup-alkali-oxygen-followup.json).
  The current three-stage classification supersedes these development epochs;
  they remain indexed for provenance and regression accounting.

- [Terminal N/O and zero-model TPSA](2026-10-18-terminal-tpsa.md):
  isolated RDKit-default TPSA source improvement on both fixed windows; all 42
  median floors and source-output gates pass. Two broader failed candidates,
  remaining COSMolKit deficits and every negative lane are retained.

- [V2000 zero-attribute tail](2026-10-17-openbabel-zero-fields.md):
  SDF writes improve 1.039–1.040x; all 30 source median floors pass, including
  a disclosed 0.6% 3D-write retreat. Five eligible fresh Open Babel lanes win
  24/24; first-slice parse remains unranked. One-host source evidence.

- [3D tetrahedral candidate guard](2026-10-17-openbabel-3d-stereo-guard.md):
  controlled 3D parse improves 6.547–6.641x versus source and 118.775–119.730x
  versus Open Babel; all 24 median floors pass, with a disclosed 0.9% write retreat.

- [V3000 3D coordinate formatting](2026-10-16-openbabel-v3000-3d-format.md):
  controlled 3D writing improves 5.330–5.543x over the source base and
  beats Open Babel by 8.556–8.560x; all 24 source median budgets pass.
- [Scalar descriptor model reuse candidate](2026-10-16-python-scalar-model-reuse.md):
  TPSA→ASA first-use improves 1.013–1.045x, but five parse-inclusive
  medians miss the retention floor; manual integration requested after disclosure,
  with the failed performance decision and all 48 source lanes retained.
- [V3000 requested stereo-loss collection](2026-10-16-openbabel-v3000-loss-work.md):
  prepared writing improves on both diverse slices; all 15 source median
  budgets pass, with a disclosed additional-slice parse retreat.
- [SDF reads and requested stereo-loss collection](2026-10-16-openbabel-sdf-loss-work.md):
  final source passes both-slice read/write/roundtrip improvement conditions
  and all 15 source median budgets; V3000 negative intervals remain disclosed.
- [Direct SD-field candidate](2026-10-16-openbabel-sdf-direct-fields.md):
  first-slice writing and diverse V3000 parsing miss the median floor;
  rejected, with every observation retained.
- [Lazy framing-buffer candidate](2026-10-16-openbabel-sdf-buffer-lazy.md):
  parse/roundtrip gains remain, prepared writing retreats; rejected.
- [Integration after the geometry update](2026-10-16-openbabel-sdf-buffered-main.md):
  initial source fails prepared-write budgets on the newer main baseline.
- [Borrowed SDF record parsing](2026-10-16-openbabel-sdf-buffered.md):
  historical b9f1e3cd baseline: exact retention and 1.9–3.0% parse / 2.4–2.8% roundtrip gains on two
  slices; all 15 source medians within the predeclared 1% budget, including
  a disclosed 0.5% V3000 roundtrip retreat. Source-built, one-host evidence.
- [Ordered-ring and V3000 parser speed follow-ups](2026-10-16-chemistry-speed-followups.md):
  exact output-retention gates and paired source timings; the NCI ring lane
  remains inconclusive and this is not a published-package claim.
- [PDBQT field preservation](2026-10-15-openbabel-pdbqt-fields.md)
  ([raw JSON](2026-10-15-openbabel-pdbqt-fields.json)): two bounded Open Babel
  3.2.1 cases covering an over-width charge and a missing charge; this is not
  a general PDBQT compatibility claim.
- [v1.1.3 evaluation fixes](2026-10-13-v113-evaluation-fixes.md): source
  compatibility, unusual-input parity, and ordinary-input timing.
- [QED stress parity](2026-10-13-qed-stress-parity.md): 4,071 comparable
  unusual-input rows against RDKit 2026.03.6.
- [Batch conformer RMSD](2026-10-14-rdkit-conformer-batch-rmsd.md): focused
  Rust/Python/WASM API parity.
- [Persistent force fields](2026-10-13-cosmolkit-persistent-force-field.md):
  Python reuse API and exploratory timing.
- [Fingerprint NumPy output](2026-10-14-cosmolkit-fingerprint-numpy.md):
  packed-bit equivalence for dense arrays.
- [CDK/Indigo bounded baseline](2026-10-10-cdk-indigo-bounded-baseline.md):
  small common matrix and smoke timing.

Older artifact audits, similarity-search work, streaming/format matrices,
3D experiments, and superseded candidates remain in
[`archive-index.md`](archive-index.md). They are retained for provenance and do
not establish current behavior.

## What counts as a speed win

A lane is called an equivalent-output speed win only when:

1. both engines perform the same declared work;
2. every output satisfies the lane's exact or numeric contract;
3. failures, refusals, unsupported rows, and skipped work are accounted for;
4. the paired 95% confidence lower bound is above `1.0`;
5. source/package, corpus, runtime, host, and raw samples are recorded.

Parse-inclusive, first-use, and prepared-input timings answer different
questions. A source-built result is never presented as a registry-package
result.

## Reproduction and maintenance

Each dated Markdown record contains its exact command, revisions, artifact
hashes, corpus hashes, and validation script. Use that record rather than a
moving "current" command in this index.

- Add each new claim as a Markdown summary plus machine-readable raw data.
- Keep superseded records and label them historical.
- Link every root-level record from this file or `archive-index.md`.
- Run `python3 scripts/check_benchmark_index.py` after changes.
- Do not infer universal accuracy or speed from one corpus, host, or API lane.

- [2026-10-17: Open Babel substituent stack-buffer source experiment](2026-10-17-openbabel-substituent-pair.md) ([numeric packet](2026-10-17-openbabel-substituent-pair.json))

- [2026-10-19: RDKit Cleanup acid-anion transfer charge audit](2026-10-19-rdkit-cleanup-acid-anion-charge-audit.json)
  and [follow-up packet](2026-10-19-rdkit-cleanup-acid-anion-followup.json);
  [acid-priority charge audit](2026-10-19-rdkit-cleanup-acid-priority-charge-audit.json)
  and [follow-up packet](2026-10-19-rdkit-cleanup-acid-priority-followup.json).
