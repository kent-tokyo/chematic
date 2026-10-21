# chematic benchmark records

This directory stores dated evidence. Every result is scoped to the named
source or package, comparator, corpus, host, runtime, options, and timed
boundary. Methodology is documented in [`docs/benchmark.md`](../docs/benchmark.md);
the [complete inventory](archive-index.md) links every retained record.

## Start here

Read the boundary column before using a number. `Published` means the named
registry artifact was measured. `Source` means a pinned revision was built
locally. `Historical` means the result remains reproducible evidence but is
not a statement about v1.1.8.

| Question | Current record | Boundary |
|---|---|---|
| Where does source CIP still trail RDKit? | [CIP lazy-preparation audit](2026-10-19-cip-lazy-preparation.md) | 24 balanced blocks; source CIP first-use 1.74× and parse-inclusive 1.62× faster. Standalone CIP/SMILES/SMARTS gaps remain; stress mismatches stay unranked. |
| What is the broadest cross-platform public-wheel chemistry packet? | [v1.1.3 published artifact packet](2026-10-14-v113-published-artifact-packet.md) | Four exact PyPI wheels; npm/WASM, crates.io, speed, and newer releases are separate. |
| What is the broadest published RDKit performance comparison? | [v1.0.42 / RDKit 2026.09.1](2026-10-11-v1042-published-rdkit-packet.md) | Exact PyPI/npm artifacts on hosted Linux x64 and macOS arm64. |
| What are the latest Python descriptor and Morgan source results? | [terminal TPSA](2026-10-18-terminal-tpsa.md) · [Morgan view/ring scratch](2026-10-18-morgan-borrowed-view-ring-scratch.md) | Source wheels on one macOS arm64 host; exact retention gates, not package or toolkit-wide superiority. |
| What do the current parse/first-use/reuse/RSS lanes show? | [published PyPI 1.1.8 Labute ASA](2026-10-19-v118-published-labute-asa-cosmolkit.md) · [lazy-preparation boundary](2026-10-19-python-lazy-preparation.md) · [10k source follow-up](2026-10-19-python-labute-asa-cosmolkit.md) | Published v1.1.8 matches all 10,000 ASA values on one host. The later source candidate shortens parse and parse-inclusive work while moving preparation cost into first use; prepared reuse is approximately flat. Clean-wheel and second-host evidence remain open. |
| How do reaction, standardization, and labelled R-group workflows compare with RDKit? | [current R-group auto-H check](2026-10-18-rdkit-rgroup-auto-hydrogen-followup.md) · [previous source check](2026-10-18-rdkit-workflows-followup.md) · [83 reaction rows](2026-10-18-rdkit-reaction-83-followup.json) · [stereo safety](2026-10-18-rdkit-reaction-stereo-safety-followup.json) · [R-group hydrogen rows](2026-10-18-rdkit-rgroup-hydrogen-followup.json) · [older source check](2026-10-18-rdkit-reaction-cleanup-rgroup-source.md) · [previous reaction rows](2026-10-18-rdkit-reaction-83-source-gate.json) · [previous stereo safety](2026-10-18-rdkit-reaction-stereo-safety.json) · [previous R-group rows](2026-10-18-rdkit-rgroup-hydrogen-source-gate.json) · [expanded workflow packet](2026-10-12-rdkit-workflows-expanded-source.md) | Bounded source results: 80/83 reaction cases have exact graph/origin/map agreement, three inputs are invalid in both engines; the prior 8/8 small R-group Python result remains, while the new automatic-H case has Rust graph/label agreement and awaits Python-extension verification. Broad accuracy remains open. |
| How did earlier Cleanup work progress? | [2026-10-19 source check](2026-10-19-rdkit-cleanup-acid-metal-followup.md) · [4,069-row source result](2026-10-19-rdkit-cleanup-acid-metal-followup.json) · [previous source check](2026-10-18-rdkit-cleanup-alkali-oxygen-followup.md) · [25-row coordination-tag fix](2026-10-18-rdkit-cleanup-nontetra-followup.md) · [earlier source sequence](2026-10-18-rdkit-cleanup-source.md) | Historical source candidate: 3,933/4,069 exact. Use the current standardization row below for the latest Cleanup count. |
| How closely do Normalize, Reionize, Cleanup and CTfile readback match RDKit? | [current standardization source](2026-10-21-rdkit-reionize-catalog-source.md) · [machine record](../validation/results/rdkit-standardization-source-2026-10-21-reionize-catalog.json) · [older source summary](2026-10-18-trust-gates-and-residuals.md) · [CTfile classification](../validation/results/rdkit-molblock-readback-current-2026-10-18.json) | RDKit 2026.09.1, 4,069 comparable rows: Normalize 4,069, Cleanup 4,053, Reionize 4,064 exact. Five Reionize stereo strings remain different. CTfile is a separate lane. |
| Does Cleanup preserve total formal charge? | [historical v1.1.8 charge audit](2026-10-18-rdkit-cleanup-charge-invariant-audit.json) · [current RDKit-contract change](2026-10-19-rdkit-cleanup-acid-metal-followup.md) | The historical audit is superseded for three rows: current source follows pinned RDKit metal representations even where total formal charge changes. This is a versioned compatibility choice, not a general charge model. |
| How do bounded workflows compare with Indigo? | [Indigo workflow packet](2026-10-12-indigo-workflows-source.md) | One-host source lanes with 15 reactions, 39 parent decisions, and 16 labelled R-group rows. |
| What is the current 3D/MMFF boundary? | [frozen holdout score](2026-10-19-trust-holdout-v1-score.md) · [RDKit 2026.09 bounds](2026-10-12-a6-rdkit-2026-09-bounds-followup.md) · [gradient rescue](2026-10-13-a6-mmff-gradient-rescue-development.md) | The one-shot holdout separates embed/stereo success from conformer-count and strict-coordinate mismatches. Two-host timing, memory and unused quality data remain open. |
| How does file interchange compare with Open Babel? | [comparison program](../docs/openbabel-file-io.md) · [SDF follow-up](2026-10-14-openbabel-sdf-read-followup.md) · [V3000 packet](2026-10-14-openbabel-v3000-speed.md) · [stereo preservation](2026-10-15-openbabel-stereo-paired.md) · [historical v1.0.42 same-process baseline](2026-10-10-openbabel-file-io-same-process-v1.0.42.md) | Open Babel 3.2.1 is the file-I/O comparator. Semantic equality is required before speed ranking; diverse residual rows stay unranked. Newer records are one-host source evidence, while the v1.0.42 baseline is historical and must not be presented as current-package performance. |
| How does COSMolKit compare? | [0.5 comparison](2026-10-10-cosmolkit-050.md) · [seeded ETKDG](2026-10-09-cosmolkit-parity-2.md) | Source-built operation and unusual-SMILES evidence; not a published-package or universal claim. |
| How does browser/WASM Morgan compare? | [three-browser packet](2026-10-03-v1030-published-linux-three-browser-paired20.md) | Published npm artifacts, 250 exact rows, and 20 blocks in Chromium, Firefox, and WebKit. |
| What is the broad historical Python operation matrix? | [v1.0.30 63-operation matrix](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | Historical published-wheel result against RDKit 2026.03.6. |

## Recent focused follow-ups

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
