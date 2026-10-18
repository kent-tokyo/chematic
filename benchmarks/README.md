# chematic benchmark records

This directory stores dated evidence. Every result is scoped to the named
source or package, comparator, corpus, host, runtime, options, and timed
boundary. Methodology is documented in [`docs/benchmark.md`](../docs/benchmark.md);
the [complete inventory](archive-index.md) links every retained record.

## Start here

| Question | Current record | Boundary |
|---|---|---|
| Do hosted Python packages reproduce the chemistry contract? | [v1.1.3 published artifact packet](2026-10-14-v113-published-artifact-packet.md) | Four exact PyPI wheels; npm/WASM, crates.io, and speed are separate. A complete v1.1.4 package packet is still open. |
| What is the broadest published RDKit performance comparison? | [v1.0.42 / RDKit 2026.09.1](2026-10-11-v1042-published-rdkit-packet.md) | Exact PyPI/npm artifacts on hosted Linux x64 and macOS arm64. |
| What are the latest Python descriptor and Morgan source results? | [sanitized-view reuse](2026-10-14-python-sanitized-view-reuse.md) · [speed frontier](2026-10-14-python-speed-frontier.md) | Source wheels on one macOS arm64 host; exact retention gates, not package or toolkit-wide superiority. |
| How do reaction, standardization, and labelled R-group workflows compare with RDKit? | [expanded workflow packet](2026-10-12-rdkit-workflows-expanded-source.md) | Source-built, bounded fixtures, 21 alternating blocks; broader MolStandardize/R-group parity remains open. |
| How closely does Cleanup match RDKit? | [Cleanup comparison](2026-10-18-rdkit-cleanup-source.md) · [254-row result](2026-10-18-rdkit-cleanup-source.json) · [full-corpus result](2026-10-18-rdkit-cleanup-full-corpus-source.json) | Source candidate: 240/254 and 3,826/4,059 exact; metal and reionization gaps remain. |
| How do bounded workflows compare with Indigo? | [Indigo workflow packet](2026-10-12-indigo-workflows-source.md) | One-host source lanes with 15 reactions, 39 parent decisions, and 16 labelled R-group rows. |
| What is the current 3D/MMFF boundary? | [RDKit 2026.09 bounds](2026-10-12-a6-rdkit-2026-09-bounds-followup.md) · [gradient rescue](2026-10-13-a6-mmff-gradient-rescue-development.md) | Frozen development corpora; quality tails, fresh holdout, two-host timing, memory, and exact packages remain open. |
| How does file interchange compare with Open Babel? | [V3000 packet](2026-10-14-openbabel-v3000-speed.md) · [SDF follow-up](2026-10-14-openbabel-sdf-read-followup.md) · [program](../docs/openbabel-file-io.md) | Semantic equality is required before speed ranking; diverse residual rows stay unranked. |
| How does COSMolKit compare? | [0.5 comparison](2026-10-10-cosmolkit-050.md) · [seeded ETKDG](2026-10-09-cosmolkit-parity-2.md) | Source-built operation and unusual-SMILES evidence; not a published-package or universal claim. |
| How does browser/WASM Morgan compare? | [three-browser packet](2026-10-03-v1030-published-linux-three-browser-paired20.md) | Published npm artifacts, 250 exact rows, and 20 blocks in Chromium, Firefox, and WebKit. |
| What is the broad historical Python operation matrix? | [v1.0.30 63-operation matrix](2026-10-02-v1.0.30-vs-rdkit-python-63op-paired20.json) | Historical published-wheel result against RDKit 2026.03.6. |

## Recent focused follow-ups

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
