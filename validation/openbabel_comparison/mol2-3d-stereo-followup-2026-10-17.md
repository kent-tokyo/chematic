# Open Babel 3.2.1: explicit-H MOL2 radicals and paired 3D stereo

## Measured result

This source-built development comparison uses the same 3D MOL2 bytes for
CheMatic and Open Babel 3.2.1. RDKit 2025.09.3 embeds each source SMILES,
gates its own V2000 read-back, and judges the two MOL2 reader outputs by
exact canonical isomeric SMILES. Open Babel writes the MOL2; its writer may
discard information before either reader sees it. See the [earlier method
and baseline](mol2-3d-stereo-2026-10-17.md).

| Selection | CheMatic before | CheMatic after | Open Babel |
| --- | ---: | ---: | ---: |
| First 160 distinct stereo identities, 125 comparable | 75/125 | **87/125** | 87/125 |
| All 177 eligible distinct stereo identities, 138 comparable | — | **95/138** | 95/138 |

The first selection gained 12 exact identities with zero regressions. Its
MOL2 hashes, source identities, and Open Babel outputs match the earlier
candidate report row for row. The broader selection adds 17 source identities,
of which 13 pass the source gate; all eight exact identities in those 13 are
shared by both readers. No Open Babel-only row remains in this bounded MOL2
set. The result is **parity**, not a measured MOL2 lead.

Eight gains preserve radicals from explicit-H MOL2 records. If a MOL2 record
contains H atoms, both CheMatic readers now use its atom/bond inventory
instead of adding valence-inferred H to heavy atoms. Records with no explicit
H retain normal valence inference. MOL2 does not declare whether a mixed
explicit/implicit-H record is complete, so this is an interpretation policy,
not proof that every such source intended radicals.

Four gains recover two mixed bridgeheads, three mutually dependent methylated
ring centers, a ring center distinguished by remote stereo, and two P=N
phosphorus centers from nonplanar coordinates. The extra recognition runs
only in the MOL2 reader. A subsequent narrow CTAB refinement is measured in
the [browser-corpus SDF holdout](unmarked-3d-browser-holdout-2026-10-17.md).
Single-methyl cyclohexane, symmetric phosphorus with duplicate methyl
ligands, and adamantane remain untagged in the new tests.

As a cross-format guard, the pinned unmarked 3D SDF set remains at 300/300
exact mapless stereo identities for CheMatic on both V2000 and V3000. Open
Babel scores 290/300 and 288/300 respectively on the same original SDF
bytes. This is an SDF reader lead on this set; it does not change the MOL2
parity result.

## Provenance and reproduction

- Main baseline is the earlier `mol2-3d-stress-candidate-2026-10-17.json.gz`
  report, built before this follow-up. Current source starts from main
  `a4ca27e63bdf689f5f52e36bbddc02d6862f65ab`.
- Current CheMatic CLI SHA-256:
  `6e83bb966b7266d061e29a13ee0c27779b4116361faa773923e1a871fdc67f94`.
  Open Babel binary SHA-256:
  `733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
- Comparison script SHA-256:
  `d7ccce2a70a5894f28e6e11a3ea9e8ca9c299232c26b94a8a0723916ae61669e`.
  Stress corpus SHA-256:
  `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
- [All-stereo MOL2 per-row report](results/mol2-3d-all-stereo-candidate-2026-10-17.json.gz)
  SHA-256:
  `86c66f0b2feadee18cdf7fa451e0703b352cf6d105570f590b8e24335e019f08`.
- [Original SDF per-row report](results/unmarked-3d-sdf-current-2026-10-17.json.gz)
  SHA-256:
  `a2bfedb04082325b61aaa37260f2c53fc3e40427afb5b90c4dda2846e43be827`.
  Its embedded provenance includes the four SDF hashes, reader script hashes,
  Open Babel hash, and CheMatic example hash.
- The four repaired MOL2 fixtures have the same SHA-256 values as source
  lines 651, 658, 2601, and 3599 in the all-stereo report. The new negative
  fixtures were written with RDKit 2025.09.3 and Open Babel 3.2.1.

From the repository root, with the pinned tools installed:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline -j1 -p chematic-cli --bin chematic
python3 scripts/compare_openbabel_mol2_3d_stress.py \
  --all-stereo --cli target/debug/chematic --report /tmp/mol2-all-stereo.json
```

The exposed development inputs and RDKit source gate do not establish
independent chemical truth, published-package behavior, or toolkit-wide
accuracy superiority over Open Babel.
