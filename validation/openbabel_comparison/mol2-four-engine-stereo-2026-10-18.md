# Four-reader 3D MOL2 stereo audit (2026-10-18)

## Method and result

RDKit 2026.09.1 parses each source SMILES, embeds an explicit-H 3D conformer,
and checks that its own V2000 read-back retains the source identity. Open
Babel 3.2.1 writes MOL2 from that V2000 record. CheMatic, Open Babel, RDKit,
and COSMolKit 0.5.0 read **the same MOL2 bytes**. The table counts exact
RDKit canonical isomeric SMILES matches to the *source text's molecular
identity*. Failed reads count as misses; failed embeddings, source read-backs,
and MOL2 writes are excluded.

| Development selection | Selected | Comparable | CheMatic | Open Babel | RDKit MOL2 | COSMolKit MOL2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Browser 10k stereo-eligible | 140 | 121 | 118 | 121 | 101 | 101 |
| ChEMBL 5k stereo-eligible, first 160 | 158 | 139 | 139 | 139 | 108 | 108 |
| `stress-v1`, all eligible | 177 | 141 | 95 | 95 | 79 | 79 |

**Exact source identity is not a 3D stereo accuracy oracle.** A source SMILES
can leave E/Z unspecified while its embedded conformer has a particular
geometry. On Browser rows `stress-line-2550`, `2551`, and `3813`, the source
leaves terminal N-H imine E/Z unspecified. Open Babel does too, while the
other readers infer geometry in some or all of those records. Suppressing
CheMatic's geometric inference would turn its 118/121 into 121/121 under
this metric, but it would also erase real stereo from a specified imine.
That candidate change was rejected; the molecule reader remains unchanged.

For a configured control, RDKit embedded opposite `C/C=N/[H]` and
`C/C=N\[H]` stereoisomers and Open Babel wrote their MOL2 files. CheMatic's
direct and source-preserving readers, RDKit, and COSMolKit recover the two
distinct source identities. Open Babel returns `CC=N` for both. The fixture
SHA-256 values are
`8083b8014ecb7c519c86ac520f39a52f312da27b035009cd6fb7d8e762c8478a`
and `c30f953656f038b86700ad018bacd9fe67cea16d6a00c52feae51a6c69e98894`.
The fixture test keeps this configured-stereo behavior from regressing.

On `stress-v1`, RDKit and COSMolKit alone retain the source's non-tetrahedral
tag on four rows (`stress-line-469`, `471`, `472`, and `985`). Each center has
repeated fluorine ligands: five F on trigonal-bipyramidal P, five F on
octahedral S (one case also has Cl), or four F on trigonal-bipyramidal P
(one case also has C). CheMatic and Open Babel leave them unspecified. Their
strict identity difference does not by itself establish a chemically
distinguishable stereoisomer. No RDKit-only or COSMolKit-only strict pass
occurs in the Browser or ChEMBL selections.

The `stress-v1` comparable count differs from an earlier 138-row Open Babel
report because that run used RDKit 2025.09.3 for embedding and source gating;
this run uses RDKit 2026.09.1. Scores across those two runs are not a
before/after comparison.

## Reproduction and evidence

The source-built CheMatic CLI has SHA-256
`2ecc605e29dfb70e6429dc192ea8f2af61bd771f7ecd2030a1c30cdd95087ee0`.
The comparison script SHA-256 is
`6ce8bcbf1ad1c7d2eacb69b4d7f4523cacd67f5d532f5100f86c30cf11e3a71a`.
Corpus SHA-256 values are Browser
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`,
ChEMBL `1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`,
and stress `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.

The compressed per-row reports are [Browser](results/mol2-browser-140-four-engine-2026-10-18.json.gz),
[ChEMBL](results/mol2-chembl-158-four-engine-2026-10-18.json.gz), and
[stress](results/mol2-stress-177-four-engine-2026-10-18.json.gz).
Their SHA-256 values, in that order, are
`c9ce8c8d4661ac663957e5c731c0362125405a8ec2bc70f1662b2f189a0c6534`,
`171e8ff684248f3551231711352ffa03348c4aa53f4b64fbfb9d854848e9b386`,
and `585ca2b5acf4763ce63a22a6b9b298d1775c09cbad3c2a096c6acf2d6e4f236e`.

```sh
python scripts/compare_openbabel_mol2_3d_stress.py --all-stereo \
  --include-rdkit-cosmolkit --cli target/debug/chematic \
  --report /tmp/mol2-stress-four-reader.json
python scripts/compare_openbabel_mol2_3d_stress.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --all-stereo --include-rdkit-cosmolkit --cli target/debug/chematic \
  --report /tmp/mol2-browser-four-reader.json
python scripts/compare_openbabel_mol2_3d_stress.py \
  --corpus scripts/chembl_accuracy_corpus_4999.smi \
  --include-rdkit-cosmolkit --cli target/debug/chematic \
  --report /tmp/mol2-chembl-four-reader.json
```

These exposed development selections are useful for finding discrepancies.
RDKit supplies both source preparation and the judge, Open Babel supplies
the MOL2 writer, and source text may omit stereo visible in 3D coordinates.
The scores do not prove accuracy on unseen molecules, published packages,
other formats, or toolkit-wide superiority over RDKit or COSMolKit.
