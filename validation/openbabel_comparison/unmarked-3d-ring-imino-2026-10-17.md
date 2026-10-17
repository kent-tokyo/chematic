# Open Babel 3.2.1: unmarked 3D ring and P=N stereo follow-up

This development comparison measures molecular identity after reading the
**same** unmarked 3D SDF bytes. RDKit 2025.09.3 generated each conformer,
removed CTAB stereo marks, and admitted a row only when RDKit read both
V2000 and V3000 blocks back as the source isomeric SMILES. RDKit also
normalizes the CheMatic and Open Babel output SMILES for scoring. These are
source-relative matches, not independent experimental truth.

The 140-row browser selection uses seed `20261021`, 70 tetrahedral and 70
alkene identities, and excludes all 400 source indices used in the first two
browser 3D comparisons. The ChEMBL selection uses seed `20261022`, 100 of
each category, and excludes 32 source indices whose normalized identity
occurs in the combined 540-row browser selection. All 740 selected identities
are distinct. Both source corpora and selection rules are exposed; these are
development checks, not sealed holdouts.

| Input and stage | CheMatic | Open Babel |
| --- | ---: | ---: |
| New browser 140, V2000/V3000, base `ad0283a8` | 139/140 | 140/140 |
| New browser 140, V2000/V3000, after ring fix | **140/140** | 140/140 |
| Disjoint ChEMBL 200, V2000/V3000, after ring fix, before P=N fix | 199/200 | 198/200 |
| Disjoint ChEMBL 200, V2000/V3000, after P=N fix | **200/200** | 198/200 |

The browser miss was a cyclohexane with one substituted CH and one
quaternary carbon carrying two different exocyclic groups. Their ring arms
must be considered together to recover the relative configuration. The
ChEMBL miss was a phosphazene with P=N as one tetrahedral slot; the CTAB
reader had excluded P=N even though its MOL2 path accepted it. Small
coordinate-pinned V2000 tests cover both cases and a ring with three CH
partners. A singly substituted ring and a carbon with duplicate exocyclic
groups remain untagged.

On the two preceding browser selections, this candidate still reads
200/200 for each of V2000 and V3000; Open Babel reads 197/200 and 198/200,
respectively. Across the 740 distinct selected molecules, CheMatic reads
**740/740** and Open Babel **733/740** from original V2000 and V3000
records. The CheMatic pass set contains Open Babel's on these fixed inputs.

An Open Babel rewrite changes the bytes and can lose source stereo. After
Open Babel rewrites V2000, both readers score 733/740. After its V3000
rewrite, CheMatic scores 740/740 and Open Babel 733/740. All 740 molecular
graphs match for both readers on the original records. These results do not
establish toolkit-wide superiority, performance, package-release behavior,
or an independent stereochemical ground truth.

## Pinned evidence

The [selection manifests, baseline and final misses, and provenance](results/unmarked-3d-ring-imino-2026-10-17.json.gz)
have SHA-256 `cbe80283000b438d9d28744f5a46a820969fd0bbdc7176f386d6cfb8acc99d78`.
The exact input SDF files are:

| Selection | V2000 gzip SHA-256 | V3000 gzip SHA-256 |
| --- | --- | --- |
| Browser 140 | [2b44a186](results/unmarked-3d-ring-imino-browser-v2000-2026-10-17.sdf.gz) `2b44a18685053af7f436d5e73d486963372023ab330e375212000e65cc42792f` | [c57a8fea](results/unmarked-3d-ring-imino-browser-v3000-2026-10-17.sdf.gz) `c57a8feaeb8ef02ed1f463a57efa80519e885c92fd842bc81cc3795e0e395368` |
| ChEMBL 200 | [ef8e313f](results/unmarked-3d-ring-imino-chembl-v2000-2026-10-17.sdf.gz) `ef8e313fdeeccf0d6cba6f0a744520b35dc30fca33d1e7863c9a41fa80e4c133` | [4e514a45](results/unmarked-3d-ring-imino-chembl-v3000-2026-10-17.sdf.gz) `4e514a4528f3bfd19f4f0ebd40fea0fe339bfe9c536a2c1696e08da17f274f86` |

The final CheMatic example binary SHA-256 is
`38e3789cb337f93098c93b0afc4cc8c43f2ae0a78f4d936de6dea5ff60353dcf`.
The package also records the exact Open Babel binary and comparator script
hashes. Decompress the SDFs, extract each `manifest` from the JSON packet,
and run `check_unmarked_3d_stereo.py --case MANIFEST SDF` with the built
`sdf_v3000_accuracy_dump` example. The source-gate script is
`validation/cosmolkit_comparison/check_unmarked_3d_stereo.py` (`prepare`),
and the comparator is `validation/openbabel_comparison/check_unmarked_3d_stereo.py`.
