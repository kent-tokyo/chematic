# Unmarked 3D SDF stereochemistry comparison (2026-10-16)

## Scope and method

`check_unmarked_3d_stereo.py` selects stereochemical rows from
`corpora/stress-v1.smi` (SHA-256
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
RDKit 2025.09.3 embeds each graph with ETKDGv3, removes hydrogens and all
CTAB atom-parity, wedge and CFG markers, then writes both V2000 and V3000.
Any record whose *actual* output format differs from the requested format, or
whose RDKit reread changes canonical isomeric SMILES, is excluded. The two
seeds select disjoint source rows: 50 tetrahedral and 50 alkene records at
seed `20261016`, then 100 of each at seed `20261017` after excluding the
first manifest's source indices. All engine outputs are normalized with the
same RDKit SMILES parser before exact comparison.

CheMatic is the source-built `sdf_v3000_accuracy_dump` example; COSMolKit is
Python 0.5.0. The tested CheMatic binary SHA-256 is
`85ff5a126b3daf9eb89649d42293ded4b5f720a30336c385a94e43b6ed6d69da`.
The compared source is this change to the MOL readers, before integration to
main; this is not a published-wheel or crate result.

| Actual format | Engine | Tetrahedral | Alkene | Total |
| --- | --- | ---: | ---: | ---: |
| V2000 | CheMatic | 150/150 | 150/150 | 300/300 |
| V2000 | COSMolKit | 150/150 | 148/150 | 298/300 |
| V3000 | CheMatic | 150/150 | 150/150 | 300/300 |
| V3000 | COSMolKit | 150/150 | 150/150 | 300/300 |

The two COSMolKit V2000 misses (source rows 716 and 1705) return
`expected SDF molecule payload, found query_graph` for wildcard-containing
records. Thus all *shared molecular payloads* agree in this sample. RDKit
is the source and reread gate, so this comparison cannot establish that
CheMatic exceeds RDKit's recognition accuracy. The first 100-row set gave
CheMatic 83/100 per format before these reader changes; it now gives
100/100. The 200-row set was selected separately and gives 200/200.

## Reproduction

From this directory, with RDKit 2025.09.3 and COSMolKit 0.5.0 installed:

```sh
cargo build -p chematic-mol --example sdf_v3000_accuracy_dump --offline
python check_unmarked_3d_stereo.py --chematic ../../target/debug/examples/sdf_v3000_accuracy_dump --output-prefix /tmp/chematic-3d-a --seed 20261016 --per-category 50
python check_unmarked_3d_stereo.py --chematic ../../target/debug/examples/sdf_v3000_accuracy_dump --output-prefix /tmp/chematic-3d-b --seed 20261017 --per-category 100 --exclude-manifest /tmp/chematic-3d-a-manifest.json
```

The generated SDF SHA-256 hashes are:

| Seed | V2000 | V3000 |
| --- | --- | --- |
| 20261016 | `e2373a90a00290d15175f03230ff296ccfa4ae4f1cf0e7c08658948473c81b69` | `cd6d727b5c101de6e792afbcc9c090496b627971d6616d8e8151410ea581c767` |
| 20261017 | `16ac9819a48e924eaa59d00999c26e84b47b49ecc8414487020d7bfb4fe2a50f` | `42c55c1df42c73884219cdfca704a931e2ed1140dfd31dc339e3406b1aa705f9` |

This test covers 300 generated, source-gated molecules and two SDF dialects.
It does not measure all chemistries, coordinate quality, speed, or every
competing API behavior.
