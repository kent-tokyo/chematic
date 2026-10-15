# V3000 SDF recognition and stereochemistry check — 2026-10-15

Scope: RDKit 2025.09.3 creates 2D V3000 MOL blocks from each source SMILES.
A record enters the comparison only when RDKit reads its own block back to
the same canonical isomeric SMILES. The same eligible blocks go to CheMatic's
`SdfDataset`, RDKit, and COSMolKit 0.5.0. CheMatic and COSMolKit output is
parsed by RDKit with `removeHs=False`, then compared to the source; atom counts
must also agree. A parser refusal or invalid output is a miss.

| Source corpus | Input | Eligible | RDKit source changed | RDKit write/read error | CheMatic | COSMolKit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `validation/cosmolkit_comparison/corpora/stress-v1.smi` | 4,072 | 4,004 | 52 | 16 | 4,004 / 4,004 | 4,004 / 4,004 |
| `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi` | 10,000 | 9,992 | 8 | 0 | 9,992 / 9,992 | 9,992 / 9,992 |

The corrections covered three independent gaps seen during diagnosis:

- CTAB wedges were applied at the incoming end as well as their owning atom,
  causing adjacent stereocentres in a cyclic sugar to lose parity.
- Wedges on pyramidal S and ring-constrained three-coordinate N were rejected
  despite source-preserving RDKit V3000 blocks.
- `VAL=` on a V3000 dative donor counted the coordinate bond against donor
  valence, removing one hydrogen from four stress records.

Reproduce with the current source and a Python environment containing the
versions above:

```sh
cargo build -p chematic-mol --example sdf_v3000_accuracy_dump --offline
python validation/cosmolkit_comparison/check_sdf_v3000_accuracy.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --chematic target/debug/examples/sdf_v3000_accuracy_dump
python validation/cosmolkit_comparison/check_sdf_v3000_accuracy.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --chematic target/debug/examples/sdf_v3000_accuracy_dump
```

The comparison program prints corpus, generated SDF, and binary SHA-256
hashes with versions, skip counts, and up to 20 mismatches per engine. On this
run the binary hash was
`7a7458898ed90f977de34a5e7ad1726fa9b0d4d7ab6f005545079d2cbf71dd26`.
The generated SDF hashes were
`a1f1b2d429590e3efeb25a9f9041a15f69592eb89ba918f765352d748fc993fe`
for stress and
`d9cf7fa54300d829e5fca58abb710d0483d8f22ed98ee36367e504d556b0db64`
for RDKit.js 10k.

This is a source-preserving SDF read test, not evidence of general molecular
recognition superiority. The 60 RDKit-source-changed records and 16 RDKit
write/read errors are excluded rather than assigned to either competitor.
