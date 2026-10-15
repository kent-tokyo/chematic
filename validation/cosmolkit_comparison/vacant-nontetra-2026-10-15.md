# Partial non-tetrahedral SMILES read-back, 2026-10-15

This comparison checks molecular recognition and stereochemistry for centers
with exactly one unoccupied coordination site. It does not establish overall
accuracy superiority over RDKit or COSMolKit.

## Fixed inputs and oracle

- Generate 4,080 rows with `corpora/make_vacant_nontetra_corpus.py`: every
  permutation of four distinct TB ligands (F, Cl, Br, I) and all 20 `@TB`
  tags; every permutation of five distinct OH ligands (F, Cl, Br, I, N) and
  all 30 `@OH` tags. The center has one vacant site in every row.
- Corpus SHA-256:
  `979670d2977b47fb6bc7d97f3952f4fff22e46dd110825705f92dd1bc265167e`.
- RDKit 2025.09.3 parses all 4,080 sources. Its canonical isomeric SMILES
  is the expected value. Parse each engine's native SMILES with the same
  RDKit build and compare the resulting canonical isomeric SMILES.
- COSMolKit version: 0.5.0. Comparison script:
  `check_native_smiles_roundtrip.py`. No source rows were skipped.

| Engine | Before | With vacant-site fix |
| --- | ---: | ---: |
| CheMatic | 136 / 4,080 | 4,080 / 4,080 |
| COSMolKit | 4,080 / 4,080 | 4,080 / 4,080 |

The before binary was built from main commit
`235df7a14ea2d11724ccf329fef6675f67a76423` and has SHA-256
`14f9dc54ccc293e7990e3ce153a791ce8764b95ff536908e77454db9fbb68b2d`.
The candidate binary SHA-256 was
`f6a22c6cbc22a0e2dc6c4323981e01dad302dc7576ef75c9c83caf5cd2032b2a`.

The independent 4,072-row `stress-v1.smi` comparison stayed at 4,050
CheMatic matches and 4,064 COSMolKit matches. Its strict string oracle may
count resonance or symmetry-equivalent spellings as differences, so those
22 CheMatic residuals require classification before being called recognition
errors. This fix closes a verified vacant-site defect; a broader superiority
claim is not supported by these results.

## Reproduction

```sh
python3 validation/cosmolkit_comparison/corpora/make_vacant_nontetra_corpus.py \
  --output /private/tmp/chematic-vacant-nontetra.smi
cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline
/private/tmp/chematic-cosmolkit-gap-venv/bin/python \
  validation/cosmolkit_comparison/check_native_smiles_roundtrip.py \
  --corpus /private/tmp/chematic-vacant-nontetra.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump
```
