# Native SMILES read-back: stereo-aware aromatic equivalence

## Finding

This source-built audit uses CheMatic main
`c84c68f297f95562d1d29c75234f2190cd5f9d75`, RDKit 2026.09.1 as the
reference parser, and COSMolKit 0.5.0 as a native read-back comparator. Its
source is the 4,072-row `stress-v1.smi` corpus (SHA-256
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
The CheMatic `canonical_roundtrip_dump` binary SHA-256 is
`2df718671c765f2d4c5f05280ae3a14d19cc48b567b297e1cc4689be9be614f9`.
The binary is identical before and after this audit change.

| Engine | Strict RDKit canonical match | Stereo-isomorphic | Raw aromatic graph | Verified equivalent | Unresolved | Invalid |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CheMatic, old oracle | 4,049 | 16 | 6 | 4,071/4,072 | 1 | 0 |
| CheMatic, corrected oracle | 4,049 | 16 | 7 | **4,072/4,072** | 0 | 0 |
| COSMolKit, either oracle | 4,064 | 0 | 0 | 4,064/4,072 | 1 | 7 |

The one changed row is zero-based index 3743. Its source and CheMatic output
use different atom orders for a fused aromatic arsenic ring. Under RDKit
2026.09.1, sanitization assigns different explicit bond orders; strict
canonical SMILES differ. The unsanitized aromatic input graphs, full
stereo-aware canonical SMILES, and molecular formulas agree. The molecule
also has specified tetrahedral stereo outside that ring. The old comparison
gate rejected every aromatic graph check when *any* stereo was present,
which made this row an unverified miss.

The corrected gate requires equal formulas and identical RDKit canonical
isomeric SMILES of the **unsanitized full input graphs**. That equality
includes atom and bond stereochemistry, isotope, atom-map, and connectivity
information. A new positive fixture pins index 3743; a corresponding
inverted stereocenter must remain unresolved. Existing controls for alkene
and tetrahedral inversions, different non-tetrahedral permutations,
connectivity, isotope, atom-map, and directional hydrogen differences pass.

The compressed [old](results/native-smiles-stress-rdkit2026-aromatic-stereo-baseline-2026-10-17.json.gz)
and [corrected](results/native-smiles-stress-rdkit2026-aromatic-stereo-2026-10-17.json.gz)
summaries have SHA-256 values
`b804fec041e48fa9863f6503321ba88e904ddaa02f8fad21c01a09e89563b7f3`
and `a2383478c5183fa72535255180056e397d95783fbe6de7b8e1a3635a5262f3ae`,
respectively. The old-oracle summary was generated from the parent script
with only its RDKit version guard adapted in memory to 2026.09.1. The new
script accepts that version explicitly:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline -j1 \
  -p chematic-smiles --example canonical_roundtrip_dump
python -m unittest discover -s validation/cosmolkit_comparison \
  -p test_native_smiles_equivalence.py
python validation/cosmolkit_comparison/check_native_smiles_roundtrip.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump \
  --expected-rdkit 2026.09.1
```

All five focused oracle tests passed under RDKit 2026.09.1. This corrects
the comparison method; it does not change CheMatic's chemistry or establish
independent chemical truth. RDKit is the reference in this lane, so its own
accuracy is not compared. The fixed development corpus does not establish
toolkit-wide superiority over RDKit or COSMolKit.
