# Native SMILES read-back on stress-v1, 2026-10-16

This is a pinned molecular-recognition comparison, not a claim of overall
superiority over RDKit or COSMolKit. RDKit supplies the reference molecule and
parses both engines' output; it is not a competing engine in this experiment.

## Inputs and method

- Source: `corpora/stress-v1.smi`, 4,072 RDKit-parseable rows, SHA-256
  `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
- CheMatic: source-built from main `ecda5c5364b1a011431b544ce5d399363d600a6f`;
  `canonical_roundtrip_dump` binary SHA-256
  `f74efb9f30cb4c882dc3e5cc3e5c0af2bc6604b03148230bb42d6ee39948af1a`.
- Comparison runtime: RDKit 2025.09.3, COSMolKit 0.5.0.
- Strict means RDKit canonical isomeric SMILES agrees after parsing native
  output. For strict differences, enumerate up to 10,000 graph mappings and
  renumber atoms before comparing noncanonical stereo SMILES. This checks non-tetrahedral
  stereochemistry as well as tetrahedral and alkene stereochemistry; RDKit's
  `useChirality` substructure match alone accepted a distinct `@SP1`/`@SP2`
  pair in a control. For rows without specified stereochemistry, an intersecting
  RDKit resonance structure set (at most 256 structures per side) also counts
  as verified equivalent. Everything else stays unresolved or invalid.

The following table records this initial method. The follow-up below adds a
check for identical raw aromatic graphs.

| Engine | Strict | Stereo-isomorphic | Resonance | Verified equivalent | Unresolved | Invalid output |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| CheMatic | 4,050 | 16 | 4 | 4,070 / 4,072 | 2 | 0 |
| COSMolKit | 4,064 | 0 | 0 | 4,064 / 4,072 | 1 | 7 |

The 16 CheMatic stereo-isomorphic rows involve symmetry-equivalent ligands at
square-planar, trigonal-bipyramidal, and octahedral centers. The four resonance
rows are indices 1795, 2010, 2138, and 3182. CheMatic indices 1686 and 3599
have different bond-order patterns in aromatic systems;
the bounded resonance search did not establish equivalence in this initial
run, so they were **not** counted as successes. COSMolKit index 304 drops an
explicitly directional H representation and is unresolved; indices 1225, 1326,
1979, 2298, 2572, 3178,
and 3394 produced no RDKit-parseable output. No source rows were skipped.

This corpus supports a measured CheMatic lead over COSMolKit for native SMILES
read-back under this conservative oracle. It does not measure RDKit as a native
round-trip competitor, cover other chemistry tasks, or establish universal
accuracy. The two unresolved CheMatic rows were investigated in the follow-up below.

## Raw aromatic graph follow-up

The source and CheMatic output at indices 1686 and 3599 have the same
canonical **unsanitized aromatic input graph** in RDKit. Their sanitized
formulas also agree, and neither molecule specifies stereochemistry. RDKit's
sanitization chooses different explicit bond orders for these two atom orders
in rings containing an unusual aromatic sulfur or arsenic atom. The earlier
resonance supplier did not cover these forms. The paired gate now checks the
unsanitized aromatic graph only after both strings have parsed successfully,
their formulas agree, and both sanitized molecules lack specified stereo.
Changed connectivity, E/Z and tetrahedral inversions, isotope and atom-map
loss, and explicit-H loss remain negative controls.

| Engine | Strict | Stereo-isomorphic | Raw aromatic graph | Resonance | Verified equivalent | Unresolved | Invalid output |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CheMatic | 4,050 | 16 | 6 | 0 | 4,072 / 4,072 | 0 | 0 |
| COSMolKit | 4,064 | 0 | 0 | 0 | 4,064 / 4,072 | 1 | 7 |

The raw-graph category includes the four rows previously counted as
resonance-equivalent, plus the two previously unresolved rows. Only those two
rows change the total. This is a correction to the comparison method; the
CheMatic parser and writer binary are unchanged (SHA-256
`f74efb9f30cb4c882dc3e5cc3e5c0af2bc6604b03148230bb42d6ee39948af1a`).
The corpus SHA-256 is
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
The comparison used RDKit 2025.09.3, COSMolKit 0.5.0, and the comparison
script based on main `8a1ab7f25fc228d33ecc9f84e5c22cb22f3cd6bb`.
RDKit remains the reference parser in this experiment, not a measured
competitor. This fixed corpus does not establish general accuracy superiority.

## Reproduction

```sh
CARGO_TARGET_DIR=/Users/k_tanabe/Documents/Documents/oss_rust/chematic/target \
  CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 \
  cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline -j 1
/private/tmp/chematic-cosmolkit-gap-venv/bin/python \
  -m unittest discover -s validation/cosmolkit_comparison \
  -p test_native_smiles_equivalence.py
/private/tmp/chematic-cosmolkit-gap-venv/bin/python \
  validation/cosmolkit_comparison/check_native_smiles_roundtrip.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --chematic /Users/k_tanabe/Documents/Documents/oss_rust/chematic/target/debug/examples/canonical_roundtrip_dump
```
