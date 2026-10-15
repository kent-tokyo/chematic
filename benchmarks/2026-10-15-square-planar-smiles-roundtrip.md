# Three-ligand square-planar SMILES round trip — 2026-10-15

This is a source-candidate check of molecular recognition and non-tetrahedral
stereo. The comparison uses RDKit 2025.09.3 as the read-back oracle and
COSMolKit 0.5.0 as the competitor. The input is the committed `stress-v1.smi`
(SHA-256 `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
For every RDKit-valid source, each toolkit writes native SMILES. RDKit parses
that output and canonicalizes it; exact canonical isomeric SMILES equality with
the source counts as a match. Refusals count as misses.

| Source | RDKit-valid rows | CheMatic matches | COSMolKit matches |
| --- | ---: | ---: | ---: |
| Before fix, `ee18d91c`, local diagnostic | 4,072 | 4,049 | 4,064 |
| Candidate with vacant-site remapping | 4,072 | 4,050 | 4,064 |

The pre-fix diagnostic output SHA-256 is
`231665a2a525d40901a859d91c8bfb0e48bc7dd1e13d5f2ede57fdb258d0fe74`;
the candidate output SHA-256 is
`75afa343a66861f9574235a0204bbf52737d05ec2efb0cf14e5118e87fbe6654`.
The candidate diagnostic binary SHA-256 is
`b69265032a4f546505468ae4f99841a24d451ad4dae59b5ed73d9b18cca24d7b`.
The one changed row is index 475: `F[Pt@SP1](Cl)Br`. Previously the native
writer dropped its square-planar tag (`Cl[Pt](F)Br`); it now emits
`Cl[Pt@SP1](F)Br`, which RDKit reads back with the source stereo.

A separate committed three-ligand fixture enumerates 6 ligand orders and all
3 square-planar tags (18 rows; SHA-256
`5d975b17ca15564728f046b3a05bc766c2e9a5a8a72a10bb7052ac7e8c72423b`).
Both CheMatic and COSMolKit match the RDKit read-back oracle on 18/18.

Reproduce from the current source with a Python environment containing the
two pinned comparators:

```sh
cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline
python validation/cosmolkit_comparison/check_native_smiles_roundtrip.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump
python validation/cosmolkit_comparison/check_native_smiles_roundtrip.py \
  --corpus validation/cosmolkit_comparison/corpora/sp3-permutations-v1.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump
```

The remaining 22 stress mismatches include non-tetrahedral centers with
repeated ligands and charged aromatic rings. This test compares RDKit canonical
strings after read-back, so a symmetry-equivalent ligand permutation can count
as a mismatch. The result establishes one corrected stereo-loss case and
does not establish general superiority over either comparator.
