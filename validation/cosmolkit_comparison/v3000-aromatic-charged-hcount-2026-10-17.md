# Charged aromatic carbon HCOUNT in concrete V3000 output

## Result

This source-built comparison starts from main `c6d10dc0218f190f88fb01f1e43fee0631575c22`.
RDKit 2026.09.1 normalizes each source SMILES and reads each writer's V3000
block back. COSMolKit 0.5.0 separately reads CheMatic's block. Exact canonical
isomeric SMILES after removing atom maps is the pass condition. The fixed
`stress-v1.smi` corpus has 4,072 rows, of which 4,070 are comparable.

| Path | Before | After |
| --- | ---: | ---: |
| CheMatic V3000 → RDKit | 4,025 / 4,070 | 4,025 / 4,070 |
| CheMatic V3000 → COSMolKit | 3,504 / 4,070 | **3,523 / 4,070** |
| RDKit V3000 → RDKit | 3,995 / 4,070 | 3,995 / 4,070 |
| COSMolKit V3000 → RDKit | 3,990 / 4,070 | 3,990 / 4,070 |

Exactly 19 CheMatic blocks changed. Each is a charged aromatic carbon with one
hydrogen and two ring bonds, including cyclopentadienyl anions and tropylium
cations. Removing the redundant `HCOUNT=1` lets COSMolKit read all 19 as
concrete molecules; RDKit still reads them as the original molecules. There
were no changed identities or block hashes on the other 4,051 comparable rows.
The ChEMBL 5,000, browser 10,000, and focused stereo 316 comparison packets
contain no RDKit-normalized `[cH-]` or `[cH+]` inputs, so this rule does not
affect their measured cases.

The writer retains `HCOUNT` for cases where it carries information, including
`[CH2-]` and `[SH2+]`. Focused Rust tests cover both changed ions and these
negative controls. The CheMatic MOL reader still renders the charged rings in
Kekulé notation; it does so with both the old and new blocks. That existing
aromaticity-recognition limit is separate from this writer interoperability
change.

## Reproduction

The source corpus SHA-256 is
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
The candidate CLI SHA-256 is
`005338829bdeb3c9bce852d51883498bbe96bf450df0676b01db964efc4e890f`;
the comparison script SHA-256 is
`fe323c8bd00a775a1e0e640814db90c5e9bcc401afce75d0e7cc3cafac64ca9c`.
The [baseline](results/v3000-stress-4070-charged-hcount-2026-10-17.json.gz)
and [candidate](results/v3000-stress-4070-aromatic-charged-hcount-2026-10-17.json.gz)
retain per-row inputs, identities, block hashes, versions, and the CLI hash.
Their compressed SHA-256 values are
`5d8c3c5906ee2ba917dd2a2e317413431bf6f7213c1a1d8b7d1d74dbf1ef6dbd`
and `f702961f49d06676880cc0ada6090fff2f4cda2fc5b830ad13df1eb76b9955dc`.

From the repository root, with RDKit 2026.09.1 and COSMolKit 0.5.0:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline -j1 \
  -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/v3000-aromatic-charged-hcount.json
```

The local `chematic-mol` library suite passed 808/808 tests. This is a
cross-reader gain on an exposed development corpus. RDKit supplies both the
input normalizer and read-back judge; this result does not establish
independent chemical truth or toolkit-wide superiority.
