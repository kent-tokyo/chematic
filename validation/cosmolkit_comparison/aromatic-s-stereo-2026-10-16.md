# Aromatic sulfur stereo SMILES preservation (2026-10-16)

RDKit 2025.09.3 enumerates two SMILES for each of five inputs containing an
aromatic `[s+]([O-])` center. The emitted variants differ by `[s@+]` and
`[s@@+]`. On reading those ten emitted SMILES back, CheMatic keeps each pair
distinct and its native canonical writer is idempotent. RDKit and COSMolKit
0.5.0 both collapse every pair to the same untagged SMILES.

| Reader/writer | Distinct pairs after read/write | Stable second CheMatic read/write |
| --- | ---: | --- |
| CheMatic native SMILES | 5 / 5 | 10 / 10 outputs |
| RDKit | 0 / 5 | — |
| COSMolKit | 0 / 5 | — |

The inputs are rows 4385, 4471, 4498, 4503 and 4528 of
`validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi` (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`).
They share the same sulfur motif in five different surrounding structures.
RDKit's own `EnumerateStereoisomers` supplies the two tagged variants for each
input; the test then reads and writes those variants through each toolkit.
The source-built CheMatic binary was from main `025fa2f6` and has SHA-256
`c98af0f6f5492ce45e9233b65cd24fb3fb9443152750e57c50402029b9dd6d85`.

This demonstrates preservation of the stereo distinction encoded in RDKit's
output on these five cases. It does not independently establish the physical
stability or absolute configuration of the proposed sulfur stereoisomers, nor
general recognition superiority across chemistries or formats. RDKit's
enumerator and its SMILES reader disagree here, so a score that normalizes
all emitted SMILES by rereading them with RDKit hides this difference.

## Reproduction

```sh
cargo build -p chematic-smiles --example canonical_roundtrip_dump --offline
python validation/cosmolkit_comparison/check_aromatic_s_stereo.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --chematic target/debug/examples/canonical_roundtrip_dump
```

Use Python with RDKit 2025.09.3 and COSMolKit 0.5.0. The script checks the
corpus hash, tool versions, all five pairs, CheMatic's second-pass stability,
and records the CheMatic binary hash in its JSON output.
