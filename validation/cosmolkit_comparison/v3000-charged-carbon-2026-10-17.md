# V3000 charged carbon valence: RDKit and COSMolKit comparison

## Measured result

This source-built development comparison uses RDKit 2025.09.3-normalized
SMILES for CheMatic, RDKit, and COSMolKit 0.5.0. RDKit rereads each writer's
V3000 block with `removeHs=False`. Exact canonical isomeric SMILES after
removing atom maps is the pass condition. The [previous packet](v3000-sulfoxide-lone-pair-2026-10-17.md)
records the same input gate and scoring method.

| Corpus | Comparable | CheMatic before → after | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 10,000 → 10,000 | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 | 4,024 → **4,025** | 4,001 | 3,990 |
| Focused stereo set | 316 | 316 → 316 | 313 | 312 |
| Charged-carbon valence cases | 4 | 4 | 4 | 4 |

Stress line 2469 is the sole gain. Its `[C@H+3]` ring center had been
written with `CHG=3 HCOUNT=1` alone. RDKit's V3000 reader then dropped
that H and the center's stereochemistry, which also reversed an adjacent
stereo label. Writing `RAD=2 VAL=4` preserves the source identity. The
writer applies RDKit's observed `RAD`/`VAL` encoding to +3 carbon with two,
three, or four drawn valence units; the [four-case fixture](corpora/charged-carbon-valence-4-v1.smi)
checks those states, including a chiral center.

The final CheMatic pass set is a strict superset of each other engine's pass
set on **each of the three fixed development corpora**. In particular, all
45 CheMatic stress misses are also misses for RDKit and COSMolKit. Compared
with the previous packet, only CheMatic's line 2469 identity changes; the
other 4,069 comparable stress rows, all browser and focused-stereo rows,
and every RDKit/COSMolKit row are identical.

## Evidence and reproduction

- CLI SHA-256: `47886d31cfbebbbafa76ecad887d7e1ac9663dc3e33068d243bc29a672f0d21d`.
  Script SHA-256: `ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
  Each per-row report embeds its corpus hash, exact CLI hash, engine versions,
  expected identity, and each engine's read-back identity.
- [Stress per-row report](results/v3000-stress-v1-charged-carbon-2026-10-17.json.gz),
  compressed SHA-256 `2e7553bd28044dfb63ec450be22365bc484a9e97504084a3d03a8cfefd9dbddb`.
- [Browser per-row report](results/v3000-browser-10k-charged-carbon-2026-10-17.json.gz),
  compressed SHA-256 `685a53413c2d8f6c40a8d9fd5649f98fdc7a8a7fe4b0ab9e9e1597e14f771868`.
- [Focused-stereo per-row report](results/v3000-stereo-316-charged-carbon-2026-10-17.json.gz),
  compressed SHA-256 `de070f4caafac2f1617694ce973c352bd5a105108f79c60d3dcbc21a2c330b43`.
- [Four-case corpus](corpora/charged-carbon-valence-4-v1.smi), SHA-256
  `0a500de8bad8040a1679a5cb53ae9740ef0fadd153b3d322e6ff74dda0f5070d`;
  [per-row report](results/v3000-charged-carbon-valence-4-2026-10-17.json.gz),
  compressed SHA-256 `1c032feaaaae8ab23d7332fc8e3033185a70ecb353351a68489e9fc5c5c6a4a5`.
- `cargo test --offline -p chematic-mol --lib --quiet`: 802 passed.
  Clippy with `-D warnings`, format, and diff checks passed.

From the repository root, with RDKit 2025.09.3 and COSMolKit 0.5.0:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline \
  -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/v3000-charged-carbon.json
```

Run the same script with the browser or four-case corpus for those rows;
omit `--corpus` for the focused stereo set. The pass-set relationship is
limited to these exposed development corpora and this V3000 interchange
path, with RDKit as input normalizer and read-back judge. It does not
establish independent chemical truth, published-package accuracy, or
toolkit-wide superiority over RDKit or COSMolKit.
