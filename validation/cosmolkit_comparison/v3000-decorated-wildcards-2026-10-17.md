# V3000 decorated wildcards: RDKit and COSMolKit comparison

## Measured result

This source-built development comparison uses the same RDKit 2025.09.3-
normalized SMILES for CheMatic, RDKit, and COSMolKit 0.5.0. RDKit rereads
each writer's V3000 block with `removeHs=False`; exact canonical isomeric
SMILES after removing atom maps is the pass condition. The
[charged-bracket comparison](v3000-charged-radicals-2026-10-17.md) describes
the input gate and the previous candidate.

| Corpus | Comparable | CheMatic before → after | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 10,000 → 10,000 | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 | 4,008 → **4,015** | 4,001 | 3,990 |
| Focused stereo set | 316 | 316 → 316 | 313 | 312 |

The seven stress gains are lines 704, 705, 711–714, and 721. They contain
isotopic, charged, or explicit-hydrogen wildcard atoms. The CheMatic writer
previously used V3000 `*` with `MASS`, `CHG`, or `HCOUNT`; RDKit's V3000
reader dropped those attributes from that atom spelling. The writer now uses
`R` only for decorated wildcards and adds `VAL` when an explicit wildcard
hydrogen count needs it. Ordinary undecorated `*` stays unchanged. The
RDKit and COSMolKit row identities are identical to the previous packet,
and CheMatic has zero changed rows outside the seven gains across all three
corpora.

Ten rows still have correct RDKit and COSMolKit read-back but incorrect
CheMatic read-back: four sulfoxide configurations, five aromatic-wildcard
structures that CheMatic does not write, and one charged ring center.

The `R` CTfile spelling can also mean an unnumbered R group. CheMatic's
V3000 reader records that label when it rereads this output. This result
establishes exact molecular SMILES identity on the stated path; it does not
establish lossless R-group annotation round trips.

## Evidence and reproduction

- Source base: PR #970 at `db93de45751a436da7f30a4e856c07e74a3caa48`.
  CheMatic CLI SHA-256:
  `7cb43c4ffb3c6dda2579750333e4401cdbdf59f040c7dae15819a8c5b261c006`.
- Comparison script SHA-256:
  `ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
  The corpus hashes and engine versions are embedded in the per-row report.
- [Stress per-row report](results/v3000-stress-v1-decorated-wildcards-2026-10-17.json.gz),
  compressed SHA-256:
  `f91b633a5586183196b3340b30a9c7db5349cb3045fd56d0de65c9bf1ae943b8`.
  The earlier [charged-bracket packet](results/v3000-stress-v1-charged-radicals-2026-10-17.json.gz)
  is the row-by-row baseline. The browser and focused-stereo results retain
  the same identities as that packet.
- `cargo test --offline -p chematic-mol --lib`: 799 passed, zero failed.
  The six source cases plus bonded `C[*H]` and ordinary `*C` were also
  checked against RDKit's direct V3000 read-back. Clippy and format checks
  passed.

From the repository root, with RDKit 2025.09.3 and COSMolKit 0.5.0:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline \
  -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/v3000-wildcards.json
```

The exposed development corpora and RDKit read-back judge establish a lead
on this V3000 interchange path. They do not establish independent chemical
truth, published-package accuracy, or toolkit-wide superiority over RDKit
or COSMolKit.
