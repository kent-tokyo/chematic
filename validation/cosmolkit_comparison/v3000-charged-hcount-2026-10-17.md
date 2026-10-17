# Charged-atom HCOUNT in concrete V3000 output

## Result

This source-built development comparison uses RDKit 2026.09.1 to normalize
each input and read each writer's V3000 block back. An exact canonical
isomeric SMILES match is required. COSMolKit 0.5.0 also reads CheMatic's
block; its result is normalized by RDKit for a separate cross-reader score.
The candidate omits `HCOUNT` on charged organic-subset atoms only when
CheMatic's bond and valence inference gives the same hydrogen count. It keeps
the property when it carries information, including `[CH2-]`, `[SH2+]`, and
aromatic `[cH-]` in the tested cases.

| Corpus | Comparable | CheMatic writer → RDKit | RDKit writer → RDKit | COSMolKit writer → RDKit | CheMatic writer → COSMolKit, before → after |
| --- | ---: | ---: | ---: | ---: | ---: |
| Stress v1 | 4,070 | 4,025 | 3,995 | 3,990 | **3,345 → 3,504** |
| Browser 10k | 10,000 | 10,000 | 9,991 | 9,992 | **9,952 → 9,957** |
| Focused stereo | 316 | 316 | 311 | 312 | **258 → 276** |
| ChEMBL | 5,000 | 5,000 | 4,981 | 4,993 | 5,000 → 5,000 |

The candidate changed CheMatic blocks on 622, 880, 45, and 309 rows,
respectively. Cross-reader gains were 159, 5, 18, and 0; no cross-reader
losses or CheMatic-to-RDKit identity changes occurred on these rows. On each
corpus, the CheMatic-to-RDKit pass set contains both other writers' pass sets.
The ChEMBL and browser/focused baselines are earlier tracked reports on the
same corpus hashes and RDKit version; the stress baseline was freshly built
from main `c84c68f297f95562d1d29c75234f2190cd5f9d75` with CLI SHA-256
`5352fbe866eefff5f5006b8d96e0c1388d14e725426fe8c6befdd84cf9f913d0`.

## Evidence

The candidate CLI SHA-256 is
`98731fd8d23d32c8f7946b387103c19069bc3a0a6237cabd223d1869ea93c8a5`.
Each packet records the corpus hash, versions, CLI hash, expected identity,
read-back identities, and block hashes per row.

| Packet | Compressed SHA-256 |
| --- | --- |
| [Stress baseline](results/v3000-stress-4070-charged-hcount-baseline-2026-10-17.json.gz) | `8d4a2afed353e19307f503c03cfb4ec3452c366db8e4a578e025f57a451802e5` |
| [Stress candidate](results/v3000-stress-4070-charged-hcount-2026-10-17.json.gz) | `5d8c3c5906ee2ba917dd2a2e317413431bf6f7213c1a1d8b7d1d74dbf1ef6dbd` |
| [Browser candidate](results/v3000-browser-10k-charged-hcount-2026-10-17.json.gz) | `320e8209b4066bfb82d56ab9eb536bb172537244a6aaf57e8251ef4cbbbcf5a7` |
| [Focused stereo candidate](results/v3000-stereo-316-charged-hcount-2026-10-17.json.gz) | `d6ee4ae30d2d8ae605ab8654de9a904e5d6e421b4a5d2e9c07c9ee5701b09610` |
| [ChEMBL candidate](results/v3000-chembl-5000-charged-hcount-2026-10-17.json.gz) | `c721b6db778f2022cb4ca4278addac6afdba74635d95bec05714e73286543afa` |

The comparison script is `scripts/compare_rdkit_cosmolkit_v3000.py`.
From the repository root, with RDKit 2026.09.1 and COSMolKit 0.5.0 installed:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline -j1 \
  -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/v3000-charged-hcount.json
```

Use `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`,
`validation/openbabel_comparison/v3000-stereo-writer-stress-316.jsonl`, or
`scripts/chembl_accuracy_corpus_4999.smi` for the other packets. Local checks:
`cargo test --offline -j1 -p chematic-mol --lib --quiet` (808 passed),
`cargo clippy --offline -j1 -p chematic-mol --all-targets -- -D warnings`,
`cargo fmt --all -- --check`, the refactor-budget script, and
`git diff --check` passed.

This result is about V3000 interchange on exposed development corpora.
RDKit supplies both the input normalizer and read-back judge, so these
counts are not independent chemical truth or evidence of toolkit-wide
superiority over RDKit or COSMolKit. COSMolKit still misses 566 CheMatic
blocks in the stress cross-reader lane and 40 in the focused-stereo lane.
