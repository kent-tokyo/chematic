# V3000 aromatic-ring wildcards: RDKit and COSMolKit comparison

## Measured result

This source-built development comparison uses RDKit 2025.09.3-normalized
SMILES for CheMatic, RDKit, and COSMolKit 0.5.0. RDKit rereads each writer's
V3000 block with `removeHs=False`. Exact canonical isomeric SMILES after
removing atom maps is the pass condition. The corpus, input gate, and scoring
script are the same as in the [previous packet](v3000-decorated-wildcards-2026-10-17.md).

| Corpus | Comparable | CheMatic before → after | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 10,000 → 10,000 | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 | 4,015 → **4,020** | 4,001 | 3,990 |
| Focused stereo set | 316 | 316 → 316 | 313 | 312 |

The five gains are stress lines 1966, 2712, 2776, 3171, and 3211. They
contain a wildcard inside an aromatic ring. The SMILES parser had treated
that wildcard as nonaromatic and split its ring bonds into single bonds;
the V3000 writer then emitted a CTfile that RDKit could not kekulize. The
parser now recognizes a wildcard with at least two aromatic ring neighbors
in a 4n+2 ring, and the writer emits deterministic single/double Kekulé
bonds for that aromatic wildcard. Pendant wildcards and inter-ring bridges
stay nonaromatic. The ring-size gate avoids promoting seven- and
eight-membered rings that RDKit de-aromatizes.

No CheMatic identity regressed among the 4,070 comparable stress rows. The
browser and focused-stereo rows are identical to the previous packet, and
RDKit and COSMolKit row identities are unchanged. Five stress rows remain
where both other engines pass and CheMatic fails: four sulfoxide stereo
configurations (367, 368, 959, 960) and one charged ring center (2469).

This targeted parser rule does not cover every aromatic wildcard pattern.
In particular, isolated five-membered rings with a wildcard and rings with
adjacent wildcards need separate aromaticity and Kekulé handling.

## Evidence and reproduction

- CLI SHA-256: `fdf5749619572aabd60b8da6966f3cbd8353954b990afe84a24b8532e8e217de`.
  Script SHA-256: `ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
  Each per-row report embeds the corpus hash, exact CLI hash, engine versions,
  expected identity, and each engine's read-back identity.
- [Stress per-row report](results/v3000-stress-v1-ring-wildcards-2026-10-17.json.gz),
  compressed SHA-256 `39e155460dc882557526a8aeac58a510fc692228688850fa2b2c1ef8c56b2fff`.
- [Browser per-row report](results/v3000-browser-10k-ring-wildcards-2026-10-17.json.gz),
  compressed SHA-256 `dc5db54a5b74a7c02c655fdb1ec4f1042ee70caaa1ca26a0a4b02f6d209e15a3`.
- [Focused-stereo per-row report](results/v3000-stereo-316-ring-wildcards-2026-10-17.json.gz),
  compressed SHA-256 `173ac22525132bdae842d58809cbb4ded894aac71f6b2c1e2b2cb2091f62b6d8`.
- `cargo test --offline -p chematic-smiles --lib --quiet`: 400 passed, 1 ignored.
  `cargo test --offline -p chematic-mol --lib --quiet`: 800 passed.
  Clippy for both libraries with `-D warnings`, format check, and diff check passed.

From the repository root, with RDKit 2025.09.3 and COSMolKit 0.5.0:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline \
  -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/v3000-ring-wildcards.json
```

Run the same script with `--corpus
validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi` for Browser 10k;
omit `--corpus` for the focused stereo set. The exposed development corpora
and RDKit read-back judge establish a lead on this V3000 interchange path.
They do not establish independent chemical truth, published-package accuracy,
or toolkit-wide superiority over RDKit or COSMolKit.
