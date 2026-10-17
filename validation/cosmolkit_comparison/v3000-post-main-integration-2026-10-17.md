# V3000 identity after integrating main

The branch merged `origin/main` at `de41dadd8f123a4b1368937c7dd86df823434588`
into `b349b80750250d9c3807a967e185cb3e3cbe86a8`. This check rebuilds
the CLI from that merge commit and repeats the three-engine comparison in the
[charged-carbon packet](v3000-charged-carbon-2026-10-17.md). The input
normalizer and V3000 read-back judge are RDKit 2025.09.3; COSMolKit is 0.5.0.

| Corpus | Comparable | CheMatic | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 10,000 | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 | 4,025 | 4,001 | 3,990 |
| Focused stereo | 316 | 316 | 313 | 312 |

All expected and read-back identities in the three new reports match the
pre-merge reports row by row. CheMatic's exact-pass set strictly contains
the other engines' exact-pass sets on each corpus. The 45 CheMatic stress
misses are shared with both comparison engines; many involve non-tetrahedral
stereochemistry outside this V3000 writer's current encoding.

The rebuilt CLI SHA-256 is
`b37a77333e5d134678145d516e55571d216c2ae551f806bebe19228cc3b29f07`.
The comparison script SHA-256 is
`ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
Each report embeds the corpus hash, CLI hash, engine versions, input identity,
and all three read-back identities:

- [Stress rows](results/v3000-stress-v1-post-main-2026-10-17.json.gz),
  compressed SHA-256 `48453ed41022ba5efe3366544fc4944a1162103dfe317211d8f27499c7aec05f`.
- [Browser rows](results/v3000-browser-10k-post-main-2026-10-17.json.gz),
  compressed SHA-256 `35900ff452b3e283b190b58bba443335d26b83f9fce4f7b629968789fe77bd8b`.
- [Focused-stereo rows](results/v3000-stereo-316-post-main-2026-10-17.json.gz),
  compressed SHA-256 `cbd3f8f7353ffec479d5523b2fc34fad03d2265b67e7c45f9d688ca1745683d7`.

Post-merge local checks passed: `chematic-mol` 802 tests, `chematic-perception`
272 tests with one ignored, `chematic-smiles` 400 tests with one ignored,
Clippy with `-D warnings` for these three libraries, `cargo fmt --all -- --check`,
and `git diff --check`.

This is a source-built development result for these fixed, exposed corpora
and one SMILES-to-V3000 interchange path. RDKit supplies both normalized
input and the read-back identity, so it is not independent chemical truth
or evidence of toolkit-wide superiority.
