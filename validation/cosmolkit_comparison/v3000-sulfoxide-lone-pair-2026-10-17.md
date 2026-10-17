# V3000 pyramidal sulfur stereo: RDKit and COSMolKit comparison

## Measured result

This source-built development comparison uses RDKit 2025.09.3-normalized
SMILES for CheMatic, RDKit, and COSMolKit 0.5.0. RDKit rereads each writer's
V3000 block with `removeHs=False`. Exact canonical isomeric SMILES after
removing atom maps is the pass condition. The [previous packet](v3000-ring-wildcards-2026-10-17.md)
records the same corpus, input gate, and scoring method.

| Corpus | Comparable | CheMatic before → after | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 10,000 → 10,000 | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 | 4,020 → **4,024** | 4,001 | 3,990 |
| Focused stereo set | 316 | 316 → 316 | 313 | 312 |
| Sulfoxide variants | 12 | 12 | 12 | 12 |

The four stress gains are lines 367, 368, 959, and 960: the two isomeric
forms of methyl phenyl sulfoxide represented as neutral `S=O` and as charged
`[S+]([O-])`. CheMatic had chosen the opposite V3000 wedge from the one
RDKit reads as the source configuration. The local three-bond parity routine
used an implicit-H-like fourth site for pyramidal sulfur, whereas its
stereogenic lone pair belongs before the three real neighbors for wedge
selection. The strict CTAB reader now stores the same three-real-neighbor
order as the SMILES parser. Twelve additional neutral and charged sulfur
variants, both isomers in each pair, pass both RDKit read-back and CheMatic's
own SMILES → V3000 → SMILES round trip.

No CheMatic identity regressed among the 4,070 comparable stress rows. The
browser and focused-stereo rows are identical to the previous packet; the
RDKit and COSMolKit rows are unchanged in all three corpora. Stress line
2469, a charged ring stereocenter, remains the one row where both other
engines pass and CheMatic fails.

The aromatic-wildcard parser now performs its ring check on a detached
copy, keeping perception's type-erased ring caches off the returned parsed
graph. The `contiguous_component_facts_match_original` regression test
passes after that repair.

## Evidence and reproduction

- CLI SHA-256: `c14dc9e94f52fc2cbaa3874ecc85df6e27ca2fa852a25bec706f540ac9ed41c5`.
  Script SHA-256: `ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
  Each per-row report embeds its corpus hash, exact CLI hash, engine versions,
  expected identity, and each engine's read-back identity.
- [Stress per-row report](results/v3000-stress-v1-sulfoxide-lone-pair-2026-10-17.json.gz),
  compressed SHA-256 `08ee6f315ba3292b0dfbeb7cf5bff976fb7dc622f7b6e286c08b8f60856e037b`.
- [Browser per-row report](results/v3000-browser-10k-sulfoxide-lone-pair-2026-10-17.json.gz),
  compressed SHA-256 `3b8226b480062fd6b4e253a669d8d529c7a2bb70223866cb542ed628098f7ce7`.
- [Focused-stereo per-row report](results/v3000-stereo-316-sulfoxide-lone-pair-2026-10-17.json.gz),
  compressed SHA-256 `9327433ec976c4fc3da4e0b37c795236ad0e9c6fa936f6a2d175c6398fef1b14`.
- [12-variant corpus](corpora/sulfoxide-stereo-12-v1.smi), SHA-256
  `10adb2a6127ec9434b39d8acbd2369de685b8d08a1c87a902b452b0016c55b29`;
  [per-row report](results/v3000-sulfoxide-12-2026-10-17.json.gz), compressed SHA-256
  `974fa36490f2d30db261f7612f477d14d84ccc1d61a84053b14b4f47f9bd12bb`.
- `cargo test --offline -p chematic-perception --lib --quiet`: 272 passed,
  1 ignored. `cargo test --offline -p chematic-mol --lib --quiet`:
  801 passed. `cargo test --offline -p chematic-smiles --lib --quiet`:
  400 passed, 1 ignored. Clippy with `-D warnings`, format, and diff checks passed.

From the repository root, with RDKit 2025.09.3 and COSMolKit 0.5.0:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline \
  -p chematic-cli --bin chematic
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/v3000-sulfoxide.json
```

Run the same script with the browser or 12-variant corpus for those rows;
omit `--corpus` for the focused stereo set. The exposed development corpora
and RDKit read-back judge establish a lead on this V3000 interchange path.
They do not establish independent chemical truth, published-package accuracy,
or toolkit-wide superiority over RDKit or COSMolKit.
