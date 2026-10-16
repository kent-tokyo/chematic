# V3000 molecular identity: bracket valence and stereochemistry

## Scope

This development comparison measures SMILES → V3000 → RDKit read-back.
RDKit 2025.09.3 first parses each source and produces one canonical isomeric
SMILES, with atom maps removed. **All three writers receive that same
normalized string.** RDKit then reads each V3000 block with
`removeHs=False`; exact canonical isomeric SMILES is the pass condition.
Writer refusal, unreadable output, or changed graph, charge, isotope, or
stereochemistry is a miss. A source is excluded only if RDKit cannot parse it
or its normalized SMILES fails RDKit's own identity read-back.

The previous [316-row report](v3000-three-engine-2026-10-17.md) passed raw
source text to CheMatic but an RDKit-normalized molecule to RDKit's writer.
At `stress-line-3524`, RDKit changes the unusual `[H]=C` spelling into a
different graph before writing; CheMatic preserved the literal double bond
and RDKit refused the resulting CTAB. The comparison script now supplies the
same normalized molecule to every writer. This changes that row's comparison
status, **not** CheMatic's parser or writer behavior on raw `[H]=C` input.

## Measured results

| Pinned corpus | Comparable | CheMatic before → after | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k, full corpus | 10,000 | 9,992 → **10,000** | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 of 4,072 | 3,975 → **4,004** | 4,001 | 3,990 |
| Focused V3000 stereo set | 316 | 316 → **316** | 313 | 312 |

The before/after inputs and RDKit/COSMolKit outputs are identical per row.
CheMatic gained eight exact identities on Browser 10k and 29 on stress-v1,
with zero regressions in either corpus. Among stereo-token sources, the
after scores are 1,687/1,687 versus 1,679/1,687 for each comparator on
Browser 10k; 1,090/1,138 versus 1,084 and 1,078 on stress-v1; and 263/263
versus 261 and 260 on the focused set. The stress corpus excludes source
lines 305 and 1327 because RDKit's normalized SMILES does not preserve its
own first-pass identity on a second read.

The eight Browser 10k gains are one bracketed sulfonyl sulfur and seven
two-bond iodine cases. The broader stress gains include isolated neutral or
charged bracket atoms and bonded B, Si, P, S, and cationic N radicals. The
V3000 writer now emits `VAL=-1` for isolated bracket atoms and uses
`RAD`/`VAL` for these ordinary-bond valence deficits. Aromatic `[nH]` is a
negative control; an initial broad rule caused regressions and was narrowed
before these final measurements.

## Provenance

- Source-built CheMatic CLI baseline from main `458cc19948f3e4a2c6e4d2df876e1495ebeb871b`:
  SHA-256 `e130addd8f720f3b515f505809b268cf10f20da0ea9a5cc6f72f371e333972b7`.
- Source-built candidate CLI: SHA-256
  `e2493adba57f134d791b149fe9a4592d6f54672a3d0f2c3ea7a1efb021671515`.
- RDKit 2025.09.3; COSMolKit 0.5.0. This is a source comparison, not a
  published-package comparison.
- Comparison script SHA-256:
  `ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
- Corpus SHA-256 values: Browser 10k
  `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`,
  stress-v1 `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`,
  focused set `f6728abce73d744260967520b5bdd07d25469c435fb5a875d7b98aaf2b82004c`.
- The deterministic compressed reports contain source text, normalized input,
  expected and observed identities, per-writer V3000 block hashes, exclusions,
  corpus and binary hashes: [Browser baseline](results/v3000-browser-10k-baseline-2026-10-17.json.gz),
  [Browser candidate](results/v3000-browser-10k-candidate-2026-10-17.json.gz),
  [stress baseline](results/v3000-stress-v1-baseline-2026-10-17.json.gz),
  [stress candidate](results/v3000-stress-v1-candidate-2026-10-17.json.gz),
  [focused baseline](results/v3000-stereo-316-baseline-2026-10-17.json.gz),
  and [focused candidate](results/v3000-stereo-316-candidate-2026-10-17.json.gz).
- Compressed report SHA-256 values, in the link order above:
  `5a73ec063b274da0de75a713f0d5d6525804cfccb436b33c46db496e6dd5f0a6`,
  `734ebb4b9151c583a0064a4f2acbca803499fe9b100ee19ec153f661d9480d86`,
  `744326eb88124babf5ea6c8c02c287149d7f23048b2715f94cac0e4ca9fd9e32`,
  `19c66c9c9da8b724f992a45699b3673d0c7425b5be7378025b8ac38407938f96`,
  `4bdb44af79f5c86c94d1564376b378030a21156bb4c6d1be03afa8013bcc818e`,
  `a74df9c55c8e0372942beecc3def686396c35754d77d0942f6d2daec39ca82b2`.

## Reproduction and limit

Build `chematic-cli` from the named source with `cargo build --offline -p
chematic-cli --bin chematic`. With RDKit 2025.09.3 and COSMolKit 0.5.0 in
the same Python environment, run:

```sh
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --cli target/debug/chematic --report /tmp/browser-v3000.json
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --cli target/debug/chematic --report /tmp/stress-v3000.json
python scripts/compare_rdkit_cosmolkit_v3000.py \
  --cli target/debug/chematic --report /tmp/focused-v3000.json
```

These exposed development corpora establish a lead on this one normalized
V3000 interchange path. RDKit is also the source normalizer and read-back
judge; the measurements do not establish independent chemical truth,
published-package accuracy, other format paths, or toolkit-wide superiority.
