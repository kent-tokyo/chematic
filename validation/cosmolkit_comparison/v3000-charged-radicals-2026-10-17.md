# V3000 charged bracket atoms: RDKit and COSMolKit comparison

## Method

This is a source-built development comparison of SMILES → V3000 → RDKit
read-back. `scripts/compare_rdkit_cosmolkit_v3000.py` gives the same
RDKit-normalized canonical isomeric SMILES to all three writers. RDKit
2025.09.3 rereads every V3000 block with `removeHs=False`; exact canonical
isomeric SMILES, with atom maps removed, is the pass condition. COSMolKit is
0.5.0. The normalizer excludes two `stress-v1` sources because its own second
parse changes identity.

| Corpus | Comparable | CheMatic before → after | RDKit | COSMolKit |
| --- | ---: | ---: | ---: | ---: |
| Browser 10k | 10,000 | 10,000 → 10,000 | 9,992 | 9,992 |
| `stress-v1.smi` | 4,070 | 4,004 → **4,008** | 4,001 | 3,990 |
| Focused stereo set | 316 | 316 → 316 | 313 | 312 |

The four stress gains are `line-142` (`[B-]`), `line-1839` (`[O+3]`),
`line-2822` (`[O+]`), and `line-3951` (`[F+3]`). The previous and new reports
have identical expected identities and RDKit/COSMolKit outputs in every row.
CheMatic has zero regressions across the three corpora. Among stereo-token
sources, CheMatic remains at 1,687/1,687, 1,090/1,138, and 263/263,
respectively.

Bracket B−, O+, and F+3 atoms with fewer ordinary bonds than their charged
valence now receive CTfile `RAD` and `VAL` attributes. Without them, RDKit's
V3000 reader invents hydrogens or aromatizes the boron ring. Full-valence
borate and isolated bracket atoms remain negative controls. Source-preserved
opaque `RAD`/`VAL` attributes still take precedence.

## Reproduction and evidence

- Source base: `a4ca27e63bdf689f5f52e36bbddc02d6862f65ab` (the parent of
  this change). Source-built CheMatic CLI SHA-256:
  `2e4a2263c2faf8ae1a7782876b3c46519de3b025e4540f16f0c41c871e4c7562`.
- Comparison script SHA-256:
  `ab8eb6eb3b21ccb7c017d372c0fb05ada78456f2ea4ae26704ba6545f5d35243`.
  Corpus hashes, in table order:
  `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`,
  `879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`,
  `f6728abce73d744260967520b5bdd07d25469c435fb5a875d7b98aaf2b82004c`.
- Per-row compressed reports:
  [Browser](results/v3000-browser-10k-charged-radicals-2026-10-17.json.gz),
  [stress](results/v3000-stress-v1-charged-radicals-2026-10-17.json.gz),
  [stereo](results/v3000-stereo-316-charged-radicals-2026-10-17.json.gz).
  Their SHA-256 values, in the same order, are
  `f83ca35fde8fd622f6723dabc39518432a004da5ec44b7162cd4e779eb16eed3`,
  `2de093646d38b8a9e8f644fbb196960b4ff182c9d26e1330c7e8b3d13b15943b`,
  `f0716114cf23cd72367c312ef18a954fe75867d383b6e8c7e171fa894bd485d9`.
- `cargo test --offline -p chematic-mol --lib`: 798 passed, zero failed.
  `cargo clippy --offline -p chematic-mol --lib -- -D warnings` and
  `cargo fmt --all -- --check` passed.

Run the comparison script with `--cli` pointing to a build of this source and
`--corpus` pointing to each named corpus; omit `--corpus` for the focused set.
The [previous packet](v3000-radical-stereo-2026-10-17.md) contains the
before reports. These exposed corpora establish a lead on one V3000
interchange path. RDKit is also the input normalizer and read-back judge; this
does not establish independent chemical truth, published-package accuracy, or
toolkit-wide superiority.
