# Open Babel 3.2.1: extended molecular-identity and stereo checks

These source-built development checks use CheMatic commit
`fc269008d86cead88493da4b3db3d0e4cb1edabe`, Open Babel 3.2.1,
and RDKit 2025.09.3. They extend the [paired 2D MOL comparison](../openbabel_stereo_2d_2026-10-16.md)
and use a new, nonoverlapping selection from the [browser 3D SDF holdout](unmarked-3d-browser-holdout-2026-10-17.md).

## Paired 2D MOL recognition

The first 1,000 distinct source identities with `@`, `/`, or `\` syntax in
`stress-v1.smi`, excluding non-tetrahedral `@SP`, `@TB`, and `@OH` syntax,
were sent through RDKit's 2D MOL writer and reader. An attempt was scored
only if RDKit emitted the requested MOL version and its read-back retained
the source identity. CheMatic and Open Babel then read the **same MOL bytes**.
RDKit normalized both output SMILES. The existing comparator accepts an
ordinary explicit-H representation difference only after proving that neither
molecule has specified stereo or a labeled/isotopic H.

| Input | Comparable | CheMatic equivalent | Open Babel equivalent | Open Babel only |
| --- | ---: | ---: | ---: | ---: |
| V2000 MOL | 958 | **958** | 806 | 0 |
| V3000 MOL | 984 | **984** | 774 | 0 |

The pass set for CheMatic strictly contains Open Babel's pass set for both
formats. RDKit excluded 42 V2000 attempts (33 version fallbacks, eight
read-back identity changes, one writer failure) and 16 V3000 attempts
(15 identity changes, one writer failure). Exclusions are not wins.

The [per-row report](results/openbabel-stereo-2d-stress-1000-2026-10-17.json.gz)
contains source identities, input MOL hashes, both converter outputs,
classification, binary hashes, and script hash. Its compressed SHA-256 is
`8657c6b0572a7cf3a1b8172a80fb6e344574811623ef25c64403ec134bb29c82`.
The CheMatic CLI SHA-256 is
`b37a77333e5d134678145d516e55571d216c2ae551f806bebe19228cc3b29f07`;
the Open Babel binary SHA-256 is
`733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.

## Independent browser 3D SDF selection

A second seed (`20261020`) selected 100 tetrahedral and 100 alkene-stereo
identities from the browser 10k corpus, excluding all 200 source indices in
the preceding browser holdout. RDKit embedded each molecule in 3D, removed
CTAB stereo markers, and admitted only V2000/V3000 records that RDKit read
back with the source identity. All engines read the **same SDF bytes**.

| Input | CheMatic | Open Babel | COSMolKit |
| --- | ---: | ---: | ---: |
| V2000 SDF, mapless stereo | 198/200 | 198/200 | 200/200 |
| V3000 SDF, mapless stereo | 198/200 | 198/200 | 200/200 |

The same counts hold after Open Babel rewrites each SDF. Both CheMatic and
Open Babel add a tetrahedral assignment at the same two centers left
unspecified by the source (browser rows 2724 and 2687). Both preserve the
mapless graph on all 200 rows. Coordinates alone do not establish whether
those source-unspecified centers should remain unspecified, so this check
records parity rather than calling either reader chemically wrong.

The [manifest and per-row report](results/unmarked-3d-browser-holdout2-2026-10-17.json.gz)
has compressed SHA-256
`7b1185c7003bbc1d548aa80cf1bdea3630404cf955e4fd9e425ef46185cef9e5`.
It records both script hashes, the reader binaries, source indices,
RDKit/COSMolKit versions, and input SDF hashes. Exact inputs are
[V2000 SDF](results/unmarked-3d-browser-holdout2-v2000-2026-10-17.sdf.gz)
(`138db6c0ea38617995e883b7ccddd8d3c854549de91efa1f66c3621d9e0b778a`)
and [V3000 SDF](results/unmarked-3d-browser-holdout2-v3000-2026-10-17.sdf.gz)
(`9bd04c3ce1ea3481af3513b891f933ca410e0702caf08bef4a3da6b0f54cc3d5`),
with compressed SHA-256 values shown in parentheses.

## Reproduce

From the repository root, with the pinned dependencies installed:

```sh
cargo build --offline -p chematic-cli --bin chematic
python3 scripts/compare_openbabel_stereo_formats.py \
  --cli target/debug/chematic \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --stress-first-unique 1000 --output /tmp/openbabel-stereo-1000.json

cargo build --offline -p chematic-mol --example sdf_v3000_accuracy_dump
python3 validation/cosmolkit_comparison/check_unmarked_3d_stereo.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --chematic target/debug/examples/sdf_v3000_accuracy_dump \
  --output-prefix /tmp/browser-holdout2 --seed 20261020 --per-category 100 \
  --exclude-manifest /tmp/previous-browser-manifest.json
python3 validation/openbabel_comparison/check_unmarked_3d_stereo.py \
  --case /tmp/browser-holdout2-manifest.json /tmp/browser-holdout2-v2000.sdf \
  --case /tmp/browser-holdout2-manifest.json /tmp/browser-holdout2-v3000.sdf \
  --chematic target/debug/examples/sdf_v3000_accuracy_dump \
  --output-prefix /tmp/browser-holdout2-openbabel --rewrite
```

Extract `/tmp/previous-browser-manifest.json` from the `manifest` field of
the preceding holdout's report before regenerating the new selection.
These results establish a lead on the paired 2D MOL path and parity on this
3D SDF path. They do not establish independent chemical ground truth,
published-package behavior, or toolkit-wide superiority over Open Babel.
