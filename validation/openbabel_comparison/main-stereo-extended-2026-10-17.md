# Open Babel 3.2.1: expanded molecular identity and stereo audit

## Measured result

This source-built audit uses CheMatic main
`c84c68f297f95562d1d29c75234f2190cd5f9d75`, Open Babel 3.2.1,
RDKit 2026.09.1, and, in the 3D holdout, COSMolKit 0.5.0. The CheMatic CLI
SHA-256 is `ff4a44ccc2617bd3a9e291c879f9053a2b13cea744c621dd5acf1d5c2cf6dc39`;
the Open Babel executable SHA-256 is
`733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
These are development source builds, not published-package measurements.

### Paired 2D MOL recognition

The [comparison script](../../scripts/compare_openbabel_stereo_formats.py)
selected all 1,216 distinct RDKit identities with tetrahedral or E/Z syntax
from `stress-v1.smi`, excluding non-tetrahedral `@SP`, `@TB`, and `@OH`
syntax. RDKit wrote each 2D MOL and reread it. A row was eligible only if
the requested MOL version and source identity survived this gate. CheMatic
and Open Babel then read the **same MOL bytes** and RDKit normalized both
output SMILES. Plain explicit-H representation differences were accepted
only after checking that neither structure carried specified stereo or a
labeled/isotopic H.

| Input | Comparable | CheMatic exact | Open Babel exact | Open Babel only | Source exclusions |
| --- | ---: | ---: | ---: | ---: | ---: |
| V2000 MOL | 1,157 | **1,157** | 974 | 0 | 59 |
| V3000 MOL | 1,192 | **1,192** | 936 | 0 | 24 |

CheMatic's pass set strictly contains Open Babel's on both fixed paths.
The [row-level packet](results/openbabel-stereo-2d-stress-1216-main-2026-10-17.json.gz)
has compressed SHA-256
`8f7a6b9c8b3ff3ad68426f71ae991fd1266d41aa714eb380ed47f085b45eecb2`.
It records source exclusions, input MOL hashes, each converter's output,
classifications, tool versions, source revision, binary hashes, and script hash.

### Disjoint unmarked 3D SDF holdout

The [generator](../cosmolkit_comparison/check_unmarked_3d_stereo.py) selected
70 tetrahedral and 70 alkene identities from the browser 10k corpus with
seed `20261021`. Its 140 source indices do not overlap the 400 rows in the
[first](unmarked-3d-browser-holdout-2026-10-17.md) or
[second](stereo-extended-holdouts-2026-10-17.md) browser 3D holdout. RDKit
embedded each source, removed CTAB stereo marks, and admitted it only when
its own reread of both V2000 and V3000 preserved the source identity. All
readers received the **same SDF bytes**.

| Input | CheMatic source identity | Open Babel source identity | COSMolKit source identity |
| --- | ---: | ---: | ---: |
| V2000 SDF | **140/140** | 139/140 | 140/140 |
| V3000 SDF | **140/140** | 139/140 | 140/140 |

Open Babel adds a tetrahedral tag at source index 2844, which has an
unspecified tetrahedral center in the source SMILES. Its mapless graph is
unchanged. CheMatic retains that unspecified center in both source records.
After Open Babel rewrites the SDF, CheMatic and Open Babel score 139/140 on
V2000, while CheMatic scores 140/140 and Open Babel 139/140 on V3000.
Coordinates alone do not prove whether an unspecified source center should
receive a tag; these are source-identity scores.

The [manifest and reader report](results/unmarked-3d-browser-holdout3-2026-10-17.json.gz)
have compressed SHA-256
`eac676c4843c50bc96f7436ee3dd29a9fb5a7e966bbd7e8cc356ad6619ddbbfc`.
The exact inputs are [V2000 SDF](results/unmarked-3d-browser-holdout3-v2000-2026-10-17.sdf.gz)
(`4a9849a2aa79189e471cedbe15b1685a9842e6e70862b0ae3f9e43b60dfab782`)
and [V3000 SDF](results/unmarked-3d-browser-holdout3-v3000-2026-10-17.sdf.gz)
(`9d8ed37b4d51ad84bf92e4246f2d7ac3543dd7878729b0197f6649294e0d8af1`),
with compressed SHA-256 values in parentheses. The CheMatic example binary
SHA-256 is `482bdc248e9b381c5c4813b91fd7a6b6c79e61cede71c7f8f3829afeb349bf69`.

### V3000 writer check

On the separate 316-case selected stereo set, RDKit reread the V3000 blocks
made by each writer. CheMatic retained 315/316 source identities, Open Babel
46/316 by default, and Open Babel 74/316 with `--gen2D`. The only CheMatic
miss, `stress-line-3524`, begins with an unusual `[H]=C` syntax; RDKit could
not read any of the three generated blocks. The [row-level
packet](results/openbabel-v3000-writer-316-main-2026-10-17.json.gz) has
compressed SHA-256
`e0e5a88663df59cfba420032c3d3c203987f79c60adfed5db1c64537a3dec30f`.
This writer comparison uses each engine's own output block, unlike the two
paired-reader comparisons above.

## Reproduce and interpret

From this revision with the pinned tools installed:

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build --offline -j1 \
  -p chematic-cli --bin chematic
python scripts/compare_openbabel_stereo_formats.py \
  --cli target/debug/chematic \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --stress-first-unique 1216 --output /tmp/openbabel-2d-1216.json
python scripts/compare_openbabel_v3000_stereo.py \
  --cli target/debug/chematic \
  --corpus validation/openbabel_comparison/v3000-stereo-writer-stress-316.jsonl \
  --report /tmp/openbabel-v3000-writer-316.json
```

For the 3D lane, decompress both stored SDFs, extract the `manifest` object
from the row-level packet, and run
`validation/openbabel_comparison/check_unmarked_3d_stereo.py` with both
`--case MANIFEST SDF` pairs and the source-built
`target/debug/examples/sdf_v3000_accuracy_dump`. The packet records the
generator/comparator hashes, corpus hash, excluded 400 indices, and per-row
outcomes.

These fixed development sets establish a measured lead or parity on their
declared molecular-identity and stereo workflows. They do not establish
independent chemical truth or toolkit-wide accuracy superiority over Open
Babel. The 3D source identity is an operational contract, and the V3000
writer uses RDKit as its read-back judge.
