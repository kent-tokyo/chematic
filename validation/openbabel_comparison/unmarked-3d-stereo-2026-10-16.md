# Open Babel 3.2.1: unmarked 3D SDF stereochemistry

## Input and scoring

The 300 molecules are the two disjoint RDKit 2025.09.3 source-gated sets
described in [the COSMolKit comparison](../cosmolkit_comparison/unmarked-3d-stereo-2026-10-16.md).
They contain 150 tetrahedral and 150 alkene records. Each embedded molecule
was written without parity, wedges or CFG and retained only if RDKit reread
the actual V2000 **and** V3000 blocks with the source's canonical isomeric
SMILES. The original SDF SHA-256 values and generation commands are in that
document. Open Babel was `/opt/homebrew/bin/obabel`, version 3.2.1; CheMatic
was the source-built `sdf_v3000_accuracy_dump` example, binary SHA-256
`8fad4ebdb83f060cab07af11b9061e25ce9f9f090815253470c910f201602901`.

`check_unmarked_3d_stereo.py` runs Open Babel's `-ocan -xn` reader and the
CheMatic example on the **same** SDF. It parses both returned SMILES with
RDKit and compares to the gated source. The primary score removes atom-map
numbers but preserves isotopes, charges, graph and stereochemistry. The
strict score also includes atom maps; the graph score removes stereochemistry
as well. These scores are exact-identity checks, not perceptual similarity.

| Input | Engine | Strict | Mapless stereo | Mapless graph |
| --- | --- | ---: | ---: | ---: |
| Original 3D V2000 | CheMatic | 300/300 | 300/300 | 300/300 |
| Original 3D V2000 | Open Babel | 252/300 | 290/300 | 299/300 |
| Original 3D V3000 | CheMatic | 300/300 | 300/300 | 300/300 |
| Original 3D V3000 | Open Babel | 250/300 | 288/300 | 297/300 |

The strict gap includes Open Babel dropping atom-map numbers. Its mapless
stereo misses include lost ring-alkene or isotope-bearing stereo and an extra
bridgehead N tag. Graph misses include changes in radical or charge
placement; exact identity does not establish which resonance depiction is
chemically preferable.

## Reciprocal writer check

Open Babel also rewrote each SDF with `-osdf` (and `-x3` to request V3000).
Both engines reread those new files and were compared to the *original*
RDKit-gated source. This tests another writer's notation as well as the
reader. Open Babel's V2000 output contains some wedge/3D conflicts; CheMatic
now reports them and uses genuine 3D geometry. V3000 `RAD=` attributes now
constrain CheMatic's inferred H count on neutral C/N/O radicals.

| Open Babel rewrite | Engine | Strict | Mapless stereo | Mapless graph |
| --- | --- | ---: | ---: | ---: |
| V2000 | CheMatic | 260/300 | 298/300 | 300/300 |
| V2000 | Open Babel | 252/300 | 290/300 | 299/300 |
| V3000 | CheMatic | 260/300 | 298/300 | 298/300 |
| V3000 | Open Babel | 250/300 | 288/300 | 297/300 |

The V2000 CheMatic misses are also Open Babel misses (source rows 3022 and
1228). The V3000 CheMatic misses are also Open Babel misses (rows 1505 and
1392): Open Babel's rewrite changes terminal `[C]` to `[CH2]`. CheMatic
correctly retains the rewrite's `RAD=2` interpretation; the rewritten record
no longer has the source's radical state. Neither engine's reciprocal score
can exceed the information retained by that writer.

## Reproduction

First generate the two input sets using the commands in the linked COSMolKit
comparison. Then, from the repository root:

```sh
cargo build -p chematic-mol --example sdf_v3000_accuracy_dump --offline
python validation/openbabel_comparison/check_unmarked_3d_stereo.py \
  --case /tmp/chematic-3d-a-manifest.json /tmp/chematic-3d-a-v2000.sdf \
  --case /tmp/chematic-3d-b-manifest.json /tmp/chematic-3d-b-v2000.sdf \
  --case /tmp/chematic-3d-a-manifest.json /tmp/chematic-3d-a-v3000.sdf \
  --case /tmp/chematic-3d-b-manifest.json /tmp/chematic-3d-b-v3000.sdf \
  --chematic target/debug/examples/sdf_v3000_accuracy_dump \
  --output-prefix /tmp/chematic-openbabel-3d --rewrite
```

The script writes a JSON report with individual misses. This is a bounded
comparison on generated, RDKit-gated 3D SDF records. It does not establish
universal superiority over Open Babel, RDKit or COSMolKit across all formats,
chemistries, or tasks.
