# Open Babel 3D MOL2 tetrahedral recognition, 2026-10-16

## Paired input and score

RDKit 2025.09.3 embedded twelve specified tetrahedral structures and one
inversion of each, giving 24 inputs. Open Babel 3.2.1 wrote each 3D MOL2
from an RDKit V2000 MOL block with explicit H atoms. RDKit then reread the **actual MOL2 bytes**;
all 24 retained the source's canonical isomeric identity and entered the
paired comparison. The same MOL2 file was passed to both converters, and
RDKit compared their output SMILES after removing ordinary explicit hydrogen
nodes. The report records every input hash, output and tool hash.

| Converter | Source identity preserved |
| --- | ---: |
| CheMatic main `b9f1e3cd`, before this change | 0/24 |
| CheMatic with 3D MOL2 tetrahedral perception | 24/24 |
| Open Babel 3.2.1 | 24/24 |

The original CheMatic reader retained MOL2 coordinates but did not use them
to recognize tetrahedral configuration. Both MOL2 graph readers now reuse the
MOL/SDF noncoplanar 3D perception routine. Flat or degenerate coordinates do
not create a stereotag. The native regression fixtures include a center with
an explicit H, an implicit H, an amino acid, and a four-heavy-atom center; coordinate
reflection flips the result. Nonfinite coordinates are rejected.

This closes the measured 3D MOL2 tetrahedral gap against Open Babel on this
generated panel. It does not measure alkene, axial, coordination or document
stereochemistry in MOL2. RDKit's source/re-read identity is a pinned
interoperability oracle, not independent experimental stereochemical truth.
No toolkit-wide superiority follows from this panel.

## Reproduce

Build the CLI from the candidate source, then run:

```sh
python3 scripts/compare_openbabel_mol2_3d_stereo.py \
  --chematic target/debug/chematic \
  --obabel /opt/homebrew/bin/obabel \
  --output /tmp/openbabel-mol2-3d-stereo.json
```

To reproduce the before/after column, add `--before` with a CLI binary built
from `b9f1e3cd`. The complete result for this run is in
[`openbabel-mol2-3d-stereo-2026-10-16.json`](openbabel-mol2-3d-stereo-2026-10-16.json).
