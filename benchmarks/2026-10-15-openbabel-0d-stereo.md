# Open Babel 0D MOL stereo preservation (2026-10-15)

The V2000 reader now interprets atom parity and Open Babel's 0D alkene bond
markers when the coordinate block is entirely zero. Open Babel 3.2.1 labels its
0D output `2D` in the header, so the decision uses the coordinates. Existing
2D wedge interpretation takes precedence, and nonzero-coordinate input does
not use the new 0D path. Identical terminal substituents suppress a parity
tag that cannot define a stereocenter, matching Open Babel's reread of the
adversarial two-fluorine case.

## Paired source check

On one macOS arm64 host, Open Babel 3.2.1 generated 0D V2000 MOL files from
the 122 unique SMILES in the [CIP corpus](../validation/cip_label_corpus.jsonl).
RDKit 2025.09.3 canonicalized each source SMILES. A case was scored only if
Open Babel reread its own MOL file to the same RDKit isomeric identity. The
pre-fix CheMatic CLI and the modified CLI then read identical MOL bytes; RDKit
canonicalized their output SMILES before identity comparison.

| Stage | Source identity preserved |
| --- | ---: |
| Open Babel self reread, source gate | 122 / 122 |
| CheMatic before the fix | 0 / 122 |
| CheMatic after the fix | 122 / 122 |

All 122 files contain at least one nonzero atom parity field. No case was
excluded by the source gate. The [JSON record](2026-10-15-openbabel-0d-stereo.json)
contains every source identity, MOL input hash, and both CheMatic outcomes.
The before binary has the same product code as base revision `cac480fb`;
the intermediate benchmark-only commits used to build it changed scripts,
not the MOL reader. The after binary was built from this branch's working tree.

Reproduce with Open Babel 3.2.1, RDKit 2025.09.3 and CheMatic CLI binaries
built before and after this change:

```bash
python3 scripts/compare_openbabel_0d_stereo.py \
  --before /path/to/chematic-before \
  --after /path/to/chematic-after \
  --report /tmp/openbabel-0d-stereo.json
```

This establishes source-identity preservation for the stated 0D MOL input
lane. It does not adjudicate absolute stereochemistry independently or show
toolkit-wide superiority over Open Babel.
