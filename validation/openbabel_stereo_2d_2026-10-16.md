# Open Babel 2D MOL stereo comparison, 2026-10-16

This is a paired MOL-to-SMILES comparison for molecular recognition and
stereochemistry. RDKit 2025.09.3 writes each 2D MOL input and reads it back.
CheMatic 1.1.4 and Open Babel 3.2.1 receive identical bytes. RDKit then reads
each output SMILES and compares its canonical isomeric SMILES with the source.
Writer failures, changed source identities, and a writer silently changing the
requested MOL version are excluded before either converter is scored.

The tested CheMatic binary was built from main `2c3ff6ff` plus the sulfur
lone-pair SMILES writer fix in this change. Its SHA-256 was
`ae81a945017fe7d273f72f18c39e6490de1490d13a815e370c7e8275005350f6`.
The Open Babel binary SHA-256 was
`733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.

| Corpus and format | Comparable | CheMatic strict | Open Babel strict | Open Babel only |
| --- | ---: | ---: | ---: | ---: |
| CIP corpus, V2000 | 155 | 155 | 143 | 0 |
| CIP corpus, V3000 | 155 | 155 | 146 | 0 |
| First 300 unique stereo stress molecules, V2000 | 291 | 290 | 247 | 0 |
| First 300 unique stereo stress molecules, V3000 | 297 | 296 | 236 | 0 |

The CIP corpus SHA-256 is
`14ce78b96a24bc66b65fd963ffc7e236e1153f72c006933555335d0bdd61730e`.
It has no source exclusions. The stress corpus SHA-256 is
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`.
The stress selection excludes `@SP`, `@TB`, and `@OH` sources, then takes the
first 300 distinct RDKit identities containing `@`, `/`, or `\` syntax. Of
600 format attempts, 12 source inputs were excluded: six were written as
V3000 despite a V2000 request, four changed identity on RDKit re-read, and
two could not be written by RDKit. The same source molecule accounts for the
one case where neither converter matched in each format.

## Nonstereogenic explicit-H adjudication on latest main

The strict table above treats stress source line 305 (`[H]/N=C(/N)N`) as a
miss for both converters. RDKit 2025.09.3 finds no specified stereo in this
source: the imine carbon has two identical N substituents. CheMatic writes
`C(N)(N)=N[H]`; Open Babel writes `N=C(N)N`. RDKit removes the ordinary
explicit H from each representation and obtains the same molecular identity.
The paired gate now accepts this H representation difference only when both
parsed molecules have no specified stereo and no labeled or isotopic explicit
H. A stereogenic H-defined imine, E/Z reversal, and isotope loss remain
negative controls. The strict result remains available per row.

| Corpus and format | Comparable | CheMatic verified equivalent | Open Babel verified equivalent | Open Babel only |
| --- | ---: | ---: | ---: | ---: |
| CIP corpus, V2000 | 155 | 155 | 143 | 0 |
| CIP corpus, V3000 | 155 | 155 | 146 | 0 |
| First 300 unique stereo stress molecules, V2000 | 291 | 291 | 248 | 0 |
| First 300 unique stereo stress molecules, V3000 | 297 | 297 | 237 | 0 |

Only line 305 changed classification for both engines and both MOL versions.
The current-main source revision was `1a7188cd6d6ea6b566892cd841fe3481f08f94b6`;
the source-built CheMatic CLI SHA-256 was
`326274218779f255c284c72f943de3c9584847ff6f343d0301754b9cdb9eea6e`.
Open Babel remained 3.2.1 with binary SHA-256
`733b2c925d60b17a5c9d8ce9ec841e918cc20c7634301e23f9008a187d1d9bc5`.
This correction removes a false negative in the *comparison*. It does not
change CheMatic's molecular parser or its SMILES writer.

The fixed defect was a three-coordinate sulfoxide. MOL wedge perception stores
its stereogenic lone pair as a virtual terminal in the neighbor order, but
both SMILES writers omitted that terminal when sulfur had zero implicit
hydrogens. The omitted terminal inverted the emitted `@`/`@@` assignment.
RDKit-written V2000 and V3000 fixtures now test both wedge directions and
both SMILES writers.

These results establish an advantage for this specific 2D MOL-to-SMILES
workflow on these corpora. The RDKit source/re-read identity is the comparison
oracle; it is not independent experimental stereochemical truth. They do not
establish toolkit-wide accuracy superiority, including other formats, 3D,
force fields, or all chemical classes.

## Reproduce

```sh
CARGO_INCREMENTAL=0 CARGO_PROFILE_DEV_DEBUG=0 cargo build -p chematic-cli --offline -j 1
python3 scripts/compare_openbabel_stereo_formats.py \
  --output /private/tmp/chematic-openbabel-cip-2d.json
python3 scripts/compare_openbabel_stereo_formats.py \
  --corpus validation/cosmolkit_comparison/corpora/stress-v1.smi \
  --stress-first-unique 300 \
  --output /private/tmp/chematic-openbabel-stress-2d.json
```
