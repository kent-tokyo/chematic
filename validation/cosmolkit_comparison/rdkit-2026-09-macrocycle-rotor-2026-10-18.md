# Large O/S rings and explicit-H rotors versus COSMolKit 0.5.0

Development comparison on macOS arm64, CPython 3.13, with RDKit 2026.09.1
as the operation-level reference. COSMolKit 0.5.0 was the installed wheel;
CheMatic 1.1.7 was built from source, starting at `f8e9d7dee5a8` and then
applying the changes described here. These are fixed-corpus results, not a
general accuracy ranking or a published CheMatic package result.

## Deficits and changes

The 256-row selection takes row `round(i * (4072 - 1) / 255)` for each
`i = 0..255` from `corpora/stress-v1.smi` (SHA-256
`879166ab5481cac805ecd2d77f37d628383ec124085b90e2ee1659d55ce76277`).
The selection SHA-256 is
`b444418f3809ae5fbc4f1a72b1142d7f9b42416ecfb1a5e316377b254d5a5bfa`.
`run_corpus.py` and `compare_corpus.py` used identical input and operation
definitions for all three engines. Expensive 3D, reaction, and tautomer-set
operations were skipped; the listed operations were not skipped.

| Operation, strict RDKit agreement | CheMatic before | CheMatic after | COSMolKit |
| --- | ---: | ---: | ---: |
| Canonical SMILES | 254/256 | 256/256 | 256/256 |
| Stereoisomer count | 254/256 | 256/256 | 256/256 |
| MolBlock read-back | 250/256 | 252/256 | 252/256 |
| Rotatable bonds | 255/256 | 256/256 | 256/256 |

Two inputs, `s1cccccccccccc1` and `c1ccccccsccccccccc1`, caused the
SMILES and stereoisomer differences. CheMatic's RDKit-compatible SMILES
aromaticity port lacked the newer RDKit rule that makes a neutral, divalent
O/S two-electron donor an aliphatic bridge when its first SSSR ring has at
least nine atoms and no incident double/triple bond. The separate
`chematic-perception` parity path already had this rule. The SMILES port now
applies it at the same first-ring typing point. Five-membered O/S rings remain
aromatic. See [RDKit's Aromaticity.cpp](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Aromaticity.cpp)
for the upstream rule.

The rotor difference was `[H]OC(=O)C([H])([H])N([H])[H]`: CheMatic counted
two bonds, versus one in RDKit and COSMolKit. Ordinary removable explicit H
atoms made a terminal OH/NH endpoint pass the strict SMARTS degree test.
The result now rejects a matched bond when either endpoint has fewer than
two neighbours after those H atoms are omitted. Isotope-labelled H remains
an explicit neighbour, matching the existing hydrogen-removal contract.

## Larger focused check

On all 4,072 `stress-v1.smi` rows, the source-built candidate matched
RDKit 2026.09.1 for **canonical SMILES, rotatable-bond count, and
stereoisomer count: 4,072/4,072 each**. COSMolKit matched 4,066/4,072 for
each because it rejected six RDKit-accepted inputs containing aromatic
`[siH]` or `[si]`; on the 4,066 inputs it parsed, its values matched.
This larger check covers only those three operations.

The 256-row comparison still has COSMolKit ahead on strict canonical-SMILES
read-back (256 versus 255), fragment SMILES (256 versus 253), and the three
`[R1]`/`[R2]`/`[R3]` SMARTS queries (256 versus 255 each). The ring-query
differences occur on one fused-ring input. These residuals need separate
investigation. Native CheMatic SMILES spelling is also a separate API from
the RDKit-compatible `rdkit_smiles` used here; an earlier full-corpus native
round-trip check found 4,049/4,072 strict RDKit strings for CheMatic versus
4,064/4,072 for COSMolKit, while all 4,072 CheMatic structures were
classified equivalent by that harness. This change does not relabel native
serialization as strict RDKit output.

Regression checks: `cargo test -p chematic-smiles --lib --offline` (408
passed, one ignored), `cargo test -p chematic-chem --lib --offline` (970
passed, one ignored), and the repository refactor-budget check. The Python
comparison used a source-built extension with `native-inchi` enabled, as in
the pre-change 256-row run; InChI is outside the three-operation full-corpus
check.
