# Chemistry and conversion issue audit (2026-10-25)

This is a small regression packet, not a claim of general superiority. The
expected structures come from the authored SMILES and stereo relationships;
another toolkit's output is never treated as the chemical oracle.

| Upstream report | Local reproduction | CheMatic check |
| --- | --- | --- |
| [RDKit #9738](https://github.com/rdkit/rdkit/issues/9738), multiple conformers with seed 0 | RDKit 2026.09.1: all four conformers after the first have RMS 0.0; seed 1 gives nonzero RMS values | Existing `crates/chematic-3d/tests/zero_seed_ensemble.rs`: five tests pass; seed 0 gives a reproducible, nonidentical ensemble |
| [RDKit #9732](https://github.com/rdkit/rdkit/issues/9732), reaction CX stereo-group index after reactant reordering | RDKit 2026.09.1 writes `CC[C@H](C)O.O>>CCC(C)=O \|&1:3\|`; the stereocentre is output atom 2 | Existing `crates/chematic-rxn/tests/reaction_cxsmiles.rs`: four tests pass, including group membership after canonical reordering |
| [Open Babel #3035](https://github.com/openbabel/openbabel/issues/3035), 0D MOL mislabelled 2D | Open Babel 3.2.1: direct InChI has `/b`, `/t`, `/m`, `/s`; SMILES → SDF → InChI loses those layers | New `crates/chematic-mol/tests/openbabel_3035_3036.rs`: import of the exact Open Babel 0D MOL output recovers both stereo features; CheMatic's laid-out MOL output survives a read-back and has nonzero depiction coordinates. CheMatic CLI MOL → Open Babel InChI also retained all four stereo layers locally. |
| [Open Babel #3036](https://github.com/openbabel/openbabel/issues/3036), variable 0D MDL stereo flags | Open Babel 3.2.1: 40 separate CLI conversions yielded the correct imine stereo 34 times and the opposite stereo 6 times | New MOL test: 40 fresh CheMatic graphs write identical bytes and retain the source stereo on read-back; the test distinguishes the wrong stereoisomer |
| [CDK #1353](https://github.com/cdk/cdk/issues/1353), sulfur-oxide aromaticity | CDK maintainer describes the result as model-dependent. RDKit 2026.09.1 marks the neutral `S=O` form nonaromatic and the charged `S+–O-` form aromatic | New `crates/chematic-perception/tests/cdk_1353_oxide.rs` pins both CheMatic policies: native counts 0/0, RDKit-parity counts 0/5, with thiophene at 5/5 as a control |
| [OpenChemLib #180](https://github.com/Actelion/openchemlib/issues/180), alternating isomeric SMILES with stable IDCode | Upstream report only; OpenChemLib was not run locally | New `crates/chematic-smiles/tests/openchemlib_180_canonical.rs`: all three reported spellings converge to one CheMatic canonical SMILES and a parse/write fixed point |

Run the focused checks with `cargo test -p chematic-mol --test openbabel_3035_3036`,
`cargo test -p chematic-perception --test cdk_1353_oxide`, and
`cargo test -p chematic-smiles --test openchemlib_180_canonical`.

The CDK case is a model-policy distinction. These fixtures do not establish a
model-independent aromaticity truth or a corpus-wide conversion rate. The
ChEMBL percentages in upstream reports were not independently reproduced here.
