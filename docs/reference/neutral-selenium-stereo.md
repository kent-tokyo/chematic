# Neutral selenium stereochemistry in 3D

The declared-stereo verifier supports a limited three-coordinate neutral
selenium class with no implicit hydrogen: exactly one double bond to nitrogen,
oxygen, sulfur or selenium, plus two ordinary single bonds. For these inputs,
it evaluates the declared `@` or `@@` parity from the three real ligand
positions and a virtual lone-pair direction. The same check is used by the
opt-in square-planar checked embedding and other existing callers of
`verify_stereo`.

This is a narrowly scoped parity check. A degenerate arrangement is reported
as unevaluable. Other selenium charges, numbers of neighbors, double-bond
partners, explicit or implicit hydrogens, and non-tetrahedral stereo retain
their existing behavior and may be unevaluable. A reported stereo match does
not certify bond lengths, electronic structure, energy, or whether a molecule
is stable. The existing square-planar angle and planarity checks remain
separate.

The checked-in Rust oracle fixture contains 24 RDKit 2026.09.2 ETKDG
coordinate sets covering both tags, multiple SMILES atom orders, seeds 7 and
42, and double bonds to N, O, S and Se. The verifier is tested on the original
coordinates and their reflected images. Additional authored structures and
seeds are compared against RDKit's independent 3D stereo reassignment in the
local development probe. These examples are small and do not establish general
selenium stereochemistry or a population-wide accuracy rate. RDKit printed
UFF typing warnings for some generated selenium examples; these coordinates
are used to test parity convention, not as validated geometry or energy.
