# Fixed fragments and bounded numerical contracts

This follows #812, whose processed main Codecov result is 95.43%. Coverage scope and exclusions are unchanged. Production algorithms and public APIs are unchanged.

- Partial drawing fragments: fix a carbon centre and zero through three of its F/Cl/Br/I ligands, rotate the fixed coordinates, then expand the remaining fragment. Check that fixed positions remain exact, all atoms are placed, attachment points are exhausted, and all four bonds have length 1.5 Å. Check the occupied/remaining angle at the two-/three-ligand boundary. Coincident fixed vectors report the normalization failure. These are private fragment-helper contracts; the public depiction API does not accept arbitrary coordinate maps.
- Distance bounds: append dimethyl/diethyl disulfide, hydrogen disulfide, dimethyl trisulfide, a cyclic disulfide, and cystine. All six new molecules have explicit hydrogens, as required by `rdkit_bounds_matrix`. Compare all eight set15/smoothing/macrocycle combinations with RDKit 2026.03.1. The original 728 reference rows remain intact; there are now 776. Hydrogen-free `SS` is outside this API contract: its UFF valence flags produce a different bound, so these tests do not claim parity for that input.
- MANCUDE enumeration: benzene has two full Kekulé matchings at the exact atom/matching limit. Starved atom, matching-count, and search-step budgets return their distinct typed errors with the exhausted limit in the message, retaining the aromatic graph.
- Workflows: exact molecule/atom limits accept small valid comparisons. Rejected comparisons identify the original input index and SMILES. Screening retains every input record and its index, with separate atom-limit and list-limit errors; selection/clustering only reference accepted records.
- Input diagnostics: exercise actual empty/unsupported embedding inputs, malformed SMILES labels and positions, POSCAR lines/values, and all three square-planar geometry rejection reasons. Collinear MOL layouts report the actual unexpressed atom/bond identity, while ordinary unstereogenic export reports no loss.
- Cartesian distance: the independent 3–4–12 triangle has distance 13, under translation, axis rotation, and scales from 1e-6 to 1e6. Check symmetry and zero self-distance.
- Alignment: recover a known rigid or reflected tetrahedron with uniform and nonuniform positive weights. Independently sum weighted residuals and check the transform determinant (1 for rotation, -1 for reflection). Reject point/weight count mismatches and nonpositive weights. RDKit's Jacobi stopping tolerance gives maximum coordinate errors about 2.85e-7 (uniform weights) and 1.21e-6 (nonuniform weights) on these exact cases, so coordinate assertions allow 2e-6 rather than assuming an exact fit; the SSR checks use 1e-9.
- MMFF bond fallback: six real interhalogens (ClBr, ClI, BrI, FCl, FBr, FI) have no direct bond-parameter row. Assigned atom types and Herschbach–Laurie empirical `kb/r0` match RDKit 2026.03.1, including reversed endpoint order and the reported empirical resolution mechanism.

Regenerate the bounds reference with `scripts/generate_bounds_boundary_fixtures.py` using exactly RDKit 2026.03.1. Its SHA-256 is `c65bbda00872d54b5af0b58a8804ce141ef938b6175ca66b8eebb24dad4e1471`.

The inline interhalogen values are reproduced with `AllChem.MMFFGetMoleculeProperties(Chem.AddHs(Chem.MolFromSmiles(smiles)))`, then `GetMMFFAtomType(i)` and `GetMMFFBondStretchParams(mol, 0, 1)`.

The alignment reference points are `(0,0,0)`, `(2,0,0)`, `(0,3,0)`, `(0,0,4)`. Probe points are `(5-s*y, -2+s*x, 7+s*z)`, with `s=1` for rotation and `s=-1` for reflection. Use `rdAlignment.GetAlignmentTransform(reference, probe, weights=[1,1,1,1]` or `[1,2,3,4], reflect=..., maxIterations=100)` from `rdkit.Numerics` to reproduce the numerical check.

Coverage is execution evidence for these contracts, not a general guarantee of chemistry correctness or RDKit equivalence.

Validation: 197 boundary tests passed with zero failures or ignored tests. Workspace all-target Clippy, formatting, and diff checks passed. The complete clean coverage command and processed Linux service result are recorded on the accompanying PR.
