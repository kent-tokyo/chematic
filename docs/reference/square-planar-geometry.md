# Square-planar angle diagnostics

`Mol.square_planar_geometry_diagnostics(coords)` measures ligand-center-ligand
angles for explicitly numbered `@SP1`, `@SP2`, and `@SP3` declarations. It
uses the SMILES-recorded neighbor order and returns the pairs declared trans
and cis by that permutation, with atom indices and angles in degrees.

```python
import chematic

mol = chematic.from_smiles("[Pt@SP1](F)(Cl)(Br)I")
coords = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]]
report = mol.square_planar_geometry_diagnostics(coords)[0]
print(report["trans"])  # [(1, 3, 180.0), (2, 4, 180.0)]
print(report["cis"])    # the other four pairs, at 90 degrees
```

`center` is the atom index, `permutation` is `SP1`, `SP2` or `SP3`, and
`neighbors` preserves the four recorded positions. One `None` neighbor means
an explicitly recorded vacant site: for a three-ligand center the report
contains only the three measurable real-ligand pairs. The two lists contain
`(atom_index, atom_index, angle_degrees)` triples, not a classification of
the measured geometry. A numbered four-ligand center has two trans and four
cis pairs; a supported one-vacancy center has one trans and two cis pairs.
The recorded stereo order, including vacancies, is required: the API never
infers it from graph iteration order.

An `unavailable` reason reports an unnumbered permutation, missing or
unsupported neighbor order, or a zero/invalid ligand vector; its angle lists
are empty. Coordinate count mismatch or nonfinite coordinates raise
`ValueError`. Molecules without square-planar declarations return an empty
list. Rust exposes the same measurements through
`chematic_3d::square_planar_geometry::square_planar_angles`.

The names cis and trans describe the *declaration*, not an assertion that a
particular input geometry is square planar or stereochemically valid. There is
no threshold, fitting, coordinate correction, embedding constraint, or
acceptance decision here. Existing checked conformer generation continues to
mark this stereo class unevaluable. In particular, a vacancy does not acquire
an inferred coordinate, and distortions, duplicate ligands, and unusual
coordination may require additional chemical interpretation.

The Python regression tests compare atom order and bond connectivity first,
then compare every reported angle with RDKit's `GetAngleDeg` on six fixed Pt
examples spanning SP1–SP3 and two SMILES atom orders. They also check three
one-vacancy examples, including two previously exposed iodine and selenium
development cases. These tests use RDKit when installed; they do not prove
general square-planar assignment or conformer accuracy. A local run used
RDKit 2026.09.2 and a source-built macOS arm64 wheel; the wheel was not a
published release.
