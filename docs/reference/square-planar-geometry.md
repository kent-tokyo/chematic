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
particular input geometry is square planar or stereochemically valid. The
diagnostic does not apply a threshold, fit, or correct coordinates. Existing
checked conformer generation continues to
mark this stereo class unevaluable. In particular, a vacancy does not acquire
an inferred coordinate, and distortions, duplicate ligands, and unusual
coordination may require additional chemical interpretation.

## Experimental constrained embedding

`mol.rdkit_embed_square_planar()` runs one ETKDG embedding attempt with the
declared ligand-pair angles added to the distance bounds. It then measures
the final coordinates and raises `RuntimeError` if embedding fails or any
measured angle falls outside the chosen policy. No coordinates are returned on
such a failure. A caller may try another seed explicitly:

```python
mol = chematic.from_smiles("[Pt@SP3](F)(Cl)(Br)I").add_hydrogens()
coords = mol.rdkit_embed_square_planar(random_seed=1)
print(mol.square_planar_geometry_diagnostics(coords))
```

By default, each trans angle must be at least 150° and each cis angle must be
within 30° of 90°. `trans_min_degrees` and `cis_tolerance_degrees` change
these limits; they must be finite, nonoverlapping and within the angular
range. These defaults are experimental geometric thresholds, not universal
chemistry criteria. Unsupported neighbor orders, an unnumbered permutation,
missing declarations and invalid thresholds raise `ValueError`. The method
requires three or four real ligands per declared center, with at most one
recorded vacancy. Multiple centers are constrained together, so ring or
other distance constraints may make the bounds unsatisfiable.

Final checking is necessary: embedding and later energy refinement can move
angles outside their initial bounds. The method validates only the measured
square-planar angles of the returned embedding; it does not validate other
stereochemistry, bond lengths, energy, planarity beyond those angle limits,
or a later optimization. It does not change `rdkit_embed`, the timed ensemble
API, or the existing checked-conformer stereo acceptance policy. Rust exposes
the same operation in `square_planar_embedding`.

`rdkit_embed_square_planar_with_retry()` offers a separate bounded retry:

```python
result = mol.rdkit_embed_square_planar_with_retry(random_seed=1)
coords = result["coordinates"]
print(result["attempts"], result["used_tightened_bounds"])
```

It returns the first result unchanged if that result passes the original angle
policy. `initial_rejection` is `None` in this case, and otherwise records why
the first coordinates failed the angle policy. Only if a measured final angle
is outside that policy does it make at most one more
attempt, with the same seed and tighter embedding bounds: the trans minimum
is raised by 20° (capped at 180°) and the cis tolerance is reduced by 10°
(floored at 0°). Both attempts are checked against the original policy, not
the tighter bounds. Setup, invalid-input and embedding failures are not
retried. If tightening cannot change the bounds, no retry is made. When both
attempts fail, the error includes both reasons. A retry repeats embedding and
can take substantially longer; it is not a time limit or a guarantee of a
result. `attempts` counts these whole embedding calls; each call can itself
perform multiple internal initialization attempts. The numeric adjustments
are an experimental heuristic and may be unsuitable for particular complexes.

The Python regression tests compare atom order and bond connectivity first,
then compare every reported angle with RDKit's `GetAngleDeg` on six fixed Pt
examples spanning SP1–SP3 and two SMILES atom orders. They also check three
one-vacancy examples, including two previously exposed iodine and selenium
development cases. These tests use RDKit when installed; they do not prove
general square-planar assignment or conformer accuracy. A local run used
RDKit 2026.09.2 and a source-built macOS arm64 wheel; the wheel was not a
published release.
