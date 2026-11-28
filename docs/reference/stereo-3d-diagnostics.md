# 3D declared-stereo diagnostics

`Mol.stereo_3d_diagnostics(coords)` reports the current verifier's result for
individual declared atom and double-bond stereochemical elements:

```python
mol = chematic.from_smiles("F/C=C/F")
coords = [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0]]
report = mol.stereo_3d_diagnostics(coords)
print(report["double_bond_stereo"])
```

`atom_stereo` contains dictionaries with `atom`, `status` and `reason`;
`double_bond_stereo` contains the same information with a `bond` index.
Statuses are `satisfied`, `violated` or `unevaluable`. An unevaluable element
has a diagnostic reason naming the verifier's reason; the other statuses have
`None`. Top-level `reported`, `satisfied`, `violated` and `unevaluable` are
counts over the reported elements. Coordinates must be finite and have one
XYZ entry for every atom in the molecule's current order, including any
explicit hydrogens. A wrong atom count or nonfinite coordinate raises `ValueError`.

These are declared-stereo checks, not a complete chemical validation. Atom
entries may include unsupported non-tetrahedral declarations such as
square-planar centers, which are reported as unevaluable. The verifier does
not cover every possible stereo class, and a zero `reported` count does not
establish that the molecule or its coordinates are correct. Bond lengths,
clashes, energies and stereochemistry not reported by this verifier need
separate assessment. The method is read-only; it neither repairs nor generates
coordinates. Square-planar ligand angles can be inspected separately with
[`square_planar_geometry_diagnostics`](square-planar-geometry.md).
On grossly distorted tetrahedral geometries, the signed-volume result can
disagree with RDKit's assignment from the same coordinates; a regression
fixture for development record 719 demonstrates this known limitation.

Set `require_regular_geometry=True` to apply a conservative tetrahedral
geometry screen before treating a reported atom as satisfied or violated:

```python
report = mol.stereo_3d_diagnostics(coords, require_regular_geometry=True)
```

Centers that fail the screen are reported as `unevaluable` with reason
`DistortedTetrahedralGeometry`; already unevaluable centers retain their prior
reason, and double-bond results are unchanged. The screen uses the same
ligand-weight and angle thresholds described for alternative repairs below.
It also rejects many centers where the original result agrees with RDKit, so it
is an opt-in uncertainty screen, not an independent stereo oracle or a proof
of validity. On the exposed development sample, the current screen made 1,100 of
2,055 otherwise comparable atom-center decisions unevaluable, including 910
that had agreed with RDKit. The fixed-sample protocol and follow-up are recorded
under `validation/results/2026-11-27-tetrahedral-geometry-screen/` and
`validation/results/2026-11-27-implicit-h-planarity/` in the repository.
The repair API can also require this screen on its output, as described below.

## Repairing declared stereo

`Mol.repair_stereo_3d(coords)` attempts to repair violated declared stereo by
moving bridge-connected substituents, and checks the result again. It returns
`ok`, `fully_verified`, `coordinates`, `repaired`, `unevaluable`, and
`failures`, `clash_check`, `regular_geometry_check`, and `retry`. Each repaired entry identifies an atom or bond and reports the
number of moved atoms and the largest displacement in Å. Unevaluable entries
identify an atom or bond. Failure entries also give a reason such as
`NoBridgeEligibleSubstituent`.

On failure, `coordinates` is `None`; partial coordinates are not returned.
`fully_verified` is true only when repair succeeds and no reported stereo was
unevaluable before repair; when `require_regular_geometry=True`, the output
must also pass that geometry screen. It is also true if the verifier reports no stereo,
so inspect `stereo_3d_diagnostics` to determine the covered elements. Unsupported
classes and unreported specifications remain outside this check. The operation
can move atoms substantially and does not minimize energy or screen for clashes;
check those separately. It does not change the input molecule or coordinates.

Set `require_regular_geometry=True` to reject an otherwise successful repair
when a reported tetrahedral center fails the geometry screen. The result then
has `ok=False`, no coordinates, and one `DistortedTetrahedralGeometry` failure
per rejected atom. `regular_geometry_check` gives `passed`, `rejected`, or
`not_evaluated_repair_failed`, along with `distorted_atoms` when evaluated.
It is `None` when the option is not selected. A `passed` geometry screen does
not mean unsupported or originally unevaluable stereo was verified; check
`fully_verified` and `unevaluable` separately.

Set `clash_cutoff` to a finite positive distance in Å to screen all nonbonded
atom pairs before and after repair. The screen includes 1–3 pairs and separate
fragments. If any pair newly falls strictly below the cutoff, the repair is
rejected with `NewNonbondedClashes` and no coordinates are returned. The
`clash_check` report gives the cutoff, status, counts, and the atom indices and
distances for new and resolved pairs. If stereo repair itself fails, its status
is `not_evaluated_repair_failed`. With no cutoff, `clash_check` is `None`.
Existing close pairs that become closer are not classified as new. A passing
check does not establish that the geometry is clash-free, chemically valid, or
free of unsupported stereochemistry. The screen has the same 2,048-atom and
coordinate-range limits as `check_geometry` and raises `ValueError` for an
invalid cutoff or input outside those limits.

With `max_alternative_attempts` set to a value from 1 through 16, a repair
rejected for new close pairs, or by the requested regular-geometry check, can try deterministic alternative choices of
bridge-connected substituents. `clash_cutoff` is required. The original repair
is attempted first, and an already accepted result is unchanged. The `retry`
report gives the number attempted, the selected variant (zero for the original,
`None` if none was accepted), and counts for stereo failures, geometry errors,
new-clash rejections, and tetrahedral geometry rejections. A repair that fails
before the clash screen is not retried. Each variant applies the same choice
index across stereo elements; this is a bounded search, not an exhaustive
enumeration of combinations.

Alternative candidates also undergo a conservative screen of every reported
satisfied tetrahedral center. Each of the center's barycentric weights in the
tetrahedron of its ligand positions must exceed `1e-6`; an implicit ligand is
reconstructed where needed. The ligand directions must form pair angles from
45° to 165°. This screen is applied to alternatives and, when explicitly
requested with `require_regular_geometry=True`, to the original repair. It rejects grossly
distorted candidates but is not an independent stereochemistry oracle and does
not guarantee valid geometry or remove pre-existing clashes. Candidates that
leave reported stereo unevaluable still return `fully_verified=False`.

For carbon with three graph neighbors and one implicit hydrogen, the screen
additionally requires the absolute scalar triple product of the three
normalized real bond directions to be at least 0.15. The same check applies to
carbon with three non-hydrogen neighbors and one explicit hydrogen when its
direction matches the reconstructed implicit-H direction within a cosine
tolerance of `1e-8`. Independently placed explicit hydrogen is not subjected
to this additional check. Zero volume indicates coplanar directions; the 0.15
margin is an experimental conservative threshold. It also applies to
`require_regular_geometry=True` and does not change the default stereo verifier.

```python
mol = chematic.from_smiles("C[N@H+](CC)CCC")
coords = mol.generate_3d_etkdg()
result = mol.repair_stereo_3d(coords, clash_cutoff=0.8, require_regular_geometry=True)
if result["ok"]:
    coords = result["coordinates"]
    print(result["repaired"], result["unevaluable"])
else:
    print(result["failures"])
```
