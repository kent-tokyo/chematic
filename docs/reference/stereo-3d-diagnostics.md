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

## Repairing declared stereo

`Mol.repair_stereo_3d(coords)` attempts to repair violated declared stereo by
moving bridge-connected substituents, and checks the result again. It returns
`ok`, `fully_verified`, `coordinates`, `repaired`, `unevaluable`, and
`failures`. Each repaired entry identifies an atom or bond and reports the
number of moved atoms and the largest displacement in Å. Unevaluable entries
identify an atom or bond. Failure entries also give a reason such as
`NoBridgeEligibleSubstituent`.

On failure, `coordinates` is `None`; partial coordinates are not returned.
`fully_verified` is true only when repair succeeds and no reported stereo was
unevaluable before repair. It is also true if the verifier reports no stereo,
so inspect `stereo_3d_diagnostics` to determine the covered elements. Unsupported
classes and unreported specifications remain outside this check. The operation
can move atoms substantially and does not minimize energy or screen for clashes;
check those separately. It does not change the input molecule or coordinates.

```python
mol = chematic.from_smiles("C[N@H+](CC)CCC")
coords = mol.generate_3d_etkdg()
result = mol.repair_stereo_3d(coords)
if result["ok"]:
    coords = result["coordinates"]
    print(result["repaired"], result["unevaluable"])
else:
    print(result["failures"])
```
