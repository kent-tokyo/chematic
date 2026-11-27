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
