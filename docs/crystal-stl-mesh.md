# Crystal atom spheres in ASCII STL

`write_crystal_stl_checked` exports ordered native periodic sites as finite atom
spheres in ASCII STL. This covers a bounded CIF-to-STL workflow requested in
[Open Babel #2679](https://github.com/openbabel/openbabel/issues/2679).
The Python API returns geometry in memory and writes no files itself.

```python
import chematic as ch

source = ch.PeriodicStructure(
    ch.Lattice.cubic(5.64),
    [ch.Site([("Na", 1.0)], (0.0, 0.0, 0.0)),
     ch.Site([("Cl", 1.0)], (0.5, 0.5, 0.5))],
)
result = ch.write_crystal_stl_checked(
    source, [0.4, 0.7], millimeters_per_angstrom=3.0, subdivisions=2
)
assert result.sphere_count == 2
assert result.triangle_count == 256
assert result.units == "millimeter"
assert result.text.startswith("solid chematic_atom_spheres_mm\n")
```

For a CIF file, construct the source with
`ch.PeriodicStructure.from_cif(path.read_text(encoding="utf-8"))`. Save the result
with `output_path.write_text(result.text, encoding="ascii")` if desired.
Radii are caller-chosen display sizes in angstroms, one per site. The required
`millimeters_per_angstrom` applies equally to Cartesian site positions and radii.
At scale 3.0, an input separation of 2 Å becomes 6 mm and radius 0.4 Å becomes
1.2 mm. This is a chosen visualization scale, not the physical size of an atom.

STL has no standard unit field. Configure the consuming application to interpret
coordinates as millimeters; the `solid` name and report record that convention
but cannot enforce it in other software. Coordinates are not wrapped, shifted
to the origin, or made positive. No lattice reduction or symmetry expansion is
performed during export.

## Geometry and preserved scope

The encoder reuses the checked sites and octahedral sphere tessellation from
[the OBJ writer](crystal-obj-mesh.md). `subdivisions` is an exact integer 0..4
(default 2); each sphere has `8 * 4**subdivisions` outward-wound triangles.
Every facet carries a finite unit normal computed from its serialized triangle
coordinates. Vertices use 17 significant digits. Faces approximate the curved
sphere surface; their vertices lie at the supplied radius.

The frozen `StlMeshExport` contains `text`, `sphere_count`, `triangle_count`,
`millimeters_per_angstrom`, `units="millimeter"`, `format="ascii_stl"` and
`omitted_features`. It omits lattice/PBC information, bonds, atom/site identity,
occupancy/chemical semantics, component names, labels, colors and materials.
STL is a triangle geometry format, not a chemical round-trip representation.

**The spheres can overlap or be disconnected.** No surface union, bond geometry,
overlap repair or printable-solid validation is performed. Exporting a valid STL
does not establish that it can be manufactured as a single solid. PDB input,
binary STL, Rust/WASM/CLI entry points, lattice outlines, periodic images and
automatic chemical radii are outside this addition. Call `make_supercell` first
if an explicit repeated structure is needed, within the same site/mesh limits.

## Checks and validation

The writer refuses empty structures, more than 1000 sites, partial/mixed
occupancy, known unexpanded CIF symmetry, invalid radii or scale, non-integer
subdivisions, overflowing/underflowing scaled geometry and numerically collapsed
triangles. The shared vertex budget is 100000, with additional STL limits of
50000 triangles and 16 MiB ASCII output. Lower resolution may be required even
below the site cap. Invalid types raise `TypeError`; unsupported values/limits
raise `ValueError`. A refused export returns no partial result.

Independent validation uses trimesh 4.9.0 to read ASCII STL, checking triangle
counts, winding, watertight components, positive volume, Euler number, radii,
centers, scale, bounds and facet normals. ASE 3.26.0 independently checks the CIF
coordinates modulo lattice-image translations. The output retains native site
positions; the selected ASE image is recorded separately. Original OBJ bytes
are compared with hashes captured before shared preparation was extracted.

Reproduction and sample CIF/STL files are under
`validation/results/checked-crystal-stl-2026-10-18/` in the repository. No slicer,
3D printer or Open Babel upstream build is claimed to have been validated.
