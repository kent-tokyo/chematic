# Crystal atom meshes in OBJ

`write_crystal_obj_checked` exports a native `PeriodicStructure` as a finite
Wavefront OBJ atom display mesh. This addresses the CIF-to-OBJ workflow requested
in [OpenBabel #2680](https://github.com/openbabel/openbabel/issues/2680).
The mesh contains one separately named, closed triangular sphere for each site.

```python
from pathlib import Path
import chematic

source = chematic.PeriodicStructure.from_cif(
    Path("structure.cif").read_text(encoding="utf-8")
)
# Caller-selected display radii, one per site, in angstroms.
radii = [0.5] * source.site_count()
result = chematic.write_crystal_obj_checked(source, radii, subdivisions=2)
Path("structure.obj").write_text(result.text, encoding="ascii")
print(result.sphere_count, result.vertex_count, result.face_count)
```

For a directory of CIF files, use the same operation for each structure:

```python
for path in sorted(Path("structures").glob("*.cif")):
    source = chematic.PeriodicStructure.from_cif(path.read_text(encoding="utf-8"))
    mesh = chematic.write_crystal_obj_checked(source, [0.5] * source.site_count())
    path.with_suffix(".obj").write_text(mesh.text, encoding="ascii")
```

The function itself returns the complete result in memory and performs no file
I/O. It uses existing Cartesian positions, retaining site order without wrapping
or lattice reduction. The radius list may vary by site; radii are display sizes
chosen by the caller, with no automatic covalent or van der Waals interpretation.
Use a list, tuple or another real-valued Python sequence; convert a NumPy radius
array with `.tolist()`.

## Mesh contract

Coordinates and radii are in angstroms. OBJ has no standard unit declaration;
configure the consuming application so one coordinate unit has that scale. The
file also records the convention in an ASCII comment. It uses standard `v`,
`vn`, `f v//vn` and `o` records with positive 1-based indices and globally shared
radial unit normals. Object names are generated as `site_1_Na`, `site_2_Cl`, etc.;
source labels cannot inject commands or asset references.

`subdivisions` is an integer from 0 through 4 (default 2). The initial octahedron
has 6 vertices and 8 triangles. Each subdivision splits every triangle into four
and projects new edge vertices onto the unit sphere. Per site, subdivision `s`
therefore gives `4**(s+1)+2` vertices and `8*4**s` triangles. Faces wind
counterclockwise viewed from outside. These are polygonal approximations, with
vertices on the specified sphere; the triangular faces lie inside its smooth
surface.

The frozen `ObjMeshExport` contains `text`, `sphere_count`, `vertex_count`,
`face_count`, `units="angstrom"` and `omitted_features`. It reports omissions:
lattice vectors, periodic boundary conditions, bonds, site labels/properties, chemical semantics,
materials and colors. Element symbols appear in object names, but OBJ geometry
is not a chemical interchange or round-trip representation.

Spheres can be disconnected or overlap. They are not merged into a union, and
the output is not a guarantee of a printable solid. Unit-cell outlines, periodic
images across boundaries, bonds, isosurfaces, material files and colors are
outside this API. `make_supercell` can construct an explicit larger structure
before export, subject to the same mesh limits.

## Validation and limits

The API rejects empty structures, more than 1000 sites, partial/mixed occupancy,
nonpositive/nonfinite radii, and known unexpanded CIF symmetry. CIF parsing uses
the native explicit-operation expansion; no space-group database is consulted
and completeness of a supplied operation list is not independently established.

Bool subdivisions/radii and wrong argument types raise `TypeError`; unsupported
values raise `ValueError`. Generation is bounded by 100000 vertices, 200000
triangles and 16 MiB of ASCII output. Overflow or floating-point coordinate
resolution that distorts a sphere's vertex directions by more than `1e-8` is
rejected. No partial result is returned on failure, and input structures are
unchanged.

The implementation adds 51 regression cases; together with 19 existing crystal
and 162 compatibility cases, 232 tests pass. The independent
`validation/validate_crystal_obj.py` runner reads outputs using trimesh 4.9.0 and
checks closed surfaces, consistent outward winding, positive volumes, normals,
centers and radii. It covers 73 conditions, 532 spheres, all five subdivisions,
three cell bases, unwrapped fractional coordinates and all 118 element symbols.
Twelve cases also pass ASE-authored CIF through chematic's native reader and
compare with ASE 3.26.0 readback. Maximum center and radius discrepancies were
`1.18e-13 Å` and `2.89e-15 Å` on those cases. No OpenBabel OBJ parity is claimed.

Evidence and a CIF/OBJ/rendered example are under
`validation/results/crystal-obj-mesh-2026-10-17/`. The example preview adds colors
in its renderer; the OBJ itself contains geometry and no materials. Tests use a
recorded existing native extension: this change is Python-only, and the crystal
Rust crate and binding are unchanged from that extension's source commit. It is
not a fresh wheel validation. See the evidence README for hashes and commands.

The ASCII record conventions follow the
[Wavefront OBJ format description](https://www.loc.gov/preservation/digital/formats/fdd/fdd000507.shtml).
No OpenBabel implementation code or external mesh package is required at runtime.
