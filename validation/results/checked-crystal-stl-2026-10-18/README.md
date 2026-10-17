# Explicit-scale ASCII STL crystal atom meshes

Base `0e7a6d9f8d50fd75a9b1bfada00336f8e320f8ec` after integrating main; initial validation base
`abcc3569e7f49d15e911945509a5083ffff8a63d`. Python, typing, stub and docs
checks were repeated after integration.
[Open Babel #2679](https://github.com/openbabel/openbabel/issues/2679) remains OPEN
and asks for crystal-to-STL output. The base has finite atom-sphere OBJ export but
no STL API. This adds a bounded Python `write_crystal_stl_checked` and frozen
`StlMeshExport`, with a required mm-per-angstrom scale and caller display radii.
It extracts existing OBJ preparation for reuse; OBJ serialization stays identical.

## Validation

- **266 Python tests passed**: 34 new STL tests, 51 existing OBJ tests, 19 native
  crystal binding tests and 162 compatibility tests. Covers ASCII facet grammar,
  uniform scaling, serialized normals/winding, source immutability, frozen metadata,
  type/value errors, symmetry/occupancy and numeric/resource boundaries.
- **75 independent mesh conditions** read by **trimesh 4.9.0**: 210 spheres,
  114,576 triangles and 57,708 merged vertices. Checks positive volume and area,
  unit/serialized facet normals, closed surfaces, Euler number, component count,
  vertex radii, center, bounds and triangle/vertex counts at subdivision levels
  0..4 and scales 0.25, 3 and 12.5 mm/angstrom. Corpus spheres do not overlap.
- **ASE 3.26.0** checks two authored CIF fixtures modulo integer cell translations.
  For C2/c it selects one image shifted by one lattice-a vector relative to the
  native result; P1 has zero image shifts. This is recorded in `oracle.json`.
  STL retains native Cartesian sites without wrapping or image changes.
- **20 pre-refactor OBJ hashes match byte-for-byte**, across four coordinate/site
  controls and all five subdivision levels. `obj-baseline.json` was captured from
  the prior fixture before helper extraction; the oracle reconstructs every case.
- Strict mypy accepts the public consumer and both mesh modules; six negative
  consumer calls/assignments produce exactly six diagnostics.
- Stub, binding, cross-binding, unsafe, refactor and product gates; strict MkDocs;
  one executed guide block and two public API anchors.

The original sample CIF/STL are `nacl-p1.cif` / `nacl-atom-spheres.stl`, using radii
0.2/0.24 Å, scale 3 mm/Å and subdivision 2. These are display spheres, not inferred
chemical radii or a bonded NaCl model.

## Open Babel control

The installed **Open Babel 3.2.1** registers a write-only STL format; its help lists
probe/radius scaling, colors and bonds. A bounded direct P1 CIF-to-STL control
timed out after 30 seconds with no captured output (`openbabel-control.json`).
This is not proof of absent STL support or an upstream format defect. An earlier
unbounded diagnostic produced binary data, which a text-only capture rejected;
that output was not retained. No output-geometry parity or timing comparison is
claimed. The installed plugin supports a broader/different rendering profile;
this addition fills the missing chematic API motivated by the open issue.

## Scope and limits

- Required finite positive `millimeters_per_angstrom`; equal scaling of centers
  and radii, no implicit unit conversion. STL itself has no standard unit field;
  consumers must use the declared mm convention.
- Caller radii are angstrom display sizes, not automatically inferred chemistry.
  Ordered fully occupied sites only; known unexpanded symmetry is refused.
- Shared OBJ site/vertex limits: 1..1000 sites, <=100000 vertices. Additional STL
  limits: <=50000 triangles, <=16 MiB. Integer subdivisions 0..4. Nonfinite,
  overflowing/underflowing, numerically collapsed geometry is refused.
- One fixed `solid` name, ASCII facet records, outward winding and finite unit
  normals; `.17g` vertex/normal tokens. Positions may be negative and are not
  shifted or wrapped. Binary STL is outside this API.
- Individual spheres may overlap or be disconnected. No surface union, bonds,
  materials/colors, lattice/PBC, site/component names or chemical round trip.
  No printable-single-solid, slicer/printer, PDB-input, Rust/WASM/CLI or performance
  claim. All output stays in memory; no output files are written by the API.

## Native fixture and reproduction

No Rust source or native binding changes are needed. The tests use a new Python
package copy with the existing **ca7e14e8** local native extension from the checked
core-removal task. That native artifact's SHA-256 is checked against its recorded
provenance. `chematic-core`, `chematic-crystal`, `chematic-mol` and the native
Python crystal binding source are byte-identical between that fixture and this
base. Exact paths/hashes are in `provenance.json`; this is **not** a fresh native
build or a published wheel validation. Python 3.13.6 / pytest 9.1.1 were used.

Set `CANDIDATE_PY` to the candidate package directory, `TRIMESH_PY` to trimesh 4.9.0
and `ASE_PY` to ASE 3.26.0. The native fixture was previously linked on this macOS
27 host with native InChI, installed SDK 26.5 and bundled LLD. No new native or
older-macOS compatibility claim is made.

```sh
PYTHONPATH="$CANDIDATE_PY" python3.13 -m pytest \
  crates/chematic-py/tests/test_crystal_stl.py \
  crates/chematic-py/tests/test_crystal_obj.py \
  crates/chematic-py/tests/test_crystal.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONPATH="$CANDIDATE_PY:$TRIMESH_PY:$ASE_PY" python3.13 \
  validation/results/checked-crystal-stl-2026-10-18/oracle.py
MYPYPATH="$CANDIDATE_PY" python3.13 -m mypy --no-incremental \
  --config-file validation/results/checked-crystal-stl-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_crystal_obj.py \
  crates/chematic-py/python/chematic/_crystal_stl.py \
  validation/results/checked-crystal-stl-2026-10-18/typing_valid.py
# typing_invalid.py: expect six errors, exit 1.
python3.13 -m mkdocs build --strict --site-dir /tmp/chematic-stl-site
```

ASCII profile references: [Library of Congress format description](https://www.loc.gov/preservation/digital/formats/fdd/fdd000506.shtml),
[trimesh STL API](https://trimesh.org/trimesh.exchange.stl.html).
