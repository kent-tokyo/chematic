# CIF to atom-sphere OBJ — 2026-10-17

OpenBabel [#2680](https://github.com/openbabel/openbabel/issues/2680) asks for
batch CIF-to-OBJ conversion. The issue is open (live GitHub API checked this
date). chematic has native periodic structures but no existing OBJ export.

`write_crystal_obj_checked` adds a Python atom display mesh, with caller-selected
radii, typed/frozen results, bounded tessellation, globally indexed vertex
normals and independent named spheres. It does not export bonds, cell topology,
materials or a merged solid. The docs include a per-file batch recipe. This is
not an OpenBabel OBJ output comparison, chemical interchange format, general
surface mesher or printable-solid claim.

Evidence:

- `pytest.log`: 232 passed = 51 new, 19 existing crystal, 162 compatibility tests.
  Tests cover radial geometry, outward winding, every edge having two incident
  faces, Euler characteristic, normals, deterministic names, unwrapped sites,
  native source immutability, explicit CIF expansion/rejection, disorder,
  numeric overflow/collapse and site/mesh/output resource boundaries.
- `mesh-oracle.json`: trimesh 4.9.0 + ASE 3.26.0; 73 conditions, 532 spheres,
  107148 vertices, 212168 triangles. Every sphere has Euler characteristic 2,
  closed faces, consistent outward winding, positive volume and expected radial
  normals. All five tessellations, three cell bases, unwrapped inputs and all
  118 elements are covered. Twelve CIF cases are independently authored and
  read with ASE. Maximum center error is 1.173e-13 Å and maximum radius error is
  2.887e-15 Å. Loader settings are `process=False`, `group_material=False` and
  `skip_materials=True`; ordinary default material grouping can merge objects.
- `nacl.cif`, `nacl.obj`, `nacl-preview.png`: full eight-site NaCl cell; 528
  vertices, 1024 faces. The preview is rendered from independently reloaded OBJ
  triangles and was visually inspected. Blue/green colors are assigned only by
  the example renderer; the OBJ contains no material assets. Reproduce using
  `validation/render_crystal_obj_example.py`.
- `typing-valid.log`: strict checks of the new helper and positive public
  consumer. `typing-invalid.log`: five expected errors. These are scoped checks,
  not a claim of full-package strict mypy cleanliness.
- `mkdocs.log`: strict docs build passed.
- Existing source gates pass: Python stubs, cross-binding manifest, binding
  surface, refactor budget, unsafe surface, product naming and diff whitespace.

## Native fixture and scope

This is a Python-only implementation over the existing periodic API. Tests use
a recorded existing native extension, not a fresh native build or wheel:

- Native source commit `71d9c12b0c16a8e6b23e4078479d5754230c831b`, version 1.1.6.
- SHA-256 `c6db880b6f5ca4ea0388f238c24e8a5b55751daab1084cf0d1d7b56e2427257b`.
- `git diff --quiet 71d9c12b0c16a8e6b23e4078479d5754230c831b HEAD -- crates/chematic-crystal crates/chematic-py/src/crystal.rs`
  succeeds at the candidate base: all crystal Rust sources and their Python
  binding are unchanged from the fixture source. There is no new native method.
- The isolated package combines candidate Python files with that extension.
  Native Cartesian getters use the existing NumPy dependency. trimesh, ASE and
  Matplotlib are development-only oracles/renderers, not runtime dependencies.
- The used `CifSymmetryStatus` attributes and `from_cif(expand_symmetry=...)`
  keyword already exist natively; this change corrects their missing stubs/docs.
- `provenance.json` records the base, source hashes, tools and extension hash.

Limits: 1..1000 ordered sites, subdivisions 0..4, ≤100000 vertices, ≤200000
faces, ≤16 MiB ASCII. A CIF with known unexpanded symmetry or partial/mixed
occupancy is refused. Native explicit-operation expansion is reused; completeness
of the listed symmetry operations is not independently verified. Coordinates
are not wrapped. Overflow or relative vertex-direction distortion over 1e-8 is
refused. Each sphere is separately closed; overlaps or coincident sites are not
unioned. The output retains coordinate scale in angstroms by convention; OBJ
has no standard machine-enforced unit field.

## Reproduction

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-crystal-obj-python python3.13 -m pytest crates/chematic-py/tests/test_crystal_obj.py crates/chematic-py/tests/test_crystal.py crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/chematic-crystal-obj-python:/private/tmp/chematic-obj-oracle:/private/tmp/chematic-xyz-atom-data-oracle python3.13 validation/validate_crystal_obj.py --output validation/results/crystal-obj-mesh-2026-10-17/mesh-oracle.json
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/chematic-crystal-obj-python:/private/tmp/chematic-obj-oracle MPLCONFIGDIR=/private/tmp/chematic-obj-matplotlib python3.13 validation/render_crystal_obj_example.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-crystal-obj-site
```

Strict mypy 1.8.0 configuration: `strict=True`, `follow_imports=silent`,
`follow_imports_for_stubs=True`, `ignore_missing_imports=True`; explicitly check
`crates/chematic-py/python/chematic/_crystal_obj.py` and
`validation/typing/crystal_obj_valid.py` with `--no-incremental`. Negative check:
set `MYPYPATH=crates/chematic-py/python` and use `--follow-imports=silent
--ignore-missing-imports --no-incremental` on
`validation/typing/crystal_obj_invalid.py`; expect five errors and exit 1.

Format reference:
[Library of Congress Wavefront OBJ description](https://www.loc.gov/preservation/digital/formats/fdd/fdd000507.shtml).
No OpenBabel implementation code or tables were copied.
