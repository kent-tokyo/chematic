# Periodic TURBOMOLE/POSCAR conversion

Base: `00596fb1ea5397524f98e3789c41f17ee8616fdc`. Checked 2026-10-19.
Open Babel [#2822](https://github.com/openbabel/openbabel/issues/2822) is open.
Chematic already reads/writes 3D TURBOMOLE lattice data and POSCAR separately;
the two new typed Python functions supply the missing checked workflow adapter.
No native parser or chemistry algorithm is replaced.

## Results and scope

- **162 Python tests passed**, including 28 new conversion tests plus existing
  TURBOMOLE, crystal and pymatgen JSON tests.
- ASE **3.26.0** independently cross-reads **288** produced POSCAR documents,
  checking the full oriented lattice, Cartesian coordinates, species and fixed
  atom indices against authored geometry. **288** reverse conversions check the
  native-parsed TURBOMOLE result against those ASE-read structures. No mismatches.
- Four cells include orthogonal, skewed, rotated and left-handed bases. Seed
  282219 selects unwrapped positions; 1/3/8 atoms cover repeated species and
  none/all/alternating frozen masks. Maximum position error is
  `3.552713678800501e-15` angstrom; maximum lattice-component error is
  `8.881784197001252e-16` angstrom. These are bounded geometry comparisons,
  not universal interoperability or simulation execution evidence.
- Open Babel **3.2.1** reproduces zero VASP lattice rows on an original Si
  TURBOMOLE input. This diagnostic is separate from adapter comparison counts.
  Its reverse VASP-to-TURBOMOLE writer already works, as recorded in the earlier
  native lattice packet; no new upstream reverse-direction bug is asserted.
- Scoped strict mypy passes the helper and valid root consumer; two invalid
  argument errors are expected. Imported legacy declarations are silent, so
  this does not claim package-wide strict typing.

Tests exercise stable species permutations, source/CRLF retention, scale modes,
cell/bohr/angstrom units, full frozen constraints, refusing partial constraints
and velocities, invalid/singular/unknown groups, UTF-8 and byte/line caps, and
oversized/malformed count headers rejected before native allocation.

The profile has no topology, charge, spin, POTCAR or simulation settings.
Positions are not wrapped. Original formatting/comments remain in source_text;
output is canonical geometry. Unknown/unrepresentable data refuses. Native
floating-point conversion and existing geometry validation still apply.

## Reproduction and runtime boundary

Local candidate Python is staged in
`/private/tmp/chematic-periodic-conversion-python`, using the existing v1.1.7
native extension from `/private/tmp/chematic-blog-v117-venv`. Its SHA256 is in
`oracle-summary.json`. TURBOMOLE core/binding and POSCAR core/crystal binding
files have no diff between tag v1.1.7 and the base. This is not a fresh native
build or a new published-wheel test. No Rust source changed.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PYTHONPATH=/private/tmp/chematic-periodic-conversion-python \
python3.13 -m pytest crates/chematic-py/tests/test_periodic_geometry.py \
  crates/chematic-py/tests/test_turbomole.py \
  crates/chematic-py/tests/test_crystal.py \
  crates/chematic-py/tests/test_pymatgen_json.py -q -p no:cacheprovider

PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/private/tmp/chematic-periodic-conversion-python:/private/tmp/chematic-xyz-atom-data-oracle \
python3.13 validation/results/periodic-geometry-conversion-2026-10-19/oracle.py

MYPYPATH=/private/tmp/chematic-periodic-conversion-python \
python3.13 -m mypy \
  --config-file validation/results/periodic-geometry-conversion-2026-10-19/mypy.ini \
  crates/chematic-py/python/chematic/_periodic_geometry.py \
  validation/results/periodic-geometry-conversion-2026-10-19/typing_valid.py
```
