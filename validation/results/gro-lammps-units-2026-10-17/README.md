# GRO to LAMMPS units and periodic-cell validation

The pinned base and native artifact/source hashes are in `summary.json` and
`source-sha256.json`. Open Babel #2669 was checked live on 2026-10-17 Japan.
The original three-atom fixture uses a 2×3×4 nm orthogonal cell and known
positions. Open Babel 3.2.1 scales its positions to Angstroms but replaces that
cell with -0.5..9.5 Angstrom bounds on all three axes. The expected bounds are
0..20, 0..30, 0..40 Angstroms. Raw output and the original probe are retained.

The native baseline has separate typed GRO and LAMMPS readers/writers, without
an explicit unit/geometry adapter. Three new tests fail with AttributeError
against a retained older native module; current base source independently
lacks this API. This is an API-absence baseline, not a full rebuild of current
main or evidence of an existing chematic cell-replacement defect.

The candidate converts to atomic data with caller-supplied types/masses and
mandatory real/metal units. Positions, cell edges and tilts scale by ten;
velocities scale by 0.01 or ten respectively. The source cell is not fitted to
atom extent. Generated IDs are unique by source row, with original IDs/names
and exact source text retained in sidecars. No topology, charges, elements,
force-field coefficients or image flags are inferred. General oriented cells,
mixed velocities, invalid assignments and unrepresentable scaling refuse.

Local checks:

- Rust: 863 existing library tests with crystal enabled + five new integration
  tests passed (`cargo test -p chematic-mol --features crystal --lib --test gro_lammps --offline`).
- Python 3.13.6: 75 passed, including 32 new adapter cases and 43 existing GRO
  cases; zero failures/skips.
- Independent NumPy geometry: 96 cases covering four orthogonal/triclinic cells,
  real/metal units, present/absent velocities and all six row permutations.
  Fractional-coordinate reconstruction and determinant ratios verify geometry;
  velocities, ID/type order and exact source retention are checked separately.
- Native parse/write round trips preserve typed LAMMPS data.
- Strict Clippy for chematic-mol/chematic-py, all targets; format, stub, refactor,
  binding inventory/manifest, naming, unsafe-surface and strict MkDocs: passed.
- Static inventory records the added Python function (214); WASM remains 317.

No LAMMPS engine executable was available, so no simulation run or engine
read_data validation is claimed. LAMMPS input scripts must separately specify
the report's unit and atom style. Outside-cell coordinates are retained in the
data text; a simulation can subsequently remap them according to its boundary
settings. Molecular force-field setup is outside this atomic adapter.

Native builds used one job, no incremental compilation/debug info, opt-level 1,
strip=none and macOS dynamic-lookup flags. Compressed logs retain actual local
checks; no published-wheel, completed hosted-CI or main-integration claim is made.

Reproduce with the built package on PYTHONPATH:

```sh
python3 validation/openbabel_issues/gro-lammps/probe.py --compare
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q \
  crates/chematic-py/tests/test_gro_lammps.py \
  crates/chematic-py/tests/test_gro.py
```
