# TURBOMOLE geometry with preserved 3D lattice

Open Babel [#2822](https://github.com/openbabel/openbabel/issues/2822) reports
loss of unit-cell information during TURBOMOLE conversion. In the locally
pinned Open Babel 3.2.1, reading a standalone `$periodic 3` / `$lattice`
document and writing VASP still produces three zero lattice vectors.
The reverse writer direction already retains lattice information in that
version; this candidate does not claim both directions remain broken.

The new Rust/Python APIs read and write standalone TURBOMOLE coordinate
geometry while retaining the original 3D lattice basis, atom order and
fixed-atom flags. They do not infer bonds or implicit hydrogen atoms.

## Checked POSCAR conversion

The Python conversion functions connect the existing native format readers and
writers without discarding the lattice:

```python
from chematic import turbomole_to_poscar_checked, poscar_to_turbomole_checked

tmol_text = """$periodic 3
$lattice angs
4 0.5 0
0 5 0.75
0.25 0 6
$coord frac
-0.25 0.5 1.25 o f
0.25 0 0.125 si
0.5 0.25 -0.5 o
$end
"""
converted = turbomole_to_poscar_checked(tmol_text)
assert converted.source_indices == (0, 2, 1)  # O, O, Si in output order
assert converted.source_text == tmol_text
assert "Direct\n" in converted.text
restored = poscar_to_turbomole_checked(converted.text)
assert "$periodic 3" in restored.text
assert "$lattice" in restored.text
```

Both return a frozen `PeriodicGeometryConversion` with original text, converted
`text`, source/target format names and `source_indices`. Entry `i` maps output
atom `i` to its zero-based source coordinate row. Species are grouped in
first-appearance order, stably within each species, for POSCAR. The reverse
conversion uses the native POSCAR writer's same grouping; duplicate species
groups therefore also have explicit correspondence.

The lattice basis, handedness, unwrapped positions and fully fixed/free atoms
survive within floating-point conversion error. TURBOMOLE `f` maps to POSCAR
`F F F`; free atoms map to `T T T`. POSCAR scale factors and coordinate modes
follow the native parser. Output is Direct-coordinate VASP 5 POSCAR or Cartesian
bohr TURBOMOLE. Comments, original units and number spelling remain in
`source_text`, rather than being replayed in another format.

This geometry-only profile requires explicit 3D periodicity and explicit
POSCAR species. It rejects unknown TURBOMOLE groups, external-file references,
partial directional constraints and POSCAR velocity/predictor-corrector data.
No cell, bonds, charges, spin, POTCAR choice or simulation settings are inferred.
Each input is limited to 16 MiB UTF-8, 50000 atoms and 1024 bytes per line.
POSCAR count headers are bounded before native site allocation. No files are
written and no partial conversion is returned on failure.

Validation for this adapter is recorded in
`validation/results/periodic-geometry-conversion-2026-10-19/`: 288 ASE POSCAR
cross-reads plus 288 reverse geometry/constraint checks. The existing native
read/write APIs below remain independently usable.

```python
import chematic

geometry = chematic.parse_turbomole(tmol_text)
assert geometry["periodicity"] == 3
positions = geometry["positions_angstrom"]
lattice = geometry["lattice_angstrom"]  # three Cartesian row vectors

output = chematic.write_turbomole(
    geometry["elements"],
    positions,
    lattice_angstrom=lattice,
    frozen=geometry["frozen"],
)
restored = chematic.parse_turbomole(output)
original_text = geometry["source_text"]
```

`parse_turbomole` returns a dictionary containing `elements`,
`positions_angstrom`, `frozen`, `lattice_angstrom`, `periodicity` and
`source_text`. Nonperiodic input has `lattice_angstrom=None` and periodicity
zero. Writer element symbols must be canonical, for example `Si` rather
than `si`; the reader accepts ASCII case variations and the writer emits
lowercase TURBOMOLE symbols. Each coordinate row has exactly three components.
Optional frozen flags must have the same length as the atom arrays.

In Rust, `parse_turbomole` / `parse_turbomole_with_limits` return a
`TurbomoleGeometry` with ordered `TurbomoleAtom` records. Its public fields
are validated by `write_turbomole`, which returns a `Result` instead of
replaying a potentially stale original string. `TurbomoleError` includes a
one-based line where available; zero denotes document-level validation.

## Units and cell semantics

- Default coordinates and lattice lengths are in bohr. `angs` selects
  angstroms independently for the coordinate and cell/lattice groups.
  Explicit `bohr` is also accepted. `$coord angs` / `$coord bohr` are
  MCTC-style extensions, rather than a claim that every TURBOMOLE version
  accepts those coordinate modifiers.
- The conversion constant is the CODATA 2022 central value
  **1 bohr = 0.529177210544 angstrom**. This differs slightly from the
  older constant used by the pinned Open Babel binary.
- `$coord frac` uses the three supplied lattice row vectors directly,
  even if the cell groups occur later in the document. Coordinates outside
  the first cell are retained without wrapping.
- `$lattice` requires three rows of three components. Orientation and
  handedness survive; vectors are never converted to lengths/angles and
  reconstructed in another orientation.
- Alternatively, `$cell` accepts one row `a b c alpha beta gamma` for a
  3D cell. Lengths follow that group's units and angles are in degrees.
  Conversion reuses chematic's validated `UnitCell` implementation and its
  conventional crystallographic orientation. `$cell` together with
  `$lattice` is rejected as conflicting input.
- The writer emits canonical Cartesian bohr coordinates and a bohr
  `$lattice`, using scientific notation with 17 decimal places. It preserves
  the numerical geometry within floating-point conversion/rounding error,
  rather than preserving the source number spelling or coordinate mode.
- The `f` atom flag is preserved. No freeze flag means that atom is not fixed.

The primary layout sources are the
[TURBOMOLE 8.0 periodic input manual](https://manual.turbomole.org/v8.0/dft-calculations-for-periodic-systems/how-to-perform-a-calculation/),
its [fixed-atom definition](https://manual.turbomole.org/v8.0/calculation-of-vibrational-frequencies-and-vibrational-spectra/vibrational-frequencies-with-fixed-atoms-using-numforce/)
and the [MCTC format documentation](https://grimme-lab.github.io/mctc-lib/page/format-tmol.html),
retrieved on 2026-10-15. The Bohr conversion follows
[NIST CODATA](https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0).
The implementation was written independently from these field definitions.

## Checked boundary

Exactly one nonempty coordinate group and a terminating `$end` are required.
`$periodic` may be omitted for nonperiodic input, or explicitly specify 0 or
3 on its header. A 3D cell requires explicit `$periodic 3` and exactly one
cell definition. Recognized groups may appear in any order.

Duplicate groups, unknown groups/options, nonblank data after `$end`,
external-file references, incomplete rows, unsupported atom flags and
non-finite/overflowing numbers fail with an error. Fortran `D`/`d` exponents
are accepted. Comments beginning with `#` and blank lines are accepted.
Input text including comments and line endings is retained in `source_text`;
canonical writing emits only the interpreted geometry and fixed-atom flags.

Cell vectors must have finite positive lengths and finite positive absolute
volume. The magnitude of the normalized determinant must exceed `1e-6`,
equivalent to a normalized Gram determinant above `1e-12`. Singular and
numerically near-singular cells fail, including on writing mutated Rust data.
This is the coordinate-reader boundary; it does not claim that the separate
crystal computation API accepts every cell above this threshold.

Default parsing limits are 64 MiB input, 1,024 bytes per physical line and
1,000,000 atom rows. All three can be lowered or raised explicitly through
Rust `TurbomoleParseLimits` or the Python keyword arguments
`max_input_bytes`, `max_line_len` and `max_atoms`. Atom-vector reservation is
fallible and counts are checked while reading rows.

1D/2D periodic documents, complete calculation control files, charge/spin
settings, isotope/ghost labels and internal-coordinate definitions are
outside this API and are rejected when encountered. No TURBOMOLE execution,
new POSCAR adapter, graph-format registry entry, topology perception or
bond-order assignment is added. This is a geometry API for callers to
connect explicitly to their existing coordinate/structure workflow.

## Validation and local artifacts

Original synthetic fixtures under MIT OR Apache-2.0 are recorded in
`validation/openbabel_issues/turbomole/lattice.json` with provenance,
minimization, expected units, tolerance and payload SHA-256. No upstream
attachment or coordinate example was copied.

The independent expected-geometry oracle uses Python `Decimal` arithmetic
at 40-digit precision, rather than the new Rust parser or the lossy Open
Babel reader. Four atoms (24 permutations), six group orders, three
coordinate modes and two cell handedness cases produce **864 comparisons**.
Every case checks input interpretation, atom order, fixed flags, source text,
all nine lattice components and read/write/read restoration. A separate
Python decoder inspects the writer's canonical layout directly. Absolute
coordinate/lattice tolerance is `1e-12` angstrom, with zero mismatches.
Triclinic cell-parameter conversion, independent group units and malformed
input are covered by Rust/Python tests.

Open Babel 3.2.1 is used only for the upstream defect diagnostic and the
separate POSCAR-to-TURBOMOLE writer control. Its VASP reader needs an absolute
file path; the control uses a task-owned temporary POSCAR file. Different
Bohr conversion constants require `2e-9` angstrom coordinate tolerance for
the reader diagnostic and `1e-8` for the writer control. The TURBOMOLE
executable was not run.

Local checks:

- `chematic-mol`: 753 library tests and 14 new integration contracts passed.
- Python: 121 tests passed, including the new geometry tests and previous
  GRO, checked mmCIF-to-PDB and PQR-to-PDBQT regressions.
- Strict Clippy for `chematic-mol` / `chematic-py`, formatting, product naming,
  static binding inventory and the refactor budget passed; no allowances added.
- Strict MkDocs build passed.

Regenerate the independent record with the candidate extension on `PYTHONPATH`:

```sh
python scripts/validate_turbomole_lattice.py --obabel /path/to/obabel
```

The result JSON records fixture/source/extension hashes and oracle settings.
This is a local, uncommitted, unpublished source candidate in
`/private/tmp/chematic-turbomole`, preserving earlier candidates through
stereo-group inversion. The portable patch is
`validation/openbabel_issues/turbomole-lattice-2026-10-15.patch`; it isolates
this feature from the preceding source snapshot and includes source, bindings,
tests, fixture, validation, documentation and CI/inventory updates.
