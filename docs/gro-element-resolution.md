# GRO records with explicit element resolution

Open Babel [#2817](https://github.com/openbabel/openbabel/issues/2817) reports
GRO atom names `BE`, `CE` and `FE` being read as elements B, C and F.
GRO stores atom names but has no element column: a name such as `CA` can refer
to an alpha carbon or calcium, depending on the source topology. A reader
cannot resolve every such name reliably from the coordinate file alone.

The new Rust/Python reader preserves names and defaults to unresolved elements.
Callers can supply canonical element symbols from their topology or other
trusted source, in the same order as the input atom rows.

```python
import chematic

raw = chematic.parse_gro(gro_text)
assert raw["atoms"][0]["element"] is None

resolved = chematic.parse_gro(gro_text, elements=["Be", "Ce", "Fe"])
assert [a["element"] for a in resolved["atoms"]] == ["Be", "Ce", "Fe"]
positions_nm = [a["position_nm"] for a in resolved["atoms"]]
original_gro = resolved["source_text"]
```

The explicit symbols are caller assignments, not element perception or a
check against a GRO element field. Names remain unchanged after assignment.
For `CA`, supplying either `C` or `Ca` is supported; omitting elements leaves
it unresolved. Invalid/noncanonical symbols and wrong element counts fail.

In Rust, `parse_gro`/`parse_gro_with_limits` return a `GroFrame`. Call
`frame.with_elements(&[Element::BE, Element::CE, Element::FE])` to obtain an
independent resolved copy. The original frame remains unchanged.

## Stored fields and units

Each atom record contains the literal five-column residue/atom numbers,
trimmed residue and atom names, `position_nm`, optional `velocity_nm_per_ps`,
and an optional element. The original text preserves field padding, numeric
spelling and line endings. Wrapped/repeated numbers are kept as supplied;
row order, rather than a supposedly unique serial number, identifies atoms.

The frame keeps the title verbatim, including any `t=` text. It also returns
`box_values_nm` in file order and `box_vectors_nm` as three row vectors.
Three-value boxes receive the documented zero off-diagonals. Zero boxes and
nine-value boxes are retained without inventing a periodic unit cell or
checking simulation-cell feasibility, invertibility or supported orientation.

Coordinates and boxes remain in **nm**, and velocities in **nm/ps**. No implicit
conversion to angstroms occurs. Consumers requiring angstroms must explicitly
multiply coordinate values by 10.

## Reader boundary

The layout follows the [GROMACS GRO documentation](https://manual.gromacs.org/current/reference-manual/file-formats.html#gro),
identified as GROMACS 2026.4 when inspected on 2026-10-15.

- Exactly one frame: additional nonblank frames/data are rejected.
- Fixed five-character identifier/name fields; atom rows use ASCII fields.
  The free title may contain Unicode.
- Decimal-point spacing determines numeric width `n+5`, for positive decimal
  precision `n>=1`. Optional velocity triples use that same width and `n+1`
  decimal places. Partial velocity triples are rejected.
- Finite numeric values, exact atom counts and complete rows are required.
  Box records contain exactly three or nine whitespace-separated real values.
- Default limits are 64 MiB input, 2,000,000 atoms and 1,024 bytes per line.
  Custom limits are available in Rust and Python. Truncated declared atom
  counts fail before reserving the atom vector; reservation is fallible.

This is a typed frame reader. It introduces no bonds, molecular graph,
implicit hydrogens, charges, isotope interpretation, general writer,
trajectory iterator or automatic element inference. It does not register
GRO as a new graph-conversion format.

## Validation and provenance

`validation/openbabel_issues/gro/two-letter-elements.gro` is an original
three-row fixture under MIT OR Apache-2.0. It uses widely separated synthetic
coordinates and copies no upstream structure, attachment or POSCAR file.
The validation JSON records its provenance, minimization, normalized hash,
expected values/units and the pinned Open Babel 3.2.1 diagnostic.

On Open Babel 3.2.1, this fixture becomes B/C/F in XYZ output. Explicit
Be/Ce/Fe assignments produce those exact elements in the native reader. A
separate explicit-element XYZ control reads back correctly through Open Babel,
isolating the GRO name interpretation from coordinate conversion.

The validator also generates all 720 row permutations of six original names
(BE/CE/FE/CA/CL/NA). Known generator widths give an independent field oracle
across three precisions, optional velocities and both box shapes. Element
assignments, names, numeric fields and complete source text have zero
mismatches. These checks do not validate force-field typing or connectivity.

```sh
cargo test -p chematic-mol --lib --test gro_contract \
  --test mmcif_pqr_integration --offline --quiet
python -m pytest crates/chematic-py/tests/test_gro.py \
  crates/chematic-py/tests/test_new_formats.py \
  crates/chematic-py/tests/test_input_errors.py -q
python scripts/validation/validate_gro_elements.py --obabel /opt/homebrew/bin/obabel
```

Local verification passed 772 Rust tests (15 new GRO contracts), 114 Python
tests (43 new GRO API tests), and all 720 permutation cases. Clippy with
warnings denied, formatting, refactor budget, binding inventory, product naming
and strict MkDocs passed.

This is a local source candidate built on the previously preserved candidates,
not an upstream Open Babel change or a published chematic package.
