# Checked XYZ per-atom data

Python `parse_xyz_atom_data` and `XyzAtomDataFrame` read a single XYZ frame
with **one explicitly sized real payload per atom**. Default `components=3`
reads the headerless `xyz+vec` layout requested in
[Open Babel #1703](https://github.com/openbabel/openbabel/issues/1703).
Use `components=1` for scalars or `components=9` for flat tensor data.

```python
from chematic import parse_xyz_atom_data, from_extxyz

source = "3\nOriginal vector fixture\nC 0 0 0 .125 -.25 .5\nCl 1.5 .1 0 -1 .75 .2\nO -.3 2 .2 .001 .002 .003\n"
frame = parse_xyz_atom_data(source)
assert frame.elements == ("C", "Cl", "O")
assert frame.values[0] == (0.125, -0.25, 0.5)

# Give the unnamed payload a caller-selected meaning for an extxyz consumer.
text = frame.to_extxyz("displacements")
native = from_extxyz(text)
assert native["properties"]["displacements"][0] == [0.125, -0.25, 0.5]
```

The record exposes immutable tuples `elements`, `positions` and `values`,
plus `components`, `comment` and exact decoded `source_text`. Atom order,
all payload columns, negative zero and representable finite binary64 values
are retained. CRLF, spacing and original number spellings remain in
`source_text`; `to_xyz()` normalizes them to LF and 17 significant digits.
The API returns data/text and writes no file.

Coordinates and payload values remain in their input units. A vector is not
automatically a force, velocity or displacement. Nine components stay flat
in their original order: no row/column-major tensor convention is imposed.
No bonds, formal charges, isotope identity or physical properties are inferred.
The optional downstream `from_extxyz` graph projection retains its existing
chemistry behavior; this feature does not establish topology equivalence.

## Explicit extxyz conversion

The original input has no `Properties=` descriptor. `to_extxyz(name)` writes
`Properties=species:S:1:pos:R:3:name:R:N`, with N equal to `components`.
The default name is `atom_data`. Names are 1–32 character ASCII identifiers
starting with a letter;
`species`, `pos`, `positions`, `numbers` and `Z` are reserved. Some consumers
assign special meanings to names such as `forces`; choose the name explicitly.

The original free-form comment is JSON-tagged quoted metadata under
`comment`, so even comments like `T`, `1 2` or `_JSON {}` remain strings in
ASE. ASE 3.26.0 decodes this automatically. The current native extxyz reader
retains the raw marker: recover it with
`json.loads(native["info"]["comment"][6:])`. The input record's `comment`
and `source_text` are always directly available.

Already named `Properties=` or periodic `Lattice=` headers are refused by
this adapter; use the existing `from_extxyz` reader for those inputs. The
new parser never guesses a column width from the first row.

## Checked profile

- One frame, 1..50,000 atoms, at most 16 MiB of UTF-8 input. Input is `str`
  or immutable `bytes`; invalid encoding and lone Unicode surrogates fail.
- Count/comment/atom lines are at most 4,096 UTF-8 bytes excluding LF/CRLF.
  LF and CRLF are accepted; bare CR and unsupported control/line-separator
  characters are refused. Trailing ASCII whitespace is allowed.
- Exactly one canonical element symbol H..Og, three coordinates and N
  real payload values per atom. Dummy/isotope aliases and atomic numbers
  are refused. The count line uses unsigned decimal digits.
- Numeric tokens use decimal notation with optional `e`/`E` exponent.
  NaN, infinity, overflow, nonzero values that round to zero, underscores,
  Fortran `D` notation and mismatched/extra columns are refused.
- Extra frames or nonblank trailing data are errors. Empty frames are
  refused. Count limits are checked before building atom rows.
- Wrong Python types raise `TypeError`; encoding, format, numerical and
  limit violations raise `ValueError`. The record is frozen and fully
  owns its parsed data.

This is an in-memory standard-library Python adapter, not a streaming
trajectory reader. Lattices, multiple independent payload fields, tensor
reshaping and native Rust/WASM bindings are outside this addition.

## Sources and local evidence — October 17, 2026

The unnamed layout follows the issue's request. Named output follows the
[ASE extxyz documentation](https://docs.ase-lib.org/_modules/ase/io/extxyz.html),
also retained in the installed ASE 3.26.0 reader's docstring. Implementation
and fixtures are original project work; no external reader code was copied.

Baseline main `35b932f9` has no `parse_xyz_atom_data` / `XyzAtomDataFrame`.
The existing named extxyz reader rejects an undeclared seven-column row;
plain XYZ reading returns coordinates without the extra payload. On the
reused native artifact these behaviors were measured. Open Babel **3.2.1**
rewrote the original seven-column fixture as four-column XYZ, dropping the
vector payload. No claim is made about every other Open Babel format/API.

Validation passed **234 tests**, including 72 new regressions and 162 existing
compatibility tests. A seeded comparison cross-read **150 records** through
both **ASE 3.26.0** and the existing native extxyz reader: 2,670 atom rows and
11,570 payload scalars matched binary64 bytes, including negative zero,
subnormal values and maximum finite payloads. An additional ASE check covered
all 118 element symbols. These are synthetic data-interchange checks.

Strict helper mypy, public valid/invalid consumers, source/stub/binding/unsafe/
product-name checks and strict MkDocs passed. Consumer typing follows imports
silently; strict cleanliness of the entire legacy package is not claimed.
Evidence and reproduction commands are in
`validation/results/xyz-atom-data-2026-10-17/`:

```sh
python3.13 validation/validate_xyz_atom_data.py --output /tmp/xyz-data-crossread.json
```

The native extension was reused for import, old tests and cross-reading;
the new parser/serializer has no native calls. This is not a fresh native
build, package release, topology validation or performance result.
