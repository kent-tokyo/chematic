# Explicit CUBE dataset selection

Open [Open Babel #2579](https://github.com/openbabel/openbabel/issues/2579)
requests the ordinary single-field CUBE dialect rather than the negative-atom
orbital dialect. Its proposed upstream [PR #2580](https://github.com/openbabel/openbabel/pull/2580)
is also open as checked on 2026-10-18. The local Open Babel 3.2.1 control still
writes negative NATOMS plus a dataset-ID record for a single scalar field.

chematic's native `VolumetricGrid.to_cube()` already writes single-dataset CUBE.
The missing capability is reading multiple datasets and explicitly selecting
one for that output. Python now supplies `read_cube_datasets`, with immutable
dataset tuples and an export sidecar retaining the exact original source.
The existing native Rust/Python/WASM scalar readers remain unchanged.

```python
from chematic import read_cube_datasets

source = """Synthetic orbital data
Two scalar fields, dataset index fastest
-1 0 0 0
1 1 0 0
1 0 1 0
2 0 0 1
6 4.0 .1 .2 .3
2 7 9
1 10 2 20
"""
document = read_cube_datasets(source)
assert document.dataset_ids == (7, 9)
assert document.datasets == ((1.0, 2.0), (10.0, 20.0))
selection = document.select(1)  # Explicit zero-based position, not orbital ID.
assert selection.dataset_id == 9
assert selection.omitted_indices == (0,)
assert selection.source_text == source
assert selection.to_grid().values.tolist() == [10.0, 20.0]
assert selection.text.splitlines()[2].split()[0] == "1"
```

## Data and output contract

The parser independently follows the
[h5cube dataset layout](https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html).
It accepts positive NATOMS with optional NVAL, or negative NATOMS followed by
a count and unique signed 64-bit dataset IDs. Negative NATOMS allows omitted
NVAL or NVAL=1. IDs may wrap across lines; their order is retained. Positive
NATOMS has `dataset_ids=None`, with no invented orbital labels. A zero-atom
positive-form scalar grid is accepted, matching the native API; negative zero
is refused because its sign has ambiguous dataset semantics.

Input values vary fastest over dataset, then the third grid axis, second axis
and first axis. `datasets[d]` is a flat scalar field with the third grid axis
fastest. Values must be finite decimal/e-notation floats. Fortran D notation,
extra values, missing values and duplicate IDs refuse. Float64 conversion is
explicit; lexical number spelling remains available in `source_text`.

`select(index)` requires an explicit valid nonnegative position. It returns
nonnegative NATOMS, no NVAL/dataset-ID record, and only the chosen scalar data,
rendered with 17 significant digits. Comments and axis/atom lines keep their
contents; the count/origin line is normalized and line endings become LF.
The exact original text, selected position/ID and omitted positions remain in
the sidecar. No implicit first-dataset selection or default file writing occurs.
Keep the sidecar when the omitted fields or original IDs are needed.

`to_grid()` creates a fresh native grid. Native element, effective nuclear
charge and finite geometry checks apply. CUBE atom charge is an ECP nuclear
charge, not an assigned ionic or partial charge. Only positive voxel counts
are accepted: all lengths remain in bohr. Ambiguous negative-axis unit flags
and ghost atoms refuse. No resampling, scalar-unit conversion, physical orbital
interpretation or GUI compatibility claim is made.

Hard caps are 8 MiB UTF-8 input, 50,000 atoms, 64 datasets, 262,144 grid points,
262,144 total scalar values and 16,384 characters per input line. The first,
atom, point and total-value limits can only be lowered. Scalar export is capped
at 16 MiB. Header counts are checked before reading/allocating the scalar data.
This API materializes a bounded document; it is not a file-streaming reader.

## Validation

106 tests pass, including 56 new dataset cases and the 50 existing format
regressions. They cover both dialects, wrapped IDs, CRLF, ECP nuclear charges,
nonorthogonal asymmetric grids, signed IDs, zero atoms, float precision,
ownership and actual 64-dataset/262144-value/50000-atom boundaries.

An independent ASE 3.26.0 oracle checks 192 authored inputs. It reads 240 fields
from positive-NVAL originals and all 480 selected scalar exports, containing
12,724 scalar values; values and geometry agree exactly. ASE accepts 24 of the
negative-form originals and refuses 72 multi-field originals. Those refusals
are recorded separately from successful cross-reads; authored exact values
also check the dataset/voxel mapping for both input forms.

A local Open Babel 3.2.1 two-field rewrite changes IDs 7/9 to 1/2 and ECP nuclear
charge 4 to 6. The new adapter retains the original IDs and charge; this is an
observed control, not a claim about every Open Babel version or API path.
VESTA, VisIt and quantum-chemistry programs were not run.

Reproduction scripts, logs and hashes are in
`validation/results/cube-datasets-2026-10-18/`. Tests reuse a hash-verified native
extension with unchanged relevant Cube/grid source. No fresh current-main
native build, published wheel or performance measurement is claimed.

## Python reference

::: chematic._cube_datasets
