# Checked XYZ output precision

`write_xyz_checked` writes the supplied elements and Cartesian coordinates in
angstroms with explicit decimal formatting. It addresses the configurable XYZ
output request in [Open Babel #1693](https://github.com/openbabel/openbabel/issues/1693).
The default is equivalent to a minimum-width `20.12f` coordinate field. Existing
`Mol.to_xyz`, `chematic_3d::write_xyz`, `XyzWriter` and `write_xyz` defaults retain
their existing byte format.

```python
import chematic

xyz = chematic.write_xyz_checked(
    ["O", "H", "H"],
    [[0.1234567890123, 0.0, 0.0], [1.08, 0.0, 0.0], [-0.2, 0.9, 0.0]],
    comment="caller geometry",
    decimal_places=12,
    min_width=20,
)

# 16 fractional digits in scientific notation = 17 significant digits:
# finite binary64 values, including signed zero/subnormals, round-trip exactly.
precise = chematic.write_xyz_checked(
    ["C"], [[1.0 / 3.0, 5e-324, -0.0]],
    notation="scientific", decimal_places=16, min_width=0,
)
```

Rust exports `XyzWriteOptions`, `XyzNotation`, `XyzWriteError` and
`write_xyz_checked(&[XyzAtom], comment, &options)` from `chematic_mol`.
The options' `validate(comment)` method permits preflight checks before gathering
geometry. Errors carry an optional zero-based atom index and a detail string.

## Formatting contract

| Option | Default | Accepted values |
| --- | --- | --- |
| `notation` | `"fixed"` | `"fixed"` or `"scientific"`; Rust uses `XyzNotation` |
| `decimal_places` | 12 | 0–17 digits **after** the decimal point in either notation |
| `min_width` | 20 | 0–64 characters per coordinate, excluding a separating space |
| `max_output_bytes` | 67,108,864 | Positive integer up to 64 MiB, including headers and padding |

Fields are right-aligned to their minimum width. Larger values expand the field;
width never truncates a coordinate. A separate space always separates fields,
including when all coordinates exceed the minimum width. Output contains the atom
count, one comment line and one canonical element/x/y/z row per supplied atom,
with LF and a final newline. Negative zero retains its sign. No locale-specific
decimal separator or caller-provided format string is used.

`decimal_places` is a fractional-digit count, **not** a significant-digit count.
Scientific notation with 16 or 17 fractional digits preserves finite f64 values
on parsing. Smaller scientific precision rounds the mantissa. Fixed notation
rounds to the requested absolute decimal grid; small nonzero values can become
zero and precision is intentionally lost. A requested rounding that cannot be
read back as a finite f64 is refused, for example the maximum f64 rounded to
`2e308`. Binary inputs near decimal half-way points round according to their
actual f64 value, not the decimal spelling originally used to construct them.

RDKit 2026.09.1 `MolFromXYZBlock` accepts the tested fixed-point output but
refuses scientific XYZ tokens, including conventional signed/two-digit exponents.
Use fixed notation for that reader. Open Babel 3.2.1 accepts both profiles;
its default XYZ re-writer rounds to five decimal places. These external reader
limits do not change the native scientific binary64 round-trip contract.

## Checked geometry and resource bounds

The writer emits precisely the caller-supplied element order and coordinates,
including explicit hydrogen rows. It performs no bond perception, atom reordering,
H addition, charge/isotope projection, unit conversion or electronic-state
inference. Plain XYZ does not represent those graph annotations. Elements in the
Python API must be canonical periodic-table symbols (including `H`, `Cl`, `Og`);
`D`, `T`, aliases, atomic numbers and dummy/query atoms are refused.

Python accepts bounded lazy iterables for elements, positions and each three-value
position. It reads at most 100,001 element rows and at most four values per
position to detect a limit/count error. Element/position counts must agree.
Iterator exceptions propagate. Integer options require ordinary Python ints and
exclude bool; coordinate booleans are refused. Coordinates must all be finite.
Empty geometry is supported as a zero-count XYZ frame.

At most 100,000 atoms and 64 MiB of aggregate output are accepted. Allocation of
the growing Rust output buffer is fallible. No partial document is returned and
there is no file I/O. Comments are limited to 4,096 printable ASCII bytes without
`=`; this refuses frame injection and unsupported extxyz `Properties=`/`Lattice=`
headers. Use the dedicated [extxyz APIs](format-capabilities.md) for structured
per-atom properties or a cell. This addition controls plain XYZ output; it does
not add an ADF output reader, custom printf directives or a streaming writer.

## Validation

Dated evidence lives under `validation/results/xyz-precision-2026-10-17/`, with
original synthetic fixtures, source/native hashes, commands and test logs.
`validation/validate_xyz_precision.py` compares the emitted decimals with Python's
independent formatter, checks scientific binary64 round trips and cross-reads
ordinary geometry with RDKit/Open Babel (12 successful cross-reads and four
explicit RDKit scientific refusals). Extreme numeric cases use numeric-only
parsing and do not assert that external chemistry engines support extreme scales.

Local validation passed 815 Rust tests (805 existing plus 10 new), 168 Python
tests (55 new), 12,288 independent decimal-format comparisons and 12,288 exact
scientific f64 numeric checks. The native Rust reader additionally round-trips
4,096 finite bit patterns across all three axes.

Public typing consumers and the exact new signature in isolation pass type checks.
Full-package strict mypy still reports 81 pre-existing missing generic type
arguments in unrelated stubs; package-wide strict typing is not claimed.
