# GAMESS-US explicit internal coordinates

`write_gamess_zmat` writes a GAMESS-US geometry template with Cartesian
`$DATA` and explicit simple internal-coordinate definitions in `$ZMAT`.
This addresses the requested `$ZMAT` conversion capability in
[Open Babel #2523](https://github.com/openbabel/openbabel/issues/2523),
verified open on 2026-10-17.

```python
import chematic

text = chematic.write_gamess_zmat(
    ["O", "H", "H"],
    [[0., 0., 0.], [0.96, 0., 0.], [-0.24, 0.93, 0.]],
    [[], [0], [0, 1]],
    charge=0,
    multiplicity=1,
)
```

Rust exposes the same function with element, position and reference slices
plus `GamessZmatWriteOptions::new(charge, multiplicity)`. Python accepts
canonical element symbols. Coordinates are in Å; references are zero-based,
distinct earlier atom indices. Each atom has `min(row_index, 3)` references.
References need not be chemically bonded; the caller defines the coordinate
system. No connectivity is inferred from XYZ.

`xyz_to_gamess_zmat(text, references, charge=..., multiplicity=...)` is the
single-frame plain XYZ adapter. It retains every atom, including hydrogen,
and accepts at most 1 MiB of input with 4,096-byte lines. It requires exactly
four fields per atom row and refuses key/value comment metadata (`=`),
additional frames and trailing non-XYZ content. Plain comments are not copied
into the output title; use the explicit `title` option when needed. This
conservative profile avoids silently discarding extxyz properties or cells.

The template uses `COORD=UNIQUE`, `UNITS=ANGS`, C1 and the source atom order
and Cartesian frame. It does not use the reorienting `COORD=ZMT` `$DATA`
profile. `$ZMAT IZMAT(1)` lists, for each supplied reference chain, a distance
(code 1), angle (code 2), and torsion (code 3), using one-based GAMESS atom
numbers. `NZVAR` is the number of emitted definitions: 3N−6 for non-linear
inputs with at least three atoms, one for two atoms, and zero for a single
atom. A single-atom template omits `$ZMAT`. C1 has no blank orientation card.
The choice of Cartesian `$DATA` keeps the explicit IZMAT indices tied to
input atom numbers.

The writer reuses the checked Cartesian-to-Z-matrix converter to reject
coincident atoms, invalid reference lists and linear/near-linear reference
angles or planes. It validates the Cartesian values again after decimal
rounding, so requested precision cannot collapse a required reference frame.
It supports 1–17 significant digits, default 17; that default round-trips
finite source f64 values, with negative zero normalized. Lower precision is
explicit quantization and does not promise exact geometry or chirality.

The profile accepts 1–1,024 atoms. Every card fits 80 columns, and output is
bounded to 256 KiB. If 17-digit extreme coordinates exceed a card's width,
the writer returns an error; callers can explicitly choose fewer digits.
The title is 1–80 printable ASCII bytes without `$` or `!`. Header injection,
zero multiplicity, nonfinite coordinates, invalid precision or undefined
rounded internal coordinates return an error without a partial document.

Charge and multiplicity are explicit metadata. The writer does not validate
their electronic-state plausibility or select a basis, method, SCF type or
run type. The result is geometry input to incorporate into a caller-configured
GAMESS job; it is not a complete calculation setup. Isotopes, atom-specific
basis/ECP data, dummy atoms, periodicity, linear bends, delocalized/redundant
coordinates and geometry constraints are outside this profile.

The independent implementation follows the official
[GAMESS input manual](https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/input.pdf),
`$CONTRL`, `$DATA` and `$ZMAT` sections. See the
[GAMESS documentation index](https://www.msg.chem.iastate.edu/gamess/documentation.html)
for the primary manuals. Local tests and interoperability evidence are in
`validation/results/gamess-zmatrix-2026-10-17/`. No GAMESS job was executed;
Open Babel coordinate cross-reading does not validate execution of `$ZMAT`
or the suitability of a calculation setup.

Local validation passed 805 existing Rust library tests, 9 new GAMESS tests,
13 existing Z-matrix integration tests and 133 Python tests (20 new). A
640-record original corpus matched the direct XYZ adapter; RDKit 2026.09.1
checked 3,584 internal values with maximum difference 3.98e-13 (Å or degrees,
as appropriate). Open Babel 3.2.1 retained Cartesian atom order/coordinates
in 80 independent cross-reads. These checks do not establish GAMESS runtime
acceptance, optimization convergence, chemistry perception or basis coverage.
