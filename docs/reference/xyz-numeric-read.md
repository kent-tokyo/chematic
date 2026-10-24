# XYZ numeric read contract

XYZ and ASE-style ExtXYZ real fields now share a checked binary64 parser in
`chematic_3d::xyz::parse_xyz_real`. A nonzero decimal such as `1e-9999` previously
became zero in both native readers. The 3D XYZ reader, used by Python
`from_xyz` and WASM `mol_from_xyz`, also accepted NaN, infinity and overflowing
coordinates. These inputs now fail explicitly.

```python
import struct
import chematic

mol, coords = chematic.from_xyz("1\nexample\nC 5e-324 -0e-999999 0\n")
assert struct.pack(">d", coords[0][0]).hex() == "0000000000000001"
assert struct.pack(">d", coords[0][1]).hex() == "8000000000000000"
output = chematic.write_xyz_checked(
    ["C"], coords, notation="scientific", decimal_places=17,
)
assert chematic.from_xyz(output)[1] == coords
```

## Accepted values and refusals

- Finite normal and subnormal values retain the result of decimal-to-binary64
  rounding. This includes `5e-324` and the largest finite binary64 value.
- Literal positive and negative zero remain valid, even with extreme exponents
  such as `-0e-999999` and `0e+999999`. Only the mantissa determines whether the
  source decimal was zero; exponent digits must not cause a false refusal.
- Nonzero decimals that round to either sign of zero are refused. This includes
  scientific and fixed-point spelling and the boundary value `2.4e-324`.
  `2.5e-324` rounds to the smallest subnormal and remains accepted.
- NaN, infinities, overflowing values and malformed real tokens are refused.
  Existing decimal grammar is retained; Fortran `D` exponents, underscores and
  hexadecimal syntax are not added.

This check prevents nonzero-to-zero loss; it does not provide arbitrary-precision
decimal storage or reject every inexact conversion. For example, `0.1` still
rounds normally to its binary64 value. Numeric source spelling is not retained.

The contract applies to all three XYZ coordinate axes, all nine ExtXYZ lattice
entries and `R` atom-property components. `S` columns and raw frame info remain
strings: a label or `energy=1e-9999` is retained without numeric interpretation.
No unit conversion, bonding, charge assignment or coordinate regeneration is
introduced. Other formats and JSON number decoding are outside this change.

## Entry points and trajectory behavior

Native `chematic_3d::parse_xyz` reports the existing `InvalidLine` error. The
shared scalar helper distinguishes invalid syntax, nonfinite values and
underflow through `XyzRealError`. Native `chematic_mol` retains the existing
`NonFiniteCoordinate` variant for nonfinite coordinates; coordinate underflow
uses `InvalidAtomLine` with axis, source line, token and reason. Cell and atom
data errors retain `InvalidLattice` and `InvalidPropertyValue` context.

Native single-frame, all-frame, in-memory, buffered-file and batch readers use
the same checks. Native trajectory iterators stop after an error, as before.
Python `from_xyz`, `from_extxyz` and `from_extxyz_all` raise `ValueError`.
Python file batches retain their recovery policy: a bad frame is counted as
rejected and later valid frames remain available. WASM topology/JSON reads throw;
the existing bounded batch JSON keeps a rejected entry at the original index.

Writer precision is unchanged. Legacy fixed-decimal XYZ output can still round
small coordinates. Use checked scientific XYZ output with sufficient precision
or the round-tripping ExtXYZ writer when retaining finite binary64 values is
required. Python and real-WASM tests check those finite read/write paths with
explicit handling of negative zero at the JavaScript JSON boundary.

## Bounded local evidence, 2026-10-25

`validation/xyz_numeric_read/cases.json` contains 14 authored finite/signed-zero
controls with explicit hexadecimal binary64 oracles and 16 refusal tokens:
five underflows, six nonfinite/overflow inputs and five malformed decimals.
Python's independent float conversion and `Decimal` zero classification check
the authored oracles. Rust additionally checks 4,096 deterministic bit patterns,
skipping nonfinite patterns, through 17-place scientific text.

| Python coordinate read path | Finite controls with exact expected bits | Explicit refusals on 16 refusal inputs | Accepted refusal inputs |
|---|---:|---:|---:|
| Earlier local XYZ | 14/14 | 5/16 | 11/16 |
| Earlier local ExtXYZ | 14/14 | 11/16 | 5/16 |
| Candidate XYZ | 14/14 | 16/16 | 0/16 |
| Candidate ExtXYZ | 14/14 | 16/16 | 0/16 |

Refusals are reported separately and do not count as preserved values. The
baseline is the preserved local `d6de9d6e` build from PR #1221, not a release.
The four relevant native reader/Python wrapper source files are verified
byte-identical to comparison parent main `6f069119`. Receipts retain baseline
reference/source hashes, candidate working-source hashes and extension hashes.

Open Babel 3.2.1 was observed through `-ixyz -oxyz -ab` on the same one-atom
inputs. Its emitted XYZ retains the expected coordinate bits for 8/14 finite
controls. On the refusal inputs it emits an observable coordinate record for
15/16; one input yields no observable record. These output observations do not
isolate its reader from its writer and do not establish 16 parser acceptance or
refusal results. The [documented XYZ writer](https://openbabel.org/docs/FileFormats/XYZ_cartesian_coordinates_format.html)
uses five decimal places, which itself loses small finite values. The separate
ORCA-AICCM `exyz` format is not measured. No toolkit-wide, whole-format, chemical
graph, speed or unseen-holdout claim follows from this packet.

The observer separately tests signed-zero/subnormal bits, emitted nonfinite
values, missing/changed atom layouts, refusal-versus-preservation accounting,
and the distinction between real zero and nonzero underflow. Missing output or
an unobservable layout is not credited as an explicit refusal. Processes have a
10-second and 1-MiB combined-output bound.

The complete receipts are:

- `validation/xyz_numeric_read/baseline-2026-10-25.json`
- `validation/xyz_numeric_read/candidate-2026-10-25.json`

Reproduce with a fresh Python extension on `PYTHONPATH`:

```sh
python scripts/validation/check_xyz_numeric_read.py --output /tmp/xyz-numeric-read.json
python -m unittest scripts.tests.test_xyz_numeric_read
```

Binding tests consume the same packet, including every coordinate axis, lattice
and `R` property fields, plus native stop-on-error and Python/WASM recovery
behavior. The original source tokens and all comparator outputs/errors remain
in the receipts; no third-party molecular fixture is copied.
