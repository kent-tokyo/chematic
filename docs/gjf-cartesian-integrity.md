# Gaussian Cartesian geometry integrity

The Cartesian GJF reader now separates an optional integer freeze code from
`x y z`. Previously, `C -1 1.25 2.5 3.75` became `(-1, 1.25, 2.5)` and the
last coordinate disappeared. Short rows could also be skipped while a partial
molecule was returned. Both behaviors are corrected in the native parser and
its existing Python `from_gjf`/common-conversion callers.

This source candidate adds geometry reading and checked Cartesian writing APIs.
The reader exposes values previously hidden by Python's atoms-only reader. It is a geometry projection, not a
complete Gaussian job representation or a released-package comparison.

## APIs

Rust exports `GjfGeometry`, `parse_gjf_geometry` and
`parse_gjf_geometry_with_limits` from `chematic_mol`. The result contains:

- `mol`: atoms in input order, without inferred bonds or assigned formal charges.
- `coords`: three finite Cartesian values per atom, in input numeric units.
- `freeze_codes`: one `Option<i32>` per atom; absent, `0` and `-1` remain distinct.
- `charge` and `multiplicity`: the job-level charge and positive spin multiplicity.

The existing `parse_gjf` tuple API shares the corrected parser. Its tuple omits
freeze codes; use the new result type when they must be retained.

Python exposes the same fields in a dictionary, plus an ordered `elements` list
for passing geometry to the checked writer:

```python
import chematic

text = "# HF/STO-3G\n\ngeometry\n\n-1 2\nC -1 1.25 2.5 3.75\n\n"
g = chematic.parse_gjf_geometry(text)
assert g["coords"] == [[1.25, 2.5, 3.75]]
assert g["freeze_codes"] == [-1]
assert (g["charge"], g["multiplicity"]) == (-1, 2)
```

A job charge does not identify which atoms carry formal charge. The reader
therefore does not distribute it over the atoms. Freeze codes are retained
as signed integers; optimization semantics and electron-count/spin compatibility
are not evaluated.

## Accepted geometry and refusals

Atom rows have exactly four fields (`element x y z`) or five fields
(`element integer-freeze-code x y z`). Whitespace and comma delimiters, plain
symbols, numeric suffix labels such as `Cl12`, atomic numbers, LF/CRLF, and
`!` comments are supported. Numeric values use finite decimal/Rust float syntax.
This implements a limited Cartesian profile of the
[Gaussian molecule specification](https://www.conflex.co.jp/gaussian_support/molspec.php).

Incomplete rows, extra atom fields, invalid or overflowing freeze codes,
non-finite coordinates and non-positive multiplicity return an error.
Z-matrices, atom isotope/fragment/MM annotations, ONIOM suffixes, fragment
charge/spin pairs and multiple `--Link1--` jobs are outside this profile.
They are refused rather than partially interpreted or silently selected.

Coordinates retain their input values. Route keywords such as `Units=Bohr`
are not interpreted or converted. The previous unconditional Ångström
label in the API documentation was inaccurate. Gaussian log standard-orientation
coordinates remain Ångström values.

`GaussianParseLimits` applies to both native entry points: defaults are 64 MiB
input, 1 MiB per line, one million physical lines, 100,000 sections and one
million atoms. Python's new geometry API defaults to 16 MiB and 100,000 atoms,
with `max_input_bytes`/`max_atoms` overrides and the remaining native limits.

Route options, title, basis sets, connectivity sections, original labels,
comments and other job sections are not returned. No coordinates or bonds are
generated. The existing GJF writer still omits freeze codes, rounds coordinates
to six decimals and fills missing coordinates with zero. The checked writer described below is a separate API; the legacy writer
keeps its existing behavior. The Python common conversion
bridge still uses the atoms-only GJF input reader: use `parse_gjf_geometry`
explicitly to obtain coordinates and job metadata. No new CLI/WASM API is added.

## Checked Cartesian writing

Rust exports `write_gjf_cartesian(elements, coords, freeze_codes, options)`,
`GjfWriteOptions`, `GjfCoordinateUnits::{Angstrom, Bohr}` and `GjfWriteError`.
The input consists of explicit geometry arrays, not a molecular graph.
`GjfWriteOptions::new(charge, multiplicity, units)` requires the electronic
state and units. `title` and `route` can be set on the options. The default
route is a template (`#`), not a validated calculation setup.

Python provides the equivalent bounded interface:

```python
# The visible input above uses angstroms. The reader does not infer route units.
output = chematic.write_gjf_cartesian(
    g["elements"], g["coords"],
    charge=g["charge"], multiplicity=g["multiplicity"],
    units="angstrom", freeze_codes=g["freeze_codes"],
    route="# HF/STO-3G", title="Explicit Cartesian geometry",
)
back = chematic.parse_gjf_geometry(output)
assert back["freeze_codes"] == g["freeze_codes"]
```

Every coordinate is written in scientific notation with 17 significant digits.
The native writer checks that parsing each scalar reproduces its f64 bits.
Signed zero and finite subnormal values retain their bits. Counts must match
exactly, and missing coordinates are never replaced by zeros. Native callers
provide one optional freeze code per atom; Python's omitted `freeze_codes`
means no atom has an explicit flag. Explicit `None`, `0` and `-1` remain distinct.

The selected units describe the supplied numbers. They are emitted as exactly
one `Units=Ang` or `Units=AU` keyword, without numerical rescaling, following the
[Gaussian Units convention](https://www.conflex.co.jp/gaussian_support/units.php).
The reader still returns input numeric units without interpreting the route:
callers must choose the matching writer unit. No source-unit inference or
arbitrary Gaussian job conversion is claimed.

The writer refuses empty/over-limit element arrays, non-elemental symbols,
missing/extra coordinate or flag rows, non-finite coordinates, zero multiplicity,
and headers containing controls or Unicode line/paragraph separators. Headers
are capped at 1024 UTF-8 bytes. Routes must be ASCII, start with one `#` and
contain no `!` comments. A conservative lexical check rejects Units keywords
and their alphabetic prefixes so the route cannot supply a second ordinary
unit declaration. This is not a Gaussian route parser: it does not establish
that arbitrary keywords describe a valid calculation or that the coordinates
will be used by the requested calculation.

The hard limit is 100,000 atoms; the resulting text is under 16 MiB. Python
bounds element/coordinate/flag iterators before materializing them, checks
exactly three coordinate values per row, and rejects bool for numeric
coordinates and integer state/flag fields. Input arrays are not modified.
Route/title/source formatting, basis-set sections, connectivity, isotopes,
fragment and ONIOM metadata are outside this explicit geometry representation.
No default conversion route is changed, and no new CLI/WASM writer is added.

### Paired writer comparison

`validation/gjf-cartesian-2026-10-22/write-cases.json` defines 12 original
geometry controls: four base geometries with three deterministic transforms.
The common lane supplies the same ordered species, coordinate numbers and
angstrom units to the checked API, unchanged legacy `write_gjf` API, and
Open Babel 3.2.1's XYZ → GJF route with bonding disabled. Python's standard
float parser independently observes the decimal output; chematic is not the
comparison's semantic oracle.

| Common writer lane | Bitwise coordinate matches | Maximum absolute coordinate error (Å) |
| --- | ---: | ---: |
| Checked API | 12/12 | 0 |
| Legacy API | 3/12 | 3.456790111933117e-7 |
| Open Babel 3.2.1 | 3/12 | 4.679135802376777e-6 |

Ordered elements match in all common outputs. This measures the selected
writing workflow, not a graph or full-job round trip. The legacy control uses
the current candidate runtime; the report verifies that its native/Python
writer function bodies equal main `654993b832f1431ce4137675102368340b049981`.
It is not described as a separately rebuilt main executable.

Separate candidate-only lanes pass eight numerical boundary controls, eight
freeze-code/state/unit controls with independent observation and public API
readback, and ten output-free refusals. These lanes are not counted as Open
Babel source-retention failures. Native property checks also cover 4096
deterministic bit patterns (excluding non-finite values) plus explicit numeric
boundaries, and the exact 100,000/100,001 atom boundary.

The writer report and raw outputs are
`validation/results/gjf-write-2026-10-22.json` and `.packet.json.gz`.
They include the fixture, source, extension and binary hashes. Observer controls
reject truncated/extra atom rows, non-finite output, duplicate/unsupported units,
zero multiplicity and multiple documents, and distinguish signed zero. Run with:

```sh
python3 scripts/validation/validate_gjf_write.py --obabel /path/to/obabel
cargo test -p chematic-mol --test gjf_cartesian_write
python3 -m pytest crates/chematic-py/tests/test_gjf_cartesian_write.py
```

Gaussian execution, arbitrary job validity, and toolkit-wide superiority have
not been tested. Numerical extremes and transformed geometries are developer
controls, not an independent chemistry holdout.

## Reader validation and Open Babel scope

The original CC0 fixtures are in
`validation/gjf-cartesian-2026-10-22/cases.json`. Forty-five cases each run with
LF and CRLF. Expectations are authored numeric/species values, not output from
one of the engines. Main at `7c34941f52c5e795fd3ca93b075aa31bc9fa748f` was built
before changing the reader, and its native probe binary was retained.

| Lane | Before | Candidate Rust | Candidate Python |
| --- | ---: | ---: | ---: |
| Cartesian values, atom count, charge/spin match | 2/20 | 20/20 | 20/20 |
| Numeric input-unit/signed-code controls match | 0/4 | 4/4 | 4/4 |
| Malformed input refused | 20/54 | 54/54 | 54/54 |
| Outside-profile input refused | 2/12 | 12/12 | 12/12 |

Python additionally checks exact freeze-code values and ordered species via XYZ
serialization. The legacy native tuple probe checks coordinates, atom count,
charge and multiplicity; it has no freeze-code field. The two numeric control
fixtures do not claim that their route/code combinations are valid Gaussian jobs.

Open Babel **3.2.1** advertises GJF as write-only, and the local binary rejects
`-igjf` with `Invalid input format`. Its stdout contains usage text. This is a
read-direction coverage difference, not a malformed-input validation score.
The pinned [Open Babel format implementation](https://github.com/openbabel/openbabel/blob/openbabel-3-2-1/src/formats/gaussformat.cpp)
also declares this direction.

Ten separate writer controls use XYZ → Open Babel GJF → candidate reader.
They retain ordered species and coordinates within the writer's five-decimal
rounding envelope; maximum coordinate error is `3.2109876535457715e-6`.
These controls do not test preservation of source freeze codes or job metadata.
Gaussian execution and Open Babel-wide superiority have not been tested.

The report and complete compressed process outputs are
`validation/results/gjf-cartesian-2026-10-22.json` and its `.packet.json.gz`.
They pin the fixture, source, extension, probe and Open Babel binary hashes.
External commands have a ten-second and combined 1 MiB pipe budget; process
memory is not capped. Crashes and usage text are not counted as safe refusals.
A pilot assumption that unsupported input produced empty stdout was corrected;
the final report records that measurement correction and both output pipes.

## Reproduction

Build the committed `validation/gjf-cartesian-2026-10-22/probe.rs` against
`chematic-mol` at the baseline commit and again against the candidate, using
`rustc --edition 2024 --extern chematic_mol=<target>/debug/libchematic_mol.rlib
-L dependency=<target>/debug/deps`. The report records build settings. Then use
a Python package freshly built from the candidate:

```sh
python3 scripts/validation/validate_gjf_cartesian.py \
  --baseline-probe /path/to/baseline-probe \
  --candidate-probe /path/to/candidate-probe \
  --obabel /path/to/obabel
cargo test -p chematic-mol --test gjf_cartesian
python3 -m pytest crates/chematic-py/tests/test_gjf_geometry.py
python3 -m unittest discover -s scripts/tests -p 'test_gjf_cartesian_validation.py'
```

The Python fixture packet and observer controls run in CI alongside the native
and binding tests. This is a development packet, not an independent holdout.
