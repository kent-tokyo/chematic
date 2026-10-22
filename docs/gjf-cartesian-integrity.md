# Gaussian Cartesian geometry integrity

The Cartesian GJF reader now separates an optional integer freeze code from
`x y z`. Previously, `C -1 1.25 2.5 3.75` became `(-1, 1.25, 2.5)` and the
last coordinate disappeared. Short rows could also be skipped while a partial
molecule was returned. Both behaviors are corrected in the native parser and
its existing Python `from_gjf`/common-conversion callers.

This source candidate adds a geometry API that exposes values previously
hidden by Python's atoms-only reader. It is a geometry projection, not a
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

Python exposes the same fields in a dictionary:

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
to six decimals and fills missing coordinates with zero. This change does not
claim source round trips or checked writing. The Python common conversion
bridge still uses the atoms-only GJF input reader: use `parse_gjf_geometry`
explicitly to obtain coordinates and job metadata. No new CLI/WASM API is added.

## Validation and Open Babel scope

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
