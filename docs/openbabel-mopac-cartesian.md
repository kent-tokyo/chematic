# Checked MOPAC Cartesian geometry

The isolated candidate adds MOPAC Cartesian input/output to `chematic-mol`
and Python. This fills a format gap identified against Open Babel's
[`mop` / `mopcrt` / `mpc` format](https://openbabel.org/docs/FileFormats/MOPAC_Cartesian_format.html).
It is a geometry and source-metadata API: it does not run MOPAC, infer bonds,
interpret charge/spin keywords or validate a quantum-chemistry calculation.
The common-format converter and WASM surface are unchanged.

The implementation was independently written from public documentation.
No Open Babel or MOPAC implementation code or upstream test data was copied.

## Rust and Python interfaces

Rust exports `parse_mopac_cartesian`, `parse_mopac_cartesian_with_limits`
and `write_mopac_cartesian`. The structured types are `MopacCartesianGeometry`,
`MopacAtom`, `MopacCartesianRow`, `MopacOptimizationFlag`, `MopacParseLimits`
and `MopacError`. Geometry retains a physical keyword block, two description lines, atom rows, ordered
translation vectors and the original `source_text`.

```python
import chematic

text = chematic.write_mopac_cartesian(
    ["Si", "Cl"],
    [[0.123456789, -1.25, 2.5], [-1.5, 0.125, 0.5]],
    keywords="PM7 XYZ CHARGE=1",
    title="original synthetic",
    comment="Cartesian coordinates in Angstrom",
    labels=["site alpha", None],
    optimization_flags=[["0", "1", "1"], None],
    translation_vectors_angstrom=[[4, 1, 0], [-1, 5, 2], [1, 0, -6]],
    translation_flags=[["0", "0", "0"]] * 3,
)
geometry = chematic.parse_mopac_cartesian(text, max_atoms=1000)
assert geometry["periodicity"] == 3
assert geometry["optimization_flags"][1] is None
assert geometry["positions_angstrom"][0][0] == 0.123456789
```

`keywords` supports bounded standalone `+` continuations as described in the
[continuation contract](mopac-keyword-continuation.md).

The result dictionary has `keywords`, `title`, `comment`, `elements`,
`positions_angstrom`, `labels`, `optimization_flags`,
`translation_vectors_angstrom`, `translation_flags`, `periodicity` and
`source_text`. Atom-array lengths must agree. Translation flag count must
match the vector count. Missing optional arrays supply absent labels/flags
and no translation vectors. To make input for Open Babel 3.2.1, explicitly
supply a three-flag row for every atom and every translation vector.

## Retained source semantics

| Field | Contract |
|---|---|
| Header | [One to 64 keyword lines](mopac-keyword-continuation.md), then title and comment; spacing/empty lines retained |
| Identity | Canonical element symbol or unsigned atomic number; writer emits symbols |
| Label | Optional immediate `(text)`, at most 38 ASCII characters; empty label differs from absent label |
| Coordinates | Three finite Cartesian values in Angstrom; `E/e` and Fortran `D/d` exponents accepted |
| Flags | Either no flags or three interleaved string codes `-2`, `-1`, `0`, `1`, `2`, `T` |
| Signed flags | `+0`, `+1`, `+2` accepted and serialized as `0`, `1`, `2` |
| Periodicity | Zero to three Cartesian `Tv` rows at the end; source order, orientation and flags retained |
| Writer | Current validated fields, 17 fractional digits in scientific notation, LF and final blank line |

Absence of flags is retained as `None`, without assigning optimizer defaults.
Flag meanings follow the public
[optimization flag documentation](https://openmopac.net/Manual/Optimization_flags.html)
and [MINI documentation](https://openmopac.net/Manual/mini.html).
The API preserves each code without checking whether it is appropriate for
the keyword line. Labels, including labels beginning with `+` or `-`, are
opaque text; no formal charge is inferred from a label.

`Tv` vectors are read directly as Cartesian row vectors, following the
[translation-vector convention](https://openmopac.net/Manual/tv.html).
One/two/three vectors describe one/two/three periodic directions; 3D is
required only for the Open Babel cell comparison below. Vector norms must
be finite and positive. Normalized cross-product magnitude (2D) or absolute
determinant (3D) must exceed `1e-6`, rejecting dependent and nearly dependent
directions. Left-handed cells are retained.

The writer validates mutable data and never returns a stale `source_text`.
It preserves binary64 coordinates through read/write, including negative
zero, subnormal values and the largest finite value. It does not preserve
numeric spelling or line endings. The original input remains available as
`source_text` for callers that need an exact source record.

## Explicit boundaries

The accepted grammar is the bounded Cartesian subset described by the
[basic usage guide](https://openmopac.net/guides/basics/),
[Cartesian coordinate reference](https://openmopac.net/Manual/cart_coord.html)
and [label reference](https://openmopac.net/Manual/Labels.html).
The parser refuses:

- Internal/hybrid coordinates, connectivity indices, isotope-mass suffixes,
  partial-charge tails and other extra geometry fields.
- Dummy atoms, sparkle encodings and the reserved `At`, `Fr`, `Es` / atomic
  numbers `85`, `87`, `99` described in the
  [legacy element definitions](https://openmopac.net/Manual/elements_in_PM3.html).
  These spellings are refused in both reader and writer rather than interpreted
  as ordinary elements.
- Atomic-unit keyword `A0`, keyword continuations `++` / `&` outside quoted values, embedded
  header controls, multiline labels and non-ASCII labels.
- More than three translation vectors, labels on vectors, or atoms after
  the first vector.
- Additional sections or jobs after the first blank geometry terminator.
  A terminator is optional; trailing blank lines are accepted.
- Non-finite numbers, zero/dependent translation vectors and unsupported
  optimization codes.

Other keywords are retained as opaque text. Accepting a keyword line does
not establish Hamiltonian availability, element parameterization, optimizer
validity, electronic state or scientific correctness. Reading MOPAC output,
archive preambles and invoking calculations are outside this feature.

Default parsing limits are 64 MiB input, 16,384 bytes per line and one million
atoms. Callers can tighten these through Rust limits or Python keyword
arguments. The structured geometry also caps headers at 16,384 bytes and
atom count at one million; serialized output is capped at 64 MiB. Errors
include source line numbers when local to a record, and line zero for
whole-geometry validation. Python maps format/validation errors to `ValueError`.

## Measured interoperability, 2026-10-15

The oracle is **Open Babel 3.2.1**, installed build
`Open Babel 3.2.1 -- Jul 11 2026 -- 19:21:27`.
The generated corpus uses original synthetic recipes: eight seeds, zero to
three vectors, flags present/absent, labels present/absent and LF/CRLF.
Inputs are authored independently of the production writer.

| Check | Result |
|---|---|
| Generated exact structured roundtrips | 256 / 256 |
| Original manifest fixture with all six flags, labels, `D/d`, skewed cell | 1 / 1 |
| Open Babel read comparison: flagged 0D/3D, optional nonempty labels | 64 / 64 |
| Open Babel write comparison: same cases, `-xu` enabled | 64 / 64 |
| Unexpected mismatches in the stated comparison | 0 |
| Maximum read position / cell error | 0 / `8.881784197001252e-16` Angstrom |
| Maximum write position / cell error | `4.761904761840796e-6` / 0 Angstrom |

Read comparison uses `obabel -imop INPUT -ocjson -ab`, with a `2e-12`
Angstrom tolerance. `-ab` disables bond inference. Write comparison uses
`obabel -imop INPUT -omop -ab -xu`, with a `5.00001e-6` Angstrom tolerance
because this Open Babel writer emits five decimal places.

Measured differences are recorded separately:

- Open Babel replaces all optimization flags with `1`, omits labels,
  replaces the keyword line with `PUT KEYWORDS HERE`, uses the filename as
  the title and clears the comment. These fields are exact native-roundtrip
  checks, not Open Babel metadata parity claims.
- Open Babel's CJSON projection reports total charge zero for the synthetic
  `CHARGE=1` keyword. chematic retains the keyword and exposes no inferred
  electronic state from this format.
- Open Babel rejects the documented unflagged Cartesian rows and an empty
  parenthesized label in the isolated diagnostics. chematic preserves both.
- Open Babel reads `1.25D-1` as `1.25` in the isolated diagnostic;
  chematic reads `0.125`. Fortran exponents are outside the common oracle lane.
- Partial periodicity is an exact native contract; no Open Babel cell parity
  claim is made for one or two vectors.

These results establish the listed geometry contracts on synthetic cases.
They do not establish complete MOPAC support, a library-wide replacement
claim or successful execution of a MOPAC calculation.

## Reproduction and delivery

The original fixture manifest is
`validation/openbabel_issues/mopac/contract-fixtures.json`. It records
provenance, license, redistribution, minimization, oracle settings, expected
values and a normalized UTF-8/LF SHA-256. Generated recipes and their input
hashes are recorded by `scripts/validate_mopac_cartesian.py`.

```sh
cargo test -p chematic-mol --lib --test mopac_cartesian_contract --test fhiaims_contract
python -m pytest crates/chematic-py/tests/test_mopac_cartesian.py
python scripts/validate_mopac_cartesian.py --obabel /opt/homebrew/bin/obabel
```

Build/install the candidate Python extension before the Python commands.
The completed local validation used Python 3.13.6 / pytest 9.1.1: Rust
753 unit tests plus 15 integration tests passed; Python 148 tests passed,
including FHI-aims, EXYZ and atom-environment regressions. Formatting, scoped
Clippy with warnings denied, the unchanged refactor budget, binding inventory
and strict MkDocs are included in the evidence report.

The source candidate is `/private/tmp/chematic-mopac`, based on
`/private/tmp/chematic-atom-environment`. The incremental patch is
`validation/openbabel_issues/mopac-cartesian-2026-10-15.patch`, and the evidence
is `validation/results/openbabel-mopac-cartesian-2026-10-15.json`.
Only this document, patch and evidence are copied to the primary checkout;
unrelated source work is preserved. The feature has not been committed,
pushed or published.
