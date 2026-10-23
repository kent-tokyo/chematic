# ORCA Cartesian integrity

`parse_orca_cartesian` reads a single embedded Cartesian geometry and returns
its original numeric values with a coordinate-unit tag. It does not convert
Bohr to Angstrom. The structured parser also refuses repeated coordinate
specifications, empty xyz blocks and multiplicity zero, instead of selecting
the final geometry or returning an invalid state.

The [ORCA coordinate documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/coordinates.html)
specifies Angstrom by default, `Angs`/`Bohrs` keywords and `%coords` unit
settings. The [input documentation](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/input.html)
describes separate mechanisms for multiple jobs. This API implements the
bounded geometry profile below; it does not validate an executable ORCA job.

## APIs and supported input

- Rust: `parse_orca_cartesian`, `parse_orca_cartesian_with_limits`,
  `OrcaCartesianGeometry`, `OrcaCoordinateUnits`, `OrcaCartesianError`.
- Python: `parse_orca_cartesian(text, *, max_input_bytes=16_777_216,
  max_atoms=100_000)` returns a dictionary.
- WASM: `orca_cartesian_geometry_json(text)` returns JSON.

```python
geometry = chematic.parse_orca_cartesian(
    "! Bohrs\n* xyz -1 2\nO 2.5 0 0\n*\n"
)
assert geometry["coordinate_units"] == "bohr"
assert geometry["coords"] == [[2.5, 0.0, 0.0]]
```

The `chematic.orca-cartesian.v1` result contains `coordinate_units`,
`units_explicit`, ordered `elements`, `coords`, `frozen` flags and `atom_extras`,
plus total `charge` and `multiplicity`. Atom extras are opaque parsed strings;
their basis or annotation semantics are not interpreted. No bonds or per-atom
formal charges are inferred. Other job settings are outside this projection.

Exact `Angs` and `Bohrs` simple keywords are case insensitive, including across
multiple `!` lines. Missing keywords mean `angstrom` with
`units_explicit=false`. Repeated identical declarations are accepted;
contradictory declarations are refused. Other keywords beginning `ang`, `bohr`
or `units` are refused as unsupported unit-like tokens. This deliberately
bounded check is not a general ORCA keyword grammar.

One nonempty `* xyz charge multiplicity` block is required. Any `%coords`
block is refused because it may override geometry or units. External xyz/gzmt
files, internal coordinates and missing geometry are refused. The structured
parser already refuses `$new_job`; this view does not resolve compound jobs,
inherited geometry or external files. Unrelated blocks remain opaque and do
not imply that the full job is valid.

## Existing APIs and preservation boundaries

`parse_orca_input` and `orca_input_to_json` already retain raw unit keywords
alongside coordinate numbers. They share the new duplicate/empty/zero-spin
refusals. Existing `OrcaInputError` variants remain unchanged.

The older WASM `orca_input_coords_json` has no unit field. It now accepts only
the checked Angstrom profile and refuses Bohr, conflicting/unsupported units
and `%coords`. Its successful JSON keys are unchanged. Call the new tagged
API for Bohr inputs. Rust `OrcaCoords::to_molecule` alone still returns raw
numbers because that enum has no enclosing keyword context.

This is a parsed geometry contract, not original-document identity. Comments,
numeric spelling and source bytes are not represented. The existing
`write_orca_input`/`write_orca_input_json` still round Cartesian coordinates to
six decimal places and normalize layout. The checked writer below preserves
represented fields but cannot restore source information already omitted.

## Checked structured output

`write_orca_input_checked` preserves the structured input's comments, ordered
keywords and opaque blocks, plus total charge/spin, ordered elements, xyz
values, freeze flags and atom extras. Coordinates use 17 significant decimal
digits in scientific notation. Every finite binary64 coordinate, including
signed zero and subnormals, must retain its bits after reparse before output
is returned. Units are neither converted nor guessed: the existing
Angs/Bohrs/default-unit Cartesian profile applies.

- Rust: `write_orca_input_checked(&OrcaInput)` and
  `write_orca_input_checked_with_limits`, with `OrcaWriteLimits` and
  `OrcaWriteError`. Strict JSON counterparts are
  `write_orca_input_json_checked[_with_limits]`.
- Python: `write_orca_input_checked(input, *, max_output_bytes=67_108_864)`.
- WASM: `write_orca_input_json_checked(input_json)`.

```python
document = chematic.parse_orca_input(text)
output = chematic.write_orca_input_checked(document)
reparsed = chematic.parse_orca_input(output)
```

Python and JSON require `comments`, `keywords`, `blocks` and `coords`. The
coordinate object requires `kind: "xyz"` (or WASM's `type` alias), `charge`,
`multiplicity` and `atoms`. Each atom requires `element`, `x`, `y`, `z` and
exactly three boolean `frozen` entries. Optional `extra` may be omitted/null.
Each block requires `name`, `raw` and boolean `has_end`. Unknown fields,
duplicate JSON keys (including aliases), boolean numeric values and missing
required fields are refused. Python dictionaries cannot represent duplicate
keys; they use the same strict adapter after JSON encoding.

The writer requires a nonempty embedded xyz geometry and positive multiplicity.
It refuses conflicting/unsupported units, any `%coords` block and internal or
external geometry. It also refuses fields that would change on structured
reparse: for example, multiline keyword injection, noncanonical extra-token
spacing, non-comment comment entries or block text the parser would trim.
Failure returns no partial output. The legacy writer remains available with
its previous defaults and rounding.

The guarantee covers fields represented in `OrcaInput`. Inline comment
placement, numeric spelling, original line endings and source bytes are not
retained. Unknown calculation settings and atom extras remain opaque; preserving
them does not validate an executable ORCA job or their chemical meaning. No
connectivity or per-atom formal charges are inferred.

| Checked writer surface | JSON input | Output text | Atom cap |
| --- | ---: | ---: | ---: |
| Rust defaults | 64 MiB | 64 MiB | 100,000 |
| Python | 16 MiB after encoding | 1 byte–64 MiB, configurable | 100,000 |
| WASM | 1,000,000 bytes | 1,000,000 bytes | 10,000 |

Rust limits are configurable. Native output uses bounded appends and applies
the structured parser's default line/keyword/block budgets on readback.
WASM also checks readback against its existing 1,024-physical-line reader
budget: an otherwise empty two-line wrapper permits 1,022 atom rows. This
earlier line cap can bind before the atom cap. Python encodes its dictionary
before checking the transport byte cap; none of these budgets is a process
memory cap. JavaScript callers should pass raw JSON when signed zero matters;
`JSON.stringify` can normalize it before the writer receives it.

### Output measurement (2026-10-24)

The authored packet `validation/orca-checked-2026-10-24/cases.json` contains
24 positive controls: 12 common geometries, six metadata/unit controls and six
numeric boundaries. Expected fields and coordinate bits come from authored
inputs, independently of either engine's output.

| Twelve common Angstrom geometries | Checked Python/WASM | Legacy Python/WASM | Open Babel 3.2.1 |
| --- | ---: | ---: | ---: |
| All coordinate bits retained | 12/12 | 2/12 | 2/12 |
| Maximum coordinate absolute error (Å) | 0 | 4.999999999e-7 | 4.567205905914307e-6 |

The common lane supplies identical ordered elements and coordinates, total
charge zero and multiplicity one. Legacy controls call the old APIs in the
same candidate runtime; their source sections are verified equal to main
`38b4c92c`. This is not a separately rebuilt main/release executable.
Open Babel runs `-ixyz -oorcainp -ab`, so its result measures the complete
XYZ-to-ORCA pipeline rather than identical internal object state. Source
comments and calculation keywords are not compared in that pipeline lane.

Both checked bindings retain all represented metadata and coordinate bits in
all 24 positive controls. Legacy bindings retain the metadata in all 24, but
their rounding changes numbers in 20. WASM/native refuse 38 invalid or
unsupported raw-JSON cases; Python refuses the 30 dictionary cases. These
candidate-only unit/freeze/numeric/refusal controls are not extra Open Babel
failures. Native tests also run 4,096 deterministic bit-pattern trials,
skipping nonfinite patterns and explicitly refusing NaN/infinities.

The independent Python observer checks emitted state, units, ordered species,
flags, extras, comments, keywords, blocks and coordinate bits. Four observer
tests detect altered metadata, signed zero, rounding, missing/extra atom rows,
nonfinite/partial output and data after the coordinate block. Results and raw
worker outputs are retained in `validation/results/orca-checked-2026-10-24.json`
and its compressed packet, including source/runtime/fixture hashes.

```sh
python3 scripts/validation/validate_orca_checked.py \
  --python-root /path/to/candidate-python \
  --wasm-module /path/to/candidate/chematic_wasm.js
```

Workers have a 30-second timeout and 1 MiB combined pipe budget, with no process
memory cap. This is a local development measurement of the declared writer
profile, not ORCA execution, a published release or toolkit-wide superiority.

## Reader resource bounds

| Surface | Input bytes | Atom cap | Physical line cap | JSON output |
| --- | ---: | ---: | ---: | ---: |
| Rust default | 64 MiB | 1,000,000 | 1,000,000 | caller controlled |
| Python default | 16 MiB | 100,000 | 1,000,000 | 16 MiB |
| WASM | 1,000,000 | 10,000 | 1,024 | 16 MiB |

Other structured-parser budgets cover line length, keywords, blocks and block
bytes. Rust limits are configurable; Python exposes byte and atom limits.
Atom limits are checked before parsing/pushing the next atom. WASM's existing
line budget binds before its atom budget: a three-line wrapper permits 1,021
atom rows. JSON caps are checked after serialization; these are operation
budgets, not process memory caps. Python/WASM JSON parsing retains signed zero;
subsequent JavaScript `JSON.stringify` can normalize it.

## Reader development measurement (2026-10-22)

The authored packet `validation/orca-cartesian-2026-10-22/cases.json` contains
16 positive controls (eight inputs each with LF and CRLF) and 29 refusals.
Positive controls independently specify coordinate binary64 bits, unit tags,
charge/spin, atom order, freeze flags and opaque extras. They include signed
zero, a subnormal value, explicit/default units and Bohr values left unscaled.

| Same 45 inputs | Candidate Python | Candidate WASM | Frozen baseline, each surface |
| --- | ---: | ---: | ---: |
| New tagged view: all fields verified | 16 | 16 | API absent |
| New tagged view: refused | 29 | 29 | API absent |
| Structured parser: accepted / refused | 29 / 16 | 29 / 16 | 39 / 6 |
| Legacy coordinate-only view: returned / refused | unavailable | 6 / 39 | WASM 29 / 16 |

Ten formerly accepted inputs are newly refused on both structured surfaces:
six repeated-coordinate variants, three zero-spin variants and an empty xyz
block. The remaining 13 tagged-view refusals are unsupported projection forms,
not additional structured-parser failures. Legacy WASM returns remain
untagged; the candidate guard restricts them to Angstrom.

The frozen runtime is retained PR candidate `891bb83f`, attested against its
PQR report. Its entire ORCA parser source is verified equal to base main
`978669e4`. It is not a rebuilt main/release executable; dependencies and runtime
identities remain separate. Raw outputs, identities and hashes are retained in
`validation/results/orca-cartesian-2026-10-22.json` and its compressed packet.
Observers have controls for missing/wrong units, scalar/boolean confusion,
signed zero, partial results and extra rows.

Local Open Babel 3.2.1 advertises `orcainp` as write-only and refuses it as an
input format. This demonstrates reader coverage only. It is not a per-field
preservation comparison, malformed-input quality ranking or toolkit-wide win.

Reproduce with separately retained baseline and candidate builds:

```sh
python3 scripts/validation/validate_orca_cartesian.py \
  --python-root /path/to/candidate-python \
  --wasm-module /path/to/candidate/chematic_wasm.js \
  --baseline-python-root /path/to/baseline-python \
  --baseline-wasm-module /path/to/baseline/chematic_wasm.js
```

Workers have a 30-second timeout and 1 MiB combined output budget, with no
process memory cap. Runtime tests also exercise native/Python atom limits,
WASM physical-line limits and exact UTF-8 input-byte boundaries.
