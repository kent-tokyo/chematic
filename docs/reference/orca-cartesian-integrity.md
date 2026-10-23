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
six decimal places and normalize layout; no new lossless writer is added.

## Resource bounds

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

## Development measurement

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
