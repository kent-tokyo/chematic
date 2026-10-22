# Checked CJSON input in common conversion

The common conversion bridges used the permissive CJSON graph reader. For two
atoms, a supplied `formalCharges: [1]` became `[1, 0]`; `order: []` became a
single bond, and `order: [7]` also became a single bond. The source's explicit
arrays were incomplete or unsupported, but conversion emitted a changed graph.
The immutable CJSON record reader already refused these inputs.

CLI `convert`, Python `convert_format` and WASM `convert_common_format` now use
that existing strict profile before projecting CJSON input. Typed validation is
shared with the record reader. Rust exposes additive
`parse_cjson_checked[_with_limits] -> CjsonProjection`; Python exposes
`from_cjson_checked(str) -> (Mol, coords)`. `CjsonProjection` is an owned graph
and source-order Cartesian coordinate vector in Angstrom. The shared helper
avoids source-text copies and opaque-path collection when only projection is
requested. Empty coordinates cover absent 3D and zero-atom projections; the
immutable record's presence flag remains available for source inspection.

```python
mol, coordinates = chematic.from_cjson_checked(source)
# Ordinary conversion now applies the same input checks.
cml = chematic.convert_format(source, "cjson", "cml")
```

## Profile and acceptance changes

Present formal-charge/isotope arrays require one correctly typed value per
atom; Cartesian coordinates require exactly 3N numeric values. Bond orders
require one supported token per declared bond. Duplicate decoded JSON keys
(including opaque subtrees and escaped spellings), self/duplicate undirected
bonds, dangling indices and unsupported order tokens refuse. The existing
record version rule (absent or integer 1), typed `properties.name` and complete
finite root `partialCharges` datasets apply too.

Only absent optional arrays use the existing local defaults: zero formal
charge, unspecified isotope, absent 3D coordinates, and single bonds if order
is absent. No supplied short array is padded. Supported order values remain
exactly 1, 1.5, 2, 3 and 4. Value 4 remains the existing **local aromatic alias**;
this is not a statement that another implementation uses that meaning.
Integral floating storage such as `2.0` is accepted by the local numeric profile.
Atomic numbers, formal charges, isotopes and indices retain the existing local
integer-token requirement. These storage rules are separate from general JSON
Schema conformance; outside-profile controls are not all invalid CJSON.
The strict reader validates representable graph fields and typed source fields;
valence, aromaticity-model equivalence, inferred H/stereo and unknown chemistry
remain outside this validation.

Common conversion acceptance intentionally tightens for inputs that the
permissive reader previously repaired or ignored. Legacy native `parse_cjson`,
Python `from_cjson` and WASM low-level graph readers keep their behavior. Record
retention/inspection and checked-conversion policies keep their existing
semantics. Source text, metadata, partial-charge datasets, stereo and coordinate
retention still follow each output route's existing contract. Input validation
does not make arbitrary output lossless. Use `convert_with_report` for audited
loss/preservation policy, or the record API for immutable source retention.

Native checked projection uses `CjsonParseLimits`: 256 MiB input, 128 configured
JSON depth, 10,000,000 items per array, 16 MiB string/key length, 1,000,000 atoms
and 2,000,000 bonds. The input byte cap precedes JSON parsing; serde_json's
recursion guard also applies. Tree limits apply after unique-key parsing,
including opaque fields; atom/bond caps apply in the native projection.
Python tightens input to 16 MiB and atoms to 100,000. WASM uses 1,000,000 bytes
and 10,000 atoms. CLI retains its 64 MiB default external byte gate and applies
10,000 atoms in checked CJSON projection; the native 256 MiB cap also applies
if a caller raises the CLI external gate. No candidate is emitted on failure.

## Development comparison

`validation/cjson_checked/cases.json` contains original MIT synthetic documents
with declared source-indexed atomic numbers, formal charges, isotopes and bond
endpoints/orders. Six common documents × LF/CRLF/compact give 18 inputs; nine
compatibility controls give 27; 35 malformed/outside-profile documents give
105. Compatibility controls include optional omissions, isotopes, integral
floating bond storage and local aromatic aliases. They are scored separately
from common inputs. This is an exposed development packet, not a holdout or an
upstream issue reproduction. Compaction strips JSON whitespace outside strings
without parsing away duplicate keys.

The final paired route is **CJSON → CML**, using the same explicit input/output
formats for all engines. An independent standard-library XML observer reads
source-ordered atom IDs/elements, formal charges/isotopes and indexed bond
endpoints/orders. Undirected endpoint orientation is normalized, but atom order
is not fitted using values. Numeric output order tokens (including zero) are
observed without declaring them chemically valid. Unknown symbols/tokens or
unreadable output remain explicit observer failures. Coordinates, hydrogen
counts, stereo, names and other metadata are outside this scorer.

| Engine | Common selected graph | Compatibility selected graph | Malformed/outside-profile: no output |
| --- | ---: | ---: | ---: |
| Baseline Python | 18/18 | 27/27 | 39/105 |
| Baseline WASM | 18/18 | 27/27 | 39/105 |
| Candidate Python | 18/18 | 27/27 | 105/105 |
| Candidate WASM | 18/18 | 27/27 | 105/105 |
| Candidate CLI | 18/18 | 27/27 | 105/105 |
| Open Babel 3.2.1 | 9/18 | 6/27 | 18/105 empty outputs, exit zero |

CLI has no attested baseline measurement in this packet. Its candidate is freshly
built and tested against the same oracles. Candidate Python/WASM maintain their
previous positive preservation results; this implementation change improves
invalid-input handling. A refusal does not count as graph preservation.

Open Babel `-icjson -ocml` emits parseable XML for 87/105 malformed/outside-profile
inputs and no output for 18/105. Empty output with exit zero is recorded
separately from exception/nonzero-exit refusal; stderr diagnostics are archived.
These controls include typed-field errors and inputs outside the local profile;
the counts are not an Open Babel-wide CJSON conformance score.

For common inputs, the nine mismatches are the charged controls: the measured
Open Babel route emits zero formal charges. In compatibility controls,
integral floating order 2.0 and aromatic 1.5 emit order zero; source isotope
extensions and order 4 have different interpretations/retention. Optional
geometry/bond controls can produce no output. These are measured workflow
observations. Some nonsensical numeric bond orders emitted from invalid input
vary between runs; raw results are retained without treating those values as
stable expected behavior. The positive graph results predate this guard change;
no new general superiority claim is inferred from the packet.

## Baselines and development revisions

Baseline Python package source equals main `b9412a06`; the native extension and
WASM module/binary match the `c713b21a` artifact attestation. Relevant CJSON/native
and cross-format common paths are checked against candidate base `9b085d59`.
The c713b21a WASM has an unrelated unmerged CJSON self-retention branch, which
the final CJSON → CML route does not execute. Previous unmerged CML/SDF/source-
retention work is not a dependency of this candidate.

The initial CJSON → CJSON pilot did execute that unmerged WASM branch, so it
cannot serve as current-main WASM before/after evidence. The final paired route
was switched to CML. The first CML observer also classified numeric order zero
as unrecognized and inspected typed graph fields before classifying outputs
from invalid input. The final observer distinguishes well-formed XML emission
from graph preservation and recognizes finite numeric order tokens without
asserting chemical validity. The product fixtures did not change. Both pilot
runners and complete baseline/candidate outputs remain in
`cjson-checked-development-pilots-2026-10-22.json.gz`.

All final outputs, commands, source/binary identities and terminal logs are
indexed in `validation/results/cjson-checked-validation-2026-10-22.json`.

## Verification and reproduction

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cjson_checked.py \
  --wasm-module /path/to/chematic_wasm.js --cli /path/to/chematic \
  --obabel /path/to/obabel --output /tmp/cjson-checked.json
```

`--baseline` records old Python/WASM without requiring candidate success and
without a CLI baseline. Subprocess time and both output pipes are bounded.
Native tests cover every fixture/layout, legacy compatibility, exact immutable
retention, numeric coordinate bits and input/tree/atom/bond limits. Python and
real Node/WASM exercise every invalid layout across seven common output formats.
CLI integration tests verify every invalid document across the same seven
outputs has nonzero exit and empty stdout; the runner checks all layouts.
WASM tests exact 10,000/10,001 atoms and the byte cap. Observer negative controls
cover compaction/duplicates, graph changes, output/refusal/resource failures and
namespaced CML references. CI discovers the native/Python/Node suites and runs
observer unit tests.

This establishes consistent strict input validation on a small local CJSON
projection profile. It does not establish full schema/chemistry validation,
metadata/coordinate/stereo preservation across arbitrary outputs, timing,
released-package behavior or toolkit-wide superiority.

References: [Open Chemistry CJSON anatomy](https://github.com/OpenChemistry/chemicaljson/blob/main/chemicaljson.md)
(atom/bond property array cardinality),
[Open Chemistry schema](https://github.com/OpenChemistry/chemicaljson/blob/main/cjson.schema).
The local optional defaults/isotope/aromatic extensions are declared separately
from the evolving external schema.
