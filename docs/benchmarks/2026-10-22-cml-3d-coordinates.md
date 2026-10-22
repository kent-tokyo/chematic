# CML Cartesian coordinates through common CJSON conversion

The Python/WASM common CML readers discarded `x3/y3/z3` before CJSON output.
A one-atom CML with coordinates `(0.123456789012345, -0.5, 2.125)` produced
CJSON with no coordinates. The existing named CML writer could already export
complete Cartesian coordinates, but there was no native reader returning them.

Rust now exposes `parse_cml_with_3d[_with_limits]`, returning `CmlWith3d`:
a source-ordered molecule and `Option<Vec<[f64; 3]>>` in Angstrom. Python exposes
`from_cml_with_3d(str_or_utf8_bytes) -> (Mol, coords_or_none)`. Coordinate rows
are Python lists. Complete source Cartesian coordinates also flow through
Python `convert_format(text, "cml", "cjson")` and WASM
`convert_common_format(text, "cml", "cjson")`.

```python
mol, coordinates = chematic.from_cml_with_3d(source)
converted = chematic.convert_format(source, "cml", "cjson")
```

No generation, scaling, axis changes or 2D-to-3D promotion occurs. When both 2D
and 3D attributes are present, the Cartesian triplet is selected. All atoms must
have finite x3/y3/z3 or all must omit them. Partial triplets, mixed presence,
malformed/nonfinite numbers and duplicate atom attributes refuse. Absence returns
`None`/no CJSON coordinate field; a present all-zero triplet remains present.
Signed zero and subnormal finite f64 values survive the tested conversion.
Python's explicit `coords` override remains available and takes precedence after
source validation.

## Profile and compatibility

This is a bounded, lightweight single-molecule reader, not a full XML/CML
schema or chemistry validator. It requires exactly one nonempty molecule and
unique nonempty atom IDs. The existing tag normalization handles compact,
LF/CRLF forms and comments/CDATA/processing instructions. Atom attributes must
fit the existing one-tag-per-line normalized profile. Namespaced element names,
array-form atom records, arbitrary XML entities, multiple/nested molecules and
packed `xyz3`, fractional coordinates and general query semantics are not
newly supported. Packed/fractional coordinate attributes refuse explicitly. The existing attribute
decoder and graph semantics remain applicable. In particular, numeric character
references in coordinate tokens are not newly decoded. Existing pending CML
character-reference/source-retention/precision PRs are not dependencies.

The graph's existing element, isotope, formal-charge and bond projections are
reused; graph-isomorphism, aromaticity-model equivalence, hydrogen/stereo
interpretation and opaque metadata are not independently verified here. CML
2D layout, names, source IDs and other format metadata are still omitted from
common CJSON output. This is interpreted coordinate projection, not raw source
retention. Existing `parse_cml`, `parse_cml_strict`, Python `from_cml` and WASM
`mol_from_cml` retain their graph/2D behavior. Other common routes are unchanged.

The CML → CJSON common route now rejects malformed XML, ambiguous IDs, multiple
molecules and incomplete Cartesian coordinates that the previous graph-only
bridge could accept. Compact XML previously could produce an empty graph;
this route normalizes tags before parsing. These are intentional acceptance
changes. Callers needing multiple molecules should split/select them explicitly.

Native defaults remain `CmlParseLimits`: 256 MiB input, 1 MiB physical line,
1,000,000 lines/atoms and 2,000,000 bonds/XML elements. Python tightens input
to 16 MiB and atoms to 100,000; WASM uses 1,000,000 bytes and 10,000 atoms.
Source byte/line checks precede structural scanning and normalization. Inserted
tag separators do not consume the source-byte budget; normalized line/element
limits still apply. Atom caps precede point accumulation, and atom/bond caps
precede graph construction. No candidate output is returned on failure.

## Exposed development comparison

`validation/cml_3d/cases.json` contains original MIT synthetic inputs and declared
source-order elements, indexed bond endpoints/orders and finite Cartesian f64
oracles. Eight complete 3D documents × LF/CRLF/compact give 24 inputs; two absent
geometry documents give six; fourteen malformed or outside-profile documents give
42. These are development fixtures, not a holdout or upstream issue reproduction.

The independent observer compares source-ordered atomic numbers and undirected
indexed bonds/orders. Coordinate equality compares IEEE-754 bits, including
signed zero. It does not fit a mapping using coordinates, infer bond chemistry,
compare full metadata or score a refusal as coordinate preservation. Original
complete stdout/stderr, exit/resource status and artifacts are archived.

| Route | Coordinate bits, complete 3D | Selected graph, complete 3D | Graph output without generated 3D, absent cases |
| --- | ---: | ---: | ---: |
| Baseline Python common | 0/24 | 16/24 | 4/6 |
| Baseline WASM common | 0/24 | 16/24 | 4/6 |
| Candidate Python common | 24/24 | 24/24 | 6/6 |
| Candidate WASM common | 24/24 | 24/24 | 6/6 |
| Open Babel 3.2.1 `-icml -ocjson` | 24/24 | 24/24 | 0/6 (no output) |

Open Babel already retains the tested complete Cartesian coordinates, including
tiny z, signed zero, subnormals and cases also carrying distinct 2D positions.
This change catches up with that coordinate behavior. Open Babel's CJSON writer
rejects absent/2D-only geometry with a diagnostic and no output; chematic can
emit a graph-only CJSON document. That is a workflow distinction, not evidence
that Open Babel manufactured coordinates on these absent cases.

| Route, malformed/outside-profile group | Valid JSON accepted outside profile | Nonempty invalid JSON | No output |
| --- | ---: | ---: | ---: |
| Baseline Python/WASM, each | 40/42 | 0/42 | 2/42 |
| Candidate Python/WASM, each | 0/42 | 0/42 | 42/42 |
| Open Babel 3.2.1 | 19/42 | 9/42 | 14/42 |

These counts concern the declared single-molecule/Cartesian profile. Multiple
molecules are legal in broader CML workflows and are outside this API's boundary;
refusing them is not a claim that the document is invalid CML. Open Babel emits
zero-filled z for partial triplets and zero x for malformed numeric x3 in these
controls. NaN, infinity and numeric overflow yield truncated, nonempty JSON in
nine controls. Duplicate IDs, missing IDs and dangling bonds can also yield
output; diagnostics and all variant-specific outcomes remain in the raw report.
No-output controls include XML errors, mixed presence and multiple molecules.
Open Babel returns exit zero on these runs even when stderr contains errors;
the report's `error_status` means an API exception or nonzero process exit only,
and does not imply absence of stderr diagnostics. No generation options are used.

Baseline Python package source matches main `b9412a06`; the native extension and
WASM artifacts match the prior `c713b21a` attestation. That WASM has an unrelated
CJSON self-route; its CML route matches main. The candidate is based on main
`b7817c64` with no earlier unmerged I/O changes. Relevant CML/native/common paths
are checked against these source identities. Full outputs, commands, binary and
source hashes, logs and summaries are indexed in
`validation/results/cml-3d-validation-2026-10-22.json`.

## Verification and reproduction

```sh
PYTHONPATH=/path/to/python-package python3 scripts/validation/check_cml_3d.py \
  --wasm-module /path/to/chematic_wasm.js --obabel /path/to/obabel \
  --output /tmp/cml-3d.json
```

`--baseline` records old bindings without requiring candidate success. Both pipes
and subprocess time are bounded. Native tests exercise all fixture layouts,
CJSON numeric read-back, named-writer read-back, missing geometry and inclusive
byte/line/atom/bond bounds, plus normalized compact line/element limits. Python
checks both public paths, UTF-8 bytes, explicit overrides and unchanged graph
charge/isotope projection. Real Node/WASM checks the same packet and exact
10,000/10,001 atom boundaries. Observer negative controls cover signed zero,
permutation, missing/zero geometry, invalid JSON, refusal classification and
graph changes. CI discovers native/Python/Node suites and runs observer tests.

The measured scope is Cartesian coordinate projection and explicit failure on
this local profile. It does not establish toolkit-wide superiority, editable
CML metadata/stereo fidelity, timing or release availability.

References: [CML molecular convention](https://www.xml-cml.org/convention/molecular)
(separate 2D/3D geometry and complete x3/y3/z3),
[Open Babel CML format](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html).
