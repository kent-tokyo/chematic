# Checked CML source retention

`convert_with_report` now supports CML → CML without projecting a document
into an editable molecular graph. It retains the exact UTF-8 source and reports
`source_text: preserved`, `source_semantics: opaque_retained`,
`verification_scope: exact_source_bytes` and `candidate_graph_verified: false`.
The existing CLI, Python and WASM report APIs delegate to the same Rust code.

```python
report = chematic.convert_with_report(
    source, "cml", "cml", policy="preserve_source"
)
assert report["output"] == source
assert report["candidate_graph_verified"] is False
```

```sh
chematic convert-with-report --input-format cml --output-format cml \
  --policy preserve_source --input source.cml --output report.json
```

The CLI file contains the JSON envelope; molecule text is its `output` field.
This route is useful when comments, property order, coordinate dimensionality
and opaque vendor payload must survive a same-format pass. It does not edit,
interpret or verify that payload. A retained stereo element, for example, is
not evidence of correct stereo interpretation. Generic graph conversion APIs
retain their existing behavior; use the checked route explicitly.

## Measured workflow

Six original two-atom methanol documents exercise: a dyadic 2D control,
coexisting 2D/3D coordinate sets, precise 3D positions, repeated properties and
unknown namespace-qualified payload, a prefixed molecule root, and comments/PIs
with CRLF. These are synthetic development inputs, not optimized geometries or
an independent/sealed holdout. Fixtures include their authorship, license,
minimization, pinned tools and normalized hashes.

The independent Python 3.13.6 ElementTree observer checks the explicit C/O
single-bond graph, neutral charge/unspecified isotope, source atom attributes,
coordinate-component presence and numbers, and original expanded
property/extension subtrees in order. It ignores prefix spelling and
indentation. Coordinate tolerance is 1e-12 Å absolute, relative tolerance zero,
with no rescaling or alignment. Added payload outside original source subtrees
is not treated as source loss. Root molecule ID/title, inferred hydrogens,
schema/chemical validity and stereo interpretation are outside this profile.
Original byte equality is a separate diagnostic, including comments and PIs.

| Specified source fields retained | Inputs |
|---|---:|
| chematic checked CML → CML, each of three policies/bindings | 6/6 |
| Open Babel 3.2.1, `-icml -ocml` defaults | 3/6 |
| Open Babel 3.2.1, plus `-xp -xm` | 3/6 |

Both Open Babel profiles retain the declared fields on the control, prefixed
root and CRLF inputs. The dual-coordinate input loses its 2D set; precise 3D
coordinates exceed the declared tolerance; the rich document loses its unknown
atom attribute and original property/extension payload. All 12 Open Babel runs
complete with exit zero. These are field differences, not refusals or crashes.
No general toolkit superiority follows from the six documents or this
source-retention workflow.

All 18 case/policy combinations return identical complete JSON reports across
source-built CLI, Python and Node/WASM. Each report retains the source exactly.
The paired record contains all 66 process outcomes: 54 binding/policy runs and
12 Open Babel runs, with commands, sources, stdout/stderr and resource status.
The runner bounds process streams; it is not a general performance benchmark.

The first observer revision incorrectly counted an additional generated
`metadataList` as loss on inputs with no source payload. Revision 2 scores the
ordered source subtrees only. A positive control now checks that added metadata
does not cause a loss, alongside negative controls for dropped coordinate sets,
unknown attributes, reordered repeated properties, changed bonds and false
graph-verification claims. Both raw runs and the original runner source are
retained; only revision 2 is used in the table.

## Validation and limits

The source profile requires one namespaced `molecule` root, optional XML 1.0
declaration, UTF-8/ASCII text and no DTD. Multiple roots, collections and nested
CML molecule hierarchies are errors. The original text is never reserialized.
Empty molecular documents and uninterpretable chemistry are retained opaquely;
numeric chemistry, atom/bond references and CML schema validity are not checked.
All three policies return the same exact source with their selected policy in
the report. No atom/bond correspondence or verified graph is fabricated.

roxmltree **0.20.0** is now a direct pure-Rust dependency (already present in the
workspace lockfile). DTD parsing is disabled, with a prolog check also rejecting
empty/external DTD declarations while permitting literal markup in comments,
PIs and CDATA. No external resources are fetched. Input bytes, physical line
bytes/count and CML atom/bond element counts use the existing shared limits.
XML parsing has a 200,000-node budget; element depth is checked after parsing
against 64. Array-form members and opaque payload are not interpreted as graph
atoms. Default CLI/Python input and output limits are 16 MiB; WASM retains its
1,000,000-byte input cap and 16 MiB output cap. Escaped report JSON has its own size
check. These are parser/profile bounds, not full CML conformance claims.

Local validation passed:

- 1,543 Rust MOL/CLI unit, integration and doc tests, including boundary and
  malformed-input tests and CLI output-file preservation on a later XML error.
- 139 source-built Python tests, including 25 new retention/error cases.
- 25 new source-built WASM retention/error cases and existing conversion-report
  contracts; 18 complete-report binding-parity conditions.
- 39 observer/bounded-runner tests; all-target MOL/CLI Clippy with warnings
  denied; formatting, stub/API docs, refactor, script-layout and naming checks.

This is source evidence, not a new wheel/npm/crates.io release. The implementation
is based on main `d288f6027d71f21130dd119f7cc7df8d58ab8daf` and is independent of
the pending character-reference and coordinate-writer PRs #1126/#1127.

## Reproduction

Build the CLI, Python extension and Node/WASM module from the same source, then:

```sh
PYTHONPATH=/path/to/source-built-python python3 -m scripts.validation.check_cml_source_conversion \
  --cli /path/to/chematic --wasm-module /path/to/chematic_wasm.js \
  --obabel /path/to/obabel --output /tmp/cml-source-new.json
python3 -m unittest scripts.tests.test_cml_source_conversion
cargo test --locked -p chematic-mol -p chematic-cli
```

The runner requires Open Babel 3.2.1 and refuses to overwrite a report.
[Official CML options](https://openbabel.org/docs/FileFormats/Chemical_Markup_Language.html)
describe the `-xp`/`-xm` settings. The validation index records the exact source,
executables, Open Babel library/plugin, runtime versions and compressed logs:

- [Fixture manifest](https://github.com/kent-tokyo/chematic/blob/feat/cml-source-conversion-20261022/validation/cml_source_preservation/cases.json)
- [Paired result](https://github.com/kent-tokyo/chematic/blob/feat/cml-source-conversion-20261022/validation/results/cml-source-conversion-paired-2026-10-22.json)
- [Validation index](https://github.com/kent-tokyo/chematic/blob/feat/cml-source-conversion-20261022/validation/results/cml-source-conversion-validation-2026-10-22.json)
