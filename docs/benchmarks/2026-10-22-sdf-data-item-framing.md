# SDF named field preservation

The shared SD data-item parser now reads named headers with preceding registry
information and treats header-shaped lines inside values as literal text.
Previously it silently lost registry-prefixed fields and split values whenever
a line started with a recognizable `> <name>` header.

An empty line ends a data item. This follows the
[MDL CTfile Formats, October 2003, chapter 6, pages 39–40](https://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf).
Header metadata may occur before or after an angle-bracketed field name.
Data values can span multiple lines. The project fixtures are original short
strings and a two-atom methanol structure; no specification figure or upstream
attachment was copied.

## Paired source result

The packet contains 12 original data cases, each combined with V2000/V3000 CTAB
and LF/CRLF line endings: **48 correlated development variants**, not 48
independent molecules or a sealed holdout. Cases cover plain/registry headers,
header-shaped first/middle/self-referencing value lines, and significant spaces
or tabs. A following `CHECK` field tests that a real empty separator still works.

The measured native workflow uses source-built Python `ForwardSDMolSupplier`
and `SDWriter(compute2d=False)`, retaining the input CTAB version. Open Babel
3.2.1 runs `-isdf -osdf` (plus `-x3` for V3000). No charge calculation or other
conversion options are requested. The installed `-L sdf` inventory is archived.

An independent Python observer splits output data items on empty lines, parses
only their first line as a header, and compares field-name/value dictionaries
with authored expectations. Leading/trailing spaces and tabs are significant;
logical multiline values use LF. The native parsed properties are also checked
before writing. The observer separately checks the simple explicit C/O graph,
neutral charge and unspecified isotope representation.

| Outcome, 48 variants | Before | Candidate | Open Babel 3.2.1 |
|---|---:|---:|---:|
| Specified field names and exact values retained | 20 | 48 | 36 |
| Field loss/change | 28 | 0 | 12 |
| Process, graph or observer failure | 0 | 0 | 0 |

The native baseline misses four registry-header cases and changes three
literal-header value cases, across four format/line-ending combinations each.
The Open Babel differences are the three padded-text/multiline/tab cases; its
conversion removes their surrounding whitespace. This local profile values
exact text preservation. It does not establish chemical significance for a
whitespace change or toolkit-wide superiority. Every baseline/candidate process
outcome, source, rewritten output, stdout/stderr and tool identity is retained.

## Reader behavior and limits

The repair is in the shared Rust field parser used by string record readers,
diagnostic and fast file readers, Python native/compatibility forward suppliers
and both WASM record JSON APIs. It does not add an API or change the writers.
Adding registry-header recognition without fixing value framing would introduce
more false headers inside values; the packet tests the two rules together.

**Compatibility change:** a header-shaped line without an intervening empty
separator now remains in the current value. The previous heuristic attempted
to recover missing separators by starting a new field. A resource-limit test
that intended two fields but omitted that separator was corrected by inserting
the empty line; its field assertions and limit checks remain. The original
failed run is archived separately from the passing final run.

The existing named-field dictionary projection remains: a genuinely repeated
key keeps its last value. Property order, registry/header annotations and
field-number-only `DTn` headers are not retained by this projection. This is not
complete SD data-item/document preservation. Coordinates, general CTAB chemistry,
stereo, property filtering and underscore-prefixed writer rules are outside the
comparison. Existing SDF input/record/line budgets remain in effect; the paired
runner also bounds process streams. No timing or generic streaming-bound claim.

## Validation

- Full MOL Rust unit/integration/doc suite: **1,478 passed**, zero failed/ignored.
  The new integration packet checks all 48 variants through string, stream and
  fast readers with 1/17/8192-byte buffers, plus both SDF writers/read-back.
  Empty values and the legacy last-value repeated-key projection are checked.
- Source-built Python: **287 passed**, including 120 new native/compatibility
  reader and writer/read-back cases, plus existing forward/I/O/compression and
  selection suites.
- Source-built debug Node/WASM: **48 variants passed** through both whole-input
  and batch record JSON APIs.
- Observer/bounded-runner tests: **38 passed**. MOL all-target Clippy with
  warnings denied, formatting, refactor/script-layout/product-name checks and
  strict MkDocs passed.

Baseline source is `a83ed2653037a7b17033e03ee22501a1ddda8bb9`; both baseline
and candidate Python extensions were built locally, with separate immutable
runtime copies. This is source evidence, not validation of a published package.
It is independent of pending CML PRs #1126–#1128.

## Reproduction

Build a Python extension from the desired source revision, then:

```sh
PYTHONPATH=/path/to/source-built-python python3 -m scripts.validation.check_sdf_data_item_framing \
  --obabel /path/to/obabel --output /tmp/sdf-fields-new.json
cargo test --locked -p chematic-mol --test sdf_data_item_framing
python3 -m unittest scripts.tests.test_sdf_data_item_framing
```

The runner requires Open Babel 3.2.1 and refuses to overwrite a result. Source
and artifact hashes, compiler/runtime versions, the installed format inventory
and compressed logs are listed in the validation index:

- [Fixture manifest](https://github.com/kent-tokyo/chematic/blob/fix/sdf-data-item-framing-20261022/validation/sdf_data_item_framing/cases.json)
- [Baseline](https://github.com/kent-tokyo/chematic/blob/fix/sdf-data-item-framing-20261022/validation/results/sdf-data-item-framing-baseline-2026-10-22.json)
- [Candidate](https://github.com/kent-tokyo/chematic/blob/fix/sdf-data-item-framing-20261022/validation/results/sdf-data-item-framing-candidate-2026-10-22.json)
- [Validation index](https://github.com/kent-tokyo/chematic/blob/fix/sdf-data-item-framing-20261022/validation/results/sdf-data-item-framing-validation-2026-10-22.json)
