# RDFile records with explicit graph projection

Open Babel [#2681](https://github.com/openbabel/openbabel/issues/2681) requests
RDFile support and batch extraction of molecular/reaction SMILES. The locally
pinned Open Babel 3.2.1 has no registered RDF/RDFile format and rejects
`-irdf`. This candidate adds native Rust/Python record reading and checked
writing, reusing chematic's MOL and RXN parsers.

```python
from pathlib import Path
import chematic

text = Path("reactions.rdf").read_text()
document = chematic.parse_rdfile(text)
smiles = chematic.rdfile_to_smiles(text)
for record, value in zip(document["records"], smiles):
    print(record["kind"], value)  # None denotes an unresolved registry reference

# Edits affect serialized records; original source text remains available.
document["records"][0]["data"].append({"name": "review", "value": "checked"})
output = chematic.write_rdfile(document)
original = document["source_text"]
```

The typed document contains `date_comment`, ordered `records`, and exact
`source_text`. Each record has:

- `kind`: `molecule` or `reaction`;
- optional `internal_registry` or `external_registry`;
- `structure`: exact embedded MOL/RXN text, or `None` for a registry-only record;
- `data`: an ordered list of `{name, value}` entries, retaining duplicate names.

Internal registry numbers remain decimal strings so leading zeros survive.
External registry values remain text; neither is looked up or resolved.
Data values remain text, including significant internal spaces, blank lines
and repeated fields. No database schema is available to interpret field types
or map arbitrary names to reaction conditions.

`rdfile_to_smiles` is an explicit graph projection, not a lossless replacement
for the RDFile document. It returns one entry per record, preserving reference
slots as `None`. Coordinates, fields, registry identifiers and source number
spelling remain in `parse_rdfile` output. Molecular/reaction graph semantics
follow the existing MOL/RXN parsers; no new bond perception is introduced.

In Rust, `parse_rdfile` / `parse_rdfile_with_limits` return an editable
`RdfileDocument`. `RdfileRecord::smiles()` returns
`Result<Option<String>, RdfileError>`. `write_rdfile` uses current records,
validates serialized output by reparsing, and checks that fields and record
boundaries survived. It never replays a potentially stale `source_text`.
`RdfileError` identifies a one-based source line where available, or zero for
constructed-document validation.

## Format profile

The implementation was written independently from the RDFile chapter of the
[BIOVIA Chemistry 2026 CTfile specification](https://www.3ds.com/assets/invest/2026-01/biovia_chemistry_2026_ctfile_formats.pdf),
retrieved on 2026-10-15, and tested against
[Indigo's RDF reader/writer API](https://lifescience.opensource.epam.com/indigo/api/index.html).
This API handles chemistry RDFile, rather than Semantic Web RDF.

Supported input includes version 1, embedded molecule/reaction records,
registry-only records, optional single registry identifiers, ordered data,
multiline values and column-81 `+` data continuation. MOL payloads can be
V2000 or V3000; embedded RXN is V2000 with V2000 components. RXN counts may
include an explicit zero agent count. Embedded blocks are consumed structurally,
so a keyword-shaped molecule title does not accidentally start another record.

Headers and RD metadata use ASCII. Header arguments fit 80 columns. Data
folding joins a physical row with a `+` in column 81 to its continuation
without adding a newline; ordinary continuation rows preserve newlines.
Plain longer data rows are accepted as an explicit common interchange
extension, still subject to the physical-line limit. The original file and
embedded structures retain their complete supplied UTF-8 bytes and line endings.

Unsupported profiles fail explicitly: V3000 RXN, nonzero agent counts,
V3000 molecules inside a V2000 RXN, nested embedded `$DATUM $MFMT` structures,
both registry identifiers on one record, continued/long headers and non-ASCII
RD metadata. Unknown keywords, missing headers/terminators/data, malformed
counts and ambiguous delimiter injection also fail. Registry-valued data
remains literal text rather than being resolved.

An empty library and empty datum are supported. Blank lines outside structures
or data are rejected; trailing blank lines inside a datum remain data.
Only a missing final LF is normalized when writing a structure block.
Canonical writing uses LF for new RD headers/data while keeping embedded block
bytes unchanged, including their original CRLF if present.

## Long-data interoperability

The default writer uses specification folding. Indigo **1.46.0** was measured
retaining the column-81 `+` as literal text instead of joining the rows.
Consequently, the default folded output is not claimed to preserve long data
when read by this pinned Indigo version.

Use the explicit unwrapped profile for this consumer:

```python
indigo_output = chematic.write_rdfile(document, wrap_data=False)
```

In Rust, use `write_rdfile_with_options(document,
RdfileWriteOptions { wrap_data: false })`. This emits plain long data lines
as the common interchange extension; it does not change the default profile.
Round-parse validation still rejects any datum whose literal position would
be confused with a continuation marker or record keyword. No data is silently
trimmed to make it fit either profile.

## Resource and edit boundaries

Default reading limits are 64 MiB input, 1 MiB per physical line, 10,000 records,
1,000 fields per record and 1,000 components per reaction. All are configurable
through Rust `RdfileParseLimits` and corresponding Python keyword arguments.
The cursor scans the source without collecting an unbounded physical-line
array. Record/field vector reservation is fallible.

The writer uses default limits, checks record/field counts and a conservative
serialized-size estimate before reserving its output buffer, then reparses the
result. It refuses malformed public Rust data, unknown Python document keys,
missing registry-only identifiers, unsupported text, conflicting registries,
newlines in headers and changes to record boundaries. A source document's
untouched original text remains available for exact replay by the caller.

This is an in-memory typed batch API. It adds no CLI format registration,
streaming file iterator, registry service, query-schema interpretation or
conversion from arbitrary rich reaction documents into RD fields.

## Validation

The original small fixture and manifest are in
`validation/openbabel_issues/rdfile/`. Three embedded records cover methanol,
carbon-13 methanol and atom-mapped oxidation; two unresolved references verify
slot retention. Duplicate names and multiline data are checked against authored
ordered vectors, since Indigo's property abstraction does not retain duplicates.
The manifest records provenance, MIT OR Apache-2.0 licensing, minimization,
pinned reference settings, expected values and the normalized payload hash.

Indigo is pinned to
`1.46.0.0-gb36d38f52-aarch64-darwin-appleclang-17000013`, RDKit to **2026.09.1**,
and Open Babel to **3.2.1**. Forty-eight dataset profiles vary record order,
LF/CRLF and empty/multiline/long/spaced data. Indigo reads 144 native records
using `wrap_data=False`; chematic reads 144 Indigo-written records. All
**288 cross-read records** match graph meaning, atom maps and data, with zero
mismatches. Reaction comparisons retain component roles; graph checks use
RDKit InChI plus mapped atom/bond counters, not lexical SMILES equality.
Forty-eight default-folded native round trips also retain all records/data.
The separate Indigo folding diagnostic is recorded as an interoperability
limitation, not counted as a passing default-output cross-read.

Local checks passed 753 Rust library tests, 18 new RDFile contracts and
80 Python tests across RDFile and previous TURBOMOLE functionality. These cover
V3000 molecule payloads, structural framing,
source/CRLF retention, duplicates, empty/reference records, mutations, limits,
continuation boundaries, malformed input and injection refusal. Strict Clippy,
formatting, static binding inventory, product naming, the refactor budget and
strict MkDocs are also checked; no allowances were added.

Regenerate the independent evidence with the candidate extension on `PYTHONPATH`:

```sh
python scripts/validate_rdfile_records.py --obabel /path/to/obabel
```

This is a local, uncommitted, unpublished source candidate in
`/private/tmp/chematic-rdfile`, preserving the previous candidates through
correlated stereo-group alternatives. The portable feature patch is
`validation/openbabel_issues/rdfile-records-2026-10-15.patch` and its evidence is
`validation/results/openbabel-rdfile-records-2026-10-15.json`. The patch includes
source, bindings, tests, fixtures, documentation and CI/inventory updates.
