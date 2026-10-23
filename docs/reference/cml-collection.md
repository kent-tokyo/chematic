# CML molecular collections

[Open Babel #2339](https://github.com/openbabel/openbabel/issues/2339), verified
open on 2026-10-17, reports malformed multi-molecule CML from Python output.
`chematic.write_cml_collection` adds complete XML framing for serialized
molecule blocks, and `read_cml_collection` returns every source-ordered molecule
without projecting the XML through a reduced molecular graph.

```python
from chematic import from_smiles, read_cml_collection, write_cml_collection

blocks = [from_smiles(s).to_cml() for s in ["CCO", "[NH4+]", "[13CH4]"]]
text = write_cml_collection(blocks)
records = read_cml_collection(text)
assert [r.molecule_id for r in records] == ["m1", "m2", "m3"]
assert [r.index for r in records] == [0, 1, 2]
assert all(r.source_text == text for r in records)
# Each record also exposes a standalone CML molecule document:
assert "hydrogenCount" in records[0].molecule_xml
```

The writer accepts serialized `str`/UTF-8 `bytes` blocks rather than `Mol`
objects. The existing `Mol.to_cml()` supplies a basic graph projection; wrapping
that output does not recover stereo, arbitrary properties or other information
that the native single-molecule writer omitted. Supply rich source CML when
those fields must survive. The collection layer does not interpret or rewrite
atom, bond, stereo, coordinate or property payloads.

## Framing and identities

Each writer input has one `molecule` root in the
`http://www.xml-cml.org/schema` namespace. The output has one XML 1.0 document
and one `cml` root containing the molecules in input order. It is not a
concatenation of XML declarations and standalone documents. No file or stream
is opened, and a malformed later block fails the whole call before returning
any output.

Existing molecule IDs are retained and must be unique. Missing IDs receive
unused `m1`, `m2`, … values; underscores avoid collisions with supplied IDs.
No existing ID is renamed. The [CML molecular convention](https://www.xml-cml.org/convention/molecular)
scopes atom and bond IDs within their eldest containing molecule. Consequently,
both sibling molecules can keep `a1`, and their local `atomRefs2`/`atomRefs4`
references remain unchanged. This reader rejects nested molecule hierarchies,
so these scopes remain explicit and independent.

Molecular XML subtrees retain their spelling, internal whitespace, comments,
CDATA, processing instructions, element/attribute order and numeric strings.
Writer changes are limited to the collection wrapper, missing molecule IDs and
a default-namespace shield when needed. A prefixed molecule with no default
namespace receives `xmlns=""`, so the new collection's default namespace does
not silently change an unqualified payload element into a CML element.
Per-block declarations and comments/PIs outside the molecule root are not
copied into the combined document.

## Reading and metadata

`read_cml_collection` accepts either one namespaced `molecule` or a namespaced
`cml` collection. The collection root contains only direct molecules plus
whitespace/comments/PIs; only optional `version="3"` is accepted as a root
attribute. Global CML properties/conventions, inherited XML language/base
attributes and nested molecular hierarchies are rejected instead of flattened.
This is a bounded XML framing profile, not full CML schema or chemistry validation.

Every frozen `CmlMoleculeRecord` provides:

- `source_text`: the exact complete decoded input, including XML preamble,
  collection-level comments/PIs and original LF/CRLF/CR line endings.
- `index`: the zero-based source molecule position.
- `molecule_id`: its existing ID, or `None` for a read input without one.
- `molecule_xml`: its molecular subtree, suitable as a standalone XML document.
  Namespace declarations inherited from the collection root are materialized on
  its opening tag. Child declarations override inherited bindings. Namespace
  aliases used in QName-valued metadata such as `dictRef` remain available.

XML byte positions come from Expat. Subtrees are sliced from the original UTF-8
bytes, so Unicode preambles do not shift record boundaries. The molecular
payload is never normalized through a serializer or chemical parser. For a
read-only consumer, the original complete document remains accessible even
when outside-root comments/PIs are absent from the standalone projection.

## Checked limits

Inputs are UTF-8 `str` or immutable `bytes`. XML must be version 1.0; explicit
encoding declarations may identify UTF-8 or ASCII. DTD/entity declarations are
rejected, and no external resource is loaded. Wrong argument/container types
raise `TypeError`; malformed XML, duplicate molecule IDs, unsupported framing
and size failures raise `ValueError`. Empty collections and a collection passed
as one writer block also raise. No bad or later molecule is silently skipped.

Limits are 1..1000 molecules, XML depth 64, 200000 elements, 16 MiB total UTF-8
input and 4 MiB per source/standalone molecular subtree. Materialized standalone
fragments are also limited to 16 MiB in total, preventing repeated inherited
namespace declarations from multiplying the returned data beyond that bound.
Writer input blocks and final output are each limited to 16 MiB; combined
output is revalidated against the total element/depth budget. These limits
bound framing, not a chemical atom count or schema-validation resource budget.

## Reproduction and validation

The current **Open Babel 3.2.1** Python binding reproduced #2339 using three
original SMILES, three `pybel.Outputfile.write` calls and an explicit `close()`.
Its output still contained multiple top-level `molecule` elements; standard XML
parsing failed with `junk after document element`. The synthetic source and
version are retained in `validation/results/cml-collection-2026-10-17/`.

The chematic regression run passed **240 tests**: 49 new collection cases,
29 existing CML hydrogen cases and 162 compatibility cases. The tests include
rich source payload, local atom IDs, namespace inheritance/overrides, prefixed
roots, precise numeric spelling, Unicode byte positions, comments, line endings,
invalid later records, explicit limits and bounded namespace materialization.

Indigo **1.46.0.0-gb36d38f52** and lxml **6.0.2** independently checked **32
collections / 272 molecular records**. Source-CML and collection-CML canonical
SMILES, atom counts and bond counts agreed. The source lane includes Indigo's
explicit tetrahedral/geometric stereo XML. Expanded element/attribute/text trees
also agreed (apart from generated root molecule IDs), and all 32 read/assemble
cycles retained the same graphs. These checks establish collection-layer
retention, not complete native CML stereo support.

Open Babel returned the expected 1/2/3/7/12 molecule counts from five generated
controls, but emitted a closing-`cml` XML diagnostic in every case despite exit
status 0. An Open Babel CLI writer/read-back control also produced an XML
diagnostic and returned only its first molecule; the new collection reader
returned all three records from that same document. Therefore **diagnostic-free
Open Babel interoperability is not established**. Validity is checked with
independent XML parsers and Indigo; Open Babel's stderr is retained separately
from its molecule counts and is not treated as clean success.

Reproduce using `validation/validate_cml_collection.py` with Indigo/lxml and
Open Babel installed. `validation/validate_cml_pybel_output.py` uses the installed
Open Babel Python binding (Python 3.14 on this host). These are validation tools,
not new project dependencies. The existing native extension was reused for
basic single-molecule fixtures/regressions; its relevant CML reader/writer sources
match the recorded candidate. No fresh native build, wheel publication, MDL/SDF
writer change or full schema/chemistry parity is claimed.
