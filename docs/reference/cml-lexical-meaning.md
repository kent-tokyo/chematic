# CML lexical layout and molecular meaning

The lightweight CML graph/2D readers now interpret compact tags and multiline
attributes consistently. Literal markup inside comments, CDATA sections and
processing instructions cannot introduce atoms, bonds or replacement attributes.

For example, splitting `elementType="N"` onto a new line previously allowed a
nitrogen atom to be read as carbon. A multiline comment containing an example
`<atom>` could add a carbon atom to the returned molecule. These cases now return
the authored graph through Rust, Python, WASM and common-format conversion.

## Reading boundary

Before reading graph fields, the lightweight scanner creates logical tag lines:

- Quoted `>` characters do not terminate tags.
- Literal CR, LF and tab inside tags become spaces; CRLF becomes one space.
  Character references are decoded afterward, so referenced characters are not
  accidentally normalized a second time.
- Comments, CDATA and processing instructions are opaque to the graph reader.
  Their contents are excluded as complete units, including multiline and Unicode
  payloads.
- Raw input byte, physical line and line-byte limits are checked first. Logical
  tag-line, tag-line-byte, XML-element, atom and bond limits also apply after
  normalization. Inserted tag separators do not count as extra source bytes.

These lexical rules follow the relevant boundaries in
[XML 1.0](https://www.w3.org/TR/xml/), including
[attribute-value normalization](https://www.w3.org/TR/xml/#AVNormalize).
This change does not establish full XML or CML schema conformance. The default
reader retains its lenient missing/empty-molecule boundary. The strict reader
adds nesting and nonempty-molecule checks; it is not a namespace-aware schema
validator. The separate Cartesian CML profile already uses namespace-aware XML
parsing and serves as a graph control in this packet.

Graph projection and source preservation have separate contracts. Reading a
molecule ignores opaque XML. `convert_with_report` with CML input/output and
`preserve_source` retains the original bytes, including comments and layout,
and reports `candidate_graph_verified=false`. A byte copy is not counted as a
verified chemical projection.

## Bounded comparison

The original MIT development packet contains 16 well-formed XML documents:
three controls and thirteen layout/opaque-content reproductions, using two
authored graph families. ElementTree independently checks source fields against
the authored oracle. This is an XML well-formedness check, not CML schema or
physical-validity validation.

The following counts cover **only the thirteen reproduction cases**. All three
controls pass before and after in every measured route.

| Route | Baseline preserved | Candidate preserved |
|---|---:|---:|
| Python default, str and UTF-8 bytes (each) | 0/13 | 13/13 |
| Python strict, str and UTF-8 bytes (each) | 2/13 | 13/13 |
| Python common CML → SMILES | 0/13 | 13/13 |
| WASM default | 0/13 | 13/13 |
| WASM strict | 2/13 | 13/13 |
| WASM common CML → SMILES | 0/13 | 13/13 |
| Python namespace-aware 3D reader, graph control | 13/13 | 13/13 |
| Open Babel 3.2.1 CML → SMILES | 13/13 | 13/13 |

Baseline default/common routes return no output for three reproductions; strict
routes return no output for two. Other baseline failures return an output with
changed observed meaning. Refusals and missing outputs are retained per case;
they never count as preservation wins. Both baseline and candidate source-copy
routes retain all 16 inputs exactly in Python and WASM, measured separately.

The baseline Python/WASM artifacts were built from
`342d1d5f9ad09aa6ade4e3079a902cba856209e2`. The affected CML reader sources are
unchanged between that commit and the candidate's main base `fcbf6d0a`.
Candidate records include production-source and runtime hashes. The recorded
HEAD is the pre-commit checkout; source hashes identify the dirty candidate.

Open Babel receives the original fixture bytes, using
`obabel -icml -osmi -xn -xa`. Its documented
[SMILES write options](https://openbabel.org/docs/FileFormats/SMILES_format.html)
omit molecule names and retain available atom classes. The emitted SMILES from
every engine is observed using RDKit 2025.09.3, without sanitization or H-node
removal. RDKit is an output observer, not the CML reader or chemical oracle.
The observer checks atomic numbers, formal charges, isotope labels, attached
hydrogen counts, atom maps (all zero in this packet) and bond orders. Unique
elements provide an unambiguous atom correspondence for these two small graphs.
There are no explicit H nodes in the packet.

Native regression tests also check input-indexed 2D coordinates exactly. The
SMILES comparison does not measure coordinates. CLI tests observe common
CML → MolJSON graph/H fields. No result here establishes stereo, aromaticity,
queries, general namespace/schema support, throughput or broad-format
superiority. On this packet the candidate reaches Open Babel's observed result.

## Reproduction and evidence

Fixture, full before/after outputs, errors and binary identities:

- [packet and reports](https://github.com/kent-tokyo/chematic/tree/main/validation/cml_lexical_meaning)
- [measurement runner](https://github.com/kent-tokyo/chematic/blob/main/scripts/validation/check_cml_lexical_meaning.py)

Build an isolated Python package and Node-target WASM module from the selected
checkout, then run:

```sh
PYTHONPATH=/absolute/python-package-root python3 scripts/validation/check_cml_lexical_meaning.py \
  --wasm-module /absolute/pkg-node/chematic_wasm.js --output /tmp/cml-lexical.json
python3 -m unittest scripts.tests.test_cml_lexical_meaning -q
cargo test -p chematic-mol --test cml_lexical_meaning
cargo test -p chematic-cli --test cml_lexical_meaning
```

The comparison requires local RDKit and Open Babel. Observer negative controls
run without RDKit in CI, alongside native, Python and Node regression suites.
