# Checked molecular ChemDoodle JSON

This isolated, unreleased candidate adds a molecular ChemDoodle JSON reader
and graph writer, motivated by
[Open Babel #2763](https://github.com/openbabel/openbabel/issues/2763).
The issue was verified OPEN on 2026-10-15 (updated 2025-01-10T08:55:05Z).
Its reported Open Babel 3.1.0 behavior loses inferred hydrogen on default
carbon and can write numeric element labels. Those examples **do not
reproduce on installed Open Babel 3.2.1**: ordinary methane reads as `C`, and
the carbonyl/alcohol writer emits oxygen symbols correctly.

Two current gaps were independently measured on Open Babel 3.2.1: explicit
`h: 0` can be replaced by inferred H, and a multi-molecule CDJSON content
document is read as its first molecule only. This candidate distinguishes
omitted/automatic H from explicit zero and returns every molecule.
It fills a missing chematic format API; it does not claim an upstream fix.

## Python API

```python
import chematic

ordinary = chematic.parse_cdjson('{"a":[{"x":0,"y":0}]}')
assert ordinary["molecules"][0]["mol"].formula == "CH4"

zero_h = chematic.parse_cdjson('{"a":[{"x":0,"y":0,"h":0}]}')
assert zero_h["molecules"][0]["mol"].formula == "C"

mol = chematic.from_smiles("C(=O)CO")
coordinates = [[0, 0, 0], [20, 0, 0], [-20, 0, 0], [-40, 0, 0]]
text = chematic.write_cdjson_graph(mol, coordinates)
assert chematic.parse_cdjson(text)["molecules"][0]["mol"].formula == mol.formula
```

`parse_cdjson(text, *, max_input_bytes=None, max_molecules=None,
max_atoms=None, max_bonds=None)` returns:

| Field | Meaning |
| --- | --- |
| `source_text` | Exact original UTF-8 source, including whitespace and CRLF |
| `molecules` | Every molecule in original content order |
| each `mol` | Native molecular graph |
| each `coordinates` | One `[x,y,z]` array per original atom; omitted z becomes zero |
| each `id` | Optional source molecule ID |
| each `atom_ids`, `bond_ids` | Optional IDs aligned with original atom/bond arrays |

Coordinates remain **unscaled drawing coordinates**. No Angstrom conversion,
axis reversal, coordinate generation or stereochemistry perception is performed.
Source IDs are separate from atom-mapping numbers. Nonempty IDs must be unique
across the content. Exact same-format preservation uses `source_text`; graph
writing generates a new document without source IDs or unrelated metadata.

`write_cdjson_graph(mol, coordinates)` writes one molecular graph under a
content `m` array. It emits element symbols, charge, optional isotope mass,
all three coordinates, ordinary bond orders and explicit native H counts.
It requires a Python list with one triple per atom. Input graph and coordinate
counts are checked before Python numeric extraction. Output coordinates use
round-trip floating-point serialization.

Rust exports `CdjsonDocument`, `CdjsonMolecule`, `CdjsonParseLimits`,
`CdjsonError`, `parse_cdjson`, `parse_cdjson_with_limits` and
`write_cdjson_graph` from `chematic_mol`. The immutable document provides
`molecules()` and `source_text()`; graph projections carry native molecules,
coordinate arrays and source IDs.

## Supported molecular profile

The profile follows the default-field and molecular-index rules in the
[official ChemDoodle JSON specification](https://web.chemdoodle.com/docs/chemdoodle-json-format/).
ChemDoodle JSON (`cdjson`) is distinct from ChemicalJSON (`cjson`).

| Input | Contract |
| --- | --- |
| Root | One molecule with `a`, optional `b`/`i`; or content with `m` and optional empty `s` |
| Atom label `l` | Omitted means C; explicit labels must be recognized case-sensitive element symbols |
| Coordinates | Finite x/y required; finite z optional, default zero |
| Charge `c` | Integer, default 0, restricted to native i8 range |
| Mass `m` | Omitted or -1 means unspecified; positive u16 isotope mass supported |
| Implicit H override `h` | Omitted/-1 uses native valence inference; 0–255 retains the explicit count, including zero |
| Radicals/lone pairs `r`/`p` | Omitted or zero accepted; nonzero values rejected |
| Bonds | Distinct in-range integer endpoints, no duplicate pair; omitted order means 1; orders 1/2/3 supported |
| Bond stereo `s` | Omitted or `none` accepted |
| IDs `i` | Optional nonempty strings; retained without interpretation |

Unknown fields, shapes, query objects, enhanced stereo, nonzero radical/lone-pair
fields, directional/ambiguous bond stereo, aromatic/half/zero bond orders and
unrecognized labels return explicit errors. The graph writer also rejects
native aromatic/query/chiral/mapped atoms, enhanced stereo groups, directional
bonds and self bonds. An explicitly Kekulé graph with ordinary bonds can be
written. The profile does not assign radical electron counts from an explicit
H override, sanitize unusual valence, or implement the full ChemDoodle drawing
model. Unsupported chemistry is not silently omitted.

The existing strict ChemicalJSON duplicate-key reader is reused internally.
Duplicate object keys are rejected recursively, including escaped equivalents.
The CJSON public contract is unchanged.

## Limits

| Resource | Hard cap |
| --- | ---: |
| Input/output JSON bytes | 16 MiB |
| Molecules per content document | 1024 |
| Aggregate atoms | 65536 |
| Aggregate bonds | 131072 |
| JSON nesting depth | 32 |
| Items per JSON array | 131072 |
| Bytes per decoded JSON string | 65536 |

Caller budgets can only tighten the first four parse limits. Input bytes are
checked before JSON decoding. Structural JSON checks and aggregate graph
counts run before molecular projection. Empty content and empty molecular
arrays are supported. Rust returns typed invalid/unsupported/resource errors;
Python maps them to `ValueError`, with unsigned conversion errors for negative
budgets. This bounds input/work; it does not introduce a general fallible
allocator contract or change other format parsers.

## Validation

- Rust: **765 passed**, comprising 753 existing library tests and 12 new
  integration tests; zero failures or ignored tests. The new suite covers H
  default versus zero, IDs, source text, coordinates, bonds, malformed values,
  unsupported fields, exact atom/molecule caps, first rejected sizes and writer
  rejection paths.
- Python: **54 passed**, comprising 31 new tests and 23 existing strict CJSON
  preservation tests.
- RDKit **2026.09.1** reference reconstruction: **544 exact molecular cases**
  and **544 native graph round trips**. Seventeen original synthetic recipes
  receive 32 atom/bond order, labeling, ID and coordinate variants. Element
  symbols, charge, isotope, H count, formula, endpoints, IDs and full-precision
  coordinates are checked separately. No unexpected mismatches occurred.
- Open Babel **3.2.1**, each input molecule read separately: **448 formula
  matches and 96 known zero-H differences** on both source reads and native
  writer cross-reads (1088 individual CLI conversions). The differing recipes
  are isolated zero-H carbon, isolated zero-H oxygen and bonded zero-H carbon.
  These differences are recorded explicitly, not counted as parity successes.
- Open Babel multi-molecule control: the same content containing **544
  molecules** returns **1** from Open Babel and **544** from chematic.
- Open Babel writer import: **12 ordinary synthetic SMILES controls** exported
  with `-ocdjson --gen2D` all import into chematic with their expected formulas.

The original fixtures, per-record provenance/license/minimization, expected
formulas and normalized payload SHA-256 are in
`validation/openbabel_issues/cdjson/fixtures.json`. No upstream attachments,
ChemDoodle code or external molecule corpus are copied. Runtime oracle versions,
normalized variant packet hash, every known external difference, source-lane
hashes and checks are saved in
`validation/results/openbabel-cdjson-2026-10-15.json`.
The reproducible oracle runner is `scripts/validation/validate_cdjson.py`.

```sh
cargo test -p chematic-mol --lib --test cdjson_contract --offline --quiet
cargo clippy -p chematic-mol -p chematic-py --all-targets --offline --no-deps -- -D warnings
python3 -m pytest crates/chematic-py/tests/test_cdjson.py crates/chematic-py/tests/test_cjson_preservation.py -q
python3 scripts/validation/validate_cdjson.py
```

Formatting, binding inventory, refactor budget, product naming, stub syntax and
strict MkDocs checks also pass. These are local debug-extension source checks,
without native InChI; package publication, WASM and full ChemDoodle application
interoperability are outside this candidate.

The isolated source is `/private/tmp/chematic-cdjson`, inherited from
`/private/tmp/chematic-fold-density`. The incremental patch is
`validation/openbabel_issues/cdjson-2026-10-15.patch`, applicable to that base.
Only this document, the patch and evidence are copied to the primary checkout.
Source changes remain uncommitted and unpushed; the sparse local baseline commit
is a diff anchor, not an integration commit. Disposable generated site/cache
files are removed after validation.
