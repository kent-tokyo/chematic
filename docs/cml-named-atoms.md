# Named atom IDs in CML output

Open Babel [#2120](https://github.com/openbabel/openbabel/issues/2120) requests
MOL2 atom names in CML instead of generated `a1`, `a2`, … identifiers.
Chematic's existing graph writer also generates those identifiers.
`write_cml_named_checked` adds an opt-in Rust/Python writer with exact
caller-supplied IDs and optional Cartesian 3D coordinates.

```python
import json
import chematic

record = json.loads(chematic.parse_mol2_record_json(mol2_text))
mol = chematic.from_mol2(mol2_text)
atom_ids = [atom["name"] for atom in record["atoms"]]
cml = chematic.write_cml_named_checked(mol, atom_ids, record["coordinates"])
```

This uses the bounded MOL2 readers' graph projection. It is a graph-export
recipe, not a full MOL2 metadata converter: partial charges, force-field atom
types, substructure records, source IDs and opaque sections are not emitted.
Retain `record`/the source text separately when that information is needed.
Coordinates and IDs must refer to the same atom insertion order as `mol`.

Rust callers can use an immutable `Mol2SourceRecord` directly:

```rust
let record = chematic_mol::parse_mol2_source_record(text)?;
let ids = record.atoms().iter().map(|a| a.name.clone()).collect::<Vec<_>>();
let xyz = record.coordinates().iter().map(|&(x,y,z)| [x,y,z]).collect::<Vec<_>>();
let cml = chematic_mol::write_cml_named_checked(record.molecule(), &ids, Some(&xyz))?;
```

Every atom and every bond reference uses the supplied ID verbatim. Names are
never repaired, truncated or deduplicated. Duplicate or invalid MOL2 names
therefore require the caller to choose an explicit, unique mapping.
Identifiers are case-sensitive and use `[A-Za-z][A-Za-z0-9_-]*`, limited to
128 bytes. The molecular convention also permits dots, but Schema 3's
`atomRefs2` list type rejects them. This checked profile uses their intersection
and refuses dotted IDs even for isolated atoms, keeping one uniform rule.
XML punctuation and non-ASCII identifiers are refused rather than escaped into invalid IDs.

The writer accepts 1..4096 atoms, at most 8192 bonds, exactly one ID per atom,
and either no coordinates or one finite `[x, y, z]` Angstrom triple per atom.
It uses scientific notation with sufficient precision for float round trips.
Python accepts iterables and consumes at most one excess item to detect a
count mismatch; source exceptions propagate. Invalid values raise `ValueError`
and wrong Python types raise `TypeError` (coordinate shape errors can also
raise `ValueError`). Rust exposes `CmlNamedError`.

The graph profile writes elements, formal charges, isotopes, total H counts,
and single/double/triple/aromatic bonds. Unsupported stereo, atom maps,
non-chemical tags, wildcard atoms, special bonds and inconsistent aromatic
flags are refused. It does not encode spin multiplicity or source-only fields.
This is not a general molecular sanitizer or a complete CML convention writer.
Existing CML APIs remain unchanged. In particular, the existing CML reader
returns graph/2D data and does not expose IDs or 3D coordinates; those are
validated from XML by independent tools. No WASM API is added.

The ID, Cartesian coordinate and bond-reference rules follow the official
[CML molecular convention](https://www.xml-cml.org/convention/molecular),
sections 9.2, 9.5–9.7 and 11.1. Independent XML checks use the official
[Schema 3](https://www.xml-cml.org/schema/schema3/).

Reproduction, schema provenance, Open Babel controls and local checks are
recorded under `validation/results/cml-named-2026-10-17/`, with synthetic
fixtures and the probe in `validation/openbabel_issues/cml-named/`.
