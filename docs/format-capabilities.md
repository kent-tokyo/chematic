# Format Capability Matrix

Read/write coverage and important limits for Rust, Python and WASM.
The matrix is a starting point; per-format notes identify coordinate,
connectivity and round-trip boundaries.

See also: [`rdkit-migration.md`](rdkit-migration.md) for a feature-by-feature
Supported/Partial/Not-supported breakdown against RDKit, and
[`language-bindings.md`](language-bindings.md) for cross-language API
mapping (function names, return shapes, `None`/`null`/`Err` divergences).

This page describes what exists today. It does not claim "full support" for
any format — each row states exactly what is implemented and names the gaps.

### Common format conversion bridge

The Python package provides `convert_format(text, input_format, output_format)`
and the WASM package provides `convert_common_format(text, input_format,
output_format)` for explicit conversion between the most interoperable molecular formats:
SMILES, MOL/SDF (V2000), MOL V3000, MOL2, CML, ChemicalJSON, MolJSON, CDXML,
PDB, XYZ, PDBQT, and Gaussian input. Extensions are accepted as aliases (for
example `.mol2`, `smi`, and `com`). PDB/XYZ/PDBQT output requires 3D
coordinates; graph-only conversions do not claim to preserve format-specific
metadata. Unsupported formats fail before conversion with a `ValueError`.

The WASM bridge intentionally exposes the topology-only subset (SMILES, MOL,
MOL V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML) so that JavaScript
and Python share a predictable core contract. Use the existing coordinate APIs
for PDB/XYZ/PDBQT in browser workflows.

This is intentionally a bounded interoperability layer, not a claim of broad
Open Babel replacement. The format-specific APIs remain available when
coordinates, metadata, multiple records, or domain-specific options must be
preserved. The narrower production-format comparison and its strict win
criteria are documented in [`openbabel-file-io.md`](openbabel-file-io.md).

MOL2 → MOL2 is a source-preserving exception in CLI, Python, and WASM:
one validated record is returned exactly, including coordinates, partial
charges, typing, comments, and opaque sections. Python rejects coordinate,
charge, name, or comment overrides on this path. Cross-format conversion
continues to use the bounded graph/coordinate bridge described above.

MOL2 → CJSON in Python/WASM now uses the shared selected-field projection:
mapped graph, source xyz, parsed molecule name and complete declared partial
charges are checked on reparse. Incomplete datasets and numeric columns without
an assigned method fail, except `NO_CHARGES` zero placeholders, which are omitted
as unassigned. Python rejects
overrides on this route. Default Rust/Python limits are 16 MiB input/output and
10,000 atoms; WASM retains its 1,000,000-byte input cap. Other MOL2 metadata and
stereo remain outside preservation. Use `convert_with_report` for the full
loss report and explicit `allow_loss` behavior. See the
[common-route contract and comparison](benchmarks/2026-10-22-mol2-cjson-common.md).

---

## Support matrix

Read (R) / Write (W) per language. A dash means that operation is not
exposed in that language. "—" in the WASM columns for CIF means CIF has
**zero** WASM exposure at all, not an oversight in this table.

| Format | Rust R | Rust W | Python R | Python W | WASM R | WASM W |
|---|---|---|---|---|---|---|
| SMILES | Y | Y | Y | Y (`Mol.smiles` property) | Y | Y |
| SMARTS | Y | — (matching, not writing) | validation only | — | Y (match) | — |
| HIN (bounded profile) | Y | exact source only | Y | exact source only | inspection JSON | exact source only |
| MOL/SDF | Y | Y | Y | Y (`Mol.to_mol_block()`/`.to_mol_block_2d()`) | Y | Y |
| PDB | Y | Y | Y | Y (`Mol.to_pdb(coords)`) | Y | Y (`ConformerHandle.get_conformer_pdb(idx)`) |
| mmCIF | Y | Y | Y | Y | Y | Y |
| CIF (small-molecule) | Y | Y | Y | — | — | — |
| PQR | Y | Y | Y | Y | Y | Y |
| XYZ | Y | Y | Y | — | Y | Y |
| DMol3 .incoor | Y | Y | Y | Y | — | — |
| MOPAC Cartesian | Y | Y | Y | Y | — | — |
| FHI-aims geometry.in | Y | Y | Y | Y | — | — |
| ORCA-AICCM EXYZ | Y | Y | Y | Y | — | — |
| Extended XYZ | Y | Y | Y | Y | Y | Y |
| QCSchema | Y | Y | Y | Y | Y | — (validation only, no writer) |
| ORCA input | Y | Y | Y | Y | Y | Y |
| ORCA output | Y | — (no writer exists in Rust) | Y | — | Y | — |
| Gaussian Cube | Y | Y | Y (via `VolumetricGrid`) | Y (via `VolumetricGrid`) | Y | Y |
| ChemDoodle JSON | Y (molecular subset) | Y (graph) | Y | Y (graph) | — | — |
| OpenDX | Y | Y (strict + `_lossy`) | Y (via `VolumetricGrid`) | Y (via `VolumetricGrid`, strict + `_lossy`) | Y | Y (strict + `_lossy`) |
| LAMMPS data | Y | Y | Y | Y | Y | Y |
| LAMMPS dump/trajectory | Y | Y | Y | Y | Y | Y |

Notes on the cells above that need qualification:

- **HIN**: one checked `mol`/`atom`/`endmol` record, reciprocal `s/d/t/a`
  connectivity, typed source metadata and exact-source retention. No edited
  graph writer or generic conversion-bridge route is added. See the
  [HIN contract and Open Babel #3000 reproduction](openbabel-hin-aromatic-bonds.md).

- **SMILES Python write / MOL-SDF Python write**: exposed as `Mol` object
  methods/properties (`Mol.smiles`, `Mol.to_mol_block()`/
  `.to_mol_block_2d()`), not as free `formats.rs` functions in the same
  shape as the read side (`from_smiles`, `from_mol_block`, etc.). Both
  directions exist; the binding shape just differs from the read side —
  see [`language-bindings.md`](language-bindings.md) for the exact names.
- **QCSchema WASM write**: `qcschema_validate_atomic_input` /
  `qcschema_validate_atomic_result` validate a JSON document; the only WASM
  writer present is `to_qcschema_molecule_json` (for `QcMolecule`, not
  `AtomicInput`/`AtomicResult`).

---

## Per-format detail

### SMILES

- **Rust**: `chematic_smiles::{parse, write, canonical_smiles, parse_cxsmiles, write_cxsmiles, parse_smi_file, write_smi_file, random_smiles}`.
- **Python**: `from_smiles`, `from_cxsmiles`, `is_valid_smiles`, `parse_smi_file`, `write_smi_file`.
- **WASM**: `parse_smiles`/`write_smiles`, `mol_block_from_smiles`, `smiles_array_to_sdf`, `parse_cxsmiles_json`, `normalize_cxsmiles`, `is_valid_smiles`.
- **Streaming**: `parse_smi_file`/`write_smi_file` materialize fully. Separately,
  `chematic-mol`'s `smiles_table` module (`SmilesRecordReader`/`SmilesRecordWriter`,
  a **different sub-format** — delimited SMILES-table files with name/property
  columns, not plain `.smi`) implements a true `BufRead`-backed streaming
  `Iterator`. It is not exposed to Python/WASM and is out of scope for the
  format matrix this page otherwise tracks; noted here because it
  bears directly on any "no format in this codebase streams" claim.
- **Streaming limits**: `TdtReaderOptions` and `SmilesReaderOptions` bound
  physical lines and retained record state; TDT also bounds tags per record.
- **Coordinate units**: N/A — SMILES encodes graph topology and stereo parity,
  not coordinates. CXSMILES can carry optional 2D coordinates as an extension.
- **Connectivity**: native — bonds are the primary content of the format.
- **Round-trip**: graph-oriented, not byte-identical (whitespace/atom order
  are not preserved). `canonical_smiles()` is not a universal dedup/cache
  key; see [compatibility scope](compatibility-scope.md).
- **Lossy operations**: none inherent to the format itself.
- **Parse limits**: `SmilesParseLimits` controls input bytes, atom count, and
  bond count; `parse` applies finite safe defaults and
  `parse_with_limits` accepts a stricter policy.
- **SMILES tables**: `parse_smi_file_with_limits` accepts
  `SmiFileParseLimits { max_input_bytes, max_line_bytes, max_records }`; the
  default `parse_smi_file` path uses finite defaults.
- **Streaming SMILES tables**: `SmilesReaderOptions` bounds line bytes,
  yielded records, and columns per row; limit violations are typed
  `SmilesTableError` results.
- **Known limitation**: coupled E/Z systems can make canonical text
  unsuitable for identity. Use the fail-closed `canonical_smiles_stable_key()`
  within its documented domain; see [compatibility scope](compatibility-scope.md).

### SMARTS

- **Rust**: `chematic_smarts::{parse_smarts, find_matches, find_matches_with_rings, ...}`.
- **Python**: `is_valid_smarts` — validation only. Full parse/substructure
  matching is exposed through `Mol.find_matches()`/`Mol.has_substructure()`
  (on the molecule object, not as a `formats.rs`-style free function) — do
  not read this as "SMARTS has no Python support," it means the shape of the
  binding differs from the other formats in this table.
- **WASM**: `smarts_match_atoms`, `smarts_match_atoms_with_chirality`, `match_smarts_smiles`, `parse_cxsmarts_json`.
- SMARTS is a query language, not a molecule storage format; coordinate and
  round-trip columns do not apply.
- **Parse limits**: `SmartsParserConfig` sets recursive `$()` depth (default
  8, clamped to 1–16). Use `parse_smarts_with_config` to choose it explicitly.

### PDBQT

- **Rust**: `chematic_mol::{parse_pdbqt, parse_pdbqt_with_limits, write_pdbqt, PdbqtError, PdbqtParseLimits}`.
- **Coordinate units**: Ångström.
- **Connectivity**: no bond perception; rigid-body PDBQT only.
- **Required atom fields**: each ATOM/HETATM record must include an explicit
  partial charge and AutoDock atom type. The parser retains charges wider than
  the usual fixed-width field when separated from the atom type by whitespace.
  The writer rounds charges to three decimal places.
- **Parse limits**: `PdbqtParseLimits` bounds input bytes, physical line
  length/count, and ATOM/HETATM records. The existing parser uses finite
  defaults.

### Checked conversion reports (initial adapter)

Rust `convert_with_report` / `convert_with_report_with_limits`, Python
`convert_with_report` (dict) / `convert_with_report_json` (JSON string), WASM
`convert_with_report_json`, and CLI `convert-with-report` share
`chematic.conversion-report.v1`. Supported routes are MOL2 → MOL2/SMILES/CJSON
and CJSON → CJSON/SMILES/MOL2 (`smi` aliases plain SMILES). Existing conversion
APIs keep their graph/metadata contracts.

```python
report = chematic.convert_with_report(source, "mol2", "smiles")
assert report["status"] == "refused"
assert report["output"] is None
lossy = chematic.convert_with_report(source, "mol2", "smiles", policy="allow_loss")
# Inspect lossy["observations"] before using lossy["output"].
```

```sh
chematic convert-with-report --input-format mol2 --output-format smiles \
  --policy allow_loss --input ligand.mol2 --output report.json
```

- `strict` is the Rust policy variant and Python/CLI default. Cross-format
  MOL2 → SMILES/CJSON and CJSON → SMILES/MOL2 refuse with `output: null` because source information is
  lost or unverified. CLI writes the refusal JSON and exits nonzero.
- `allow_loss` returns parseable candidate output with detected losses and
  uncertainty. A graph mismatch is `unverified`; it cannot claim preservation.
  A candidate that fails to parse is refused even under this policy.
- `preserve_source` requires MOL2 → MOL2 or CJSON → CJSON. All three policies
  retain exact source bytes on these routes after bounded structural validation. The scope
  is `exact_source_bytes`; `source_semantics: opaque_retained` explicitly
  excludes interpreted/editable chemical semantics. It is not a graph gate.
- Observations identify a source field group, status, stage and reason in
  source record zero. Coordinate/partial-charge/type/name/substructure/bond
  metadata losses and opaque sections are reported before graph projection.
  The aggregate `source_text` loss also covers comments and unparsed columns.
- Candidate output is reparsed. SMILES emitted atom order or CJSON source-order
  atom arrays establish correspondence. Counts, mapped elements, aromaticity, projected formal charges and bond
  orders are checked. `candidate_atom_mapping` / `candidate_bond_mapping`
  carry source IDs and zero-based source/candidate indices only if this check
  passes. These mappings refer to the candidate even when strict policy
  refuses to emit it. `verification_scope` defines the checked fields.
- MOL2 formal-charge interpretation beyond `N.4` and source implicit-hydrogen
  assignment remain unverified. Noncoplanar 3D tetrahedral stereo is perceived
  in graph reads, but checked conversion does not verify stereo; other
  coordinate/extension stereo remains unverified. SMILES hydrogen
  counts are explicitly inferred under the graph valence model; CJSON does
  not emit an inferred hydrogen count. A verified projected graph does not establish chemical equivalence to the full source record.
- MOL2 → CJSON carries source x/y/z components and projected formal charges.
  `coordinates: preserved` requires verified graph correspondence and bitwise
  equality of reparsed f64 coordinates, including signed zero. Decimal spelling,
  physical dimensionality and coordinate-derived stereo are outside this check.
  `projected_formal_charges` is checked separately from the unverified source
  formal-charge interpretation. Source atom types/names/IDs,
  substructures, comments and opaque sections remain reported losses.
  Complete declared partial-charge datasets now emit separately as
  `partialCharges.<source method>` and are verified bit for bit through mapped
  atoms. Values are not computed, rounded into formal charges or filled when
  missing. Empty/`NO_CHARGES` method declarations do not produce a dataset.
  The parsed molecule name emits as `properties.name` and is checked on
  reparse. The source method label is retained without verifying its chemistry.
  A declared method without source values is now explicitly `unverified` for
  `partial_charges`. Omitted raw zero columns under `NO_CHARGES` are reported
  as `unassigned_zero_columns: dropped`, without inventing a charge dataset.
  Atom names/types/IDs and unsupported metadata remain losses.
- Invalid input, unsupported routes/policies and limits are errors (Python
  `ValueError`, WASM throw, CLI nonzero without a report). CLI `--output` is
  always a JSON report destination, never a standalone molecule destination.

Default Rust/Python input and report/output limits are 16 MiB with 10,000
reader atoms; WASM input is 1,000,000 bytes. This new SMILES route additionally
caps projection at 512 atoms because the current writer traverses recursively.
The cap does not apply to same-format source retention or CJSON/MOL2 output. JSON
size checks include escaping and envelope overhead. The MOL crate enables serde_json's
`float_roundtrip` feature: a finite f64 emitted as decimal JSON can be read
back without the previous 1 ULP rounding error. This also applies to existing
Rust/Python/WASM CJSON readers. Integral bond orders emit JSON integer
numbers, while aromatic order 1.5 remains fractional; this avoids the observed
Open Babel 3.2.1 CJSON reader converting floating integer orders into zero.
[Paired development checks](https://github.com/kent-tokyo/chematic/blob/main/validation/interchange/README.md) record both
engines' losses and targeted cross reads. No cross-engine superiority is established.

### Immutable ChemicalJSON records

Rust `parse_cjson_record` / `parse_cjson_record_with_limits`, `CjsonRecord`,
and `write_cjson_record` retain the exact original UTF-8 source. Python/WASM
`parse_cjson_record_json` expose `chematic.cjson-record.v1` inspection;
`roundtrip_cjson_record` validates and returns the source. The record is
immutable: opaque references are retained without editing or interpretation.

The strict local profile accepts absent `chemicalJson` or integer version 1,
atomic numbers, optional formal charges/isotopes/3D coordinates, bond
indices/orders, `properties.name` and root `partialCharges` datasets. Present arrays must have matching lengths and correct types.
Duplicate decoded JSON keys, self/duplicate bonds, dangling indices and
unsupported orders are rejected. Supported order values are exactly 1, 1.5,
2, 3 and 4; 4 is an explicit local alias for aromatic bond 1.5. Atom aromaticity
and coordinate stereo are not perceived. The legacy `parse_cjson` API keeps
its permissive behavior; strict validation is additive.

Root `partialCharges` must be an object of complete finite numeric arrays, one
value per atom; `properties` must be an object and its optional `name` a string.
These fields are validated even during same-format retention. The inspection
API exposes `name` and `partial_charge_datasets`; Rust has corresponding getters.
Method keys remain source declarations, with no charge calculation or model
validation. `atoms.partialCharges` remains opaque and is not conflated with the
root dataset representation.

Unknown fields, including other properties and stereo extensions, are listed
as sorted JSON Pointer paths in `opaque_fields`. Whole opaque subtrees are
represented by their root paths. They survive same-format source retention;
CJSON → SMILES reports their loss and separately reports typed name/charge
dataset losses. That route also reports coordinate loss,
defaulted charges/isotopes/orders where absent, inferred hydrogen counts and
unverified stereo. Candidate verification checks mapped elements, aromatic
flags, formal charges, isotopes and bond orders. CJSON mapping `source_id`
values are zero-based array indices; opaque user IDs are not interpreted.
Strict CJSON → SMILES refuses even when the projected graph verifies.

```python
import json
record = json.loads(chematic.parse_cjson_record_json(source))
assert chematic.roundtrip_cjson_record(source) == source
assert chematic.convert_with_report(source, "cjson", "cjson")["output"] == source
lossy = chematic.convert_with_report(source, "cjson", "smiles", policy="allow_loss")
```


Checked CJSON → MOL2 accepts explicit 3D coordinates and at most one declared
partial-charge dataset. It writes shortest round-trip decimal coordinates and
charges, reparses the candidate and checks exact f64 bits (including signed
zero/subnormals) through verified graph correspondence. `properties.name` and
the charge-method label are checked separately on reparse. Missing datasets
emit `NO_CHARGES` with no partial-charge columns; missing coordinates are never
filled. Multiple datasets, unrepresentable method tokens and names with
boundary whitespace, control characters, empty text or section markers refuse
even under `allow_loss`. Exact CJSON retention remains available for valid
records refused by this adapter.

MOL2 atom IDs/names, element-only types and `SMALL` header are generated;
charged atom rows use placeholder substructure `1 LIG`. Hybridization and
force-field types are not inferred from chemistry. The sole nonzero formal
charge representation is +1 nitrogen via `N.4`. Other nonzero formal charges
and all specified isotopes are reported dropped; they disable projected graph
verification and candidate mappings. Partial charges are never substituted for
formal charges. Aromatic bond tokens retain the local projection, without
perceiving atom aromaticity. Opaque fields and stereo remain lost/unverified.
Strict policy refuses this cross-format route even when the candidate verifies.

```python
report = chematic.convert_with_report(source, "cjson", "mol2", policy="allow_loss")
# Inspect status and observations; valid source data may still be refused.
```

Rust/Python record APIs use `CjsonParseLimits` defaults (256 MiB input,
1,000,000 atoms); WASM record APIs cap input at 1,000,000 bytes and 10,000
atoms, with 16 MiB inspection output. Checked conversion uses the smaller
16 MiB / 10,000 atom defaults described above. JSON tree limits apply to
opaque fields too. Existing `convert_format` remains a graph bridge for CJSON;
use the record or checked API for source retention.

### MOL2 (Tripos)

Python also provides [checked explicit-H row grouping](mol2-hydrogen-order.md),
with original/output correspondence and updated bond/root/charge-attribute
references. It is an opt-in source edit with a bounded section profile.

- **Multiple records / RMSD (Rust/Python)**: `Mol2SourceReader` lazily parses
  bounded records from an in-memory string; Python `parse_mol2_records_json`
  retains each source. `mol2_reference_rmsds` / Python
  `mol2_reference_rmsds_json` compare every reference with one probe, checking
  source/indexed-graph identity and optionally applying a proper rigid alignment.
  See [contracts, limits and the upstream diagnostic](mol2-reference-rmsd.md).
Generated/editable Rust output and Python/WASM graph output support an
[opt-in adjacent ATOM/BOND section order](mol2-section-order.md). The default
retains UNITY before BOND for measured Open Babel charge compatibility;
source-preserving APIs retain the original order.

- **Rust**: `chematic_mol::{parse_mol2, parse_mol2_with_limits, write_mol2}` retains the published graph/editable-record contract. Exact-source inspection uses the additive `Mol2SourceRecord`, `parse_mol2_source_record`, `parse_mol2_source_record_with_limits`, and `write_mol2_source_record` APIs. `Mol2SourceAtomRecord` and `Mol2SourceBondRecord` describe their metadata; the shared `Mol2Error`/`Mol2ParseLimits` types remain compatible.
- **Python/WASM**: `parse_mol2_record_json` exposes the shared
  `chematic.mol2-record.v1` inspection schema; `roundtrip_mol2_record` validates
  and returns the original source. Rust exposes the same schema through
  `Mol2SourceRecord::to_json()`.
- **Coordinate units**: Ångström.
- **Immutable source-profile connectivity**: native MOL2 bond section. `am` projects to single;
  `ar` to aromatic. Unknown bond tokens reject instead of becoming single.
  Duplicate IDs/edges, self bonds, dangling references, mismatched counts,
  nonfinite numbers, and multiple records/duplicate core sections reject.
- **Immutable source-profile charges**: partial charges are validated separately from formal charge.
  The graph recognizes `N.4` as +1; it does not round partial charges or
  implement general formal-charge inference/UNITY semantics. Optional charge
  columns remain optional in the record API, with the source charge-method
  header retained. `NO_CHARGES` means numerical columns are not a claim of
  assigned partial charges.
- **Round-trip**: immutable `Mol2SourceRecord` retains the exact source and typed
  atom/bond metadata, including source IDs, names, atom types, substructure
  references, partial charges, and original bond tokens. Extra sections are
  listed as opaque; source retention does not prove their interpretation.
  Noncoplanar 3D coordinates are used to perceive tetrahedral configuration;
  planar and degenerate coordinates do not produce a tag. Other
  geometry-derived stereo remains unverified.
- **Graph-only write**: `write_mol2` emits `NO_CHARGES`, simplified atom types,
  and zero numerical charge placeholders; it does not compute Gasteiger
  charges or retain source metadata. Missing coordinates still default to
  zero on this compatibility API. The published writer encodes formal charges
  with UNITY attributes; isotope, atom-map and stereo annotations are not
  generally encoded. Use the immutable source API for exact source retention.
- **Parse limits**: `Mol2ParseLimits` bounds input bytes, physical line
  length/count, section count, atom records, and bond records. The existing
  parser uses finite defaults.

### CML

- **Rust**: `chematic_mol::{parse_cml, parse_cml_strict, parse_cml_with_limits, parse_cml_strict_with_limits, write_cml, CmlError, CmlParseLimits}`.
- **Strict boundary**: `parse_cml_strict` is opt-in; it requires a non-empty
  molecule, balanced XML-ish elements, exactly one root element, and rejects
  unterminated quoted attributes, comments, CDATA sections, and processing
  instructions. The historical `parse_cml` path remains lenient for
  compatibility.
- **Coordinate units**: Ångström in returned 2D coordinates.
- **Connectivity**: native CML bond elements.
- **Parse limits**: `CmlParseLimits` bounds input bytes, physical line
  length/count, XML element count, and atom/bond elements. The existing parser
  uses finite defaults.

### MolJSON

- **Rust**: `chematic_mol::{parse_moljson, parse_moljson_with_limits, write_moljson, MolJsonError, MolJsonParseLimits}`.
- **Parse limits**: `MolJsonParseLimits` bounds JSON input bytes/depth,
  array/string resources, and atom/bond records. The existing parser uses
  finite defaults.

### ChemicalJSON

- **Rust**: `chematic_mol::{parse_cjson, parse_cjson_with_limits, write_cjson, CjsonError, CjsonParseLimits}`.
- **WASM**: `mol_from_cjson` provides the bounded topology-handle entry point;
  coordinates and CJSON metadata remain on the format conversion path.
- **Parse limits**: `CjsonParseLimits` bounds JSON input bytes/depth,
  array/string resources, and atom/bond records. The existing parser uses
  finite defaults and rejects numeric type/range truncation.

### CDXML

- **Stereo options**: [Python and WASM opt-in inference](cdxml-stereo-options.md)
  exposes the existing Rust Bold/Hash/Dash policy. Directional wedge parsing
  uses the narrow endpoint and converts screen coordinates for parity;
  source coordinates are retained.

- **Rust**: `chematic_mol::{parse_cdxml, parse_cdxml_with_limits, parse_cdxml_all, parse_cdxml_all_with_limits, write_cdxml, CdxmlError, CdxmlParseLimits, CdxmlParseOptions}`.
- **Document API**: `CdxmlDocument` supports bounded multi-page summaries,
  exact-source writing, and loss-preserving page/object command edits. Opaque
  presentation objects are retained; arbitrary ChemDraw presentation
  interpretation/rendering is outside the compatibility claim.
- **Parse limits**: `CdxmlParseLimits` bounds input bytes, physical line
  bytes/count, attribute bytes, atoms, bonds, and fragments. Existing parsers
  use finite defaults; option and limit policies can be combined with
  `parse_cdxml_all_with_options_and_limits`.

### MOL/SDF

- **Rust read**: `parse_mol`, `parse_mol_with_coords`, `read_mol_with_diagnostics`, `parse_mol_v3000*`, `SdfReader`, `SdfFileReader`, `SdfRecordReader`, `SdfDataset`.
- **Rust write**: `write_mol`, `write_mol_with_conformer[_checked]`, `write_sdf*`, `write_mol_v3000*`.
- **Python**: `from_mol_block`, `from_mol_block_with_coords`, `from_mol_block_with_diagnostics`, `parse_sdf_with_coords`, `from_mol_v3000[_with_coords|_with_diagnostics]`, `iter_sdf`, `iter_sdf_batched`, and `SdfDataset`.
- **WASM**: `mol_from_v3000_block`, `mol_from_sdf_block`, `to_mol_block`, `to_mol_v3000_block`, `roundtrip_mol_v3000_block` (preserves opaque SGROUP metadata and COLLECTION groups), `sdf_to_smiles_json`, `sdf_to_records_json`, `sdf_from_records_json`, `mol_block_stereo_diagnostics_json`, `mol_v3000_stereo_diagnostics_json`, `mol_block_coords_json`.
- **Streaming**: `SdfFileReader<R: BufRead>` is a true I/O-streaming `Iterator`
  (does not require the whole file in memory up front). `SdfBatchReader<R>` is
  the bounded pull-based variant: each item contains at most the configured
  batch size, preserves source order, retains rejected records, and can be
  cancelled at a batch boundary. Pulling the next item is the backpressure
  contract; no background queue is created. `SdfReader`/
  `SdfRecordReader` are lazy iterators over an already-loaded `&str` (do not
  eagerly collect every record into a `Vec`, but do require the full text in
  memory). Python exposes file-backed `iter_sdf` and bounded
  `iter_sdf_batched`; WASM bindings still materialize SDF input.
- **Coordinate units**: Ångström, per the Ctab standard.
- **Connectivity**: native Ctab bond table.
- **Round-trip**: semantic, not byte-identical — `write_mol` regenerates a
  Ctab from the `Molecule`'s own perceived bonds/coordinates, not a copy of
  the original text.
- **Lossy operations**: none named beyond the stereo-write scope below.
- **Square-planar stereo (`@SP1`/`@SP2`/`@SP3`)**: read automatically on
  every MOL/SDF parse. Write is **opt-in only**, via exactly 3 functions —
  `write_mol_with_conformer_checked`, `write_mol_v3000_with_conformer_checked`,
  `write_sdf_record_with_conformer_checked`. Do not describe MOL/SDF writing
  in general as "square-planar supported" — the plain `write_mol`/`write_sdf`
  path does not perceive or emit it.
- **Parse limits**: `SdfParseLimits` bounds input bytes, physical line bytes,
  individual record bytes, and yielded record count; use
  `SdfFileReader::with_limits`, `SdfRecordReader::with_limits`, or
  `parse_sdf_with_limits` for an explicit policy. The existing `parse_sdf`,
  `SdfReader`, and `SdfRecordReader::new` use finite defaults.
- **Indexed access**: `SdfDataset::open` scans delimiter lines once and stores
  byte ranges plus titles, not parsed molecules. `len`, `metadata`, `get`,
  `record_text`, and `read_range` then provide bounded-memory random access.
  Python exposes the same boundary as `SdfDataset`; malformed records retain
  their index and fail only when requested.

### PDB

- **Rust**: `chematic_3d::{PdbAtom, parse_pdb_atoms, parse_pdb_atoms_strict, pdb_to_molecule, write_pdb}` —
  geometry-based reading lives in `chematic-3d`. `chematic_mol::parse_pdb_conect`
  provides the separate explicit-graph projection described below.
- **Python**: `from_pdb` (lenient) and `from_pdb_strict` (fixed-column
  validation) for reading; both delegate to `chematic_3d`;
  `Mol.to_pdb(coords)` for writing. `protein_from_pdb_json` exposes the shared
  amino-acid-only model/chain/residue/atom projection.
- **WASM**: `mol_from_pdb` (lenient), `mol_from_pdb_strict` (fixed-column
  validation), and `pdb_coords_json` for reading. Writing is
  exposed as a method, not a free function:
  `ConformerHandle.get_conformer_pdb(idx)` returns conformer `idx` as a PDB
  string (or `null` if `idx` is out of range), delegating to
  `chematic_3d::write_pdb` internally — this is real write capability, just
  shaped around an existing conformer handle rather than a
  `write_pdb_json(mol, coords)`-style free function.
- **Protein projection**: Rust
  `chematic_mol::project_protein_from_pdb`, Python
  `protein_from_pdb_json`, and WASM `protein_from_pdb_json` preserve model,
  chain, residue, atom, coordinate, occupancy, B-factor, alternate-location,
  and insertion-code fields and expose one-letter sequences. The view is
  intentionally amino-acid-only: unrecognized residues, water, ligands, ions,
  and nucleic acids are excluded rather than guessed.
- **Coordinate units**: Ångström (PDB standard).
- **Connectivity**: `pdb_to_molecule` perceives bonds from geometry
  (distance-based), unlike most of the other formats in this table, which
  either carry an explicit bond table or never infer one.
  [Explicit CONECT graph reading](pdb-conect-graph.md) is available through
  Rust `parse_pdb_conect`, Python `from_pdb_conect`, and WASM
  `mol_from_pdb_conect` / `pdb_conect_frame_json`. This narrow, unmodelled
  single-location profile preserves explicit atoms/charges/XYZ/serials and
  repeated-target bond multiplicities without adding proximity edges.
- **Round-trip**: not characterized in this pass — no documented round-trip
  guarantee found; treat as best-effort.
- **Lossy operations**: the explicit CONECT API projects graph fields only;
  it drops residue/name/chain and other PDB metadata and performs no stereo,
  aromaticity or valence perception. Legacy extra contact slots require opt-in.
- **Parse limits**: `PdbParseLimits` bounds input bytes, physical line length,
  ATOM/HETATM records and MODEL records. Use `parse_pdb_atoms_with_limits`
  for typed resource-limit errors. `parse_pdb_atoms_strict` additionally
  rejects invalid fixed-column serial, residue sequence and XYZ fields.

### mmCIF

- **Rust**: `chematic_mol::{parse_mmcif, parse_mmcif_with_limits, write_mmcif, MmcifAtomRecord, MmcifError, MmcifParseLimits, MmcifResult}`.
- **Python**: `parse_mmcif`, `write_mmcif`, and
  `protein_from_mmcif_json`.
- **WASM**: `mol_from_mmcif`, `mmcif_coords_json`, `mmcif_to_json`,
  `write_mmcif_json`, and `protein_from_mmcif_json`.
- **Protein projection**: Rust
  `chematic_mol::project_protein_from_mmcif` and the Python/WASM JSON
  functions use the same amino-acid-only hierarchy and residue policy as the
  PDB projection above. This is not a complete biological-assembly/entity API.
- **Streaming**: no streaming reader type exists for mmCIF.
- **Coordinate units**: always orthogonal Ångström (unlike small-molecule CIF's fractional coordinates).
- **Connectivity**: **no bond table, ever** — this module never infers bonds.
- **Round-trip**: lossy on read for the `label_*`/`auth_*` tag-pair collapse —
  `auth_*` is preferred when both are present, both are treated as one
  logical field, except `label_seq_id` (kept as its own field since it is
  legitimately `.` for non-polymer atoms even when `auth_seq_id` is normal).
  `write_mmcif` writes the single stored value into both columns, so a
  read→write round trip loses the distinction if the source file's `label_*`
  and `auth_*` genuinely disagreed. Unrecognized `_atom_site` columns are
  collected into `unhandled_columns` rather than silently dropped, but are
  not necessarily re-emitted verbatim on write.
- **Lossy operations**: the `label_*`/`auth_*` collapse above.
- **Parse limits**: `MmcifParseLimits { max_input_bytes, max_atoms, max_line_len }`.

### CIF (plain, small-molecule)

- **Rust**: `chematic_mol::{parse_cif, write_cif, write_cif_checked, write_cif_with_bonds_checked, CifError, CifBondWriteError, CifResult, UnitCell}`, plus `#[cfg(feature = "crystal")]`-gated periodic-structure variants.
- **Python**: `parse_cif`, `write_cif_checked`, `write_cif_with_bonds_checked`.
- **WASM**: **none.** Zero WASM exposure — confirmed by grep, not by omission from this table.
- **Coordinate units**: fractional cell coordinates converted to orthogonal Ångström on read (IUCr convention).
- **Connectivity**: the reader returns atoms with no bonds. The opt-in
  [molecular bond writer](cif-molecular-bonds.md) emits explicit graph edges as
  `_geom_bond_*` rows, including Cartesian distances and S/D/T/A bond types.
  It does not infer connectivity, wrap periodic images or read that bond table.
- **Round-trip**: symmetry expansion is **not** performed on the plain (non-`crystal`-feature) path — only atoms literally listed in `_atom_site_*` are returned (effectively P1 treatment). See `docs/crystal_scope.md` for the `crystal`-feature periodic-structure variants.
- **Lossy operations**: symmetry-operation information is not applied/round-tripped on the plain path.
- **Parse limits**: `CifParseLimits { max_input_bytes, max_line_bytes,
  max_tokens, max_atoms }`; use `parse_cif_with_limits` for an explicit
  policy. The existing `parse_cif` parser uses finite defaults.

### PQR

- **Rust**: `chematic_mol::{parse_pqr, parse_pqr_with_limits, write_pqr, infer_element, PqrAtomRecord, PqrError, PqrParseLimits, PqrResult}`.
- **Python**: `parse_pqr`, `write_pqr`, `infer_element`.
- **WASM**: `mol_from_pqr`, `pqr_coords_json`, `pqr_to_json`, `write_pqr_json`, `pqr_infer_element`.
- **Coordinate units**: Ångström.
- **Connectivity**: **no bond table, none inferred** (unlike PDB, which does infer from geometry).
- **Element**: PQR has no element column in the format; `Element` is
  inferred from the atom name via `infer_element`, a documented heuristic
  matching classic PDB-parser convention for blank element columns.
- **Round-trip**: not characterized beyond the element-inference note above.
- **Lossy operations**: none named beyond element inference (which is a
  read-side reconstruction, not a write-side loss).
- **Parse limits**: `PqrParseLimits { max_input_bytes, max_atoms, max_line_len }`.

### XYZ / Extended XYZ

Two distinct functions share the name `parse_xyz` in different crates —
disambiguate by crate, not by name alone:

- `chematic_3d::parse_xyz -> Result<(Molecule, Coords3D), XyzError>` —
  returns every listed atom and coordinate, including H, with no bonds.
- `chematic_mol::parse_xyz -> Result<XyzFrame, XyzError>` — no `Molecule`,
  no bonds; the Extended-XYZ-superset document type. Also has
  `parse_extxyz`/`parse_extxyz_all`/`write_extxyz`/`write_xyz`,
  `XyzReader`/`ExtxyzReader`/`ExtxyzWriter`.
- **Python**: `from_xyz` uses `chematic_3d`'s atoms/coordinates version.
  `from_extxyz`/`from_extxyz_all`/`to_extxyz` use `chematic_mol`'s
  non-bond version.
- **WASM**: `mol_from_xyz`/`to_xyz` use `chematic_3d`'s version.
  `mol_from_extxyz`/`extxyz_frame_json`/`to_extxyz_json` use
  `chematic_mol`'s version.
- **Element symbols**: both native paths accept canonical symbols and
  [uppercase two-letter XYZ aliases](xyz-element-symbols.md). This rule is
  local to the species column; core lookup, PDB names and SMILES are unchanged.
- **Streaming**: `XyzFileReader<R: BufRead>` is a true file-backed streaming
  iterator, and `XyzBatchReader<R: BufRead>` adds bounded pull-based batches
  with input order, cancellation, and a versioned progress manifest.
  `XyzReader`/`ExtxyzReader` remain lazy iterators over an already-loaded
  `&str`. `ExtxyzFileReader<R: BufRead>` provides the same file-backed frame
  boundary for Extended XYZ, and `ExtxyzBatchReader<R: BufRead>` provides the
  corresponding bounded batch contract. Python/WASM still materialize XYZ
  frames.
- **Coordinate units**: Ångström (standard XYZ/extended-XYZ convention).
- **Connectivity**: `chematic_3d::parse_xyz` infers bonds by distance;
  `chematic_mol::parse_xyz`/`parse_extxyz` never do (no `Molecule` is even
  produced on that path).
- **Round-trip**: not characterized beyond the two-crate split above.
- **Lossy operations**: none named.
- **Parse limits**: both XYZ implementations have finite defaults. The
  `chematic_3d::XyzParseLimits` API bounds input bytes, atom count, and physical
  line length for the bond-inferring parser; `chematic_mol::XyzParseLimits`
  additionally bounds frame count for multi-frame input. Their respective
  `parse_xyz_with_limits`,
  `parse_xyz_all_with_limits`, and `parse_extxyz_with_limits` accept an
  explicit policy; default single/all-frame parsers use finite defaults.

### KET (Ketcher JSON)

- **Rust**: `chematic_mol::{parse_ket, parse_ket_with_limits, parse_ket_3d,
  parse_ket_3d_with_limits, write_ket, write_ket_3d, KetError,
  KetParseLimits}`.
- **Parse limits**: `KetParseLimits { max_input_bytes, max_atoms,
  max_bonds }`; default parsing uses finite defaults and the explicit APIs
  return typed `KetError::ResourceLimit` failures.
- **Python/WASM**: existing KET entrypoints use the bounded default parser
  path; structured Rust errors are surfaced as the binding's normal value/JS
  error representation.

### QCSchema

- **Rust**: `chematic_mol::{QcMolecule, AtomicInput, AtomicResult, parse_qcschema_molecule, write_qcschema_molecule, parse_atomic_input, write_atomic_input, parse_atomic_result, write_atomic_result, chematic_to_qc_molecule, qc_molecule_to_chematic}`.
- **Parse limits**: `QcSchemaParseLimits` bounds input bytes, JSON nesting depth,
  array entries, and string bytes. `parse_qcschema_molecule_with_limits`,
  `parse_atomic_input_with_limits`, and `parse_atomic_result_with_limits` accept
  an explicit policy; the existing parsers use finite defaults.
- **Python**: same names, routed through Python's stdlib `json` module rather than a hand-mapped dict.
- **WASM**: `mol_from_qcschema_molecule`, `qcschema_molecule_coords_json`, `to_qcschema_molecule_json`, `qcschema_validate_atomic_input`, `qcschema_validate_atomic_result`.
- **Coordinate units**: `QcMolecule.geometry` is explicitly **Bohr (a0)** in Rust; Python and WASM bindings convert to Ångström for convenience.
- **Connectivity**: `connectivity: Option<Vec<(usize, usize, f64)>>` is an
  optional explicit bond table — used if present, never fabricated if absent.
- **Round-trip**: open extensibility bags (`extras`/`keywords`/`protocols`/
  `native_files`/`properties`) round-trip losslessly via `BTreeMap`, kept
  distinct from an `unknown_fields` bag for genuinely undocumented
  top-level keys.
- **Lossy operations**: none named — every numeric leaf is rejected if
  non-finite rather than silently coerced.

### ORCA input

- **Rust**: `chematic_mol::{parse_orca_input, parse_orca_input_with_limits, write_orca_input, OrcaInput, OrcaInputError, OrcaInputParseLimits, OrcaBlock, OrcaCoords, OrcaAtom}`.
- **Python**: `parse_orca_input`, `write_orca_input`.
- **WASM**: `mol_from_orca_input`, `orca_input_coords_json`, `orca_input_to_json`, `write_orca_input_json`.
- **Coordinate units**: Ångström.
- **Connectivity**: no bond perception anywhere in this module.
- **Round-trip**: lossless preservation of unknown `%block ... end` blocks, including nested sub-blocks.
- **Lossy operations**: none named.
- **Parse limits**: `OrcaInputParseLimits` bounds input bytes, physical line
  length/count, keywords, blocks, block bytes, and coordinate atoms. The
  existing `parse_orca_input` parser uses finite defaults.

### ORCA output

- **Rust**: `chematic_mol::{parse_orca_output, parse_orca_output_with_limits, OrcaOutput, OrcaOutputError, OrcaOutputParseLimits, OrcaTermination, OrcaOptConvergence}` — **read-only, no writer, in all 3 languages.**
- **Python**: `parse_orca_output`.
- **WASM**: `orca_output_to_json`.
- **Coordinate units**: Ångström.
- **Connectivity**: no bond perception.
- **Semantics**: `termination` and `optimization_convergence` are two
  independently-reported fields — `ORCA TERMINATED NORMALLY` alone does
  **not** imply a requested geometry optimization converged; check both.
- **Parse limits**: `OrcaOutputParseLimits` bounds input bytes, physical line
  bytes/count, geometry frames/atoms, and frequency values. The existing parser
  uses finite defaults.

### Gaussian input/log

- **Rust**: `chematic_mol::{parse_gjf, parse_gjf_with_limits, write_gjf, parse_gaussian_log, parse_gaussian_log_with_limits, GaussianLogResult, GaussianParseLimits, GaussianError}`.
- **Parse limits**: `GaussianParseLimits` bounds input bytes, physical line
  length/count, section count, and atom count for both GJF and log parsing.
- **Coordinate units**: Ångström.
- **Connectivity**: no bond table is inferred.

### ChemDoodle JSON

Rust exposes `parse_cdjson`, `parse_cdjson_with_limits`, `CdjsonDocument` and
`write_cdjson_graph`. Python exposes `parse_cdjson` and `write_cdjson_graph`.
The reader returns all molecules, source IDs, exact input text and unscaled
drawing coordinates. Explicit zero H overrides are retained. Unknown fields,
shapes, query/stereo/radical/lone-pair content and unsupported bonds are rejected.
The graph writer emits element symbols and explicit native H counts; source
IDs and unrelated metadata require same-format source preservation.
This format is distinct from ChemicalJSON. See the
[checked profile and external diagnostics](chemdoodle-json.md).

### Gaussian Cube

- **Rust**: `chematic_mol::{parse_cube, parse_cube_with_limits, write_cube, CubeError, CubeFileReader, CubeParseLimits}`, plus the shared `VolumetricGrid` type.
- **Python**: scalar `VolumetricGrid.from_cube()` (staticmethod), `.to_cube()`, `.values_3d`; a separate [bounded dataset reader and explicit selection](cube-dataset-selection.md) supplies `read_cube_datasets`.
- **WASM**: `mol_from_cube`, `cube_grid_json`, `write_cube_json`, plus typed-array `cube_values_f64`, `cube_shape_u32`.
- **Streaming**: `CubeFileReader<R: BufRead>` streams the *input reading*
  only (reduces peak parse-time memory) — the returned
  `VolumetricGrid.values` is always one fully-materialized `Vec<f64>`. This
  is a different sense of "streaming" than `SdfFileReader`/`LammpsDumpReader`
  (this native reader only supports one dataset).
- **Coordinate units**: `GridUnits::{Bohr, Angstrom}` explicit tag; values
  are never silently normalized between them. The sign of the first-axis
  voxel count in the header disambiguates Bohr vs. Ångström (documented
  against a real spec source).
- **Connectivity**: no bonds, ever. `GridAtom.charge` is effective nuclear
  charge, not partial/formal charge.
- **Multi-dataset files**: native scalar readers return
  `CubeError::MultiDatasetUnsupported`. The separate Python dataset reader
  accepts positive NVAL and negative NATOMS inputs and explicitly exports
  one scalar field, retaining original text/IDs and omitted positions.
- **Round-trip**: not byte-identical (no claim of verbatim footer/comment
  preservation is made in the module docs).
- **Lossy operations**: native scalar read/write does not silently truncate.
  The Python dataset selection exports one chosen field with a source sidecar.
- **Parse limits**: `CubeParseLimits { max_input_bytes, max_atoms, max_grid_points }`.

### OpenDX

- **Rust**: `chematic_mol::{parse_opendx, parse_opendx_with_limits, write_opendx, write_opendx_lossy, OpenDxError, OpenDxParseLimits}`.
- **Python**: no free functions — only via the `VolumetricGrid` pyclass: `.from_opendx()`, `.to_opendx()`, `.to_opendx_lossy()`. No `mol_from_opendx` — OpenDX has no atom section at all.
- **WASM**: `opendx_grid_json`, `write_opendx_json`, `write_opendx_lossy_json`, plus typed-array `opendx_values_f64`, `opendx_shape_u32`.
- **Coordinate units**: OpenDX has **no in-file unit tag** — every parsed
  grid is tagged `Angstrom` by convention (an assumption about the dominant
  real-world producer, not information recovered from the file itself).
- **Connectivity**: no bonds — OpenDX carries no atom section.
- **Round-trip**: footer boilerplate is regenerated on write, not preserved byte-for-byte.
- **Lossy operations, explicitly named**:
  - `write_opendx` **fails closed** (`OpenDxError::NonAngstromUnits`) on a `Bohr`-tagged grid — it will not silently write a wrong-unit file.
  - `write_opendx_lossy` is the explicit opt-in Bohr→Ångström conversion. It rescales `origin`/`axes` only — it **never** rescales `values`.
  - Both writers refuse (`OpenDxError::AtomsNotSupported`) any grid with non-empty atoms — there is no lossy-atom-dropping path.
- **Parse limits**: `OpenDxParseLimits { max_input_bytes, max_grid_points }` — no `max_atoms` (the format has none).

### LAMMPS data

- **Rust**: `chematic_mol::{parse_lammps_data, parse_lammps_data_with_limits, write_lammps_data, LammpsData, LammpsDataError, LammpsDataParseLimits, LammpsAtom, LammpsAtomStyle, LammpsBond, LammpsBox, LammpsMass, LammpsVelocity}`.
- **Python**: `parse_lammps_data`, `write_lammps_data` (plain dict).
- **WASM**: `lammps_data_to_json`, `write_lammps_data_json`.
- **`atom_style` handling**: `atom_style` **cannot be inferred from the
  file** and must be passed explicitly. `atomic`/`charge`/`molecular`/`full`
  are supported; anything else maps to `LammpsAtomStyle::Other` and is
  rejected with `LammpsDataError::UnsupportedAtomStyle` — there is no
  best-effort guess, because `charge` and `molecular` rows are both
  genuinely ambiguous 6-field rows by column count alone.
- **Connectivity**: raw index topology (atom-type/bond indices), not
  chemically perceived — this is a standalone document type, not integrated
  with `Molecule`.
- **Section handling**: 4 sections are typed (`Masses`/`Atoms`/`Velocities`/
  `Bonds`); everything else is preserved byte-for-byte verbatim in
  `unparsed_sections` (an ordered `Vec<(String, String)>`, not a
  `HashMap`). The writer emits typed sections first, then opaque sections in
  their original relative order — interleaving between typed and opaque
  sections is **not** preserved, though each opaque section's own content is.
- **LAMMPS "Type Labels" extension**: explicitly **unsupported**, fails
  closed — any section whose name ends in `"Type Labels"` is rejected
  rather than silently mishandled.
- **Coordinate units**: none stated in-file — LAMMPS units come from the
  simulation's own `units` command, not this file type; chematic does not
  infer or convert them.
- **Round-trip**: opaque sections round-trip verbatim; typed sections
  round-trip through the typed representation (not guaranteed
  byte-identical, e.g. whitespace).
- **Parse limits**: `LammpsDataParseLimits` bounds input bytes, physical line
  length, header counts, typed-section rows, opaque-section bytes, and section
  count. Use `parse_lammps_data_with_limits`; the existing parser uses finite
  defaults.

### LAMMPS dump/trajectory

- **Rust**: `chematic_mol::{parse_lammps_dump_frame, parse_lammps_dump_frame_with_limits, write_lammps_dump_frame, write_lammps_trajectory, box_bounds_to_true, true_to_box_bounds, LammpsDumpFrame, LammpsDumpReader, LammpsDumpParseLimits, LammpsDumpError}`.
- **Python**: `parse_lammps_dump_frame`, `parse_lammps_dump_all` (**materializes the whole trajectory** — does not expose `LammpsDumpReader`'s streaming; a documented scope choice), `write_lammps_dump_frame`, `write_lammps_trajectory`, `box_bounds_to_true`, `true_to_box_bounds`; pyclass `LammpsDumpFrame` (`.column()`, `.cartesian_positions()`).
- **WASM**: `lammps_dump_frame_to_json_str`, `lammps_trajectory_to_json` (**also materializes, does not stream** — same disclosed choice), `write_lammps_dump_frame_json`, `write_lammps_trajectory_json`, `lammps_dump_cartesian_positions_json`, plus typed-array `lammps_dump_rows_f64`, `lammps_dump_cartesian_positions_f64`.
- **Streaming**: `LammpsDumpReader<R: BufRead>` is a **true streaming
  `Iterator`** over any `BufRead` — reads and yields one frame at a time
  without holding the whole trajectory in memory. Together with
  `SdfFileReader`, this is one of two `BufRead`-backed streaming readers
  among these documented Rust formats, and the only one whose Python/WASM bindings
  both deliberately materialize instead of streaming (a disclosed
  trade-off, not an oversight).
- **Box bounds**: triclinic box-bounds↔true-box conversion
  (`box_bounds_to_true`/`true_to_box_bounds`) independently verified
  against LAMMPS's own documentation. `LammpsBox` is reused from the
  `lammps_data` module; only this module performs the bound↔true conversion.
- **`cartesian_positions()`**: resolves `x y z` (passthrough) or
  `xs ys zs` (triclinic-aware scaled transform), but **deliberately does
  not** resolve `xu yu zu` ("unwrapped" coordinates) — returns `None` if
  neither of the first two triples is present, even if unwrapped
  coordinates are, since conflating "current position" with "unwrapped
  position" would mislead. Callers use `.column("xu")` directly for that case.
- **Divergence between JSON and typed-array WASM bindings**: when a frame
  has no resolvable coordinate columns,
  `lammps_dump_cartesian_positions_json` returns JSON `null`, but
  `lammps_dump_cartesian_positions_f64` returns `Err` — a `Float64Array`
  has no `null` representation. Both are doc-commented; this is disclosed,
  not accidental. See [`language-bindings.md`](language-bindings.md).
- **Coordinate units**: not stated in-file (same as LAMMPS data).
- **Parse limits**: `LammpsDumpParseLimits { max_input_bytes,
  max_line_bytes, max_atoms_per_frame, max_columns, max_frames }` bounds both
  single-frame input and the streaming reader; use
  `parse_lammps_dump_frame_with_limits` or `LammpsDumpReader::with_limits`.

---

## Formats that never fabricate bonds

mmCIF, PQR, ORCA (input and output), Gaussian Cube, OpenDX, LAMMPS data, and
plain CIF never infer or fabricate a bond table — they either have no bond
concept in the format at all, or (PQR specifically) have an element-inference
step that is documented separately from connectivity. PDB and
`chematic_3d::parse_xyz` are exceptions: both infer bonds
from 3D geometry (distance-based), which is a documented, disclosed choice,
not a default you should assume applies elsewhere.

## Cross-reference

For the exact per-language function signatures, return-value shapes, and
`None`/`null`/`Err` divergence points summarized above, see
[`language-bindings.md`](language-bindings.md). For how these facts map to
RDKit's equivalent APIs, see [`rdkit-migration.md`](rdkit-migration.md).


The immutable `Mol2SourceRecord` profile was integrated during the v1.1.0
development line without changing the published editable `Mol2Record` and
`parse_mol2_record`/`write_mol2_record` APIs. Checked conversion and JSON
inspection use the stricter immutable profile; legacy editable/graph parsing
does not claim that contract. Earlier v1.0.20 and v1.1.0 integration records
remain historical and are not current-package evidence.

The original integration validation is recorded in
`validation/results/main-integration-2026-10-11/summary.json`.
The earlier v1.0.20 development measurements remain historical evidence.

## ORCA-AICCM EXYZ

`parse_exyz` / `write_exyz` preserve explicit atom coordinates, virtual flags,
three lattice row vectors and the cell offset in angstroms. Single-document
validation rejects malformed counts, unsupported fields and singular cells.
No bonds, charge, spin or hydrogens are inferred. This dialect uses `%PBC`
and `%VIRTUAL`; the existing ASE `parse_extxyz` API uses `Lattice=` and
`Properties=` instead. See [EXYZ scope and validation](openbabel-exyz.md).

## FHI-aims geometry.in

`parse_fhiaims` / `write_fhiaims` retain coordinate kind, ordered lattice
vectors, Cartesian component constraints and optional starting-density charge
and spin moment. Cartesian geometry supports zero to three periodic directions;
fractional coordinates require three independent vectors. Custom species and
unsupported directives fail explicitly. See [scope and validation](openbabel-fhiaims.md).

MOPAC Cartesian preserves Angstrom coordinates, optional optimization flags,
atom labels, three header lines, and zero to three translation vectors. Internal
coordinates, atomic units, dummy/sparkle encodings and additional job sections
are refused; see [the checked contract](openbabel-mopac-cartesian.md).

DMol3 standalone `.incoor` retains Cartesian Bohr values and an optional 3D
cell, with explicit CODATA 2022 Angstrom projections. Output logs, extra
sections and partial periodicity are refused; see [scope](openbabel-dmol3.md).
