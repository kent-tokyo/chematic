# Format capability matrix

This page is the short map of chematic's format surfaces. A check mark means an
API exists; it does not mean that every feature of the external format is
supported or preserved. Use checked APIs and format-specific references when
coordinates, stereo, charges, properties, or opaque metadata matter.

See [language bindings](language-bindings.md) for API names,
[errors and limits](error-and-limits.md) for refusal behavior, and
[Open Babel file-I/O scope](openbabel-file-io.md) for comparison rules.

## Common format conversion bridge

Python `convert_format` accepts SMILES, MOL/SDF V2000, MOL V3000, MOL2, CML,
ChemicalJSON, MolJSON, CDXML, PDB, XYZ, PDBQT, and Gaussian input. The WASM
`convert_common_format` surface intentionally covers the topology-focused
subset: SMILES, MOL/V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML.

The bridge is for bounded interoperability. It does not promise byte identity
or preservation of every format-specific field. PDB, XYZ, and PDBQT output
requires 3D coordinates. Multiple-record SDF work should use the SDF document
or iterator APIs instead of relying on implicit first-record selection.

Exact-source routes are available for validated same-format MOL2, CJSON, and
SDF documents. Cross-format conversion reconstructs the declared graph,
coordinate, and metadata subset and reports or refuses unsupported loss.

## Core support matrix

`R` and `W` mean read and write. `Bounded` means a documented subset. A dash
means the operation is not exposed on that binding.

| Format | Rust | Python | WASM / Node | Main boundary |
|---|---|---|---|---|
| SMILES | R/W | R/W | R/W | Graph/stereo text; not an identity key |
| SMARTS | match | validate/match | match | Query evaluation, not a general writer |
| HIN | bounded R/source W | bounded R/source W | inspect/source W | No edited-graph writer |
| MOL/SDF V2000/V3000 | R/W | R/W | R/W | Use checked/source routes for stereo and SD fields |
| MOL2 | R/W | R/W | R/W | Typed/charge/metadata preservation is profile-specific |
| CML | R/W | R/W | R/W | Single-molecule checked profiles; unknown XML is not interpreted |
| ChemicalJSON / MolJSON | R/W | R/W | R/W | Strict local ChemicalJSON profile |
| CDXML | bounded R/W | bounded R/W | bounded R/W | Rich editing and all ChemDraw semantics are out of scope |
| PDB | R/W | R/W | R/W | Model/record and metadata handling is API-specific |
| mmCIF | R/W | R/W | R/W | Molecular and crystallographic profiles are distinct |
| CIF, small molecule | R/W | R | — | No Python writer or WASM surface |
| PQR | R/W | R/W | R/W | Numeric atom fields are checked |
| XYZ / Extended XYZ | R/W | R/bounded W | R/W | Plain XYZ does not carry bonds |
| PDBQT | R/W | R/W | R/W | Rigid checked-output profile |
| QCSchema | R/W | R/W | validate/read | No general WASM AtomicInput/Result writer |
| ORCA input | R/W | R/W | R/W | Explicit units; checked xyz output retains represented fields and numeric bits |
| ORCA output | R | R | R | No writer |
| Gaussian input/log | R/W or R | R/W or R | bounded | Input and output are separate APIs |
| Gaussian Cube | R/W | R/W | R/W | Volumetric grid profile |
| OpenDX | R/W | R/W | R/W | Strict and explicitly lossy writers |
| ChemDoodle JSON | bounded R/W | bounded R/W | — | Molecular graph subset |
| LAMMPS data | R/W | R/W | R/W | Units and atom styles are explicit |
| LAMMPS dump/trajectory | R/W | R/W | R/W | Streaming/materialization differs by binding |
| DMol3, MOPAC, FHI-aims, EXYZ | R/W | R/W | — | Cartesian geometry profiles |

## Checked conversion reports (initial adapter)

`convert_with_report`-style APIs separate output from observations about dropped
or unsupported meaning. Inspect the report before accepting an output. Strict
paths refuse when the requested contract cannot be preserved; lossy paths must
be selected explicitly.

This contract is narrower than successful parsing. A parser may accept a file
while a requested cross-format conversion still refuses because the target
cannot carry required coordinates, charges, stereo, properties, or metadata.

## Immutable ChemicalJSON records

The local CJSON profile rejects duplicate keys, malformed known arrays,
self/duplicate bonds, unsupported order tokens, and conflicting molecule-name
declarations. Unknown fields can be retained as opaque source data, but they are
not treated as interpreted chemistry.

Same-format checked conversion validates and returns exact source. CJSON → MOL2
requires a corresponding graph and explicit charge/name agreement; ambiguous
or incomplete mappings refuse instead of choosing a value.

Related records:

- [checked CJSON input](../benchmarks/2026-10-22-cjson-checked-input.md)
- [source conversion](../benchmarks/2026-10-22-cjson-source-conversion.md)
- [CJSON ⇄ MOL2 names](../benchmarks/2026-10-21-cjson-name-preservation.md)

## Source-preserving document profiles

- **MOL2:** a validated same-format route retains coordinates, atom typing,
  partial charges, comments, and opaque sections. Edited graph output is a
  different contract.
- **SDF:** the complete-document route retains record order, SD fields, UTF-8
  bytes, and line endings after validating every local-profile record.
- **CML:** checked source retention preserves the original XML document. The
  separate graph/coordinate projection validates the interpreted subset and
  rejects partial or conflicting coordinate declarations. Its namespace-aware
  profile excludes foreign annotations from the graph and resolves forward IDs.
- **V2000/V3000:** coordinate projection keeps source-order finite x/y/z values, but
  target formats retain only the fields they can represent.

## LAMMPS dump/trajectory

Rust, Python, and WASM expose bounded dump/trajectory readers and writers.
Coordinate conversion depends on declared box and coordinate style; callers
must not treat scaled, wrapped, and Cartesian positions as interchangeable.
See [language bindings](language-bindings.md)
for ownership and return-shape differences.

## Formats that never fabricate bonds

Plain XYZ and geometry-only quantum-chemistry profiles do not encode molecular
connectivity. Their readers preserve atoms and coordinates without inventing
bonds. Call a separate perception step only when that policy is appropriate for
the workflow.

For molecular formats that do carry connectivity, malformed or unsupported
bonds are errors or typed unsupported outcomes. They are not silently replaced
with single bonds.

## Format-specific references

Use these pages when the short matrix is not enough:

- [MOL2 hydrogen order](mol2-hydrogen-order.md) and
  [Kekulé output](mol2-kekule-bond-output.md)
- [CIF molecular bonds](cif-molecular-bonds.md) and
  [mmCIF chain records](mmcif-chain-records.md)
- [PDB atom selection](pdb-atom-selection.md),
  [model splitting](pdb-model-splitting.md), and
  [checked PDBQT output](pdbqt-checked-output.md)
  and [explicit PDBQT isotope projection](pdbqt-isotope-projection.md)
- [CML collections](cml-collection.md) and [named atoms](cml-named-atoms.md)
- [CML namespace-aware projection](cml-namespace-projection.md)
- [ChemDraw block formats](chemdraw-block.md) and
  [CDXML stereo options](cdxml-stereo-options.md)
- [ORCA Cartesian integrity](orca-cartesian-integrity.md),
  [MOPAC Cartesian](openbabel-mopac-cartesian.md),
  [DMol3](openbabel-dmol3.md), and [FHI-aims](openbabel-fhiaims.md)
- [XYZ precision](xyz-output-precision.md),
  [per-atom data](xyz-atom-data.md), and [GRO/LAMMPS units](gro-lammps-units.md)
- [V2000 source coordinate projection](v2000-coordinate-projection.md)

The generated [Python API reference](../api/chematic.md) is authoritative for
current Python signatures. Rustdoc and TypeScript declarations remain the
binding-specific sources of truth for Rust and WASM/Node.
