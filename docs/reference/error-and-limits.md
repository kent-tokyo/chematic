# Errors and resource limits

This page summarizes how malformed, unsupported, ambiguous, and oversized
inputs fail. It is an operational contract, not a broad security claim.

Use the [format matrix](format-capabilities.md) for format coverage and the
[language-binding guide](language-bindings.md) for current signatures.

## Error model

Rust format modules return typed enums such as `MolParseError`, `CmlError`,
`CifError`, and `CubeError`. Python parsers normally raise
`ChematicInputError`, a `ValueError` subclass with stable `category`, `code`,
and `format` fields. Legacy Python writers and non-parser paths may still raise
plain `ValueError`.

WASM/Node support is mixed:

- checked workflow APIs return JSON envelopes with status/reason fields;
- many parser APIs still throw a message string;
- newer domain APIs, including nucleic-acid documents, return stable tagged
  error envelopes.

Do not parse human-readable messages when a checked or envelope API is
available.

## Parser limits

The following readers expose dedicated limits. Omitted values use the Rust
defaults; language bindings may apply smaller limits before calling the core.

| Input | Limits type or options | Main bounded dimensions |
|---|---|---|
| mmCIF | `MmcifParseLimits` | bytes, atoms, line length |
| PQR | `PqrParseLimits` | bytes, atoms, line length |
| Gaussian Cube | `CubeParseLimits` | bytes, atoms, grid points |
| OpenDX | `OpenDxParseLimits` | bytes, grid points |
| MOL/SDF | `SdfParseLimits` | bytes, record/line size, record count |
| XYZ / Extended XYZ | `XyzParseLimits` | bytes, atoms/frames, lines, properties |
| KET | `KetParseLimits` | bytes, atoms, bonds |
| SMILES table | `SmiFileParseLimits` | bytes, lines, records |
| Streaming SMILES/TDT | reader options | lines, records, fields/tags |
| InChI | `InchiParseLimits` | bytes, atoms |
| MDL RXN | `RxnFileParseLimits` | bytes, reactants, products, molecules |
| ORCA output | `OrcaOutputParseLimits` | lines, frames, atoms |

SMILES, SMARTS, PDB, plain CIF, QCSchema, ORCA input, and some simulation
formats do not all expose a dedicated native limits type. Callers handling
untrusted data should enforce an outer byte/time boundary even when the
binding supplies its own cap.

## Binding-level bounds

WASM workflow and rendering APIs bound input bytes, record counts, and atoms
before materializing molecule collections. Format JSON helpers also cap output
size and nested arrays. Reaction, MCS, R-group, MMP, tautomer, library, and 3D
workflows add operation-specific limits such as match visits, products,
enumerated results, or search nodes.

The MCP server independently bounds request framing, SMILES/SMARTS and MolJSON
sizes, atom/bond counts, search embeddings, response size, upstream PubChem
bodies, timeouts, and concurrent connections. CLI output is bounded before
writing to stdout or a file.

These limits are implementation details that can tighten in a compatible
release. Code that needs a specific budget should use the checked API and pass
an explicit option rather than relying on a default.

## Checked workflow outcomes

Checked APIs distinguish successful output from refusal or exhaustion:

| Workflow | Machine-readable outcome |
|---|---|
| Reaction | profile, status, reason, detail, products, and match counts |
| MCS | timeout and resource-limit flags, plus optional search budget |
| 3D pipeline | stage/cause diagnostics and final termination |
| Nucleic-acid documents | stable `error.code`, path, and message |
| Checked conversion | output plus loss/unsupported observations |

Resource exhaustion never returns a partial result as a full success. A typed
refusal is safer than a guessed answer, but it is still not a match in
validation.

## Reject rather than guess

Important fail-closed examples include:

- unsupported LAMMPS atom styles and type-label sections;
- multi-dataset Cube input on the scalar reader;
- OpenDX output with atoms, or Bohr-tagged coordinates on the strict writer;
- invalid valence in SMILES, including over-bonded neutral O or F;
- ambiguous reaction stereo that depends on source syntax not retained by the
  molecular graph;
- cross-format conversion that cannot preserve required stereo, coordinates,
  charges, properties, or metadata.

Explicitly lossy operations use a `_lossy` name. For example,
`write_opendx_lossy` opts into Bohr-to-angstrom coordinate conversion. The
strict writer refuses instead.

## Numeric safety

- Volumetric-grid indexing uses checked bounds and arithmetic.
- QCSchema rejects non-finite numeric leaves.
- LAMMPS box conversion rejects non-finite coordinates.
- Checked geometry and force-field APIs validate shapes, finite values, and
  operation-specific size limits before evaluation.

## Connectivity and streaming

Geometry-only formats do not automatically imply a molecular graph. Plain XYZ
and selected quantum-chemistry readers preserve atoms and coordinates without
inventing bonds; call a separate perception step when appropriate. Molecular
formats that carry connectivity reject malformed or unsupported bonds rather
than replacing them silently.

Rust provides true streaming readers for SDF and LAMMPS trajectories. Python
SDF iterators preserve bounded batches and progress accounting, while several
Python/WASM convenience APIs materialize complete results. Choose the
streaming API for large or untrusted inputs and inspect its rejection counts.

## Minimal examples

Strict OpenDX writing refuses a unit conversion:

```rust
let err = write_opendx(&bohr_grid).unwrap_err();
// Use write_opendx_lossy only when conversion is intentional.
```

QCSchema connectivity is optional:

```rust
let molecule = parse_qcschema_molecule(text)?;
if let Some(bonds) = &molecule.connectivity {
    // Connectivity was present in the source document.
}
```

For exact limits and return shapes, treat Rustdoc, generated Python docs, and
TypeScript declarations as authoritative. Versioned probes and abstention
audits are indexed from the [validation report](validation.md).
