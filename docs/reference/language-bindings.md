# Language binding cross-reference

Rust, Python, and WASM share the same Rust chemistry core but intentionally use
different container and error shapes. This page lists the differences that
matter when porting code. Format coverage is in
[format capabilities](format-capabilities.md); limits are in
[errors and limits](error-and-limits.md).

## Surface summary

| Surface | Typical API shape | Errors | Large numeric data |
|---|---|---|---|
| Rust (`chematic-*`) | Typed values and `BufRead` readers where available | Typed enums | Owned Rust vectors |
| Python (`chematic-py`) | `Mol`, pyclasses, and free functions | Usually `ValueError` with Rust error text | Fresh NumPy copies |
| WASM (`chematic-wasm`) | JSON-returning free functions plus selected handles | Thrown JS error | Fresh typed-array copies where exposed |

Python arrays and WASM typed arrays are copies, not views into Rust memory.
Parsers remain format-specific; there is no `parse_any(text)` auto-detection.

## Atom provenance and output order

| Capability | Rust | Python | WASM / Node |
|---|---|---|---|
| SMILES atom output order | `write_with_atom_order`, `canonical_smiles_with_atom_order` | `Mol.smiles_with_atom_order()` | Not exposed |
| Component source indices | `Molecule::fragments_with_source_atoms()` | `Mol.connected_components_with_atom_indices()` | Not exposed |
| Reaction atom provenance | traced `PreparedReaction` / runner APIs | `run_smirks_checked()` | `run_reactants_checked()` JSON |
| Rejection diagnostics | typed detailed runner | `rejection_diagnostics` | same JSON field |
| Label-stable R-group columns | `rgroup_decompose_labeled()` | `rgroup_decompose_columns()` | `rgroup_decompose_columns_json()` |

For SMILES, `order[k]` is the input atom index of the `k`-th written atom.
New product-template atoms have no reactant origin. Python reaction arrays use
the returned `Mol` order; WASM arrays use the returned canonical-SMILES parse
order. R-group output is input-aligned and uses `None`/`null` for unmatched or
absent rows.

## V3000 either bonds and unknown alkene configuration

[Open Babel #2429](https://github.com/openbabel/openbabel/issues/2429) requests
retention of V3000 bond `CFG=2`. The
[BIOVIA CTfile specification](https://www.3ds.com/assets/invest/2026-01/biovia_chemistry_2026_ctfile_formats.pdf)
defines that bond configuration as either (printed page 13). Chematic retains
this marker on single and double bonds in `MolMetadata`. A wavy single bond
also prevents coordinate-based E/Z assignment to double bonds adjacent to
either endpoint. A remote wavy bond does not suppress unrelated alkene stereo.

| Path | Preserved meaning |
|---|---|
| Rust `read_mol_v3000_with_diagnostics` followed by `write_mol_v3000` with its metadata | Original single-bond marker and endpoint order, with normalized output IDs |
| Rust `write_mol_v3000_with_conformer` with source metadata | Original either marker retained in the 3D block |
| Python `from_mol_v3000` / `Mol.to_mol_v3000` with complete real 2D coordinates | Unknown alkene remains unknown on reparse; output may use double-bond `CFG=2` |
| WASM `mol_from_v3000_block` / `to_mol_v3000_block` | Same unknown-alkene contract through generated 2D coordinates |

Graph-only bindings do not own the original marker's source position. The
2D V3000 writer marks undeclared potential alkene stereo as either when given
complete real coordinates, so the layout cannot create an E/Z assignment on
reparse. Complete zero coordinates and missing-coordinate defaults keep their
existing behavior. Already declared E/Z stereo follows the existing depiction
path. Use Python `roundtrip_mol_v3000_source` or WASM
`roundtrip_mol_v3000_block` for exact unedited source retention; those APIs
predate this change.

This addresses alkene interpretation and V3000 rewriting. It does not add SVG
wavy-bond rendering, a core distinction between unknown and absent stereo,
general tetrahedral either-bond perception, or changes to V2000 output.
Atom-line CFG is separate and remains outside this bond configuration profile.

The fixed packet `validation/openbabel-2429-2026-10-24/` contains 16 authored
wavy-alkene variants (left/right carrier, reversed endpoints, mirrored geometry,
2D/3D), a direct double-either control and an unmarked geometry control. Reading,
writing and reparsing preserve the expected defined/unknown alkene state in
18/18 patched Python cases and 18/18 independent RDKit 2026.09.1 cases. The
prior local Chematic runtime and Open Babel 3.2.1 each match 1/18 in their tested
pipelines. Literal SMILES spelling is not compared. These are local fixed-input
regression results, not general stereo coverage or published-release evidence.

## RDKit compatibility namespaces

Python provides an explicit lightweight namespace:

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import AllChem, rdFingerprintGenerator

mol = Chem.MolFromSmiles("CCO")
legacy = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
stable = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048
).GetFingerprint(mol)
assert legacy.GetOnBits() == stable.GetOnBits()
```

Unsupported options raise `NotImplementedError`; they are not replaced by a
different algorithm. This namespace is an API migration layer, not a blanket
claim that every RDKit implementation detail is reproduced.

WASM `MolHandle` exposes separately named `rdkit_*` writers, fingerprints,
alignment, stereo, scaffold, seeded-embedding, and bounds-matrix operations.
RDKit-style block readers return a result containing the molecule, coordinate
JSON, and sanitized SMILES where applicable. Bounds include one million
alignment matches, 4,096 stereo isomers, and 512 atoms for ETKDG/bounds output.

## Copy and serialization semantics

| Boundary | Contract |
|---|---|
| Python fingerprint/grid arrays | Fresh NumPy copy; fingerprint arrays are `uint8` in packed-bit order |
| WASM numeric helpers | Fresh `Float64Array`/`Uint32Array` copy |
| WASM JSON helpers | Serialized copy; simpler but slower for large arrays |
| QCSchema unknown fields | Preserved through Rust, Python, and WASM round trips |

For volumetric grids, the third axis is fastest. Python `values_3d` reshapes to
`(nx, ny, nz)` without changing that order.

## Unit ownership

| Data | Conversion rule |
|---|---|
| Gaussian Cube | Core preserves the Bohr/angstrom tag and does not silently convert |
| OpenDX | Core tags coordinates as angstrom because the format has no unit field |
| QCSchema geometry | Rust keeps Bohr; Python/WASM expose convenience angstrom coordinates |
| mmCIF/PQR/ORCA output | Core follows the documented format units |
| ORCA input | Checked APIs report explicit units without rescaling |

## `None`, `null`, and errors

Rust `Option` normally maps to Python `None` and JSON `null`. A typed-array API
cannot return `null`, so missing numeric data may become a thrown error. For
example, LAMMPS Cartesian positions behave as follows:

| Surface | Missing position columns |
|---|---|
| Rust `cartesian_positions()` | `None` |
| Python | `None` |
| WASM JSON | JSON `null` |
| WASM typed array | thrown error |

Typed Rust variants are flattened to their `Display` text in Python and WASM;
bindings do not expose every structured error field.

## Strict and lossy operations

Strict/lossy choices remain explicit across bindings. OpenDX is representative:

| Surface | Strict | Explicit lossy |
|---|---|---|
| Rust | `write_opendx` | `write_opendx_lossy` |
| Python | `grid.to_opendx()` | `grid.to_opendx_lossy()` |
| WASM | `write_opendx_json` | `write_opendx_lossy_json` |

The strict writer rejects a Bohr-tagged grid. The lossy writer rescales origin
and axes, never scalar values. Neither writer silently drops atoms.

## Streaming versus materialization

| Format | Rust | Python | WASM |
|---|---|---|---|
| MOL/SDF | `SdfFileReader<R: BufRead>` | `iter_sdf` / `iter_sdf_batched` | bounded resumable batches over a string |
| XYZ/Extended XYZ | file-backed readers | recoverable batch readers | bounded resumable batches over a string |
| LAMMPS dump | streaming iterator | materialized list | materialized JSON |
| Cube | streamed input, materialized grid | materialized grid | materialized grid |
| Other formats | materialized unless documented otherwise | materialized | materialized |

Materialization is a disclosed binding choice, not evidence that records were
silently omitted.

## Shared fail-closed rules

The core enforces these rules for every binding:

- strict OpenDX rejects non-angstrom grids and grids containing atoms;
- scalar Cube readers reject multi-dataset input unless a bounded dataset API
  and explicit selection are used;
- unsupported LAMMPS atom styles and Type Labels are rejected;
- checked coordinate and reaction APIs reject inconsistent or unsupported
  states instead of guessing.

See [validation](validation.md) for artifact evidence and the
[benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md)
for versioned measurements.
