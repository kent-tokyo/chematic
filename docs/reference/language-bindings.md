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

### Zero-seed conformer ensembles

RDKit [#9738](https://github.com/rdkit/rdkit/issues/9738) reports that
`EmbedMultipleConfs(randomSeed=0)` repeats the same conformer. Chematic's Rust
`rdkit_embed_multiple_confs` and Python `Mol.rdkit_embed_multiple_confs` fix
this by using attempt seeds `0, 2, 3, ...`. Seed one is skipped because zero
and one alias in the underlying random generator. The first result preserves
the existing single-conformer seed-zero behavior. Positive base seeds retain
the existing `base * (index + 1)` schedule.

```python
import chematic

mol = chematic.from_smiles("CCCCCCO").add_hydrogens()
conformers = mol.rdkit_embed_multiple_confs(5, random_seed=0)
assert conformers == mol.rdkit_embed_multiple_confs(5, random_seed=0)
assert conformers[0] == mol.rdkit_embed(random_seed=0)
```

Candidate ordering remains independent of native worker count. Corrected
gradient and MMFF ensemble APIs share this schedule; the quality profile also
advances across additional seed-zero batches. Impossible signed seed ranges
are rejected before allocating the candidate ensemble. Failed embeddings are
still omitted. Different seeds do not guarantee different conformational
minima, and no RMS pruning is added. Seed-zero ensembles intentionally differ
from RDKit 2026.09.1 after the first candidate. WASM currently exposes single
RDKit-style embedding, which is unchanged.

The fixed diagnostic in `validation/rdkit-9738-2026-10-24/` records four flexible
SMILES, seeds 0/1/42 and five unoptimized conformers each. With seed zero, the
prior local Chematic runtime and RDKit 2026.09.1 show no change in heavy-atom
internal distances across candidates; the patched runtime varies by more than
0.1 angstrom for all four inputs. All 12 runs preserve repetition and prefix
stability. All 12 first conformers and all eight positive-seed ensembles match
the prior runtime exactly on the measured host. This small regression packet
does not establish ensemble coverage, energy quality or a published release.

WASM `MolHandle` exposes separately named `rdkit_*` writers, fingerprints,
alignment, stereo, scaffold, seeded-embedding, and bounds-matrix operations.
RDKit-style block readers return a result containing the molecule, coordinate
JSON, and sanitized SMILES where applicable. Bounds include one million
alignment matches, 4,096 stereo isomers, and 512 atoms for ETKDG/bounds output.

## Isotope-aware total formula

[Open Babel #1464](https://github.com/openbabel/openbabel/issues/1464) requests
non-hydrogen isotope labels and total charge in molecular formulas. The
isotope-aware total formula is available as Rust
`chematic_chem::calc_mol_formula_with_isotopes(&mol)`, Python's read-only
`Mol.formula_with_isotopes` property and WASM `MolHandle.formula_with_isotopes()`.

```python
mol = chematic.from_smiles("[13CH3]C(=O)[O-]")
assert mol.formula == "C2H3O2-"
assert mol.formula_with_isotopes == "C[13C]H3O2-"
assert chematic.from_smiles("[2H]O[2H]").formula_with_isotopes == "[2H]2O"
```

The operation includes explicit atom counts, bracket/inferred hydrogens and
net formal charge using the existing formula descriptor's hydrogen/model
convention, including its original-graph fallback. It does not add or remove
atoms from the caller's graph. The existing isotope-blind `formula` output is
unchanged. Total charge already worked in that API and in the measured Open
Babel 3.2.1; the added capability distinguishes isotope-specific counts.

Carbon and hydrogen families come first, then other symbols alphabetically.
Within an element, unspecified atoms precede numeric isotope labels in ascending
mass order. Isotopes use `[13C]`, `[15N]`, `[2H]` and `[3H]`, with a following
count when greater than one. Explicit `[12C]` remains distinct from unspecified
carbon. Mass zero means unspecified; declared nonzero mass numbers are retained
without validating physical isotope existence or abundance. Net charge uses
`+`, `-`, `+2`, `-3`. Wildcards retain `*` or a bracketed isotope-qualified `*`.

This is a total-composition string. It does not encode isotope locations,
connectivity, stereo, atom maps or radical state and is not a molecular identity
key. Existing `parse_formula` does not preserve isotope labels; it is not a
round-trip reader for this output. The core `Molecule::formula_with_isotopes`
helper remains a separate explicit-atom-only view; use the chemistry API above
for hydrogens and net charge. Existing graph/parser budgets still apply.

Authored controls and the open issue snapshot are retained under
`validation/openbabel-1464-2026-10-24/`; the report and raw outputs are in
`validation/results/openbabel-1464-2026-10-24*`. The packet covers 47 compositions,
including isotope order/multiplicity, charge cancellation, explicit/inferred
hydrogens, aromatic heteroatoms and separate wildcard/mass-label boundaries.
Native tests also verify formula invariance when adding ordinary H graph atoms.
Python/WASM tests verify that formula calls leave the graph and legacy output
unchanged.

```sh
python3 scripts/validation/validate_isotope_formula.py \
  --python-root /path/to/candidate-python \
  --rdkit-python /path/to/rdkit-python \
  --wasm-module /path/to/candidate/chematic_wasm.js
```

The Open Babel lane uses the same SMILES with `-ismi -omolreport`, normalized by
an independent formula observer. D/T aliases count as equivalent to `[2H]`/`[3H]`.
Repeated charge signs (`++`, `--`) also normalize to `+2`, `-2`.
Five wildcard/mass-label boundary controls are candidate-only. Expected fields
are authored before execution; RDKit's isotope-separated output provides an
additional cross-check. Source/runtime/fixture hashes are retained. This is
local development evidence, not a published release or a toolkit-wide ranking.

The candidate Python/WASM APIs match all 47 authored compositions. The old
isotope-blind APIs in the same runtime match 12; their complete composition
implementation is verified equal to main `38b4c92c`. This is not a separately
rebuilt main/release binary. In the 42-input common lane, Open Babel 3.2.1
matches 16 compositions: 23 differ only in isotope-specific counts, and three
also undercount ordinary hydrogen. Those hydrogen regressions are separate
from the new isotope capability. No charge differences are introduced as a
new win. RDKit 2026.09.1 matches composition on 46 accepted controls, with one
equivalent ordering difference, and refuses the chemically invalid H–H chain
control. Workers have 30-second/1 MiB combined pipe budgets, with no process
memory cap.

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
