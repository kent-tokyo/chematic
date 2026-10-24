# Compatibility scope

This page defines the public 1.x boundary. “Compatible” means an operation is
implemented under its stated contract. It does not mean chematic is a drop-in
replacement for RDKit, ChemDraw, or another toolkit.

For exact counts and versions, use the [validation report](validation.md).
Rich-document issue boundaries are listed in the
[document matrix](issue-473-rich-document-boundary.md).

## 1.x contract

- `canonical_smiles()` is a representation, not a universal identity key.
  `canonical_smiles_stable_key()` refuses structures whose stability is not
  established.
- Aromaticity and CIP have explicit modes. Unresolved stereochemistry is not
  guessed.
- Python `RWMol`, RDKit-style APIs, CDXML editing, Markush/polymer expansion,
  and R-group decomposition are documented subsets.
- 3D generation and MMFF94 are experimental. A successful result does not
  establish universal parameter coverage, convergence, or conformer quality.
- Rust is the reference implementation. Python and WASM/Node expose selected,
  versioned operations rather than every Rust type or RDKit method.

The v1.2.1 source is qualified before publication. The latest completed hosted
cross-platform runtime and chemistry-gate evidence remains v1.2.0 until the
v1.2.1 tag workflow finishes. Older package measurements keep their original
scope and are not relabelled as v1.2.1 results.

## Main bounded surfaces

### Python `RWMol`

Supported operations are construction, `AddAtom`, `AddBond`, `RemoveAtom`,
`RemoveBond`, `GetNumAtoms`, `GetNumBonds`, and `GetMol`, with the documented
single, double, triple, and aromatic bond orders.

Arbitrary property mutation, proxy iteration, conformer/query/reaction
editing, sanitization flags, and full RDKit exception parity are not included.
Unsupported options fail explicitly.

### Wildcards and R-group labels

MOL/SDF, V3000, and bounded CDXML paths distinguish carbon, `*`, unnumbered
`R`, and `R1..R9999`. SMILES has no matching label side channel and writes
these forms as `*`.

R-group decomposition supports terminal labelled wildcards and mapped core
atoms with input-aligned columns in Rust, Python, and WASM/Node. Duplicate or
non-terminal labels, ambiguous attachments, and resource exhaustion fail
explicitly. Multi-core alignment, tautomer/core enumeration, multiple
substituents per site, and RDKit GA scoring remain outside the contract.

### CDXML and rich documents

`CdxmlDocument` preserves multiple pages, unknown attributes/objects, and
untouched XML. Bounded commands can inspect or replace selected page and
presentation fields.

chematic is not a complete ChemDraw editor. It does not interpret or
regenerate every reaction, template, query, polymer, Markush, or presentation
object from an edited molecular graph.

### Markush, polymer, and nucleic-acid documents

The `chematic.semantic.v1` contract validates explicit Markush selections and
bounded repeat expansion. Repeat endpoints must be explicit; branching,
cross-linking, stochastic/nested expansion, and unsafe endpoint inference are
rejected.

The separate `chematic.nucleic-acid.v1` contract preserves bounded DNA/RNA
document metadata and adjacent linear phosphodiester linkages. It does not
infer residues from a graph, create a molecular graph, or accept branching,
cross-strand/cyclic linkages, and unknown modification semantics.

Generic molecular graph APIs can carry glycans, but WURCS notation and residue
or linkage semantics are not implemented.

## RDKit-oriented compatibility

`chematic.rdkit_compat` is a selected 2D compatibility layer. Matching names
do not imply identical defaults or every RDKit option.

MOL/SDF readers follow the [MDL property precedence](v2000-property-precedence.md)
for sparse `M  CHG`/`M  ISO` records. In the bounded mixed atom-block/property
packet, RDKit 2025.09.3 retains some old isotope labels that chematic and
Open Babel treat as unspecified. This difference also applies to Python's
`SDMolSupplier`; it is not full RDKit input-interpretation parity.

- Native ECFP/Morgan uses a different hash and must not be mixed with RDKit
  bit positions. Use the named RDKit-compatible Morgan API when bit identity is
  required under its documented options.
- The default and RDKit-like aromaticity models remain distinct.
- Accurate CIP is opt-in. Unsupported, tied, budget-limited, wildcard, and
  oracle-sensitive cases may return a typed abstention.
- Canonical SMILES, aromaticity, CIP, and stereo each retain separate residuals
  and must be validated independently.

See [RDKit migration](rdkit-migration.md) for API mappings and
[RDKit comparison](rdkit-comparison.md) for dated evidence.

## 3D, force fields, and depiction

The legacy 3D path is retained for compatibility. The named fixed-seed
RDKit-compatible ETKDGv3 path has bounded corpus evidence, not a promise for
all seeds, topologies, platforms, or public artifacts.

MMFF94 exposes typed parameters, energies, gradients, prepared evaluation,
and minimization for its supported domain. Missing parameters, iteration
limits, and final stereo violations stay observable. UFF or DREIDING may be a
better fit for structures outside that domain.

Conformer generation is a search over candidate basins. One successful
embedding is not evidence that the lowest-energy or experimentally relevant
conformer was found. Compare typing, energy, gradient, convergence, stereo,
clash, quality, time, and memory separately.

2D layout is rule-based rather than RDKit `Compute2DCoords`; coordinates will
differ. Checked MOL and depiction paths report or refuse unexpressed stereo
instead of treating successful serialization as full semantic preservation.

## Binding contract

Shared operations should return the same documented chemistry and failure
classification across bindings. Return shapes can differ where ownership or
serialization requires it.

Atom-output order and connected-component source indices are available in
Rust and Python. Checked reaction APIs expose product-origin and template-map
data in Python and WASM/Node, but their index conventions differ as documented
by each binding. Product SMILES alone does not encode provenance.

For operational limits and error shapes, see [Errors and limits](error-and-limits.md).
For per-format support, see the [format matrix](format-capabilities.md).
