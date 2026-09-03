# Compatibility scope

This is the final v0.89 compatibility boundary. “Compatible” means the
operation is implemented and covered by its stated contract; it does not mean
that chematic is a drop-in reimplementation of RDKit or ChemDraw.

## Python `RWMol`

Supported: construction from an empty molecule or `Mol`; `AddAtom`, `AddBond`,
`RemoveAtom`, `RemoveBond`; `GetNumAtoms`, `GetNumBonds`, and `GetMol`; and the
documented `SINGLE`, `DOUBLE`, `TRIPLE`, and `AROMATIC` bond orders.

This is an intentionally small subset. Arbitrary atom/bond property mutation,
mid-edit RDKit-style atom/bond proxy iteration, conformer/query/reaction
editing, sanitization flags, and full exception/return-value parity are not
supported. Unsupported options fail explicitly.

## CDXML

Supported: bounded molecular parsing with coordinates; loss-preserving
`CdxmlDocument` parsing and exact-source `write()`; bounded page/object
attribute edits, insertion/removal, replacement, and nested `ReplaceObjectPath`;
and preservation of unknown attributes, unknown objects, and untouched XML.

This is not a complete ChemDraw editor/writer. Arbitrary semantic editing of
every nested CDXML object, regeneration from an edited molecular graph, and
full templates, reactions, Markush/R-group, polymer, query, and presentation
semantics are not supported.

## Polymer and Markush expansion

Expansion is bounded and explicit. A repeat requires a count and either a
`repeat_smiles` fragment with exactly two `[*]` linkage markers or two explicit
`repeat_endpoint_atoms`. Expansion returns `ExpandedSemantic` with a
`source_to_expanded` provenance mapping.

Implicit endpoint guessing, unrestricted Markush selection, branching,
cross-linking, stochastic/nested/topology-changing polymer expansion, and
flattening unsupported semantic objects into `Molecule` are not supported.
Ambiguous or unsafe input returns typed `SemanticError`.

## RDKit and Morgan compatibility

`chematic.rdkit_compat` is a selected 2D compatibility layer, not a full RDKit
clone. The default ECFP/Morgan implementation is shape- and workflow-
compatible but not bit-identical: its FNV-based hashing must not be mixed with
RDKit bit positions. Named RDKit-Morgan APIs have their own documented
radius/size/options contract and do not expose every RDKit option.

Canonical SMILES, aromaticity, CIP, and stereo retain separate documented
residuals and model choices. Matching function names do not promise
algorithmic identity.

## Binding contract

Rust is the reference implementation. Python and WASM expose selected,
versioned operations; Node uses the WASM surface. Binding parity means the
same documented result/error contract for shared fixtures, not every Rust type
or every RDKit method.
