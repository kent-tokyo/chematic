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

## Aromaticity and CIP

The default aromaticity path is the compatibility-preserving per-SSSR Hückel
model. `AromaticityAlgorithm::RdkitLike` is an explicit opt-in whole-graph
model for fused and non-alternant systems. Purine and azulene are covered by
the RDKit-like regression gate, but the two models intentionally remain
distinct; neither mode claims universal RDKit aromaticity parity. Known
bridgehead-N and other fused-ring gaps remain documented residuals.

The default CIP path remains the fast legacy assignment. Accurate hierarchical
CIP is opt-in through `CipMode::Accurate` (and its named Python/WASM binding
endpoints). It improves coverage and parity for many cases but is not a
universal guarantee: symmetric cages, tied Rule 4b/pseudoasymmetric cases,
and unsupported or budget-limited structures may remain unresolved. An
unresolved result is not a guessed R/S label.

## 3D and force fields

3D generation and MMFF94 minimization are experimental. The legacy
`generate_coords`/`generate_3d` path is retained for compatibility and does
not promise ETKDGv3-quality conformers or success on every topology. The
bounded `embed_pipeline_v2` path is opt-in and reports stage/provenance
evidence with typed failure outcomes; a successful result still guarantees
sanity checks, not conformational-quality or RDKit parity.

MMFF94 exposes the implemented energy terms and minimization APIs, but
typing, parameter coverage, charged/metal/fused-ring cases, and convergence
remain incomplete. Missing parameters or failed minimization must remain
observable as failure; callers must not interpret the presence of all seven
term implementations as universal MMFF94 coverage. Use UFF or DREIDING when
the selected molecule and policy require a different supported scope.

## Binding contract

Rust is the reference implementation. Python and WASM expose selected,
versioned operations; Node uses the WASM surface. Binding parity means the
same documented result/error contract for shared fixtures, not every Rust type
or every RDKit method.
