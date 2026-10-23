# Issue #473: rich document and semantic boundary

This document is the downstream contract for issue #473. It deliberately
separates data that is preserved, data that can be edited or expanded within a
bounded model, and constructs that are rejected. A consumer must not infer
full ChemDraw, MDL RXN, polymer, or nucleic-acid semantics from a successful
parse of a simpler input.

## Capability matrix

| Requested area | Stable entry point | Preserved or supported | Explicitly unsupported / loss behavior |
| --- | --- | --- | --- |
| Rich RXN | `ReactionDocument`, `parse_rxn_document`, `write_rxn_document` | Ordered steps in JSON, reactant/agent/product roles, coefficients, conditions, provenance, deterministic component IDs, optional validated atom-map identities | RXN V2000 has one reaction boundary and no agent/condition/provenance/coefficient channels. Writing such a document returns typed `ReactionLoss` data; it never silently flattens the document. Atom-map identities are component-local map-number/index pairs and remain optional for older authored JSON. |
| CDXML presentation | `CdxmlDocument`, `CdxmlEdit`, `parse_cdxml_document_json`, `edit_cdxml_document_json` | Source XML, page IDs/attributes, nested object paths, unknown objects and attributes, bounded command edits | Arbitrary semantic ChemDraw regeneration, full styles/fonts/annotations/reaction semantics, and unrestricted nested graph editing are not interpreted. Unknown presentation XML remains opaque; malformed or over-limit edits return `CdxmlError`. |
| Markush/polymer | `SemanticModel`, `SemanticCommand`, `semantic_*_json` | Explicit R-group selection, explicit repeat units, two `[*]` linkage markers or two endpoint atom indices, bounded expansion, `source_to_expanded` mapping | No endpoint guessing, implicit alternatives, branching/cross-linking, stochastic or topology-changing expansion, or flattening of unsupported constructs. Ambiguous input returns `SemanticError`. |
| Nucleic acids | `NucleicAcidDocument`, `NucleicAcidCommand`, `nucleic_acid_*_json` | Ordered DNA/RNA strands, explicit residues and known modifications, sugar identity, one-to-one atom ownership, adjacent linear phosphodiester linkages, annotations, bounded metadata edits | No graph-to-sequence inference, graph generation, cross-strand/cyclic/branched topology, unknown modification semantics, PDB/mmCIF biopolymer interpretation, or silent flattening. Unsupported or ambiguous input returns `NucleicAcidError`. |
| Other biopolymers | No stable first-class API | Ordinary molecular fragments may use existing molecule formats when their semantics fit those formats | Protein, peptide, carbohydrate, and unrestricted biopolymer residue/linkage semantics remain upstream or are rejected rather than flattened. |

## Binding entry points

The same bounded contracts are exposed as follows:

- Rust: `chematic_rxn::ReactionDocument` and `chematic_mol::{CdxmlDocument,
  SemanticModel, NucleicAcidDocument}`.
- Python: `parse_cdxml_document_json`, `edit_cdxml_document_json`,
  `semantic_model_json`, `semantic_apply_json_command`,
  `semantic_expand_json`, `nucleic_acid_validate_json`,
  `nucleic_acid_apply_json_command`, and RXN document conversion functions.
- WASM: `rxn_document_from_rxn`, `rxn_document_to_rxn`,
  `cdxml_document_json`, `edit_cdxml_document_json`, and the three
  `semantic_*_json` functions, plus the two `nucleic_acid_*_json` functions.
- Node: the WASM-generated JSON boundary; Node must consume the same JSON
  schema and error categories rather than reimplementing the model.

All binding-facing parsers apply input limits. JSON-facing operations return a
typed or stable error string for invalid schema, unsupported richness, resource
limits, or ambiguous topology. They do not silently discard requested data.

## Evidence and regression coverage

The contract is backed by repository-local tests for:

- typed reaction JSON round trips and RXN loss rejection;
- CDXML page/object preservation, nested object-path edits, and limit errors;
- semantic JSON round trips, explicit Markush selection, bounded repeat
  expansion, source mapping, and unsafe-topology rejection;
- shared Rust/Python/WASM nucleic-acid fixtures for short DNA, short RNA, an
  explicit modification, annotation edits, malformed atom ownership, unknown
  modification, unsupported topology, and resource limits;
- Python and WASM shared semantic/CDXML/RXN contract calls.

These tests establish the bounded v1.0 contracts. They are not evidence of
complete ChemDraw, polymer, Markush, nucleic-acid, or general biopolymer
interoperability. Full semantics remain separate future work and require new
schemas, fixtures, and per-binding compatibility gates before being advertised.

## Downstream integration rule

An editor should branch on the returned contract/error, not guess from the
presence of fields:

1. use the typed document or semantic JSON when preservation is required;
2. use an explicit command or bounded expansion when the operation is within
   the documented subset;
3. keep the original source or reject the operation when the requested
construct is unsupported or would lose meaning.
