# Semantic model API

The semantic API is an explicit, loss-aware layer for structures that cannot
be represented safely as an ordinary `Molecule`.

## Stability contract

The schema marker `chematic.semantic.v1` is the public interchange contract.
`SemanticModel::validate` must succeed before edits or expansion. An R-group
alternative is never selected implicitly. A polymer repeat is expandable only
when it has an explicit repeat count and either a `repeat_smiles` value with
exactly two `[*]` linkage placeholders or two explicit `repeat_endpoint_atoms`.

`SemanticModel::apply` returns a new model for command-style editing. Expansion
returns `ExpandedSemantic`, including a `source_to_expanded` mapping for undo,
re-edit, and provenance display. Unsupported or ambiguous input returns a
typed `SemanticError`; callers must not treat it as a best-effort molecule.

## CDXML document contract

`CdxmlDocument` preserves the original XML and exposes page/object summaries.
`CdxmlEdit` applies bounded page-attribute or opaque-object replacements and
re-parses the result before returning it. `ReplaceObjectPath` addresses a
multiline nested object by its parent-to-child sibling path (for example,
`[0, 1]` is the second child of the first page-level object), while retaining
unknown attributes and all untouched objects in `write()` output. The binding JSON marker is
`chematic.cdxml-document.v1`.
