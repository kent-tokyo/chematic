# Semantic model API

The semantic API is an explicit, loss-aware layer for structures that cannot
be represented safely as an ordinary `Molecule`.

## Stability contract

The schema marker `chematic.semantic.v1` is the public interchange contract.
`SemanticModel::validate` must succeed before edits or expansion. An R-group
alternative is never selected implicitly. An R-group alternative may contain
one or more `[*]` linkage placeholders; placeholders are paired with
`attachment_atoms` in source order and are never copied into the expanded
graph. A polymer repeat is expandable only when it has an explicit repeat
count and either a `repeat_smiles` value with exactly two `[*]` linkage
placeholders or two explicit `repeat_endpoint_atoms`.

When `polymer_units[].end_groups` is non-empty it must contain exactly two
SMILES fragments, ordered `[left, right]`. Each fragment must contain exactly
one `[*]` marker with one neighbor. The marker is replaced by a bond to the
corresponding polymer attachment atom; the mapping exposes the stable keys
`<unit>.end_group_left` and `<unit>.end_group_right`. Other Markush/polymer
topologies remain explicitly rejected with `SemanticError`.

`SemanticModel::apply` returns a new model for command-style editing. Expansion
returns `ExpandedSemantic`, including a `source_to_expanded` mapping for undo,
re-edit, and provenance display. Unsupported or ambiguous input returns a
typed `SemanticError`; callers must not treat it as a best-effort molecule.
The stable JSON contract can be decoded with `SemanticModel::from_json` and
expanded through the Python and WASM/Node `semantic_expand_json` APIs.
Markush alternatives are selected explicitly via
`semantic_apply_json_command`; no alternative is inferred.

## CDXML document contract

`CdxmlDocument` preserves the original XML and exposes multi-page/page-object
summaries, including presentation-only objects. `CdxmlEdit` applies bounded
page-attribute or opaque-object replacements and re-parses the result before
returning it; `apply_json_edit` is the binding-neutral command boundary.
`ReplaceObjectPath` addresses a multiline nested object by its parent-to-child
sibling path (for example, `[0, 1]` is the second child of the first page-level
object), while retaining unknown attributes and all untouched objects in
`write()` output. Python exposes `parse_cdxml_document_json` and
`edit_cdxml_document_json`; WASM/Node expose the corresponding
`cdxml_document_json` and `edit_cdxml_document_json` functions. The binding
JSON marker is `chematic.cdxml-document.v1`.
