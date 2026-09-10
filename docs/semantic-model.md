# Semantic model API

The semantic API is an explicit, loss-aware layer for structures that cannot
be represented safely as an ordinary `Molecule`.

## Stability contract

The schema marker `chematic.semantic.v1` is the public interchange contract.
`SemanticModel::validate` must succeed before edits or expansion. An R-group
alternative is never selected implicitly. An R-group alternative may contain
one or more `[*]` linkage placeholders; placeholders are paired with
`attachment_atoms` in source order and are never copied into the expanded
graph. An R-group may contain typed `nested_groups`; nested choices are
validated recursively, selected explicitly, and expanded in stable parent
before-child order. A polymer repeat is expandable only when it has an explicit repeat
count and either a `repeat_smiles` value with exactly two `[*]` linkage
placeholders or two explicit `repeat_endpoint_atoms`. A null `repeat_count`
is a valid unselected editing state and can be set with the typed
`SetPolymerRepeatCount` command; expansion rejects it until selected.
`ClearPolymerRepeatCount` returns a selected unit to that explicit
non-expandable state without flattening it into an ordinary molecule.

When `polymer_units[].end_groups` is non-empty it must contain exactly two
SMILES fragments, ordered `[left, right]`. Each fragment must contain exactly
one `[*]` marker with one neighbor. The marker is replaced by a bond to the
corresponding polymer attachment atom; the mapping exposes the stable keys
`<unit>.end_group_left` and `<unit>.end_group_right`. Other Markush/polymer
topologies remain explicitly rejected with `SemanticError`.
For typed interchange, `polymer_units[].end_group_definitions` accepts exactly
two `PolymerEndGroup` records (`id`, `smiles`) in left/right order; their IDs
are used as stable `source_to_expanded` keys. The legacy `end_groups` and typed
form cannot be supplied together.

`SemanticModel::apply` returns a new model for command-style editing. Expansion
returns `ExpandedSemantic`, including a `source_to_expanded` mapping for undo,
re-edit, and provenance display. Unsupported or ambiguous input returns a
typed `SemanticError`; callers must not treat it as a best-effort molecule.
Expansion is bounded by `SemanticExpansionLimits`: `expand` uses the default
10,000 repeats per polymer unit and 100,000 output atoms, while trusted callers
can use `expand_with_limits` with an explicit budget. A limit failure is a
typed `SemanticError::ExpansionLimit` and is checked before each fragment is
added to the output graph.
`ExpandedSemantic::contract` returns the exact base graph captured at
expansion time; it does not attempt to infer a base graph from an arbitrary
edited expansion. Binding JSON includes the same provenance-backed graph as
`contracted_smiles`, so callers can restore the source graph without guessing.
The stable JSON contract can be decoded with `SemanticModel::from_json` and
expanded through the Python and WASM/Node `semantic_expand_json` APIs.
Those bindings also expose `semantic_expand_json_with_limits` with explicit
`max_atoms` and `max_repeat_count` arguments; the legacy function continues to
use the same finite defaults.
Markush alternatives are selected explicitly via
`semantic_apply_json_command`; no alternative is inferred.
`ClearRGroupAlternative` restores the unselected Markush state and provides a
lossless semantic contraction boundary.

## CDXML document contract

`CdxmlDocument` preserves the original XML and exposes multi-page/page-object
summaries, including presentation-only objects. `CdxmlPage::bounding_box`
provides the optional page dimensions as a typed `left, top, right, bottom`
tuple and rejects malformed values. `CdxmlDocument::diagnostics`
reports unknown presentation tags with page/object locations while retaining
their raw XML. `CdxmlObject::transform` and `CdxmlObject::z_order` provide
typed access to common `Matrix` and `ZOrder` attributes and reject malformed
values. `CdxmlObject::text_style` provides typed access to common text/caption
font, size, and alignment attributes; malformed or non-finite sizes are
rejected. Unknown attributes remain available in the object attribute map.
`CdxmlEdit` applies bounded
page-attribute or opaque-object replacements and re-parses the result before
returning it; `apply_json_edit` is the binding-neutral command boundary.
`ReplaceObjectPath` addresses a multiline nested object by its parent-to-child
sibling path (for example, `[0, 1]` is the second child of the first page-level
object), while retaining unknown attributes and all untouched objects in
`write()` output. Python exposes `parse_cdxml_document_json` and
`edit_cdxml_document_json`; WASM/Node expose the corresponding
`cdxml_document_json` and `edit_cdxml_document_json` functions. The binding
JSON marker is `chematic.cdxml-document.v1`; diagnostics are included in that
JSON summary. Full ChemDraw presentation semantics remain outside the contract
and are reported rather than silently interpreted.
