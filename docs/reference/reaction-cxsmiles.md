# Reaction CXSMILES atom groups

RDKit [#9732](https://github.com/rdkit/rdkit/issues/9732) reports an enhanced
stereo group pointing to the wrong atom after canonical reactant reordering.
Chematic provides an additive, checked reaction CXSMILES profile for enhanced
tetrahedral **atom** groups.

```python
import chematic

text = "O.C[C@H](O)CC>>CCC(C)=O |&1:2|"
output = chematic.normalize_reaction_cxsmiles(text)
assert chematic.normalize_reaction_cxsmiles(output) == output
```

The result's group refers to the emitted stereocentre, even when water changes
position. Its numeric atom reference may differ from both the input and RDKit's
output because the molecular SMILES traversal can differ.

| Surface | API |
|---|---|
| Rust parse | `parse_reaction_cxsmiles`, `parse_reaction_cxsmiles_with_limits` |
| Rust write | `write_reaction_cxsmiles(&reaction, canonical)` |
| Python | `normalize_reaction_cxsmiles(text, canonical=True)`; text accepts UTF-8 bytes |
| WASM / Node | `normalize_reaction_cxsmiles(text, canonical)` |

## Supported representation

The optional compact CX block accepts `a:...` (ABS), `&n:...` (AND) and
`on:...` (OR) groups. Atom references are unsigned, zero-based indices across
reactants, agents and products, in input component/atom order. AND/OR numbers
are positive `u32` values. A kind/number denotes one group across the whole
reaction, including across participants and roles. Parsing stores each
participant's local members in its core stereo groups; writing combines the
same kind/number again. This global-number convention is explicit for these
APIs and should be retained when editing their returned Rust `Reaction`.

With `canonical=True`, each molecule is canonicalized and participant strings
are sorted within each role. The writer maps each local group atom through
the actual returned molecular atom-output permutation and then adds its new
role/component offset. `canonical=False` keeps participant order but still
remaps through the actual molecular traversal. Group numbers are preserved.

Non-stereogenic references without a tetrahedral annotation, overlaps, repeated
group declarations, bad indices/numbers and empty groups are errors. Writers
also reject programmatic invalid members, bond groups and a grouped atom that
loses its tetrahedral annotation during molecular serialization. This is an
annotation/representation check, not a certificate of chemical stereogenicity.
Other CX fields (labels, coordinates, atom properties, radicals, wavy/bond
groups), multiple blocks and titles are refused rather than silently discarded.
Plain reaction SMILES and empty reactions are accepted. Native reaction limits
apply, with the input byte limit including CX data. WASM also applies its
existing component/atom/output limits.

This is not a globally unique enhanced-stereo reaction key, full CX metadata
exchange, reaction execution support for group correlations, or a change to
legacy reaction writers. Exact source text is not retained. Group labels and
equivalent orientations are not canonicalized at the reaction level.

## Fixed validation packet

`validation/rdkit-9732-2026-10-24/` contains 20 authored positive cases in both
writer modes (40 runs), plus 22 expected refusals. Cases cover the upstream
example, participant reordering, ABS/AND/OR groups, all three roles, shared and
independent groups, maps, isotopes, explicit H, dative arrows and group-number
boundaries. Native tests also renumber atoms while retaining their group/map
identity and CIP configuration.

An independent RDKit 2026.09.1 observer reparses captured output and checks
group kind, member partitions, role, rooted atom identity and CIP configuration.
It ignores equivalent numeric group renaming and uniform inversion of an
AND/OR group, with six observer controls. Candidate Python passes 40/40.
In the same fixed packet RDKit's writer has 29 matches, eight atom-reference
differences and three differences in cross-participant grouping. The latter
are recorded separately from the wrong-atom issue because this API declares
a global correlation convention. These are local fixed-input representation
results, not toolkit-wide stereo correctness, speed or published-release proof.
The retained prior runtime lacks the new API; it is not counted as failed
chemistry. Raw captures and source/runtime/fixture hashes accompany the report.

The native, Python and WASM regressions run in the existing integration-test
discovery paths. The standalone capture/observer can be reproduced with
`scripts/validation/validate_reaction_cxsmiles.py --help`.
