# Correlated full alternatives of stereo groups

RDKit [#9705](https://github.com/rdkit/rdkit/issues/9705) describes a mixed
relative group whose atom centers and geometric double bonds must change
together. In a group containing one tetrahedral center and one double bond,
the drawn configuration and its full alternative might be `(E, R)` and
`(Z, S)`. Flipping only the center introduces `(E, S)`, which was not one
of the two correlated configurations.

This candidate builds on [checked group construction](checked-stereo-groups.md)
and [full group inversion](stereo-group-inversion.md). It expands all AND/OR
group choices with explicit bounds, keeping every atom/bond member correlated.
ABS groups and ungrouped configurations stay fixed. It provides a separate
API from `enumerate_stereoisomers`, which assigns unspecified carbon centers.
That existing API's behavior is unchanged.

```python
import chematic

source = chematic.from_smiles("C/C=C/[C@H](F)C[C@@H](Cl)/C=C/Br")
source = source.with_stereo_groups([
    ("and", 1, [3], [1]),
    ("or", 2, [6], [8]),
])
indices, representatives = source.enumerate_stereo_group_alternatives()
assert indices == [0, 1]
assert len(representatives) == 4
assert all(m.stereo_groups == source.stereo_groups for m in representatives)
```

The return value is `(relative_group_indices, representatives)`.
Output index is an inversion mask. Bit `k` selects full inversion of original
stored group `relative_group_indices[k]`; zero keeps that group's drawn
configuration. Relative groups follow stored order, skipping ABS. Group
numbers do not determine bit positions.

| Output index | Group 0 | Group 1 |
| --- | --- | --- |
| 0 | drawn | drawn |
| 1 | fully inverted | drawn |
| 2 | drawn | fully inverted |
| 3 | fully inverted | fully inverted |

AND and OR use the same full-alternative operation, while their different
sample meanings remain encoded in the original group kinds. Each output
retains that metadata; it is a drawn representative of the original enhanced
record, rather than a conversion to a pure-isomer record. Its plain SMILES
represents the drawing without enhanced-group identity. Use the retained
groups and existing V3000 interchange when that identity is needed.

## Contract and limits

Rust `enumerate_stereo_group_alternatives(molecule, max_alternatives)` returns
`Result<StereoGroupAlternatives, StereoGroupAlternativesError>`. The result
holds `relative_group_indices` and `molecules`. Python
`Mol.enumerate_stereo_group_alternatives(max_alternatives=64)` returns detached
lists and raises `ValueError` for invalid or unsupported inputs.

- At most **six relative groups**, **64 representatives**, **1,024 atoms** and
  **4,096 bonds**. The caller may lower `max_alternatives` to any integer in
  `1..=64`. It must fit the complete `2^n` family; there is no truncation.
- Budget and graph-size checks precede generation. Group references and
  disjoint membership are validated before graph access. Output-vector
  reservation is fallible. A generation error returns no partial family.
- Every relative atom member must have specified tetrahedral chirality and
  a complete reference frame. Every relative bond member must have supported,
  specified double-bond geometry. The underlying joint direction solver
  preserves unselected geometric facts, including conjugated double bonds
  that share a directional marker.
- ABS members remain fixed. Unspecified and ungrouped centers are not assigned
  or inverted. With no relative groups there is one independent unchanged
  representative, including for an empty Rust molecular graph.
- Atom/bond indices, connectivity, tags, reference frames, raw group membership
  and Python typed molecule properties are retained. Legacy duplicate atom
  memberships remain in metadata, but a repeated member flips only once.
- Every representative, including mask zero, starts with cold derived caches
  and cleared stored atom CIP annotations. Priority labels can depend on a
  neighboring configuration and must be recomputed. Input metadata/caches
  remain unchanged.

Errors distinguish invalid limits, atom/bond/group/count bounds, invalid
stored groups, allocation failure and the original stored group index of a
failed full inversion. Axial/non-tetrahedral members, unspecified geometry,
incomplete frames, small-ring double bonds below eight atoms and conflicting
or unrepresentable directions are rejected by the reused inversion boundary.

The API intentionally retains symmetry-equivalent masks: **the number of
representatives is not a distinct stereoisomer or compound count**. It does
not perform meso reduction, canonicalization, identity registration,
coordinate generation or 3D feasibility checks. Original drawing coordinates
are not updated; regenerate appropriate coordinates before geometry export.

Same-kind/same-number atom/bond grouping follows the #9705 proposal and
chematic's existing storage contract. The issue itself calls for confirming
this convention with BIOVIA; this candidate does not claim such confirmation
or invent a bond-indexed CXSMILES extension.

## Validation

The pinned independent reference is RDKit **2026.09.1**, legacy stereo
perception, explicit H retained (`removeHs=False`). The reference directly
inverts the selected atom tags and swaps selected double-bond E/Z states,
then clears and regenerates neighbor directions. It does not use RDKit's
mixed-group enumerator or enhanced-group canonicalizer as expected truth.
Complete indexed atom attributes and undirected bond graphs are compared
before the configuration comparisons. Native/reference bond members are
matched by atom endpoints, accounting for ring-closure endpoint orientation.

Six original synthetic fixtures cover two mixed groups, conjugated bond-only
groups, atom groups with ABS, explicit H, a symmetric diol, and three
independent disconnected groups. Twenty-four atom numberings and stored/reverse
group orders produce **1,056 representative comparisons**, with zero mismatches
in drawn InChI, grouping, graph identity, source immutability or mask restoration.
Every mask remains present even where multiple masks have the same drawn InChI.
InChI is used only for the drawn configuration, not enhanced-record identity.

A separate diagnostic reproduces the upstream gap: RDKit's
`EnumerateStereoisomers(..., onlyUnassigned=True, unique=False)` produces four
outputs for the two mixed groups while leaving both grouped double bonds
fixed. Its configuration set differs from the correlated four-mask reference.
Atom-only group enumeration with an ABS center is a separate passing control.

Local validation passed:

- `chematic-chem`: 948 library tests, 14 new family contracts and 15 existing
  group-inversion contracts; one pre-existing library test remains ignored.
- Python: 107 tests across family expansion, full inversion and checked groups.
- Strict Clippy, formatting, product naming, binding inventory, refactor budget
  and strict MkDocs build. No refactor allowances were added.

The fixture manifest is
`validation/rdkit_issues/stereo/group-alternatives.json`, under
MIT OR Apache-2.0 with provenance, minimization, pinned oracle and normalized
payload hash. No upstream attachments or proprietary corpus were copied.
Regenerate the source/fixture/compiled-extension evidence with:

```sh
python scripts/validation/validate_stereo_group_alternatives.py
```

The implementation is a local, uncommitted, unpublished source candidate in
`/private/tmp/chematic-group-alternatives`, preserving prior candidates through
TURBOMOLE. The feature patch is
`validation/rdkit_issues/stereo-group-alternatives-2026-10-15.patch` and the
record is `validation/results/rdkit-stereo-group-alternatives-2026-10-15.json`.
The patch requires the preceding checked-group and full-inversion candidates.
This implements correlated candidate expansion associated with #9705; full
canonicalization, axes, meso reduction and CXSMILES remain outside its scope.
