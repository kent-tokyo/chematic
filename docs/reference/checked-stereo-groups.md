# Checked enhanced stereo-group construction

RDKit [#9705](https://github.com/rdkit/rdkit/issues/9705) describes how losing a
bond member changes the correlation between atom and double-bond
configurations. A mixed group must retain both kinds of members to represent
that correlation.

chematic already stores bond members separately from its public `StereoGroup`
struct and reads/writes V3000 STE/STEB collections. This candidate adds a
checked construction path in Rust and Python, plus Python inspection. It
reuses that storage and interchange contract.

```python
import chematic

mol = chematic.from_smiles("C/C=C/[C@H](O)CC")
coords = mol.rdkit_2d_coords()
# (kind, number, atom indices, bond indices), all native indices zero-based
mixed = mol.with_stereo_groups([("and", 1, [3], [1])])
assert mixed.stereo_groups == [("and", 1, [3], [1])]
assert mol.stereo_groups == []

block = mixed.to_mol_v3000(coords, strict=True)
restored = chematic.from_mol_v3000(block)
assert restored.stereo_groups == mixed.stereo_groups
```

The mixed group is written as `STERAC1 ATOMS=(1 4)` and
`STEBRAC1 BONDS=(1 2)`. The existing reader merges these into one group.
A bond-only group has an empty atom list and a nonempty bond list; it is never
serialized as an empty atom collection.

## API contract

`Mol.stereo_groups` returns a detached list of
`(kind, number, atom_indices, bond_indices)` tuples. Mutating the returned
lists does not modify the molecule. The index order is the stored native graph
order, available through `atom_table` and `bond_table`; it is not a promise
about a different toolkit's atom or bond order.

`Mol.with_stereo_groups(groups)` replaces the full group list on an independent
copy. Typed molecule properties and names are copied independently. The indexed
graph, existing stereo declarations/frames, atom tags and other graph metadata
remain intact. Derived caches start empty. Supplying `[]` clears group metadata
on the copy without changing the source.

- Kinds are `abs`, `and`, and `or`.
- ABS uses number 0. AND/OR require a positive 32-bit integer.
- Atom and bond indices must be nonnegative 32-bit integers within the stored
  graph's bounds. Python rejects bool, float, text and negative values.
- A group must have at least one atom or bond member.
- Each `(kind, number)` key occurs once. AND and OR may use the same number.
- Each atom or bond belongs to at most one group. Repeated members within a
  single group are deduplicated in first-seen order.
- Invalid input raises `ValueError` without changing the source. A structurally
  malformed Python tuple/list may raise `TypeError` during argument conversion.

In Rust, use:

```rust,ignore
let grouped = mol.with_checked_stereo_groups(vec![(
    chematic_core::StereoGroup::new(
        chematic_core::StereoGroupKind::And(1),
        vec![chematic_core::AtomIdx(3)],
    ),
    vec![chematic_core::BondIdx(1)],
)])?;
```

`StereoGroupError` identifies invalid references, empty groups, invalid numbers,
repeated group keys and overlapping membership. Existing unchecked metadata
setters retain their signatures and behavior.

## Interchange boundary

[BIOVIA CTfile Chemistry 2026](https://www.3ds.com/assets/invest/2026-01/biovia_chemistry_2026_ctfile_formats.pdf)
defines separate atom and bond collection tokens for absolute, AND and OR
stereochemistry. The shared-number interpretation of mixed STE/STEB records
is the existing chematic convention and the proposal in #9705. The published
examples use distinct atom/bond group numbers and do not explicitly resolve
that correlation question. This candidate does not establish interoperability
of mixed groups with BIOVIA products.

This API validates metadata structure. It does not perceive stereogenicity,
verify that a selected bond is a valid stereogenic axis/double bond, assign
new geometry, invert groups, enumerate isomers, canonicalize enhanced stereo
identity, or perform meso reduction. Atropisomer bond membership can be stored;
axis configuration is not created or validated by this API.

Plain SMILES gives the drawn representative and does not encode enhanced
group membership. The existing RDKit-compatible CXSMILES method continues to
reject enhanced groups. Use the V3000 path above to retain bond members; no
new CXSMILES group syntax is introduced here.

## Validation and provenance

`validation/rdkit_issues/stereo/checked-groups.json` is original project data
under MIT OR Apache-2.0: two small dienes and one mixed alcohol, with no copied
upstream attachment or BIOVIA example structure. It records minimization and
a pinned RDKit 2026.09.1 diagnostic/geometry oracle.

The validator checks ABS, AND and OR over all three fixtures, each with 24
deterministic atom numberings: 216 V3000 round trips. Exact memberships and
collection tokens use source-authored expectations. RDKit InChI checks the
drawn configuration independently, without treating a drawn representative
as the identity of an enhanced group.

On RDKit 2026.09.1, the separate diagnostic reproduces loss of bond members
from bond-only and mixed AND groups during V3000 export/read. Mixed-group
canonicalization also drops its bond member. Those behaviors are diagnostics,
not the membership oracle for this implementation.

```sh
cargo test -p chematic-core --lib --test checked_stereo_groups --offline --quiet
cargo test -p chematic-mol --lib --test v3000_semantic_interop --offline --quiet
python -m pytest crates/chematic-py/tests/test_stereo_groups.py \
  crates/chematic-py/tests/test_pickle.py \
  crates/chematic-py/tests/test_stereo_candidate_budget.py -q
python scripts/validation/validate_checked_stereo_groups.py
```

Tests cover graph metadata and warm/cold caches, detached inspection, atomic
validation errors, independent typed properties, pickle, permutations, bond-only
and mixed V3000 groups, and the shipped typing declarations.

Local validation passed 936 Rust tests (including nine new checked-group tests)
and 70 Python tests (including 44 new API tests), plus 216 reference cases.
Clippy with warnings denied, explicit include-file formatting, refactor budget,
binding inventory, product naming, and strict MkDocs also passed.

This is an uncommitted local source candidate built on the previous candidates.
It is not an upstream RDKit fix or a published chematic package.
