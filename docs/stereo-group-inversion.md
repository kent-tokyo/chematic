# Full alternatives of relative stereo groups

RDKit [#9705](https://github.com/rdkit/rdkit/issues/9705) proposes treating an
atom/bond group's alternative as the simultaneous inversion of all of its
members. Inverting only the tetrahedral center of a mixed group loses the
intended correlation with its geometric bond members.

This candidate adds Rust `invert_stereo_group` and Python
`Mol.invert_stereo_group(index)`. It builds on the
[checked group construction API](checked-stereo-groups.md) and reuses the
existing joint geometric-direction solver.

```python
import chematic

source = chematic.from_smiles("C/C=C/[C@H](O)CC")
source = source.with_stereo_groups([("and", 1, [3], [1])])
alternative = source.invert_stereo_group(0)
assert alternative.stereo_groups == source.stereo_groups
assert alternative.invert_stereo_group(0).smiles == source.smiles
```

Here both the tetrahedral declaration and the double-bond geometry change.
The returned `smiles` represents the drawn configuration; it does not encode
the enhanced group's identity. Membership remains available through
`stereo_groups` and the existing V3000 interchange path.

## Contract

The group index is zero-based and follows `Mol.stereo_groups` / Rust
`Molecule::stereo_groups`. Only AND/OR groups have this relative alternative.
ABS is rejected. The result is an independent molecule; failures leave the
source unchanged.

- Every selected atom must have a specified tetrahedral tag and a complete
  four-ligand reference frame matching its graph neighbors and implicit H.
  Its tag is flipped against that unchanged frame.
- Every selected bond must be a supported double bond with fully specified
  relative geometric directions. One reference side is inverted.
- All geometric facts are solved together, retaining unselected geometry.
  A conjugated single bond shared by two double bonds may need a different
  direction token so one double bond can change while the other stays fixed.
- Group kinds, numbers, membership, atom/bond indices, connectivity, atom tags
  and stereo reference frames are preserved. Legacy duplicate atom memberships
  are retained verbatim. Python typed properties and
  names are copied independently.
- Derived caches start empty. All atom CIP annotations are cleared, since
  changing a nearby stereo assignment can change priority calculations.
  No new CIP labels are assigned by this operation.

In Rust:

```rust,ignore
let alternative = chematic_chem::invert_stereo_group(&molecule, 0)?;
```

`StereoGroupInversionError` distinguishes invalid group metadata/index,
absolute groups, unsupported atom/bond members, incomplete tetrahedral frames,
unspecified geometry, small-ring constraints, conflicting directions and
representations the direction solver cannot safely rewrite.

Python raises `ValueError` for those errors, and rejects bool, float, text,
negative and oversized group indices. The group itself is also validated
using the previous checked construction contract.

## Scope

The [BIOVIA CTfile Chemistry 2026 definitions](https://www.3ds.com/assets/invest/2026-01/biovia_chemistry_2026_ctfile_formats.pdf)
describe separate atom and geometric/axial bond collections. The simultaneous
alternative and shared-number mixed-group convention here follow the #9705
proposal and existing chematic storage; this is not a verified claim about
BIOVIA product interoperability.

This operation supports tetrahedral atoms and geometric double bonds. It
rejects axial/non-double bond members, non-tetrahedral atoms, incomplete or
unspecified configurations, double bonds in rings smaller than eight atoms,
and conflicting/unsupported direction representations. Eight-membered rings
are covered, without a claim of 3D feasibility.

It does not detect stereogenicity, count distinct stereoisomers, generate
coordinates, canonicalize enhanced groups, perform meso reduction or add
CXSMILES syntax. A full alternative may be the same compound as the source;
the symmetric fixture demonstrates that case. Group kinds/membership remain
unchanged rather than being converted to ABS.

Coordinates belonging to the original drawing are not updated. Regenerate a
layout reflecting the returned configuration before coordinate-based export;
do not reuse the original drawing as if it had also been inverted.

## Validation

`validation/rdkit_issues/stereo/group-inversion.json` contains original small
fixtures under MIT OR Apache-2.0, with provenance and expected invariants. No
upstream attachments or proprietary structures were copied. The result JSON
records normalized fixture, source and compiled-extension hashes.

The independent reference is pinned to RDKit 2026.09.1, legacy stereo
perception and `removeHs=False`. It directly inverts selected atom tags and
swaps selected bond stereo, then regenerates neighbor directions. It does not
use RDKit's lossy enhanced-group canonicalizer as the reference. Atom
attributes and the complete indexed graph are checked before comparing.
Ring-closure endpoint orientation is normalized as an undirected edge, and
selected native bonds map to reference bonds by their atom endpoints.

Eleven fixtures cover atom-only, bond-only and mixed groups, shared markers,
explicit H, an unchanged center, a symmetric full alternative and an
eight-membered ring. AND/OR and 24 deterministic atom numberings per fixture
produce 528 comparisons of drawn alternatives, membership, graph identity and
two-application restoration, with zero mismatches.

RDKit's mixed-group canonicalization bond loss is reproduced separately as an
upstream diagnostic. InChI validates the drawn configurations, not the identity
of an enhanced group or the number of its physically distinct alternatives.

```sh
cargo test -p chematic-chem --lib --test stereo_group_inversion --offline --quiet
python -m pytest crates/chematic-py/tests/test_stereo_group_inversion.py \
  crates/chematic-py/tests/test_stereo_groups.py \
  crates/chematic-py/tests/test_stereo_candidate_budget.py -q
python scripts/validate_stereo_group_inversion.py
```

Local validation passed 963 Rust tests, with one pre-existing real-budget VF2
calibration ignored, and 92 Python tests. This includes 15 new Rust contracts
and 31 new Python API tests. Clippy with warnings denied, formatting, refactor
budget, product naming, binding inventory and strict MkDocs also passed.

This is a local source candidate built on the previous checked-group candidate.
It does not resolve all of RDKit #9705 or constitute package publication.
