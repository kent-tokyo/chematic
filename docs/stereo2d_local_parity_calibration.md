# P1-S1a-core: local parity sign convention, measured against RDKit

`chematic_perception::local_parity_from_wedges` (`crates/chematic-perception/src/stereo2d_local.rs`)
computes `Atom.chirality` + `stereo_neighbor_order` from wedge bonds and 2D
coordinates, with **zero CIP involvement**. This note records how its sign
convention was pinned down, so a future change doesn't have to re-derive it
by analogy (the failure mode this note exists to prevent).

## Why not just reason it out on paper

`chematic_core::Chirality`'s own doc comment already states the convention
("`@` = counterclockwise looking from the first neighbor"), but two things
are *not* obvious from that sentence alone:

1. Whether chematic's chosen "first neighbor" (`mol.neighbors(center)`
   insertion order) is even in the same frame as whatever RDKit's
   `CHI_TETRAHEDRAL_CW`/`CCW` tag is relative to.
2. What the right formula is for a 3-heavy+implicit-H center, where there's
   no explicit position for the 4th (H) neighbor.

Both were resolved by direct measurement against real RDKit output on MOL
V2000 fixtures, not by extending the existing CIP-based `assign_rs`'s
lowest-priority-behind convention by analogy — a sign bug picked to match a
handful of fixtures is indistinguishable from a correct one until it's wrong
on a fixture you didn't test.

## Frame alignment (measured, not assumed)

For every fixture, chematic's `mol.neighbors(center)` order (bond
declaration order) and RDKit's `atom.GetNeighbors()` order were dumped side
by side from the *same* MOL block text. They were identical in every case,
including after reversing the bond block. Both engines index a stereocenter's
neighbors by bond-declaration order, so raw parity values are directly
comparable without any permutation translation.

## 4 explicit neighbors

Formula: apex = `order[0]`, viewed = `order[1..3]`.
`vol = signed_volume(pos[order[1]], pos[order[2]], pos[order[3]], pos[order[0]])`.

Measured on 5 fixtures (wedge on first bond, wedge on last bond, hash
inversion, full bond-order reversal, two simultaneous wedges):

| vol sign | RDKit tag |
|---|---|
| negative | `CHI_TETRAHEDRAL_CCW` |
| positive | `CHI_TETRAHEDRAL_CW` |

Consistent across a hash-vs-wedge sign flip and a full 4-element bond-order
reversal (an even permutation — sign correctly preserved in both engines).
The two-simultaneous-wedges fixture disagreed with a naive sign read (RDKit
still emits *some* tag via its own dual-volume fallback logic) — this is the
concrete evidence behind the "contradictory wedges → no assignment" rule:
`local_parity_from_wedges` refuses rather than picks a side.

## 3 heavy neighbors + implicit H

RDKit's own `atomChiralTypeFromBondDirPseudo3D` (nNbrs==3 path, read from
source at the pinned commit below) does not synthesize a position for the
missing H at all — it computes a triple product of the three *real* bond
vectors **from the center atom**. `signed_volume(p1, p2, p3, p4)` is already
exactly `p4`'s the pivot; setting `p4 = center` reproduces RDKit's formula
verbatim:

`vol = signed_volume(pos[order[0]], pos[order[1]], pos[order[2]], center_pos)`

Measured on 4 fixtures (wedge on first bond, hash inversion, wedge on a
different bond, full 3-element bond-order reversal):

| vol sign | RDKit tag |
|---|---|
| negative | `CHI_TETRAHEDRAL_CW` |
| positive | `CHI_TETRAHEDRAL_CCW` |

Note the mapping is the **opposite** polarity from the 4-neighbor case —
expected, since the pivot moved from "first neighbor" to "center", a
genuinely different geometric setup, not a bug.

This was cross-checked a second, independent way: RDKit's own non-canonical,
rooted SMILES output for the wedge fixture is `[C@H](F)(Cl)Br` (order
`[H, F, Cl, Br]`, i.e. H *first*, since the bracket atom is the traversal
root with no preceding atom — matching chematic's parser's own
`begin_stereo_if_chiral` convention for a root atom). Moving H from the front
to the back of a 4-slot list (chematic's chosen convention: H goes *last* in
`stereo_neighbor_order` for the 3-heavy case) is a 4-cycle, i.e. 3
transpositions, an odd permutation — so the symbol must flip: `@` → `@@`.
That predicts `Chirality::Clockwise` for the same wedge fixture, which is
exactly what the direct RDKit-tag comparison above also gives. Two
independent derivations agreeing is the actual basis for confidence here, not
either one alone.

## Contradictory-wedge and degenerate-geometry rejection

`local_parity_from_wedges` returns `None` when more than one neighbor bond
carries a wedge/hash from the same center (see the two-simultaneous-wedges
finding above), or when the computed volume is within `1e-6` of zero
(coplanar/degenerate, including the "no wedge at all" case).

## Reproducing this measurement

The calibration was done with a throwaway script (not checked in) that
builds MOL V2000 block text by hand, feeds the identical text to RDKit
(`.venv/bin/python`, rdkit==2026.03.3) and to hand-derived arithmetic
mirroring `signed_volume`/`wedge_z`, and compares. RDKit source referenced at
commit `8afba32ec539dcb2369bc84549d802aca3f7eb39` (same pin used by
`docs/stereo2d_reader_integration_rfc.md`), specifically
`Code/GraphMol/Chirality.cpp`'s `atomChiralTypeFromBondDirPseudo3D`.

## Scope note

This module never calls `chematic_perception::cip_priority` and never writes
`Atom.cip_code`. It is not called from any reader yet, and the SMILES writer
is unchanged — both are deliberately out of scope for this step (see
`docs/stereo2d_reader_integration_rfc.md` for the full integration-boundary
design and the reason the writer's existing standalone-wedge bug — it emits
a meaningless `/`/`\` token for a wedge bond with no adjacent double bond,
independently re-encountered while validating this module's output — must be
fixed before any reader wires this in automatically).
