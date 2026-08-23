# RFC: connectivity-ordered ring/chain placement for `dg::generate_coords` (Phase 0 — design + fixtures, no implementation)

Issue #256 (`generate_coords` places every ring in a component before any
chain atom, breaking chain-bridged ring islands like bibenzyl) and #277
(closed as a duplicate, its 17-molecule regression population absorbed
into #256) need a staged fix. This document is **Phase 0 only**: design +
test fixtures, no behavior change. Nothing under `crates/*/src/**` outside
`crates/chematic-3d/src/dg.rs`'s `#[cfg(test)]` module is touched by this
round, and nothing in that module changes production behavior — every new
test calls the existing, unmodified `generate_coords`.

## 1. The two known defects in `place_rings`

`place_component` (`crates/chematic-3d/src/dg.rs:195`) calls `place_rings`
(`:331`) for every ring in a connected component, THEN `dfs_place` (`:584`)
walks chain atoms from whatever got placed. `place_rings` branches on live
`placed` state per ring:

- `shared_atoms.len() >= 2` — fuse to a shared edge. **Issue #255**: the
  new ring's center is placed at a fixed `ring_cy = (p0.y + p1.y) / 2.0 +
  r`, i.e. always `+y` from the fusion bond's midpoint. Correct only when
  `+y` happens to be the direction away from the rest of the structure.
- `shared_atoms.len() == 1` — fuse to one shared atom (spiro). Unaffected
  by either issue below; confirmed still correct under the redesign this
  document sketches (see §4).
- `shared_atoms.is_empty()` — either a direct bond to an already-placed
  atom exists (biphenyl/terphenyl's "new island" case, fixed correctly in
  PR #253 via the specific ring's own centroid), or it doesn't, and the
  code falls back to a fixed `x_offset + 5.0` anchor with no real
  constraint on the eventual bond length. **Issue #256**: this is exactly
  what happens when a ring is reachable only through not-yet-placed chain
  atoms (bibenzyl's `-CH2CH2-` bridge) — `place_rings` runs entirely
  before `dfs_place` ever walks a chain atom, so the bridge's real
  endpoint position doesn't exist yet to anchor against.

**These are two independent defects in the same function**, not the same
bug. #255 is referenced throughout this document because a
connectivity-ordered rewrite touches the exact same code path and could
easily leave #255 unfixed or (per its own history below) regress it
further — but **fixing #255 is not authorized work for this PR**. Only
#256's fixtures and design are in scope; #255's fixtures are included
because the redesign will pass through its exact failure mode regardless,
and because of its specific regression history:

> #255's own comment thread records a reverted overnight fix attempt: a
> 2-point rigid rotation of the ring's own template onto the real
> shared-edge positions. Result: anthracene fixed exactly (all bonds at
> 1.4000 Å), but phenanthrene got two EXACT atom coincidences plus a 3.7 Å
> stretch, and pyrene got one exact atom coincidence — both worse than the
> original distortion. Diagnosed as a missing reflection/winding-handedness
> step: a pure rotation assumes the ring's own listed atom order's
> chirality matches the real-space winding needed, which held for
> anthracene's specific fusion geometry but not phenanthrene's angular or
> pyrene's multi-ring fusion.

## 2. Current characterization (measured this round, not assumed)

Every fused-ring or chain-bridged case measured against **today's
unmodified** `generate_coords` turned out to be broken — not a mix of
sound and broken cases as originally guessed before measuring:

| Molecule | Topology | Worst bonded-pair distance | Ideal | New test |
|---|---|---|---|---|
| naphthalene | 2-ring, 1 shared edge | 2.2644 Å | ~1.40 Å | `generate_coords_naphthalene_fusion_seam_known_broken` |
| quinoline | 2-ring fused heterocycle | 2.2644 Å (identical topology to naphthalene) | ~1.40 Å | `generate_coords_quinoline_fused_heterocycle_known_broken` |
| phenanthrene | 3-ring, angular fusion | 3.3856 Å | ~1.40 Å | `generate_coords_phenanthrene_angular_fusion_known_broken` |
| pyrene | 4-ring fusion | 3.8974 Å | ~1.40 Å | `generate_coords_pyrene_multi_ring_fusion_known_broken` |
| diphenylmethane | ring-CH2-ring (bridge 1) | 8.0738 Å | ~1.5 Å | `generate_coords_diphenylmethane_chain_bridge_length_1_known_broken` |
| bibenzyl | ring-(CH2)2-ring (bridge 2) | 8.7358 Å | ~1.5 Å | pre-existing (`generate_coords_bibenzyl_chain_bridged_ring_islands_known_broken`) |
| 1,3-diphenylpropane | ring-(CH2)3-ring (bridge 3) | 8.9731 Å | ~1.5 Å | `generate_coords_diphenylpropane_chain_bridge_length_3_known_broken` |
| 1,4-diphenylbutane | ring-(CH2)4-ring (bridge 4) | 8.7060 Å | ~1.5 Å | `generate_coords_diphenylbutane_chain_bridge_length_4_known_broken` |

An earlier framing of this issue (a comment on #255, 2026-08-11) described
naphthalene as "unaffected, fine" — that referred specifically to the
*reverted fix attempt's* regression pattern (naphthalene didn't get worse
under that particular broken fix), not to naphthalene's status under
today's unmodified code, which this round's direct measurement shows is
itself already distorted by the same `shared_atoms.len() >= 2` bug as
anthracene. Corrected here since it directly affects fixture scope.

anthracene (pre-existing test, `generate_coords_anthracene_terminal_rings_not_superimposed`)
and the new-island cases (biphenyl, terphenyl, 3-phenylpyridine, all
pre-existing) remain as they were — anthracene deliberately un-checked at
the bonded-pair level (see that test's own comment), new-island cases
sound (`assert_bonded_pairs_sane`). Single-ring-plus-chain
(`generate_coords_ring_with_tail_substituent_all_atoms_placed`,
ibuprofen-shaped) and pure-chain, no-ring
(`generate_coords_propane_linear`) cases were already covered before this
round — no new fixture needed for either.

## 3. The proven precedent: `chematic-depict/layout.rs`

Issue #347 (2D depiction layout) hit the same class of problem —
independently-placed pieces (rings, chain runs) need to anchor to their
real, connectivity-discovered attachment point, not a blind fixed
direction — and its fix is directly relevant here in two separable parts.

### 3a. Traversal structure (portable)

`grow_layout` (`crates/chematic-depict/src/layout.rs:701`) + `place_ring_system`
(`:328`): a single `VecDeque<AtomIdx>` worklist, seeded with whatever's
already placed (a seed ring system, or an isolated chain start). Pop an
already-placed atom; for each of its unplaced neighbors, compute an
outgoing direction and walk a zigzag chain (`dfs_zigzag`) that itself
detects and places any newly-discovered ring system along the way (via
`place_ring_system`, anchored to its real entry atom/position), then
pushes the newly-ring-placed atoms back onto the worklist. Rings and
chains are discovered and placed in true connectivity order — never "all
rings, then all chains." **This structure is what `place_component`
needs; nothing 2D-specific about the worklist/discovery logic itself.**

### 3b. Ring-fusion center + winding (verified working, not assumed)

`place_ring_anchored` (`:467`) is the 2D analog of the exact geometric
problem #255 describes, and it is demonstrably correct on every topology
that broke the reverted 3D fix attempt. Verified this round via
`mol.depict_data()`:

```
naphthalene:   n_bonds=11  all bonds == BOND_LEN (40.0), zero outliers
phenanthrene:  n_bonds=16  all bonds == BOND_LEN (40.0), zero outliers
pyrene:        n_bonds=19  all bonds == BOND_LEN (40.0), zero outliers
anthracene:    n_bonds=16  all bonds == BOND_LEN (40.0), zero outliers
```

Its technique, and why it avoids the reverted attempt's winding bug:

1. Compute **both** candidate ring centers on the shared edge's
   perpendicular bisector (`cand1`, `cand2`), and choose whichever is
   **farther from the centroid of every currently-placed atom in the
   molecule** (`centroid_of_placed`, not just this ring's own two anchor
   atoms — the two anchors' own centroid is always exactly the edge
   midpoint, an uninformative tie).
2. Derive the winding direction (clockwise vs. counterclockwise
   placement of the ring's remaining atoms) from the **real, measured
   signed angle** between the two anchor points around the *chosen*
   center (`delta = angle_to_a2 - angle_to_a1`, sign taken directly from
   that measurement) — never from the ring's own listed atom order's
   implicit chirality.

The reverted 3D attempt did a pure 2-point rigid rotation instead — no
mechanism to choose a reflection, so it silently trusted the ring
template's own winding to already match real space. That assumption held
for anthracene's specific fusion geometry and failed for phenanthrene's
angular and pyrene's multi-ring fusion. `place_ring_anchored`'s
center-then-empirical-winding approach has no such assumption: winding
falls out of measurement, not of trusting a fixed template.

## 4. What must be independently redesigned for 3D

Not resolved here — this section frames each as an open question with a
recommendation, for Phase 1 to implement.

**Ring-plane orientation.** 2D has exactly one binary choice per fusion
(which side of the shared edge/line). 3D adds a full normal-vector degree
of freedom per ring: given an entry bond direction, the ring's plane can
be tilted around that bond axis arbitrarily. Candidates:
  (a) orient the plane to contain the entry bond and the parent atom's
      other placed neighbor (mimics a "flat" extension, matches sp2/sp3
      local geometry reasonably for a single substituent bond, degenerate
      when the parent has no other placed neighbor — the common case for
      the very first ring in a component);
  (b) pick a fixed reference normal (e.g. the XY plane, matching this
      file's existing single-ring convention) and only rotate around the
      entry bond to satisfy the 2D-style center/winding constraint,
      accepting that fused ring systems end up coplanar even when real
      3D fused systems (e.g. non-planar aliphatic fusions) would not be;
  (c) derive the normal from three real anchor points once at least 3
      atoms are already placed nearby (over-constrained but self-
      correcting, degenerate with fewer than 3).
  **Recommendation**: (b) for Phase 1's first cut — matches this file's
  existing behavior for every already-working case (new-island,
  single-shared-atom), keeps the change scoped to fixing #256/#255's
  specific center/winding defects rather than introducing a new
  out-of-plane degree of freedom in the same pass. Revisit only if
  Phase 2's differential evaluation shows real fused/aliphatic systems
  need non-coplanar fusion.

**Dihedral / parent-grandparent maintenance.** `dfs_place`'s current
scheme (`:608-618`) rotates by a fixed `bend_angle` from the incoming-bond
direction, then spaces multiple unplaced neighbors of the *same* atom
120° apart around that axis — it does **not** reference a grandparent
atom's position at all. This is a real gap, not an existing mechanism to
port: state this plainly rather than imply grandparent-referenced
dihedral maintenance already exists somewhere to extend. Whether the new
engine should add real grandparent-referenced dihedral placement (a
larger change, closer to true ETKDG-style staggering) or keep the current
scheme unchanged and scope this rewrite strictly to the ring-anchor
defects (#255/#256) is an open decision for Phase 1 — recommend the
latter for a first cut, given `generate_coords`'s own doc comment already
disclaims "not a full distance-geometry solver."

**Fused/spiro/bridged ring handling.** Spiro (`shared_atoms.len() == 1`)
is untouched by both defects and should remain behaviorally identical
under the new traversal — confirm via the existing
`generate_coords_spiro_ring_adjacency_unaffected` test, which this Phase 0
round does not modify. Bridged bicyclics (two rings sharing more than one
non-adjacent atom, e.g. norbornane) are not currently covered by any
fixture in this file and are out of scope for Phase 0 — flagged as a gap
for Phase 1's own fixture additions if it touches that case.

**Interaction with declared stereo.** `dg.rs` has no stereo awareness
today (`generate_coords`/`generate_coords_etkdg` are both legacy,
pre-dating this project's stereo-aware `distance_geometry_v2`/
`embed_pipeline_v2` path). Recommend the new engine stays stereo-blind
too, matching its current scope exactly — adding stereo constraints here
would be new functionality, not a fix, and production callers needing
stereo-correct 3D already use `embed_pipeline_v2`, which never calls this
code (confirmed: `distance_geometry_v2.rs` and `pipeline_v2.rs` have no
reference to `dg::generate_coords`).

**Atom-order reproducibility.** The current code already collects
`component_set`/`shared_atoms` into `Vec`/filters over slices, not
`HashSet` iteration order — no known nondeterminism today. The new
worklist-based traversal must preserve this: seed order and per-atom
neighbor visitation order must come from `mol.neighbors()`'s own stable
order (as `grow_layout` does, sorting an explicit `Vec` rather than
draining a hash-ordered collection), not from any `HashMap`/`HashSet`
iteration.

**Re-anchoring an already-placed ring.** Recommend **no** — once a ring
is placed by the traversal, it is never rigidly moved again, even if a
later-discovered chain connection would have suggested a different
anchor. This matches `place_ring_system`'s own one-shot-placement
convention in chematic-depict (every `place_first_ring_anchored`/
`place_ring_anchored` call fills in only currently-unplaced atoms; it never
revisits a placed one) and avoids a much harder "re-derive downstream
placements after moving an anchor" problem that current `dfs_place`
doesn't have to solve either.

## 5. Staged plan (restated; only Phase 0 is this PR)

- **Phase 0 (this PR)**: design (this document) + fixtures (§2's new
  tests). No behavior change.
- **Phase 1**: new, private, additive function (e.g.
  `generate_coords_connectivity_ordered`) implementing §3a's worklist
  structure with §3b's center/winding technique adapted per §4's
  recommendations. Existing `generate_coords` untouched — this is a
  parallel implementation, never an in-place rewrite.
- **Phase 2**: differential evaluation against every fixture in §2, plus
  #277's original 17-molecule regression population
  (`chembl_tier_b_0003, _0008, _0010, _0021, _0025, _0026, _0027, _0035,
  _0065, _0066, _0095, _0143, _0146, _0147, _0153, _0156, _0157`;
  `validation/results/pipeline_v2_vs_rdkit_3point_paired_diff_summary.json`,
  key `legacy_etkdg_soundness_diff_pre_2b_to_v0_12_0`). Check that BOTH the
  gross-clash-count improvement PR #253 already achieved for that
  population AND the new elongation defect it introduced improve together
  — not just that the `sound` gate's aggregate count recovers.
- **Phase 3**: routing decision. Strict improvement, no material behavior
  change elsewhere → switch `generate_coords` to the new engine. Broad
  divergence → ship as a new API, deprecate `generate_coords`, point
  callers at `embed_pipeline_v2`.

## 6. Priority / non-goals

Legacy-path-only. Confirmed this round: `distance_geometry_v2.rs` and
`pipeline_v2.rs` (the production embedding path) never call
`dg::generate_coords`. `generate_coords_etkdg`
(`crates/chematic-3d/src/etkdg.rs:35`) calls `dg::generate_coords`
directly as its starting geometry before adding torsion noise, and the
legacy conformer ensemble path also depends on it — both are live,
non-trivial callers a half-finished rewrite would put at risk without
benefit. Not scheduled ahead of the current roadmap thread (A1/A2/A2.1
conformer-ensemble work and the best-of-N-vs-RDKit benchmark are already
done; this Phase 0 document is the next item after that, per the
project's own stated sequencing). Phase 1 implementation is separate,
future, unauthorized-by-this-PR work.
