# Generalized stereo-configuration geometry

Date: 2026-08-15. Branch: `feat/generalized-stereo-geometry`. Module:
`crates/chematic-core/src/stereo_geometry.rs`.

## 1. Why this exists

Before this PR, chematic had two independent, hand-written stereo-remapping
algorithms, both living in `chematic-smiles/src/canonical.rs`, both solving
the same underlying problem with unrelated code and unrelated proofs:

- `permutation_is_odd` -- tetrahedral `@`/`@@` parity, via classic
  cycle-counting permutation parity.
- `remap_square_planar` -- `@SP1`/`@SP2`/`@SP3` remapping, via matching a
  trans-pair-of-positions partition against 3 hand-written templates
  (`docs/rfcs/square_planar_stereo_rfc.md`).

Both answer the same question -- "given a declared stereo tag against one
neighbor ordering, what tag describes the same physical arrangement against
a *different* neighbor ordering?" -- but neither shares a line of code, a
data structure, or a proof strategy with the other. Every future coordination
geometry (trigonal-bipyramidal, octahedral) would otherwise need a third,
again-unrelated bespoke algorithm, and there would be no shared vocabulary to
even ask "is this new algorithm consistent with the old ones?"

This PR replaces both with one idea, standard in stereochemistry and
crystallography: a stereo configuration is a **coordination geometry** plus
an ordering of ligand-slot ids, and two orderings describe the **same
physical arrangement** iff one is reachable from the other by a **proper
rotation** of that geometry -- a rotation realizable by physically rotating
the rigid coordination shape in 3-space. Reflections are excluded on
purpose: a reflection of a chiral center generally produces the *other*
stereoisomer, which is exactly the distinction `@`/`@@` (or `@SPn`) exists to
preserve. The set of proper rotations of a geometry forms a group acting on
ligand-slot permutations (the geometry's "proper rotation group"); two
configurations are equivalent iff they reduce to the same representative
under that group.

This replaces both algorithms for the geometries' **common, 4-substituent
case** -- which is all of square-planar and the overwhelming majority of
tetrahedral. One real exception surfaced during integration (an allene *end*
carbon's declared stereo has only 3 real substituent positions, not 4 --
`StereoGeometry::Tetrahedral` is fixed at 4 slots by this PR's scope and
doesn't model that case): `canonical.rs` keeps the original, unmodified,
length-generic `permutation_is_odd` as an explicit fallback for that one
non-4-length shape. See §8's second delta for the concrete fixture and why
this is the correct call, not an unfinished replacement.

Scope of this PR: **Tetrahedral** and **Square-planar** only (the two
geometries chematic already models). Trigonal-bipyramidal and octahedral are
architected for (§9) but not implemented.

## 2. Slot convention per geometry

Both implemented geometries use exactly 4 ligand slots, numbered 0-3, stored
as a fixed `[u32; 4]` (raw ligand/atom ids, never a `Vec` -- see §7 for why
duplicate raw ids are still meaningful data, not necessarily an error).

- **Tetrahedral**: slot order is `Molecule::stereo_neighbor_order`'s existing
  convention -- the literal SMILES chirality-neighbor order, with
  `STEREO_H_SENTINEL` (`u32::MAX`) standing in for an implicit H. Unchanged
  from before this PR; this module doesn't touch how tetrahedral order is
  recorded, only how two recorded orders are compared.
- **Square-planar**: slot order is also the SMILES chirality-neighbor order,
  but additionally has a fixed *geometric* meaning independent of any
  particular molecule: slots 0/2 are one trans-pair, slots 1/3 are the
  other -- exactly `SquarePlanarPermutation::SP1`'s own convention
  (`trans_pairs() == [(0,2),(1,3)]`). A `(tag, neighbor_order)` pair is
  converted into this fixed base convention by `to_base_slots` (§4) before
  being handed to the rotation-group machinery at all.

## 3. Rotation-group derivation

### 3.1 Tetrahedral: A4, order 12

A rigid tetrahedron's proper rotation group is the alternating group **A4**
(all *even* permutations of its 4 vertices), order 12 -- the standard
identification from elementary group theory: every physical rotation of a
tetrahedron permutes its 4 vertices by an even permutation, and every odd
permutation of the vertices requires a reflection (not realizable by a
proper rotation). This is exactly why an odd permutation of a tetrahedral
center's neighbor order flips `@`<->`@@`: it's not realizable as a rotation,
so it necessarily describes the *other* enantiomer.

24 total orderings of 4 distinct ligands / 12 rotations = **2** orbits,
matching the existing 2-state `@`/`@@` tag.

`TETRAHEDRAL_ROTATIONS` lists: the identity, the 8 three-cycles
(`(012),(021),(013),(031),(023),(032),(123),(132)`), and the 3
double-transpositions (`(01)(23),(02)(13),(03)(12)`) -- `1 + 8 + 3 = 12`.
Cross-checked in-source two ways (both required by the task, both actually
implemented, not just asserted):

- `tetrahedral_rotations_form_a_group_of_order_12`: identity present,
  closed under composition, every element has an inverse, order exactly 12,
  no duplicate rows.
- `tetrahedral_rotations_are_independently_confirmed_even_permutations`: an
  **independent, from-scratch** brute-force cycle-counting parity function
  (distinct code from `remap_tetrahedral_parity`'s own internals, which
  don't even compute parity directly -- see §3.3) computes the even-
  permutation subset of all 24 permutations of `[0,1,2,3]` and asserts it
  equals `TETRAHEDRAL_ROTATIONS` exactly, as sets.

### 3.2 Square-planar: D4-flavored stabilizer, order 8 -- NOT the naive order-4 group

This is the crux of the whole design, and the part most likely to be gotten
wrong silently: **a square-planar center's proper rotation group has order
8, not 4.**

The naive guess -- "a square-planar complex only has 4 in-plane 90-degree
rotations about its own normal axis" -- is wrong. A physical square-planar
complex also has proper rotations by 180 degrees about an axis lying *in*
the molecular plane, through the midpoints of two opposite edges of the
coordination square. Such a rotation is a genuine rigid rotation of the
complex (not a reflection), and it **swaps which trans-pair is which** --
it maps ligands at positions `{0,2}` to positions `{1,3}` and vice versa. A
group that omits these 4 additional rotations under-counts the true
symmetry and gives the wrong orbit structure.

**Derivation (orbit-stabilizer theorem).** S4 (order 24, all permutations of
the 4 ligand slots) acts on the set of the 3 ways to partition `{0,1,2,3}`
into two unordered pairs: `{0,1}|{2,3}`, `{0,2}|{1,3}`, `{0,3}|{1,2}`. This
action is transitive (any partition can be mapped to any other by some
permutation of S4) with orbit size 3, so by orbit-stabilizer the stabilizer
of any one partition has order `|S4| / 3 = 24 / 3 = 8`.

**Explicit enumeration of the stabilizer of `{0,2}|{1,3}`** (SP1's own
partition, `trans_pairs() == [(0,2),(1,3)]`), in two cases:

- **Block-preserving** (maps `{0,2}`->`{0,2}` and `{1,3}`->`{1,3}`):
  independently permute within each 2-element block -> `2 x 2 = 4`
  elements: identity, `(02)`, `(13)`, `(02)(13)`.
- **Block-swapping** (maps `{0,2}`->`{1,3}` and `{1,3}`->`{0,2}`): a
  bijection `{0,2}->{1,3}` (2 choices) combined with a bijection
  `{1,3}->{0,2}` (2 choices), all 4 combinations giving genuine
  permutations of `{0,1,2,3}` -> `2 x 2 = 4` elements: `(01)(23)`,
  `(0123)`, `(0321)` (the 4-cycle's inverse), `(03)(12)`.

`4 + 4 = 8`, matching the orbit-stabilizer count exactly.

**Tying the table to `trans_pairs()`, not a hardcoded literal.** Per the
review that shaped this PR: the naive wrong group (the 4 in-plane rotations
only) is *also* closed, *also* has an identity, *also* has inverses for
every element -- every group-axiom self-test above would pass on it
silently, and it would still give `24 / 4 = 6` orbits, which cannot recover
the 3 real SP1/SP2/SP3 tags (the flagship orbit-count test, §4, is what
actually catches this class of error -- 6 != 3). To make the *specific*
partition this table stabilizes verifiably tied to the pre-existing,
oracle-verified `SquarePlanarPermutation::SP1.trans_pairs()` (not just a
`{0,2}|{1,3}` literal typed by hand in this RFC or in a test), the
in-source test `square_planar_rotations_stabilize_sp1_partition` derives
the expected partition from `trans_pairs()` *at runtime* via the same
`to_base_slots` helper production code uses (§4), and then asserts every
one of the 8 table rows stabilizes it. This is load-bearing, not
decorative: it is the one test in this PR that would fail if
`SQUARE_PLANAR_ROTATIONS` were accidentally built as the stabilizer of a
*different* one of the 3 partitions (e.g. `{0,1}|{2,3}`) -- a mistake the
group-axiom tests and even the orbit-count test alone cannot distinguish
from the correct table, since all three of S4's order-8 pair-partition
stabilizers are equally valid groups giving equally-3 orbits in the
abstract, just orbits that don't line up with `trans_pairs()`'s actual
SP1/SP2/SP3 semantics.

`SQUARE_PLANAR_ROTATIONS` lists, in the same order as the two cases above:
identity, `(02)`, `(13)`, `(02)(13)`, `(01)(23)`, `(0123)`, `(0321)`,
`(03)(12)`.

## 4. Configuration-class / canonicalization definition

`canonicalize_configuration(geometry, ligand_order: [u32; 4])` applies every
element of `geometry`'s rotation group to `ligand_order` (`apply(perm,
arr)[i] = arr[perm[i]]`, `const [[u8; 4]; N]` tables, hand-derived, not
runtime-generated) and keeps the lexicographically-smallest result as the
`CanonicalStereoConfiguration`'s private `representative`. Two
configurations are `equivalent_under_rotation` iff their canonical forms are
equal (currently literally `a == b`, since `representative` and `geometry`
both participate in derived equality).

**Tetrahedral bridge (`remap_tetrahedral_parity`).** Two orderings of the
same 4 distinct ids differ by an even permutation (no `@`<->`@@` flip) iff
they canonicalize to the same Tetrahedral representative -- true by
construction, since `TETRAHEDRAL_ROTATIONS` (A4) *is* the even-permutation
group. `remap_tetrahedral_parity(original, canonical)` is exactly
`canonicalize(Tetrahedral, original).representative !=
canonicalize(Tetrahedral, canonical).representative`. This is a *provable*
restatement of classic cycle-counting parity, not a different rule (see the
independent cross-check in §3.1).

**Square-planar bridge (`remap_square_planar_tag`).** `to_base_slots(tag,
order)` converts `(tag, order)` into the fixed base convention (§2): given
`tag.trans_pairs() == [(a,b),(c,d)]`, the reorder `[order[a], order[c],
order[b], order[d]]` always puts `order[a]`/`order[b]` at slots 0/2 and
`order[c]`/`order[d]` at slots 1/3 -- built directly from `trans_pairs()`
for any of the 3 tags, not a per-tag hand-written special case.
`remap_square_planar_tag(tag, original, canonical)` canonicalizes
`to_base_slots(tag, original)`, then tries all 3 candidate tags against
`to_base_slots(candidate, canonical)` and returns the (unique, when one
exists) candidate whose canonical form matches.

**Why this is provably equivalent to the removed `remap_square_planar`, not
just observed to agree on test fixtures.** A `SquarePlanar` rotation orbit
is *exactly* the set of orderings sharing the same unordered
trans-pair-of-ids partition `{{slots[0],slots[2]}, {slots[1],slots[3]}}`
(every element of `SQUARE_PLANAR_ROTATIONS` stabilizes that partition by
construction -- §3.2 -- and the group has no smaller stabilizer for 4
distinct ligand ids, since they have trivial individual stabilizers). That
partition is precisely what the old algorithm computed and matched against
the 3 templates. So `remap_square_planar_tag` and the old
`remap_square_planar` compute the *same function*, expressed two different
ways -- one directly manipulating a pair-of-pairs, the other via an
explicit rotation-group orbit. This is checked, not just argued: the
144-case oracle table (§12) and all 8 pre-existing fixture tests pass
byte-for-byte identically post-rewiring (§13).

**Idempotence.** `canonicalize(canonicalize(x).representative) ==
canonicalize(x)` for both geometries, checked exhaustively over all 24
orderings of `[0,1,2,3]` in `canonicalization_is_idempotent`.

## 5. Mapping to OpenSMILES/RDKit tags

- Tetrahedral: `@` = `Chirality::CounterClockwise`, `@@` =
  `Chirality::Clockwise`. Unchanged by this PR -- `Chirality`'s public shape
  is untouched (a hard constraint of this task); this module only backs the
  *parity computation* `canonical.rs` already performed.
- Square-planar: `@SP1`/`@SP2`/`@SP3` = `Chirality::SquarePlanar(SP1/SP2/SP3)`,
  via the pre-existing, oracle-verified `SquarePlanarPermutation::trans_pairs()`
  (`docs/rfcs/square_planar_stereo_rfc.md`) -- this PR reuses that mapping as
  ground truth (§12), it does not re-derive or re-validate it.

## 6. Atom-renumbering transformation

`StereoConfiguration::renumber(&self, id_map: impl Fn(u32) -> Option<u32>)`
remaps every slot id through `id_map`, failing closed with
`StereoGeometryError::UnknownLigandId` the instant any slot has no answer,
rather than silently dropping or zeroing it.

**Renumbering invariance**, checked exhaustively
(`renumber_preserves_rotation_equivalence`, both geometries, every rotation
in each group, a deliberately *non-monotonic* bijective `id_map`): if two
configurations are rotation-equivalent before renumbering, they remain
rotation-equivalent after applying the *same* renumbering map to both --
even when the map doesn't preserve numeric ordering (so a naive "just
re-derive which candidate was lexicographically smallest" argument would not
obviously survive relabeling; the actual proof is that the *set* of
rotation-orbit candidates commutes with an elementwise id map, which is
what's checked).

## 7. Duplicate-ligand handling and why

`canonicalize_configuration`/`equivalent_under_rotation`/`remap_*` operate
purely on raw `u32` slot ids -- **never chemical identity**. This matches
what `remap_square_planar` already did before this PR (it compared raw
neighbor atom ids, never element symbols or CIP rank). It's why cisplatin
(`@SP1`, 2xCl + 2xNH3 *cis*) and transplatin (`@SP2`, same composition,
*trans*) stay distinguishable even though both have two chemically-identical
Cl ligands and two chemically-identical NH3 ligands: the geometry layer sees
4 distinct **atom ids** regardless of element repeats. "These two slots
happen to be the same chemical species" is a separate, molecule-canonical-
ranking-layer concern (Morgan/CIP-style priority) that this module
deliberately does not need to know about -- it operates one layer below
chemical identity, on raw connectivity-graph ids only.

A genuinely **duplicate raw id in two slots** (the same atom id appearing
twice -- a data-integrity problem, not "two chemically-identical atoms",
which always have distinct ids) is the one case this module does actively
reject: `StereoGeometryError::DuplicateSlotId`, checked via a plain O(1)-
bounded pairwise scan (no `HashSet`, no heap allocation, deterministic).

## 8. Unspecified/invalid handling (fail-closed)

- `canonicalize_configuration`: `Err(DuplicateSlotId)` on a repeated slot id.
- `StereoConfiguration::renumber`: `Err(UnknownLigandId)` the moment
  `id_map` has no answer for any slot.
- `remap_square_planar_tag`: `None` (not a guessed tag) whenever no
  candidate tag's canonical form matches -- which happens exactly when
  `original`/`canonical` don't name the same 4 distinct ids (mismatched
  set, or either array has a duplicate), since every candidate then either
  errors (`DuplicateSlotId`) or lands on a representative containing an id
  `original` didn't have, which can never equal `original`'s own
  representative. No extra "same id set" guard is needed -- it falls out of
  the canonicalization machinery for free.
- `remap_tetrahedral_parity`: `Err(DuplicateSlotId)` on a repeated id in
  either array. At its `canonical.rs` call site this is treated the same as
  the pre-existing "no verifiable order" pass-through-unchanged fallback --
  the documented safe no-op for a 2-state tag (`docs/rfcs/square_planar_stereo_rfc.md`'s
  own framing: unchanged-against-unverified is still a *valid* state for a
  2-state tag, just not provably correct; it is *not* safe for a 3-state
  tag, which square-planar's `None`-on-any-doubt handling reflects).

No panics anywhere in this module on any input; no silent modulo-wrapping of
out-of-range ids (ids are opaque `u32`s compared for equality only, never
indexed or wrapped).

**Two real, benign behavioral deltas from the pre-PR code**, found during
integration testing (not by inspection) and recorded here rather than
smoothed over:

1. A **duplicate id in a 4-element tetrahedral order** now produces
   `Err(DuplicateSlotId)` -> the pass-through-unchanged fallback, where the
   old `permutation_is_odd`'s `HashMap`-based lookup would have silently
   computed *some* deterministic-but-essentially-arbitrary parity value
   (via `pos.get(v).unwrap_or(&0)`, the first duplicate "wins" the position
   0 slot). This is strictly fail-safer, not fail-worse, and unreachable
   from any real parsed molecule (a genuine tetrahedral center always has 4
   distinct substituent atoms); no existing fixture exercises it.
2. A **tetrahedral order whose length isn't 4** (see §12.1: this is not a
   hypothetical -- it's the real, legitimate shape of an allene *end*
   carbon's declared stereo, e.g. `F[C@@H]=[C]=[C@H]Cl`'s F-bearing atom
   has a 3-element `stereo_neighbor_order`, `[F, implicit-H sentinel, =C
   partner]`, since an sp2 allene-end carbon has only 3 real attachment
   points) is **not** routed through `StereoGeometry::Tetrahedral` at all --
   that geometry is fixed at 4 slots by design (matching this PR's stated
   scope). `canonical.rs`'s call site keeps the original, unmodified,
   length-generic `permutation_is_odd` as an explicit fallback for this
   case, so behavior for allene-end stereocenters is provably unchanged
   (byte-identical output, §13) rather than silently dropped or
   reinterpreted. `StereoGeometry::Tetrahedral` genuinely only covers the
   common 4-substituent case; a length-3 "tetrahedral-like" geometry for
   cumulated systems is out of scope for this PR and not modeled here.

## 9. Extending to TBP / Octahedral later (sketch, not implemented)

`StereoGeometry` is `#[non_exhaustive]` specifically so this extension
doesn't require a breaking change. Adding trigonal-bipyramidal (5 ligands,
`@TB1`-`@TB20` in OpenSMILES) or octahedral (6 ligands, `@OH1`-`@OH30`)
would need:

1. A new `StereoGeometry` variant.
2. A slot-count generalization: `StereoConfiguration`/`CanonicalStereoConfiguration`
   are currently hardcoded to `[u32; 4]` for the minimal-allocation, exactly-
   two-geometries scope of this PR; a 5/6-slot geometry would need either a
   const-generic `[u32; N]` (cleanest, but touches every signature in this
   module) or a small fixed-capacity array type wide enough for the largest
   supported geometry (e.g. `[u32; 6]` with a `len` discriminant) -- a real
   design decision to make at that time, not decided here.
3. A new rotation-group table, derived the same way §3 derives A4/D4. Both
   cross-checked below against the OpenSMILES tag count for that geometry
   (`5! / |group| = 20`, `6! / |group| = 30`) -- the same orbit-count
   arithmetic the flagship test in §3/§4 relies on, and the cheapest
   correctness anchor available for a section with no actual test behind
   it, since neither geometry is implemented in this PR.

   - **TBP** (5 positions: 3 equatorial + 2 axial; OpenSMILES `@TB1`-`@TB20`).
     The full point group of a trigonal bipyramid is D3h, order 12 (E, 2C3,
     3C2, sigma_h, 2S3, 3sigma_v); its **proper rotation subgroup is D3,
     order 6**, not the naive "independently permute the 3 equatorial
     positions AND independently swap the 2 axial positions" guess (which
     would wrongly give `3! x 2 = 12`, treating the two choices as
     unconstrained/independent when physically they are not). D3's 6
     elements: identity; 2 rotations by +-120 degrees about the axial
     axis (C3, C3^2) -- these fix both axial positions and cyclically
     permute the 3 equatorial positions; and 3 rotations by 180 degrees
     about an axis running through one equatorial vertex and perpendicular
     to the axial axis (C2) -- each such rotation **swaps the 2 axial
     positions AND simultaneously transposes the other 2 equatorial
     positions**, fixing the equatorial vertex the axis passes through (a
     C2 axis through, say, equatorial position E1 flips top<->bottom,
     which swaps the 2 axial atoms, and 180-degree rotation about that
     axis also swaps E2<->E3). An "equatorial transposition with the
     axials left fixed" is *not* in this group -- it's a reflection
     (one of D3h's 3 sigma_v planes), not a proper rotation. Orbit-count
     cross-check: `5! / 6 = 120 / 6 = 20`, matching `@TB1`-`@TB20` exactly.
   - **Octahedral** (6 positions; OpenSMILES `@OH1`-`@OH30`). The proper
     rotation group of an octahedron (equivalently, of the 6 vertices of an
     octahedral coordination sphere) has order 24, isomorphic to S4 (acting
     on the 4 body-diagonals connecting opposite vertex-pairs of the dual
     cube). Orbit-count cross-check: `6! / 24 = 720 / 24 = 30`, matching
     `@OH1`-`@OH30` exactly.

   Both derivations above are stated here as a sketch, cross-checked only
   by tag-count arithmetic; neither has the from-scratch orbit-stabilizer
   argument, explicit element enumeration, or exhaustive group-axiom/orbit-
   count self-tests §3.1/§3.2 give A4/D4 -- that work is required before
   either would be trustworthy enough to implement, exactly like this PR's
   own square-planar derivation was, not assumed by analogy to it.
4. New bridge functions analogous to `remap_square_planar_tag`, once
   OpenSMILES's `@TBn`/`@OHn` tag semantics are derived (this codebase has
   no oracle-verified TBP/octahedral tag semantics yet -- that derivation
   work, likely following the same RDKit-oracle-empirical methodology
   `docs/rfcs/square_planar_stereo_rfc.md` used for square-planar, is a
   separate, unstarted prerequisite, not part of this sketch).

None of this is implemented in this PR.

## 10. Independent derivation statement

This module was derived from group-theory fundamentals (the orbit-stabilizer
theorem, explicit permutation enumeration, standard proper-rotation-group
identification for a tetrahedron) and this codebase's own pre-existing,
independently oracle-verified `SquarePlanarPermutation::trans_pairs()` --
the generalization is "coordinate-geometry-plus-rotation-group" thinking
common in stereochemistry and crystallography, worked out from scratch for
this codebase's own data model. Zero dependency on any external
cheminformatics library for this design, and zero code, comments, tables, or
fixtures copied from any third-party source. The only external reference
used anywhere in this PR is the local RDKit installation, and only as a
**regression** anchor for the pre-existing `trans_pairs()`/`remap_square_planar`
semantics (§12) -- this PR performs no *new* RDKit oracle queries and
derives no new chemistry from RDKit.

## 11. Known gaps (not fixed by this PR)

Two pre-existing, real, unrelated bugs were identified while implementing
this PR and are deliberately **not** fixed here, per the approved scope:

- **`chematic-mol/src/mol2000.rs`'s `wedge_vs_3d_conflicts`** gates on
  `atom.chirality == Chirality::None` (an equality check) rather than
  `!atom.chirality.is_tetrahedral()` -- the same "equality vs. exhaustive/
  `is_tetrahedral()` check" gap `docs/rfcs/square_planar_stereo_rfc.md`
  already found and fixed at 3 other call sites when `SquarePlanar` was
  introduced. Confirmed **unreachable today**: chematic-mol's MOL/SDF
  readers never produce `Chirality::SquarePlanar` (that variant's producers
  are the SMILES parser and direct `MoleculeBuilder` construction only), so
  this specific gate can never actually see a `SquarePlanar` atom in
  practice. Left as a documented, latent gap, not fixed in this PR.
- **`chematic-core/src/molecule.rs`'s `Molecule::fragments()`** silently
  drops `stereo_neighbor_order` entirely when splitting a molecule into
  connected components (the per-fragment `MoleculeBuilder` never calls
  `set_stereo_neighbor_order` for any atom). Any declared tetrahedral or
  square-planar stereo on an atom that survives fragmentation loses its
  neighbor-order record, silently falling back to `corrected_chirality`'s
  "no recorded order" path on next canonicalization (pass-through-unchanged
  for tetrahedral, `Chirality::None` drop for square-planar) -- a real,
  pre-existing, unrelated data-loss bug, out of scope for and not touched by
  this PR.

## 12. RDKit differential fixtures / oracle provenance

**Local RDKit version actually checked for this PR** (not restated from an
earlier PR's citation): `.venv/bin/python -c "import rdkit;
print(rdkit.__version__)"` -> **`2026.03.4`**. This is a newer patch version
than the `2026.03.3` the original square-planar RFC
(`docs/rfcs/square_planar_stereo_rfc.md`) pinned when the
`SquarePlanarPermutation::trans_pairs()`/`remap_square_planar` semantics
were first oracle-derived (144-case table,
`scripts/square_planar_permutation_oracle.py`). Recorded honestly rather
than smoothed over, per this task's explicit instruction not to silently
pick one version's answer.

**No new RDKit oracle validation was performed by this PR.** Per the task's
own instruction, this PR does not re-derive or re-run the oracle script
against the (now `2026.03.4`) local RDKit install -- the 144-case table is
reused strictly as a **regression** fixture proving this PR's new
rotation-group-based code reproduces the *previously-validated*
`trans_pairs()`/`remap_square_planar` semantics, not as a fresh oracle
validation. If a future PR needs to confirm `2026.03.3`-vs-`2026.03.4`
parity for the underlying chemistry, that is a separate, explicit task (and
would need to handle the "RDKit shows version-dependent behavior" stop
condition this task's brief calls out, which this PR did not need to
because it made no new oracle calls).

### 12.1 Fixtures used

- **The 144-case oracle table** (24 neighbor permutations x 3 tags x 2
  molecule shapes, `scripts/square_planar_permutation_oracle.py`), reused
  two ways as a regression check:
  1. End to end, unchanged, via the pre-existing
     `oracle_verified_permutation_table_matches_end_to_end` test in
     `crates/chematic-smiles/tests/square_planar_stereo.rs` (parser +
     rewired canonical writer).
  2. **New**, direct/unit-level, added by this PR:
     `remap_square_planar_tag_matches_oracle_table_directly` in the same
     file -- calls `chematic_core::remap_square_planar_tag` directly against
     the same 24x3 sweep, sharing the file's own `TAGS`/`predict`/
     `permutations_of_4` fixture-generation code (extracted to file scope
     by this PR so both tests use one definition, not two hand-transcribed
     copies).
- **88 stereo-bearing fixture SMILES**, assembled from every `@`-containing
  string literal already present in `chematic-smiles/src/canonical.rs`,
  `crates/chematic-smiles/tests/canonical_robustness.rs`, and
  `crates/chematic-smiles/tests/square_planar_stereo.rs`, plus all 18 SMILES
  in `validation/platinum/pt_corpus.jsonl` (the pre-existing "before/after
  byte-identical" precedent from `docs/rfcs/square_planar_stereo_rfc.md`'s
  own verification section, reused here rather than re-invented), plus 2
  synthetic isotope/atom-map cases. Dumped via a throwaway
  `chematic-smiles/examples/stereo_baseline_dump.rs` (never committed --
  used locally, then removed before staging this PR's changes) to
  `input\tcanonical` pairs, once before any code change and once after,
  diffed byte-for-byte. **First pass found a real discrepancy** (see
  §8's item 2 and §13) -- the diff-based check caught a behavior change the
  pass/fail unit test suite alone did not (those tests check relative
  invariance -- "two enantiomers differ", "round-trips are stable" -- not
  exact golden output, so they kept passing throughout even while the
  literal output value briefly changed).
- The two named CIP-gating regression tests
  (`chematic-chem/src/cip.rs::square_planar_center_never_gets_a_bogus_tetrahedral_cip_code`,
  `chematic-cip/src/tests.rs::square_planar_center_is_skipped_not_assigned_a_bogus_cip_code`),
  run and confirmed passing unchanged -- proving this PR didn't disturb an
  already-correct boundary elsewhere.

## 13. Performance

Benchmark: `crates/chematic-smiles/benches/parse_bench.rs`'s pre-existing
`canonical_smiles_10mol` (criterion, includes glucose,
`OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O`, 5 tetrahedral stereocenters -- the
only stereo-dense entry in that benchmark's fixed 10-molecule set). Reused
rather than adding a new bench target, since a new `[[bench]]` entry would
require a `Cargo.toml` change (forbidden by this task).

**Honest result: inconclusive on this machine, no directional claim made.**
Six runs total, criterion 3-point estimates (low/median/high), in the order
taken:

| # | Tree | Note | time |
|---|---|---|---|
| 1 | `main` (pristine -- taken right after `git checkout -b`, before `stereo_geometry.rs` existed) | baseline, early in this session | 124.99 / **134.53** / 145.48 µs |
| 2 | this PR's code | concurrent w/ an unrelated `cargo test` -- contaminated | 168.00 / 188.87 / 212.72 µs |
| 3 | this PR's code | `ps aux` confirmed idle | 101.40 / **108.40** / 115.74 µs |
| 4 | this PR's code | `ps aux` confirmed idle, immediately after #3 | 97.65 / 104.03 / 110.60 µs |
| 5 | `main` (`git stash` to the exact pre-PR tree, again) | `ps aux` confirmed idle | 117.17 / **124.60** / 132.25 µs |
| 6 | this PR's code (`git stash pop` back) | `ps aux` confirmed idle, immediately after #5 | 127.01 / 138.85 / 152.40 µs |

Before-side (runs 1, 5): medians 134.53 and 124.60 µs, consistent with each
other. After-side (runs 3, 4, 6 -- run 2 excluded as contaminated): medians
108.40, 104.03, 138.85 µs -- a wider spread that itself overlaps the
before-side's range. This is exactly the overlap the conclusion below draws
on, now stated from a table with both sides correctly attributed.

Runs 3-6 were all taken with `ps aux` confirming zero other `cargo`/`rustc`
processes of this session running, yet still span 98-152 µs -- because
`uptime` during this window reported **load averages of 14-17** on a
machine with 2 logged-in users, i.e. substantial background activity from
outside this session that `ps aux`'s narrower "no cargo/rustc process of
mine" check cannot see or control for. Criterion's own paired significance
test, run pairwise against whatever it had stored from the *immediately
preceding* invocation (not always the run listed as "before" above -- this
was gotten wrong once during this work and corrected: criterion overwrites
its stored baseline on every run, so a "run N vs. the original before
baseline" claim is only valid if nothing else was benchmarked in between),
called some adjacent pairs "improved," others "regressed," and others "no
change" across these 6 runs -- entirely explained by which two runs
happened to land on which side of the external load's fluctuation, not by
any code difference (run 4 vs. run 3 -- same code, immediately sequential --
already shows a 97.65-115.74 µs spread with "no change" only at `p=0.88`).

**Conclusion (exact wording, do not paraphrase into a speed claim in either
direction): No clear performance degradation was observed; the measured
difference is within noise and cannot be statistically confirmed as either
an improvement or a regression on this shared, noisy machine.** The
*before* (runs 1 and 5, on `main`) and *after* (runs 3/4/6, on this PR's
tree) absolute-time envelopes overlap substantially (roughly 98-152 µs on
both sides). Per this task's stop condition ("if the generalized approach
turns out to be measurably slower ... report the numbers and your
assessment, don't hide it") the numbers above are reported in full rather
than a favorable subset -- but the honest reading of that full data is "not
measurable either way under these conditions," not "probably faster" or
"comparable" or any other directional gloss. Architecturally, not
empirically: the new code path removes a `HashMap<T, usize>` allocation+
build entirely from the common (4-element) tetrahedral/square-planar write
path -- `permutation_is_odd`'s `HashMap` construction now only runs for the
non-4-length allene-end fallback (§8, rare in practice); the common case
instead goes through `canonicalize_configuration`'s fixed-size array scan
over a `const` 8-or-12-row table (24 or 16 plain `[u32; 4]` copies-and-
compares, no heap allocation at all). This architectural difference is noted
for a reader curious *why* no degradation would be expected; it is
explicitly not offered as evidence of an improvement, since the measurement
above cannot confirm one. A re-measurement on a quieter/dedicated machine
(or a much longer criterion run to average out the observed noise) would be
needed to actually confirm or refute any directional claim.

## 14. Verification summary

See the PR body for the full command-by-command verification gate output
(`cargo test --workspace --all-features`, `--no-default-features`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`,
`cargo fmt --all -- --check`, `cargo check --target wasm32-unknown-unknown
-p chematic-core`, the zero-third-party-dependency grep described in §10
returning nothing, and the byte-identical fixture diff described in
§12.1).
