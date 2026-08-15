//! Generalized stereo-configuration geometry: a coordination geometry plus
//! the equivalence class of ligand-slot permutations under that geometry's
//! *proper rotation group*.
//!
//! # Why this exists
//!
//! Before this module, chematic had two independent, hand-written
//! stereo-remapping algorithms living in `chematic-smiles/src/canonical.rs`:
//! `permutation_is_odd` (tetrahedral `@`/`@@` parity via cycle counting) and
//! `remap_square_planar` (`@SP1`/`@SP2`/`@SP3` trans-pair-partition matching).
//! Both solve the same underlying problem -- "given a declared stereo tag
//! against one neighbor ordering, what tag describes the same physical
//! arrangement against a *different* neighbor ordering?" -- with unrelated
//! code, unrelated proofs, and no shared vocabulary. Every future geometry
//! (trigonal-bipyramidal, octahedral) would otherwise need its own bespoke
//! third algorithm.
//!
//! This module replaces both with one idea, standard in stereochemistry and
//! crystallography: a stereo configuration is [`StereoGeometry`] (a
//! coordination shape) plus an ordering of ligand-slot ids, and two orderings
//! describe the *same physical arrangement* iff one can be reached from the
//! other by a **proper rotation** of that geometry (a rotation realizable by
//! physically rotating the rigid coordination shape in 3-space -- reflections
//! excluded, since a reflection generally produces a different stereoisomer).
//! The set of proper rotations of a geometry forms a group acting on the
//! ligand-slot permutations; "canonicalizing" a configuration means picking
//! the lexicographically-smallest ordering reachable under that group, and
//! two configurations are equivalent iff they canonicalize to the same
//! representative.
//!
//! See `docs/rfcs/generalized_stereo_geometry_rfc.md` for the full
//! derivation, oracle/regression provenance, and the TBP/octahedral
//! extension sketch.
//!
//! # Scope
//!
//! [`StereoGeometry::Tetrahedral`] and [`StereoGeometry::SquarePlanar`] only.
//! `StereoGeometry` is `#[non_exhaustive]` so future geometries (TBP,
//! octahedral) can be added without a breaking change, but none are
//! implemented here.
//!
//! # Independent derivation
//!
//! This module was derived from group-theory fundamentals (orbit-stabilizer
//! theorem, explicit permutation enumeration) and this codebase's own
//! previously oracle-verified `SquarePlanarPermutation::trans_pairs()` --
//! zero dependency on, and zero code copied from, any third-party
//! cheminformatics library.

use crate::atom::SquarePlanarPermutation;

// ---------------------------------------------------------------------------
// Core types
// ---------------------------------------------------------------------------

/// A coordination geometry with a defined ligand-slot count and proper
/// rotation group. `#[non_exhaustive]`: trigonal-bipyramidal and octahedral
/// are architected for (see the RFC's extension sketch) but not implemented.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StereoGeometry {
    /// 4-coordinate tetrahedral center (`@`/`@@`). Proper rotation group:
    /// A4 (alternating group on 4 points, order 12 -- all even permutations).
    Tetrahedral,
    /// 4-coordinate square-planar center (`@SP1`/`@SP2`/`@SP3`). Proper
    /// rotation group: the order-8 stabilizer, under S4, of the trans-pair
    /// partition `{0,2}|{1,3}` (see [`SQUARE_PLANAR_ROTATIONS`]'s doc for the
    /// orbit-stabilizer derivation).
    SquarePlanar,
}

impl StereoGeometry {
    /// The proper rotation group for this geometry, as `apply(perm,
    /// arr)[i] = arr[perm[i]]`-convention permutations of the 4 ligand
    /// slots. `const`, hand-derived -- not runtime-generated.
    fn rotation_group(self) -> &'static [[u8; 4]] {
        match self {
            Self::Tetrahedral => &TETRAHEDRAL_ROTATIONS,
            Self::SquarePlanar => &SQUARE_PLANAR_ROTATIONS,
        }
    }
}

/// A declared stereo configuration: a geometry plus the raw ligand ids
/// occupying its 4 slots, in the order that geometry's tag semantics were
/// declared against (e.g. SMILES chirality-neighbor order). `slots` holds
/// raw `u32` ids only -- never chemical identity -- so callers that need to
/// treat chemically-identical-but-distinct-atom ligands specially (e.g.
/// duplicate-ligand detection for CIP-style priority) must do so at a layer
/// above this one; see the RFC's "duplicate ligands" section.
///
/// `pub(crate)`, not `pub`: this type, [`CanonicalStereoConfiguration`],
/// [`canonicalize_configuration`], and [`equivalent_under_rotation`] are all
/// hardcoded to `[u32; 4]`, which only fits the two 4-coordinate geometries
/// this PR implements. Publishing them now would commit chematic-core's
/// public API to "every geometry has exactly 4 slots" -- a claim
/// [`StereoGeometry`]'s own `#[non_exhaustive]` explicitly declines to make,
/// since trigonal-bipyramidal (5 slots) and octahedral (6 slots) are
/// architected for (see the RFC's extension sketch). Only the two bridge
/// functions actual callers need --
/// [`remap_tetrahedral_parity`]/[`remap_square_planar_tag`], both already
/// geometry-specific and arity-fixed by their own OpenSMILES tag semantics,
/// not by an assumption this module bakes in -- are `pub`. Fields are
/// private even at `pub(crate)` scope: the only way to build one is
/// [`StereoConfiguration::new`], which runs the same duplicate check
/// [`canonicalize_configuration`] does, so no code path inside this crate
/// can construct an unvalidated configuration either.
// No production caller constructs a `StereoConfiguration` today -- both
// production bridge functions (`remap_tetrahedral_parity`/
// `remap_square_planar_tag`) operate directly on raw `[u32; 4]` arrays via
// `canonicalize_configuration`, not through this wrapper type. This type,
// `new`, and `renumber` exist as tested, validated infrastructure per this
// PR's required API (atom-renumbering transformation) ahead of a production
// consumer -- same shape as `stereo_constraints.rs`'s own
// `TetrahedralConstraint`/`StereoConstraintSet::unsupported`, which carry an
// identical `#[allow(dead_code)]` for the same reason (see that module).
// Promote by removing this `allow` once a real caller constructs one.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct StereoConfiguration {
    geometry: StereoGeometry,
    slots: [u32; 4],
}

#[allow(dead_code)]
impl StereoConfiguration {
    /// Construct a configuration, rejecting a duplicate slot id up front --
    /// the same check [`canonicalize_configuration`] runs, applied here too
    /// so a `StereoConfiguration` can never exist in an already-invalid
    /// state (struct-literal construction is unavailable outside this
    /// module; private fields, no other constructor).
    pub(crate) fn new(
        geometry: StereoGeometry,
        slots: [u32; 4],
    ) -> Result<Self, StereoGeometryError> {
        if let Some(dup) = find_duplicate(slots) {
            return Err(StereoGeometryError::DuplicateSlotId(dup));
        }
        Ok(Self { geometry, slots })
    }

    /// Remap every slot id through `id_map` (e.g. an atom-renumbering table),
    /// preserving `geometry`. Fails closed two ways, both returning
    /// [`StereoGeometryError`] rather than silently producing a bad
    /// configuration: [`StereoGeometryError::UnknownLigandId`] the moment
    /// any slot's id has no answer in `id_map`; and
    /// [`StereoGeometryError::DuplicateSlotId`] if `id_map`, while total on
    /// the input slots, is not *injective* on them -- e.g. mapping two
    /// distinct input ids to the same output id -- which would otherwise
    /// silently manufacture a duplicate that could never have been accepted
    /// by [`Self::new`]/[`canonicalize_configuration`] directly. The input
    /// itself is already known duplicate-free (a `StereoConfiguration` can
    /// only exist via [`Self::new`]'s own check), so any duplicate detected
    /// here was created BY the renumbering, not carried over from it.
    pub(crate) fn renumber(
        &self,
        id_map: impl Fn(u32) -> Option<u32>,
    ) -> Result<StereoConfiguration, StereoGeometryError> {
        let mut new_slots = [0u32; 4];
        for (i, slot) in self.slots.iter().enumerate() {
            new_slots[i] = id_map(*slot).ok_or(StereoGeometryError::UnknownLigandId(*slot))?;
        }
        StereoConfiguration::new(self.geometry, new_slots)
    }
}

/// A [`StereoConfiguration`] reduced to its canonical representative under
/// its geometry's proper rotation group -- the "configuration class" /
/// equivalence-class identity. Two configurations describe the same physical
/// arrangement iff their canonical forms are equal (see
/// [`equivalent_under_rotation`]). Fields are private: the only way to
/// compare two configurations is through this type's own equality /
/// [`equivalent_under_rotation`], never by peeking at which specific
/// group-orbit member happened to sort first. `pub(crate)`: see
/// [`StereoConfiguration`]'s doc for why this whole family of types is not
/// yet public.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct CanonicalStereoConfiguration {
    geometry: StereoGeometry,
    representative: [u32; 4],
}

impl CanonicalStereoConfiguration {
    /// The lexicographically-smallest slot ordering in this configuration's
    /// rotation orbit. Exposed so callers can feed it back into
    /// [`canonicalize_configuration`] (e.g. to verify idempotence) or into
    /// [`StereoConfiguration::renumber`]; not meaningful as "the" canonical
    /// spelling of anything outside this module -- only equality between two
    /// [`CanonicalStereoConfiguration`] values is a meaningful comparison.
    pub(crate) fn representative(&self) -> [u32; 4] {
        self.representative
    }
}

/// Fail-closed errors for this module. No panics, no silent
/// modulo-wrapping, no lossy fallback for out-of-range/duplicate/unmapped
/// ligand ids anywhere in this module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoGeometryError {
    /// The same raw ligand id appeared in two (or more) slots of a
    /// configuration -- a data-integrity problem, not a valid 4-distinct-
    /// ligand arrangement. Also returned by [`StereoConfiguration::renumber`]
    /// when a non-injective `id_map` *creates* a duplicate that wasn't in
    /// the original configuration -- reusing this variant rather than adding
    /// a separate "renumbering created a duplicate" one, since the resulting
    /// state is identical either way (two slots now share one id) and every
    /// caller's fail-closed handling is the same regardless of which
    /// operation produced it.
    DuplicateSlotId(u32),
    /// [`StereoConfiguration::renumber`]'s `id_map` had no answer for a
    /// slot's id.
    UnknownLigandId(u32),
    /// [`remap_tetrahedral_parity`]'s `original` and `canonical` arrays
    /// don't name the same 4 distinct ids (as a set) -- e.g. a foreign id
    /// present in one but not the other. Computing a parity flip by
    /// comparing canonical representatives is only meaningful when both
    /// arrays are genuine permutations of the *same* 4 ids; without this
    /// check, two arrays naming different id sets would (correctly, but
    /// meaninglessly) canonicalize to unequal representatives, which the
    /// parity computation would misread as "needs a flip" -- a
    /// confident-looking wrong answer for malformed input, not the honest
    /// "can't tell" this variant exists to report instead.
    MismatchedLigandSet,
}

impl core::fmt::Display for StereoGeometryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DuplicateSlotId(id) => {
                write!(f, "duplicate ligand id {id} appears in two stereo slots")
            }
            Self::UnknownLigandId(id) => {
                write!(f, "no renumbering answer for ligand id {id}")
            }
            Self::MismatchedLigandSet => {
                write!(
                    f,
                    "original and canonical orders do not name the same ligand ids"
                )
            }
        }
    }
}

impl std::error::Error for StereoGeometryError {}

// ---------------------------------------------------------------------------
// Rotation groups
// ---------------------------------------------------------------------------

/// Apply a rotation-group permutation to a 4-slot ligand order.
/// `apply(perm, arr)[i] = arr[perm[i]]`.
fn apply(perm: &[u8; 4], arr: [u32; 4]) -> [u32; 4] {
    [
        arr[perm[0] as usize],
        arr[perm[1] as usize],
        arr[perm[2] as usize],
        arr[perm[3] as usize],
    ]
}

/// The tetrahedral proper rotation group: **A4**, the alternating group on
/// 4 points (all even permutations), order 12. A physical rotation of a
/// rigid tetrahedron permutes its 4 vertices by an even permutation only
/// (this is the standard identification of the tetrahedron's rotation
/// group with A4; an odd permutation of the 4 vertices requires a
/// reflection, which is *not* realizable by a proper rotation and is
/// exactly the operation that flips `@`<->`@@`). 24 total orderings / 12
/// rotations = 2 orbits, matching the existing `@`/`@@` (CW/CCW) 2-state
/// tag this table backs.
///
/// Listed as: identity, the 8 three-cycles, the 3 double-transpositions
/// (1 + 8 + 3 = 12). Cross-checked in this module's tests against an
/// independently-written brute-force parity function (not
/// [`remap_tetrahedral_parity`]'s own internals) over all 24 permutations
/// of `[0,1,2,3]`.
const TETRAHEDRAL_ROTATIONS: [[u8; 4]; 12] = [
    [0, 1, 2, 3], // identity
    // 3-cycles (8): (012),(021),(013),(031),(023),(032),(123),(132)
    [1, 2, 0, 3],
    [2, 0, 1, 3],
    [1, 3, 2, 0],
    [3, 0, 2, 1],
    [2, 1, 3, 0],
    [3, 1, 0, 2],
    [0, 2, 3, 1],
    [0, 3, 1, 2],
    // double-transpositions (3): (01)(23),(02)(13),(03)(12)
    [1, 0, 3, 2],
    [2, 3, 0, 1],
    [3, 2, 1, 0],
];

/// The square-planar proper rotation group: the order-**8** stabilizer,
/// under S4 (order 24), of the trans-pair partition `{0,2}|{1,3}` --
/// [`SquarePlanarPermutation::SP1`]'s own partition (`trans_pairs()`). NOT
/// the naive "4 in-plane rotations only" group (order 4, which would give
/// 24/4 = 6 orbits and cannot recover the 3 real SP1/SP2/SP3 tags) -- a
/// square-planar center's rotation symmetry includes the 4 rotations that
/// swap the two trans-pairs (a 90 degree rotation about an axis in the
/// molecular plane, through the midpoints of two opposite edges of the
/// square, is a genuine proper rotation of the physical complex and swaps
/// which pair is "trans-pair A" vs "trans-pair B") in addition to the 4
/// that preserve each trans-pair individually.
///
/// Derivation (orbit-stabilizer theorem): S4 acts transitively on the 3 ways
/// to partition `{0,1,2,3}` into two unordered pairs (`{0,1}|{2,3}`,
/// `{0,2}|{1,3}`, `{0,3}|{1,2}`) -- any of the 3 partitions can be mapped to
/// any other by some permutation, and the action is transitive with a
/// 3-element orbit, so by orbit-stabilizer the stabilizer of one partition
/// has order `|S4| / 3 = 24 / 3 = 8`.
///
/// Explicit enumeration of the 8-element stabilizer of `{0,2}|{1,3}`, in two
/// cases:
/// - **Block-preserving** (maps `{0,2}` to itself and `{1,3}` to itself):
///   independently permute within each block -> 2 x 2 = 4 elements:
///   identity, `(02)`, `(13)`, `(02)(13)`.
/// - **Block-swapping** (maps `{0,2}` to `{1,3}` and vice versa): a
///   bijection `{0,2}->{1,3}` (2 choices) combined with a bijection
///   `{1,3}->{0,2}` (2 choices) -> 4 elements: `(01)(23)`, `(0123)`,
///   `(0321)`, `(03)(12)`.
///
/// 4 + 4 = 8, matching the orbit-stabilizer count. This exact 8-element set
/// is asserted against `SquarePlanarPermutation::SP1.trans_pairs()` at
/// runtime in this module's tests (`square_planar_rotations_stabilize_sp1_partition`),
/// not just hand-verified in this comment -- the table below is
/// *load-bearing*, verifiably tied to the pre-existing, oracle-verified
/// `trans_pairs()` semantics, not an independently-guessed group.
const SQUARE_PLANAR_ROTATIONS: [[u8; 4]; 8] = [
    [0, 1, 2, 3], // identity
    [2, 1, 0, 3], // (02)
    [0, 3, 2, 1], // (13)
    [2, 3, 0, 1], // (02)(13)
    [1, 0, 3, 2], // (01)(23)
    [1, 2, 3, 0], // (0123)
    [3, 0, 1, 2], // (0321) -- inverse of (0123)
    [3, 2, 1, 0], // (03)(12)
];

// ---------------------------------------------------------------------------
// Canonicalization
// ---------------------------------------------------------------------------

/// Detect a duplicate raw ligand id among the 4 slots, without heap
/// allocation (plain O(1)-bounded pairwise scan over a fixed 4-element
/// array -- no `HashSet`, deterministic regardless of any hashing).
fn find_duplicate(slots: [u32; 4]) -> Option<u32> {
    for i in 0..4 {
        for j in (i + 1)..4 {
            if slots[i] == slots[j] {
                return Some(slots[i]);
            }
        }
    }
    None
}

/// Reduce `ligand_order` to its canonical representative under `geometry`'s
/// proper rotation group: the lexicographically-smallest ordering reachable
/// by applying any element of the rotation group.
///
/// Fails closed with [`StereoGeometryError::DuplicateSlotId`] if the same
/// raw id occupies two slots -- not a valid 4-distinct-ligand arrangement,
/// and letting it through would make "lexicographically smallest" pick an
/// arbitrary, meaningless tie-break among rotation-equivalent duplicates.
///
/// `pub(crate)`: see [`StereoConfiguration`]'s doc for why this isn't public
/// yet (hardcoded `[u32; 4]` arity, deferred until a second geometry family
/// forces the real generalization).
pub(crate) fn canonicalize_configuration(
    geometry: StereoGeometry,
    ligand_order: [u32; 4],
) -> Result<CanonicalStereoConfiguration, StereoGeometryError> {
    if let Some(dup) = find_duplicate(ligand_order) {
        return Err(StereoGeometryError::DuplicateSlotId(dup));
    }
    let group = geometry.rotation_group();
    let mut best = apply(&group[0], ligand_order);
    for perm in &group[1..] {
        let candidate = apply(perm, ligand_order);
        if candidate < best {
            best = candidate;
        }
    }
    Ok(CanonicalStereoConfiguration {
        geometry,
        representative: best,
    })
}

/// `true` iff `a` and `b` describe the same physical arrangement -- i.e. one
/// is reachable from the other by a proper rotation of their (shared)
/// geometry. Configurations of different geometries are never equivalent.
/// Currently just `a == b` (both fields, including the private
/// `representative`, participate in equality), kept as a named function so
/// call sites read as a geometric claim rather than an incidental struct
/// comparison. `pub(crate)`: see [`StereoConfiguration`]'s doc.
pub(crate) fn equivalent_under_rotation(
    a: &CanonicalStereoConfiguration,
    b: &CanonicalStereoConfiguration,
) -> bool {
    a == b
}

// ---------------------------------------------------------------------------
// Bridge functions -- replace `chematic-smiles/src/canonical.rs`'s
// `permutation_is_odd` (tetrahedral) and `remap_square_planar`
// (square-planar) at their exact call sites.
// ---------------------------------------------------------------------------

/// Whether remapping a declared tetrahedral tag from `original` neighbor
/// order to `canonical` neighbor order requires flipping `@`<->`@@`.
///
/// Two orderings of the same 4 distinct ids differ by an even permutation
/// (no flip) iff they canonicalize to the same [`StereoGeometry::Tetrahedral`]
/// representative -- by definition, since [`TETRAHEDRAL_ROTATIONS`] (A4) is
/// exactly the even-permutation group. This is a direct restatement of
/// classic cycle-counting permutation parity through the rotation-orbit
/// abstraction, not a different rule -- see this module's
/// `tetrahedral_rotations_are_independently_confirmed_even_permutations`
/// test for a from-scratch cross-check.
///
/// Fails closed with [`StereoGeometryError::DuplicateSlotId`] if either
/// array has a repeated id, or [`StereoGeometryError::MismatchedLigandSet`]
/// if `original` and `canonical` don't name the same 4 ids as a set (the
/// caller -- `canonical.rs`'s `corrected_chirality` -- treats any `Err` the
/// same as its pre-existing "no verifiable order" pass-through-unchanged
/// fallback, matching the documented safe-no-op behavior for a 2-state tag;
/// see that call site's own comment for why this is the right fallback,
/// not a workaround).
///
/// The mismatched-set check matters, not just for symmetry with
/// [`remap_square_planar_tag`]'s analogous guard: without it, two arrays
/// naming *different* id sets (e.g. `[1,2,3,4]` vs `[1,2,3,5]`) each
/// canonicalize successfully (no duplicates in either one alone) to
/// necessarily-*unequal* representatives -- purely because they contain
/// different ids, not because one is an odd permutation of the other -- and
/// the naive `orig.representative != canon.representative` test would
/// misread that as "needs a parity flip," a confident-looking wrong answer
/// for malformed input rather than the honest "can't tell."
pub fn remap_tetrahedral_parity(
    original: [u32; 4],
    canonical: [u32; 4],
) -> Result<bool, StereoGeometryError> {
    // Each array's OWN internal-duplicate check first, so a genuinely
    // duplicated id is always reported as `DuplicateSlotId` (the more
    // specific diagnosis) even when the two arrays also happen to differ as
    // sets -- `MismatchedLigandSet` below is reserved for the case where
    // *neither* array has an internal duplicate but they still don't name
    // the same 4 ids.
    if let Some(dup) = find_duplicate(original) {
        return Err(StereoGeometryError::DuplicateSlotId(dup));
    }
    if let Some(dup) = find_duplicate(canonical) {
        return Err(StereoGeometryError::DuplicateSlotId(dup));
    }
    let mut sorted_original = original;
    sorted_original.sort_unstable();
    let mut sorted_canonical = canonical;
    sorted_canonical.sort_unstable();
    if sorted_original != sorted_canonical {
        return Err(StereoGeometryError::MismatchedLigandSet);
    }
    let orig = canonicalize_configuration(StereoGeometry::Tetrahedral, original)?;
    let canon = canonicalize_configuration(StereoGeometry::Tetrahedral, canonical)?;
    Ok(orig.representative() != canon.representative())
}

/// Convert `(tag, order)` into the geometry's own base-convention slot
/// array: an ordering where position 0/2 are one of `tag`'s trans-pairs and
/// position 1/3 are the other -- [`SquarePlanarPermutation::SP1`]'s own
/// convention (`(0,2)`/`(1,3)`), which is also exactly the partition
/// [`SQUARE_PLANAR_ROTATIONS`] stabilizes. Built directly from
/// `tag.trans_pairs()`, not a per-tag hand-written special case: given
/// `trans_pairs() = [(a,b),(c,d)]`, the reorder `[a,c,b,d]` always puts `a`
/// at 0 and `b` at 2 (the pair `(a,b)`) and `c` at 1 and `d` at 3 (the pair
/// `(c,d)`).
fn to_base_slots(tag: SquarePlanarPermutation, order: [u32; 4]) -> [u32; 4] {
    let [(a, b), (c, d)] = tag.trans_pairs();
    [
        order[a as usize],
        order[c as usize],
        order[b as usize],
        order[d as usize],
    ]
}

/// Remap a declared square-planar tag from `original` neighbor order to
/// `canonical` neighbor order -- the [`StereoGeometry`]-based counterpart to
/// the removed `remap_square_planar`.
///
/// `(tag, original)` and `(candidate, canonical)` describe the same physical
/// arrangement iff [`to_base_slots`]'s outputs for each are in the same
/// [`StereoGeometry::SquarePlanar`] rotation orbit (their unordered
/// trans-pair-of-ids partition, `{{slots[0],slots[2]}, {slots[1],slots[3]}}`,
/// is exactly what a `SquarePlanar`-orbit reduces to -- see
/// [`SQUARE_PLANAR_ROTATIONS`]'s doc). Tries all 3 tags against `canonical`
/// and returns the (unique, when one exists) match.
///
/// `None` -- never a guessed tag -- whenever no candidate matches: this
/// happens exactly when `original`/`canonical` don't name the same 4
/// distinct ids (mismatched id set, or a duplicate id in either array),
/// since in that case every candidate's canonicalization either errors
/// (`DuplicateSlotId`) or lands on a representative array containing ids
/// [`original`] didn't have, which can never equal `original`'s own
/// representative.
///
/// Unlike [`remap_tetrahedral_parity`], this function does NOT need an
/// explicit mismatched-id-set guard: it compares via
/// [`equivalent_under_rotation`] (full array equality, both `geometry` and
/// `representative`), not a boolean not-equal test, so a candidate whose
/// representative merely differs *for the wrong reason* (different id set,
/// not a real rotation-inequivalence) still correctly fails the equality
/// check and falls through to `None` rather than being misread as a match
/// -- verified by `remap_square_planar_tag_none_on_mismatched_id_set` below.
pub fn remap_square_planar_tag(
    tag: SquarePlanarPermutation,
    original: [u32; 4],
    canonical: [u32; 4],
) -> Option<SquarePlanarPermutation> {
    let canon_original =
        canonicalize_configuration(StereoGeometry::SquarePlanar, to_base_slots(tag, original))
            .ok()?;
    [
        SquarePlanarPermutation::SP1,
        SquarePlanarPermutation::SP2,
        SquarePlanarPermutation::SP3,
    ]
    .into_iter()
    .find(|&candidate| {
        let slots_candidate = to_base_slots(candidate, canonical);
        match canonicalize_configuration(StereoGeometry::SquarePlanar, slots_candidate) {
            Ok(canon_candidate) => equivalent_under_rotation(&canon_original, &canon_candidate),
            Err(_) => false,
        }
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// All 24 permutations of `[0,1,2,3]`, as `[u8;4]` for group-element
    /// tests and `[u32;4]` for configuration tests.
    fn permutations_of_4_u8() -> Vec<[u8; 4]> {
        let mut out = Vec::with_capacity(24);
        for a in 0..4u8 {
            for b in 0..4u8 {
                if b == a {
                    continue;
                }
                for c in 0..4u8 {
                    if c == a || c == b {
                        continue;
                    }
                    for d in 0..4u8 {
                        if d == a || d == b || d == c {
                            continue;
                        }
                        out.push([a, b, c, d]);
                    }
                }
            }
        }
        out
    }

    fn permutations_of_4_u32() -> Vec<[u32; 4]> {
        permutations_of_4_u8()
            .into_iter()
            .map(|p| [p[0] as u32, p[1] as u32, p[2] as u32, p[3] as u32])
            .collect()
    }

    /// Compose two rotation-table permutations so that `apply(compose(g,h),
    /// arr) == apply(g, apply(h, arr))` for all `arr` -- i.e.
    /// `compose(g,h)[i] = h[g[i]]`. Used only by the group-axiom tests
    /// below, not by production code.
    fn compose(g: &[u8; 4], h: &[u8; 4]) -> [u8; 4] {
        [
            h[g[0] as usize],
            h[g[1] as usize],
            h[g[2] as usize],
            h[g[3] as usize],
        ]
    }

    const IDENTITY: [u8; 4] = [0, 1, 2, 3];

    fn assert_is_group(table: &[[u8; 4]], expected_order: usize, name: &str) {
        assert_eq!(table.len(), expected_order, "{name}: wrong group order");

        // No duplicate rows.
        for i in 0..table.len() {
            for j in (i + 1)..table.len() {
                assert_ne!(table[i], table[j], "{name}: duplicate row at {i},{j}");
            }
        }

        // Every row is actually a permutation of 0..4.
        for row in table {
            let mut sorted = *row;
            sorted.sort_unstable();
            assert_eq!(
                sorted,
                [0, 1, 2, 3],
                "{name}: row {row:?} not a permutation"
            );
        }

        // Identity present.
        assert!(
            table.contains(&IDENTITY),
            "{name}: identity element missing"
        );

        // Closure.
        for g in table {
            for h in table {
                let gh = compose(g, h);
                assert!(
                    table.contains(&gh),
                    "{name}: not closed, compose({g:?},{h:?})={gh:?} not in table"
                );
            }
        }

        // Every element has an inverse in the table.
        for g in table {
            let has_inverse = table.iter().any(|h| compose(g, h) == IDENTITY);
            assert!(has_inverse, "{name}: {g:?} has no inverse in table");
        }
    }

    #[test]
    fn tetrahedral_rotations_form_a_group_of_order_12() {
        assert_is_group(&TETRAHEDRAL_ROTATIONS, 12, "TETRAHEDRAL_ROTATIONS");
    }

    #[test]
    fn square_planar_rotations_form_a_group_of_order_8() {
        assert_is_group(&SQUARE_PLANAR_ROTATIONS, 8, "SQUARE_PLANAR_ROTATIONS");
    }

    /// Independent second derivation for the tetrahedral table: brute-force
    /// cycle-counting parity over all 24 permutations, written from scratch
    /// here (not reusing `remap_tetrahedral_parity`'s internals, which don't
    /// even compute parity directly anymore) -- must select exactly the 12
    /// rows in [`TETRAHEDRAL_ROTATIONS`].
    #[test]
    fn tetrahedral_rotations_are_independently_confirmed_even_permutations() {
        fn is_odd(p: [u8; 4]) -> bool {
            let mut visited = [false; 4];
            let mut num_cycles = 0usize;
            for start in 0..4 {
                if !visited[start] {
                    num_cycles += 1;
                    let mut j = start;
                    while !visited[j] {
                        visited[j] = true;
                        j = p[j] as usize;
                    }
                }
            }
            (4 - num_cycles) % 2 == 1
        }

        let mut brute_force_even: Vec<[u8; 4]> = permutations_of_4_u8()
            .into_iter()
            .filter(|&p| !is_odd(p))
            .collect();
        brute_force_even.sort_unstable();

        let mut table_sorted = TETRAHEDRAL_ROTATIONS.to_vec();
        table_sorted.sort_unstable();

        assert_eq!(
            brute_force_even, table_sorted,
            "TETRAHEDRAL_ROTATIONS must equal the brute-force even-permutation set"
        );
    }

    /// Load-bearing check tying [`SQUARE_PLANAR_ROTATIONS`] to the
    /// pre-existing, oracle-verified `SquarePlanarPermutation::SP1::trans_pairs()`
    /// -- not just to a partition literal written by hand in this file. Every
    /// group element, applied to the reference order `[0,1,2,3]`, must
    /// preserve SP1's own trans-pair partition as an *unordered pair of
    /// unordered pairs* (individual pairs may swap which "half" they land in
    /// -- that's exactly the block-swapping half of the group).
    #[test]
    fn square_planar_rotations_stabilize_sp1_partition() {
        let reference: [u32; 4] = [0, 1, 2, 3];

        // Derive the reference partition FROM `trans_pairs()` at runtime,
        // via the same `to_base_slots` production code uses -- not a
        // hardcoded `{0,2}|{1,3}` literal in this test. `to_base_slots` for
        // SP1 must be a no-op: SP1 *is* the base convention by definition.
        let sp1_base = to_base_slots(SquarePlanarPermutation::SP1, reference);
        assert_eq!(
            sp1_base, reference,
            "SP1.trans_pairs() must already match the (0,2)/(1,3) base convention"
        );

        let partition_of = |arr: [u32; 4]| -> [[u32; 2]; 2] {
            let mut p1 = [arr[0], arr[2]];
            let mut p2 = [arr[1], arr[3]];
            p1.sort_unstable();
            p2.sort_unstable();
            let mut both = [p1, p2];
            both.sort_unstable();
            both
        };
        let expected = partition_of(sp1_base);

        for perm in &SQUARE_PLANAR_ROTATIONS {
            let rotated = apply(perm, reference);
            assert_eq!(
                partition_of(rotated),
                expected,
                "rotation {perm:?} does not stabilize SP1's trans-pair partition"
            );
        }
    }

    /// The flagship test: exactly 2 orbits for Tetrahedral, exactly 3 for
    /// SquarePlanar, over all 24 orderings of 4 distinct ids -- would catch
    /// a wrong group (e.g. 6 orbits from the naive order-4 in-plane-only
    /// square-planar group).
    #[test]
    fn orbit_counts_are_2_for_tetrahedral_and_3_for_square_planar() {
        for (geometry, expected_orbits) in [
            (StereoGeometry::Tetrahedral, 2),
            (StereoGeometry::SquarePlanar, 3),
        ] {
            let mut representatives: Vec<[u32; 4]> = permutations_of_4_u32()
                .into_iter()
                .map(|order| {
                    canonicalize_configuration(geometry, order)
                        .expect("4 distinct ids never duplicate")
                        .representative()
                })
                .collect();
            representatives.sort_unstable();
            representatives.dedup();
            assert_eq!(
                representatives.len(),
                expected_orbits,
                "{geometry:?}: expected {expected_orbits} orbits, got {}: {representatives:?}",
                representatives.len()
            );
        }
    }

    #[test]
    fn canonicalization_is_idempotent() {
        for geometry in [StereoGeometry::Tetrahedral, StereoGeometry::SquarePlanar] {
            for order in permutations_of_4_u32() {
                let once = canonicalize_configuration(geometry, order).unwrap();
                let twice = canonicalize_configuration(geometry, once.representative()).unwrap();
                assert_eq!(once, twice, "idempotence failed for {geometry:?} {order:?}");
            }
        }
    }

    #[test]
    fn duplicate_slot_id_fails_closed() {
        let err =
            canonicalize_configuration(StereoGeometry::Tetrahedral, [1, 2, 1, 3]).unwrap_err();
        assert_eq!(err, StereoGeometryError::DuplicateSlotId(1));
    }

    #[test]
    fn equivalent_under_rotation_matches_equality_and_respects_geometry() {
        let a = canonicalize_configuration(StereoGeometry::Tetrahedral, [1, 2, 3, 4]).unwrap();
        let b = canonicalize_configuration(
            StereoGeometry::Tetrahedral,
            apply(&TETRAHEDRAL_ROTATIONS[3], [1, 2, 3, 4]),
        )
        .unwrap();
        assert!(equivalent_under_rotation(&a, &b));

        let odd = canonicalize_configuration(StereoGeometry::Tetrahedral, [2, 1, 3, 4]).unwrap();
        assert!(!equivalent_under_rotation(&a, &odd));

        let sp = canonicalize_configuration(StereoGeometry::SquarePlanar, [1, 2, 3, 4]).unwrap();
        // Different geometry, same raw slots -- never equivalent.
        let te = canonicalize_configuration(StereoGeometry::Tetrahedral, [1, 2, 3, 4]).unwrap();
        assert!(!equivalent_under_rotation(&sp, &te));
    }

    // -------------------------------------------------------------------
    // remap_tetrahedral_parity
    // -------------------------------------------------------------------

    #[test]
    fn remap_tetrahedral_parity_matches_hand_cases() {
        // Identity: never odd.
        assert!(!remap_tetrahedral_parity([1, 2, 3, 4], [1, 2, 3, 4]).unwrap());
        // Single transposition: odd.
        assert!(remap_tetrahedral_parity([1, 2, 3, 4], [2, 1, 3, 4]).unwrap());
        // Double transposition: even.
        assert!(!remap_tetrahedral_parity([1, 2, 3, 4], [2, 1, 4, 3]).unwrap());
        // 3-cycle: even.
        assert!(!remap_tetrahedral_parity([1, 2, 3, 4], [2, 3, 1, 4]).unwrap());
    }

    #[test]
    fn remap_tetrahedral_parity_fails_closed_on_duplicate() {
        assert_eq!(
            remap_tetrahedral_parity([1, 1, 3, 4], [1, 2, 3, 4]).unwrap_err(),
            StereoGeometryError::DuplicateSlotId(1)
        );
    }

    /// The bug an independent review caught: `original`/`canonical` naming
    /// *different* id sets (neither internally duplicated) must fail closed
    /// with `MismatchedLigandSet`, not silently return a confident-looking
    /// `Ok(true)`/`Ok(false)`. Before this check existed,
    /// `remap_tetrahedral_parity([1,2,3,4], [1,2,3,5])` returned `Ok(true)`:
    /// the two arrays canonicalize to necessarily-unequal representatives
    /// (they contain different ids), which the naive `!=` parity test
    /// misread as "needs a flip" -- a wrong answer for malformed input, not
    /// an honest "can't tell."
    #[test]
    fn remap_tetrahedral_parity_fails_closed_on_mismatched_ligand_set() {
        assert_eq!(
            remap_tetrahedral_parity([1, 2, 3, 4], [1, 2, 3, 5]).unwrap_err(),
            StereoGeometryError::MismatchedLigandSet
        );
        // Same multiset, different order -- must NOT be flagged as
        // mismatched (this is the ordinary, common case this function
        // exists to compute a real parity answer for).
        assert!(remap_tetrahedral_parity([1, 2, 3, 4], [1, 2, 3, 4]).is_ok());
        assert!(remap_tetrahedral_parity([1, 2, 3, 4], [4, 3, 2, 1]).is_ok());
    }

    // -------------------------------------------------------------------
    // remap_square_planar_tag
    // -------------------------------------------------------------------

    #[test]
    fn remap_square_planar_tag_identity_is_a_no_op() {
        for tag in [
            SquarePlanarPermutation::SP1,
            SquarePlanarPermutation::SP2,
            SquarePlanarPermutation::SP3,
        ] {
            assert_eq!(
                remap_square_planar_tag(tag, [10, 20, 30, 40], [10, 20, 30, 40]),
                Some(tag)
            );
        }
    }

    /// Duplicate-*chemistry*, distinct-*ids* case, shaped like
    /// cisplatin/transplatin's real coordination chemistry: `2xCl + 2xNH3`
    /// ligand composition (whether this particular slot arrangement happens
    /// to be the cis or trans isomer specifically depends on which of the
    /// two identical-composition slots the caller places where, which this
    /// test deliberately does not commit to -- see the note below). This
    /// module only ever sees 4 raw slot ids -- it has no concept of "these
    /// two are chemically the same ligand" at all (see the RFC's
    /// duplicate-ligand section) -- so two chemically-identical Cl ligands
    /// still get two distinct ids here (`CL1`/`CL2`), same for the two NH3
    /// nitrogens (`N1`/`N2`). Proves directly, at the geometry-module level
    /// (independent of `square_planar_stereo.rs`'s end-to-end
    /// `cisplatin_and_transplatin_have_distinct_canonical_identity`, which
    /// checks the same property through the full parser/writer instead):
    /// SP1/SP2/SP3, applied to this one fixed 2xCl+2xN slot assignment,
    /// canonicalize to *three different* representatives -- chemical-
    /// identity duplication never collapses any of them into one orbit.
    #[test]
    fn duplicate_chemistry_distinct_ids_keeps_cisplatin_transplatin_shaped_tags_distinct() {
        const CL1: u32 = 101;
        const CL2: u32 = 102;
        const N1: u32 = 103;
        const N2: u32 = 104;
        // Deliberately not naming which of SP1/SP2/SP3 is "the cis one" or
        // "the trans one" for this specific order -- that reading is a real
        // but easy-to-get-backwards derived fact (which pair of positions a
        // given tag makes trans depends on both the tag AND which ligand
        // sits at which position), and getting it right or wrong doesn't
        // change what this test actually checks: all 3 tags, applied to the
        // SAME slot assignment, must land in 3 distinct orbits.
        let order: [u32; 4] = [CL1, N1, CL2, N2];

        let sp1 = canonicalize_configuration(
            StereoGeometry::SquarePlanar,
            to_base_slots(SquarePlanarPermutation::SP1, order),
        )
        .unwrap();
        let sp2 = canonicalize_configuration(
            StereoGeometry::SquarePlanar,
            to_base_slots(SquarePlanarPermutation::SP2, order),
        )
        .unwrap();
        let sp3 = canonicalize_configuration(
            StereoGeometry::SquarePlanar,
            to_base_slots(SquarePlanarPermutation::SP3, order),
        )
        .unwrap();

        assert!(
            !equivalent_under_rotation(&sp1, &sp2),
            "SP1-shaped and SP2-shaped configurations must NOT collapse to one orbit even \
             though both slot assignments repeat 2xCl+2xN chemistry"
        );
        assert!(!equivalent_under_rotation(&sp1, &sp3));
        assert!(!equivalent_under_rotation(&sp2, &sp3));

        // remap_square_planar_tag itself, the actual production bridge
        // function, must also keep them distinct end to end (original order
        // unpermuted, i.e. "canonical" == "original" -- this is the
        // identity-remap case, so the returned tag must be the same tag
        // fed in, for both SP1 and SP2, distinctly).
        assert_eq!(
            remap_square_planar_tag(SquarePlanarPermutation::SP1, order, order),
            Some(SquarePlanarPermutation::SP1)
        );
        assert_eq!(
            remap_square_planar_tag(SquarePlanarPermutation::SP2, order, order),
            Some(SquarePlanarPermutation::SP2)
        );
    }

    #[test]
    fn remap_square_planar_tag_full_24_by_3_table() {
        // Direct unit-level version of the same 24-permutations x 3-tags
        // sweep `square_planar_stereo.rs`'s end-to-end oracle test performs,
        // exercised against `remap_square_planar_tag` itself rather than
        // through the SMILES parser/writer. Reference/"canonical" order is
        // fixed at the identity [0,1,2,3] (ligand ids 0..3 by value), same
        // convention `square_planar_stereo.rs`'s own `predict` helper uses.
        let tags = [
            SquarePlanarPermutation::SP1,
            SquarePlanarPermutation::SP2,
            SquarePlanarPermutation::SP3,
        ];
        let mut checked = 0;
        for order in permutations_of_4_u32() {
            for &tag in &tags {
                // `order[i]` = which ligand id sits at original position i.
                let predicted = tags
                    .into_iter()
                    .find(|&candidate| {
                        let a = to_base_slots(tag, order);
                        let b = to_base_slots(candidate, [0, 1, 2, 3]);
                        let ca =
                            canonicalize_configuration(StereoGeometry::SquarePlanar, a).unwrap();
                        let cb =
                            canonicalize_configuration(StereoGeometry::SquarePlanar, b).unwrap();
                        equivalent_under_rotation(&ca, &cb)
                    })
                    .expect("exactly one of the 3 tags must match (3 orbits, 3 tags)");
                assert_eq!(
                    remap_square_planar_tag(tag, order, [0, 1, 2, 3]),
                    Some(predicted),
                    "order={order:?} tag={tag:?}"
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 24 * 3);
    }

    #[test]
    fn remap_square_planar_tag_none_on_mismatched_id_set() {
        // `canonical` doesn't contain the same 4 ids as `original` (5
        // replaces 4) -- must fail closed, never guess a tag.
        assert_eq!(
            remap_square_planar_tag(SquarePlanarPermutation::SP1, [1, 2, 3, 4], [1, 2, 3, 5]),
            None
        );
    }

    #[test]
    fn remap_square_planar_tag_none_on_duplicate_original() {
        assert_eq!(
            remap_square_planar_tag(SquarePlanarPermutation::SP1, [1, 1, 3, 4], [1, 2, 3, 4]),
            None
        );
    }

    // -------------------------------------------------------------------
    // StereoConfiguration::new
    // -------------------------------------------------------------------

    #[test]
    fn stereo_configuration_new_rejects_duplicate_slot_id() {
        // Struct-literal construction is unavailable outside this module
        // (private fields) -- `new` is the only way in, and it must run the
        // same duplicate check `canonicalize_configuration` does, so no
        // `StereoConfiguration` can ever exist in an already-invalid state.
        assert_eq!(
            StereoConfiguration::new(StereoGeometry::Tetrahedral, [1, 2, 1, 3]).unwrap_err(),
            StereoGeometryError::DuplicateSlotId(1)
        );
    }

    #[test]
    fn stereo_configuration_new_accepts_distinct_ids() {
        assert!(StereoConfiguration::new(StereoGeometry::Tetrahedral, [1, 2, 3, 4]).is_ok());
    }

    // -------------------------------------------------------------------
    // StereoConfiguration::renumber
    // -------------------------------------------------------------------

    #[test]
    fn renumber_remaps_every_slot() {
        let cfg = StereoConfiguration::new(StereoGeometry::Tetrahedral, [1, 2, 3, 4]).unwrap();
        let renumbered = cfg
            .renumber(|id| Some(id * 10))
            .expect("total map succeeds");
        assert_eq!(renumbered.geometry, StereoGeometry::Tetrahedral);
        assert_eq!(renumbered.slots, [10, 20, 30, 40]);
    }

    #[test]
    fn renumber_fails_closed_on_unmapped_id() {
        let cfg = StereoConfiguration::new(StereoGeometry::Tetrahedral, [1, 2, 3, 4]).unwrap();
        let err = cfg
            .renumber(|id| if id == 3 { None } else { Some(id) })
            .unwrap_err();
        assert_eq!(err, StereoGeometryError::UnknownLigandId(3));
    }

    /// The bug an independent review caught: a non-injective `id_map` --
    /// total on every input slot, but mapping two *distinct* input ids to
    /// the *same* output id -- must fail closed, not silently manufacture a
    /// duplicate that could never have been accepted by
    /// [`StereoConfiguration::new`] directly. Before this check existed,
    /// mapping both `1` and `2` to `10` would succeed and return
    /// `slots: [10, 10, 3, 4]`.
    #[test]
    fn renumber_rejects_id_map_that_creates_a_duplicate() {
        let cfg = StereoConfiguration::new(StereoGeometry::Tetrahedral, [1, 2, 3, 4]).unwrap();
        // Both 1 and 2 map to 10 (non-injective); 3 and 4 map to themselves.
        let err = cfg
            .renumber(|id| Some(if id == 1 || id == 2 { 10 } else { id }))
            .unwrap_err();
        assert_eq!(err, StereoGeometryError::DuplicateSlotId(10));
    }

    /// Renumbering invariance: if two configurations describe the same
    /// physical arrangement (same rotation orbit), applying the *same*
    /// (possibly non-order-preserving) renumbering map to both must still
    /// leave them describing the same physical arrangement.
    #[test]
    fn renumber_preserves_rotation_equivalence() {
        for geometry in [StereoGeometry::Tetrahedral, StereoGeometry::SquarePlanar] {
            let a = StereoConfiguration {
                geometry,
                slots: [1, 2, 3, 4],
            };
            for perm in geometry.rotation_group() {
                let b = StereoConfiguration {
                    geometry,
                    slots: apply(perm, a.slots),
                };
                // Sanity: same orbit before renumbering.
                let ca = canonicalize_configuration(geometry, a.slots).unwrap();
                let cb = canonicalize_configuration(geometry, b.slots).unwrap();
                assert!(equivalent_under_rotation(&ca, &cb));

                // Deliberately non-monotonic bijection on {1,2,3,4}.
                let map = |id: u32| -> Option<u32> {
                    match id {
                        1 => Some(40),
                        2 => Some(5),
                        3 => Some(77),
                        4 => Some(1),
                        _ => None,
                    }
                };
                let a2 = a.renumber(map).unwrap();
                let b2 = b.renumber(map).unwrap();
                let ca2 = canonicalize_configuration(a2.geometry, a2.slots).unwrap();
                let cb2 = canonicalize_configuration(b2.geometry, b2.slots).unwrap();
                assert!(
                    equivalent_under_rotation(&ca2, &cb2),
                    "renumbering broke rotation-equivalence for {geometry:?} perm {perm:?}"
                );
            }
        }
    }
}
