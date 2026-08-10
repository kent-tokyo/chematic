//! Stereochemistry quality validation.
//!
//! Detects common stereochemistry errors that arise when reading molecular
//! files or manually editing structures:
//!
//! - [`StereoErrorKind::ImpossibleCenter`] — a chirality annotation (`@`/`@@`)
//!   is present on an atom with fewer than 4 distinct heavy-atom neighbours.
//! - [`StereoErrorKind::ConflictingWedges`] — the same atom is the base of two
//!   or more stereo bonds (Up or Down) pointing in opposite directions.
//! - [`StereoErrorKind::RedundantStereo`] — the annotated atom is topologically
//!   equivalent to a neighbour (same Morgan rank), so the stereo specification
//!   is chemically meaningless.

use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule};
use std::fmt;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Kind of stereochemistry error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoErrorKind {
    /// Chirality annotation on an atom with < 4 heavy-atom neighbours (or all
    /// neighbours identical).
    ImpossibleCenter,
    /// Two or more Up/Down bonds originate from the same atom with conflicting
    /// directions (both Up and Down from the same center).
    ConflictingWedges,
    /// Stereo annotation on a topologically symmetric atom (all neighbours
    /// have the same Morgan rank — no priority ordering possible).
    RedundantStereo,
}

/// A detected stereochemistry error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoError {
    /// 0-based atom index of the problematic center.
    pub atom_idx: usize,
    pub kind: StereoErrorKind,
}

impl fmt::Display for StereoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind_str = match &self.kind {
            StereoErrorKind::ImpossibleCenter => {
                "impossible stereocenter (< 4 distinct neighbours)"
            }
            StereoErrorKind::ConflictingWedges => "conflicting wedge directions",
            StereoErrorKind::RedundantStereo => "redundant stereo on symmetric atom",
        };
        write!(f, "atom {}: {}", self.atom_idx, kind_str)
    }
}

impl std::error::Error for StereoError {}

/// Summary of stereocenters in a molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoCompleteness {
    /// Stereocenters with an explicit `@`/`@@` annotation.
    pub specified: usize,
    /// Stereocenters with 4 distinct heavy-atom neighbours but no annotation.
    pub unspecified: usize,
    /// `specified + unspecified`
    pub total_centers: usize,
}

// ---------------------------------------------------------------------------
// Internal: lightweight Morgan ranks (avoids chematic-smiles dependency)
// ---------------------------------------------------------------------------

/// Compute simple Morgan connectivity ranks for atoms in `mol`.
/// Uses initial invariant = atomic_number * 1_000_000 + charge_term * 1000 + degree.
fn simple_morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();
    let mut ranks: Vec<u64> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let deg = mol.neighbors(idx).count() as i64;
            // `atom.charge` (i8) sign-extends on a plain `as u64` cast for
            // negative values (e.g. -1i8 as u64 == u64::MAX), which made the
            // old `atom.charge as u64 * 1000` overflow `u64` unconditionally
            // for any negatively-charged atom (issue #267). Computing the
            // whole invariant in i64 and reinterpreting as u64 only once, at
            // the end, avoids that overflow while leaving the value
            // bit-for-bit identical to before for the (already-correct)
            // charge >= 0 case. `atomic_number()` is always >= 1 (never a
            // 0/wildcard sentinel here), so `an * 1_000_000` always dominates
            // even the most extreme i8 charge magnitude (128_000 at most),
            // keeping the final sum non-negative for every realistic and
            // every representable i8 charge.
            let an = atom.element.atomic_number() as i64;
            let charge = atom.charge as i64;
            (an * 1_000_000 + charge * 1000 + deg) as u64
        })
        .collect();

    let hash_round = |r: u64, nbrs: &[u64]| -> u64 {
        let mut h: u64 = 14695981039346656037u64;
        let prime: u64 = 1099511628211u64;
        h ^= r;
        h = h.wrapping_mul(prime);
        for &nb in nbrs {
            h ^= nb;
            h = h.wrapping_mul(prime);
        }
        h
    };

    for _ in 0..(n + 2) {
        let old_distinct = {
            let mut v = ranks.clone();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        let new_ranks: Vec<u64> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let mut nb_ranks: Vec<u64> = mol
                    .neighbors(idx)
                    .map(|(nb, _)| ranks[nb.0 as usize])
                    .collect();
                nb_ranks.sort_unstable();
                hash_round(ranks[i], &nb_ranks)
            })
            .collect();
        let new_distinct = {
            let mut v = new_ranks.clone();
            v.sort_unstable();
            v.dedup();
            v.len()
        };
        ranks = new_ranks;
        if new_distinct <= old_distinct {
            break;
        }
    }

    // Normalise to consecutive ordinals.
    let mut sorted = ranks.clone();
    sorted.sort_unstable();
    sorted.dedup();
    ranks
        .iter()
        .map(|r| sorted.partition_point(|&u| u < *r) as u64)
        .collect()
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Validate the stereochemistry of `mol` and return any errors found.
///
/// An empty `Vec` means the stereo annotations are chemically consistent.
pub fn validate_stereo(mol: &Molecule) -> Vec<StereoError> {
    let ranks = simple_morgan_ranks(mol);
    let mut errors = Vec::new();

    for (idx, atom) in mol.atoms() {
        let i = idx.0 as usize;

        // Only inspect atoms with explicit chirality.
        if atom.chirality == Chirality::None {
            continue;
        }

        let heavy_neighbors: Vec<AtomIdx> = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .map(|(nb, _)| nb)
            .collect();

        // Rule 1: ImpossibleCenter — fewer than 4 distinct heavy neighbours.
        // (3 heavy + 1 implicit H is OK; < 3 heavy is definitely impossible.)
        let implicit_h = chematic_core::implicit_hcount(mol, idx);
        let total_groups = heavy_neighbors.len() + implicit_h as usize;
        if total_groups < 4 {
            errors.push(StereoError {
                atom_idx: i,
                kind: StereoErrorKind::ImpossibleCenter,
            });
            continue; // no point checking further
        }

        // Rule 2: ConflictingWedges — Up and Down bonds from same center.
        let mut has_up = false;
        let mut has_down = false;
        for (_, bid) in mol.neighbors(idx) {
            let bond = mol.bond(bid);
            if bond.atom1 == idx {
                match bond.order {
                    BondOrder::Up => has_up = true,
                    BondOrder::Down => has_down = true,
                    _ => {}
                }
            }
        }
        if has_up && has_down {
            errors.push(StereoError {
                atom_idx: i,
                kind: StereoErrorKind::ConflictingWedges,
            });
        }

        // Rule 3: RedundantStereo — all heavy neighbours have the same rank.
        if !heavy_neighbors.is_empty() {
            let first_rank = ranks[heavy_neighbors[0].0 as usize];
            let all_same = heavy_neighbors
                .iter()
                .all(|nb| ranks[nb.0 as usize] == first_rank);
            // Also check the center itself doesn't break ties via implicit H.
            if all_same && implicit_h == 0 {
                errors.push(StereoError {
                    atom_idx: i,
                    kind: StereoErrorKind::RedundantStereo,
                });
            }
        }
    }

    errors
}

/// Return the tetrahedral stereocenter candidates in `mol` that this
/// classifier recognises: an sp3 atom with 4 distinct heavy-atom-or-
/// implicit-H groups, paired with whether it carries an explicit `@`/`@@`
/// chirality annotation (`true` = specified, `false` = a valid candidate
/// left unannotated).
///
/// This is the single source of truth for stereocenter classification;
/// [`stereo_completeness`] is defined in terms of it.
pub fn stereo_centers(mol: &Molecule) -> Vec<(AtomIdx, bool)> {
    let ranks = simple_morgan_ranks(mol);
    let mut centers = Vec::new();

    for (idx, atom) in mol.atoms() {
        // Skip aromatics and obvious non-centers.
        if atom.aromatic {
            continue;
        }

        let heavy_nbs: Vec<AtomIdx> = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .map(|(nb, _)| nb)
            .collect();
        let implicit_h = chematic_core::implicit_hcount(mol, idx) as usize;
        let groups = heavy_nbs.len() + implicit_h;

        if groups != 4 {
            continue;
        } // only tetrahedral candidates

        // Check all neighbours have distinct ranks (including implicit H as a
        // sentinel rank). `simple_morgan_ranks` normalises to consecutive
        // ordinals starting at 0 (see its tail: `partition_point(|&u| u <
        // *r)`), so a real heavy-atom neighbour can legitimately carry rank
        // 0 -- it's the ordinary "lowest invariant in the molecule" rank,
        // not a reserved value. Using the literal `0` as the implicit-H
        // sentinel therefore collided with real rank-0 neighbours (issue
        // #267's follow-up bug): `dedup()` merged the two, `sorted.len()`
        // dropped below 4, and a genuine 4-distinct-group stereocenter was
        // silently skipped. The maximum normalised rank is (number of
        // distinct invariants - 1), which is always <= atom_count() - 1, so
        // `atom_count()` itself is never a reachable real rank and is safe
        // to use as the sentinel here.
        let implicit_h_rank_sentinel = mol.atom_count() as u64;
        let mut nb_ranks: Vec<u64> = heavy_nbs.iter().map(|nb| ranks[nb.0 as usize]).collect();
        if implicit_h > 0 {
            nb_ranks.push(implicit_h_rank_sentinel);
        }

        let mut sorted = nb_ranks.clone();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() < 4 {
            continue;
        } // symmetric neighbours — not a stereocenter

        centers.push((idx, atom.chirality != Chirality::None));
    }

    centers
}

/// Summarise how many stereocenters in `mol` have been specified vs left open.
///
/// A potential stereocenter is an sp3 atom with 4 distinct heavy-atom
/// neighbours (counting one implicit H as a distinct group when present).
pub fn stereo_completeness(mol: &Molecule) -> StereoCompleteness {
    let centers = stereo_centers(mol);
    let specified = centers
        .iter()
        .filter(|(_, is_specified)| *is_specified)
        .count();
    let unspecified = centers.len() - specified;

    StereoCompleteness {
        specified,
        unspecified,
        total_centers: centers.len(),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_valid_chiral_center_no_errors() {
        // L-alanine: valid R/S center with 4 distinct groups.
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let errors = validate_stereo(&mol);
        assert!(
            errors.is_empty(),
            "L-alanine should have no stereo errors: {errors:?}"
        );
    }

    #[test]
    fn test_impossible_center_explicit_h_zero() {
        // A carbon with chirality annotation, 1 heavy bond, and explicit H=0
        // gives total_groups = 1 → ImpossibleCenter.
        use chematic_core::{Atom, BondOrder, Chirality, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let mut c = Atom::new(Element::C);
        c.chirality = Chirality::CounterClockwise;
        c.hydrogen_count = Some(0); // force 0 implicit H
        let ci = b.add_atom(c);
        let cl = b.add_atom(Atom::new(Element::CL));
        b.add_bond(ci, cl, BondOrder::Single).unwrap();
        let mol = b.build();
        let errors = validate_stereo(&mol);
        assert!(
            errors
                .iter()
                .any(|e| e.atom_idx == 0 && e.kind == StereoErrorKind::ImpossibleCenter),
            "should detect ImpossibleCenter (1 group total): {errors:?}"
        );
    }

    #[test]
    fn test_stereo_completeness_alanine() {
        // L-alanine has 1 specified stereocenter, 0 unspecified.
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let sc = stereo_completeness(&mol);
        assert_eq!(sc.specified, 1);
        assert_eq!(sc.unspecified, 0);
        assert_eq!(sc.total_centers, 1);
    }

    #[test]
    fn test_stereo_completeness_unspecified() {
        // Alanine without stereo annotation: 1 unspecified center.
        let mol = parse("NC(C)C(=O)O").unwrap();
        let sc = stereo_completeness(&mol);
        assert_eq!(sc.specified, 0);
        assert!(sc.unspecified >= 1, "should detect unspecified center");
    }

    #[test]
    fn test_no_centers_in_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let sc = stereo_completeness(&mol);
        assert_eq!(sc.total_centers, 0);
    }

    // Regression tests for issue #267: `atom.charge as u64` sign-extends for
    // negative i8 charges (e.g. -1i8 as u64 == u64::MAX), which made the
    // Morgan-rank invariant's `* 1000` multiply overflow `u64` -- panicking
    // in debug builds and silently corrupting the invariant in release.

    #[test]
    fn test_stereo_completeness_negative_charge_no_panic() {
        // Acetate: the [O-] atom must not trigger an overflow panic.
        let acetate = parse("CC(=O)[O-]").unwrap();
        let sc = stereo_completeness(&acetate);
        assert_eq!(sc.total_centers, 0);
    }

    #[test]
    fn test_stereo_completeness_positive_charge_no_panic() {
        // Small positive i8 charge never overflowed, but keep it as a
        // regression fixture per the issue's exact repro.
        let cation = parse("C[NH3+]").unwrap();
        let sc = stereo_completeness(&cation);
        assert_eq!(sc.total_centers, 0);
    }

    #[test]
    fn test_stereo_completeness_mixed_salt_no_panic() {
        // Disconnected salt combining both a negative and a positive charge.
        let salt = parse("CC(=O)[O-].C[NH3+]").unwrap();
        let sc = stereo_completeness(&salt);
        assert_eq!(sc.total_centers, 0);
    }

    #[test]
    fn test_stereo_completeness_doubly_negative_charge_no_panic() {
        // A doubly-deprotonated phosphate: charge -2, more extreme than the
        // issue's -1 repro, to make sure the fix isn't a -1-only special case.
        let phosphate = parse("[O-]P(=O)([O-])OC").unwrap();
        let sc = stereo_completeness(&phosphate);
        assert_eq!(sc.total_centers, 0);
        // Also exercise the other public entry point sharing the same helper.
        let errors = validate_stereo(&phosphate);
        assert!(errors.is_empty());
    }

    // Regression tests for the implicit-H rank-0 sentinel collision (issue
    // #267 follow-up, distinct from the overflow bug above): `simple_morgan_ranks`
    // normalises ranks to consecutive ordinals starting at 0, so an ordinary
    // heavy-atom neighbour can legitimately carry rank 0. Using the literal
    // `0` as the implicit-H stand-in collided with such a neighbour, `dedup()`
    // merged them, and a genuine 4-distinct-group stereocenter was silently
    // dropped (`specified` undercounted).

    #[test]
    fn test_stereo_completeness_rank_zero_collision_chain() {
        // Atom indices: 0=C(methyl) 1=C@@H(chiral) 2=Cl 3=C(quaternary) 4=Br
        // 5=F 6=I. Verified via `simple_morgan_ranks`: ranks = [0, 6, 2, 5, 4,
        // 1, 3] -- the chiral atom's methyl neighbour (atom 0) is the
        // lowest-invariant atom in the whole molecule and normalises to rank
        // 0, colliding with the old implicit-H sentinel. Before the fix this
        // atom was dropped entirely (specified=0, unspecified=1, total=1,
        // counting only the separate, unrelated unspecified center at the
        // quaternary carbon 3, whose 4 heavy neighbours -- 1, Br, F, I -- are
        // trivially distinct and is not itself affected by this bug); after
        // the fix both stereocenters are counted.
        let mol = parse("C[C@@H](Cl)C(Br)(F)I").unwrap();
        let sc = stereo_completeness(&mol);
        assert_eq!(
            sc.specified, 1,
            "annotated chiral carbon must not be dropped by the rank-0 collision: {sc:?}"
        );
        assert_eq!(
            sc.total_centers, 2,
            "expected 2 stereocenters total (1 fixed + 1 pre-existing, bug-unrelated \
             unspecified center at the quaternary carbon): {sc:?}"
        );
    }

    #[test]
    fn test_stereo_completeness_rank_zero_collision_different_neighbor() {
        // Same collision class as the chain test above, but with a
        // structurally different colliding neighbour to confirm the fix
        // isn't specific to "a bare terminal methyl happens to be rank 0".
        // Atom indices: 0=O 1=C(CH2, bonded to O) 2=C@@H(chiral) 3=N
        // 4=C(quaternary) 5=Br 6=F 7=I. Verified via `simple_morgan_ranks`:
        // ranks = [1, 0, 5, 3, 7, 6, 2, 4] -- here it's the chiral atom's
        // *substituted* CH2-OH neighbour (atom 1, not a plain methyl) that
        // lands on rank 0 and collides with the implicit-H sentinel.
        let mol = parse("OC[C@@H](N)C(Br)(F)I").unwrap();
        let sc = stereo_completeness(&mol);
        assert_eq!(
            sc.specified, 1,
            "annotated chiral carbon must not be dropped by the rank-0 collision: {sc:?}"
        );
        assert_eq!(
            sc.total_centers, 2,
            "expected 2 stereocenters total (1 fixed + 1 pre-existing, bug-unrelated \
             unspecified center at the quaternary carbon): {sc:?}"
        );
    }

    #[test]
    fn test_stereo_centers_mixed_specified_and_unspecified() {
        // Atom 0 ([C@]) carries an explicit chirality annotation and has 4
        // distinct fully-explicit heavy neighbours (F, Cl, Br, atom 4) ->
        // specified stereocenter.
        // Atom 4 (C(I)(N)O) has no annotation but also has 4 distinct heavy
        // neighbours (atom 0, I, N, O) -> unspecified candidate.
        let mol = parse("[C@](F)(Cl)(Br)C(I)(N)O").unwrap();
        let centers = stereo_centers(&mol);

        assert_eq!(
            centers.len(),
            2,
            "expected exactly 2 stereocenter candidates: {centers:?}"
        );
        assert!(
            centers.contains(&(AtomIdx(0), true)),
            "atom 0 should be a specified stereocenter: {centers:?}"
        );
        assert!(
            centers.contains(&(AtomIdx(4), false)),
            "atom 4 should be an unspecified stereocenter candidate: {centers:?}"
        );

        // stereo_completeness must agree exactly (single source of truth).
        let sc = stereo_completeness(&mol);
        assert_eq!(sc.specified, 1);
        assert_eq!(sc.unspecified, 1);
        assert_eq!(sc.total_centers, 2);
    }

    // Regression tests, on `stereo_centers` itself, for the two bugs fixed
    // upstream of this PR's rebase onto `main` (issue #267 overflow fix,
    // commit a99fc9b; implicit-H rank-0 sentinel collision fix, commit
    // 5790bb0). Both were previously exercised only via `stereo_completeness`
    // aggregate counts; these confirm the new `stereo_centers` API itself is
    // correct now that it's the single source of truth both bugs lived in.

    #[test]
    fn test_stereo_centers_negative_formal_charge_no_panic() {
        // Issue #267's exact repro: acetate's [O-] atom must not overflow
        // u64 in simple_morgan_ranks (it used to sign-extend and panic in
        // debug builds). No stereocenters expected -- just confirm
        // stereo_centers runs to completion with a sensible, empty result.
        let acetate = parse("CC(=O)[O-]").unwrap();
        let centers = stereo_centers(&acetate);
        assert!(
            centers.is_empty(),
            "acetate has no stereocenters: {centers:?}"
        );
    }

    #[test]
    fn test_stereo_centers_rank_zero_sentinel_collision_fixed() {
        // This PR's own body cited this exact repro as a known, pre-existing
        // limitation: atom 1 is a real, `@@`-annotated stereocenter that
        // `stereo_centers` used to silently drop because its methyl
        // neighbour (atom 0) normalises to Morgan rank 0, the same sentinel
        // value `stereo_centers` used to stand in for the implicit H. Fixed
        // upstream in commit 5790bb0; confirm atom 1 is now correctly
        // reported as (AtomIdx(1), true) directly from stereo_centers.
        let mol = parse("C[C@@H](Cl)C(Br)(F)I").unwrap();
        let centers = stereo_centers(&mol);
        assert!(
            centers.contains(&(AtomIdx(1), true)),
            "atom 1 must be reported as a specified stereocenter: {centers:?}"
        );
    }
}
