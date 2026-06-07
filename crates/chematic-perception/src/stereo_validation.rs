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

use std::fmt;
use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule};

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
            StereoErrorKind::ImpossibleCenter => "impossible stereocenter (< 4 distinct neighbours)",
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
/// Uses initial invariant = atomic_number * 1000 + degree.
fn simple_morgan_ranks(mol: &Molecule) -> Vec<u64> {
    use std::collections::HashSet;
    let n = mol.atom_count();
    let mut ranks: Vec<u64> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let deg = mol.neighbors(idx).count() as u64;
            atom.element.atomic_number() as u64 * 1_000_000 + atom.charge as u64 * 1000 + deg
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
        let old_distinct = { let mut v = ranks.clone(); v.sort_unstable(); v.dedup(); v.len() };
        let new_ranks: Vec<u64> = (0..n).map(|i| {
            let idx = AtomIdx(i as u32);
            let mut nb_ranks: Vec<u64> = mol.neighbors(idx).map(|(nb, _)| ranks[nb.0 as usize]).collect();
            nb_ranks.sort_unstable();
            hash_round(ranks[i], &nb_ranks)
        }).collect();
        let new_distinct = { let mut v = new_ranks.clone(); v.sort_unstable(); v.dedup(); v.len() };
        ranks = new_ranks;
        if new_distinct <= old_distinct { break; }
    }

    // Normalise to consecutive ordinals.
    let mut sorted = ranks.clone();
    sorted.sort_unstable();
    sorted.dedup();
    ranks.iter().map(|r| sorted.partition_point(|&u| u < *r) as u64).collect()
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
            errors.push(StereoError { atom_idx: i, kind: StereoErrorKind::ImpossibleCenter });
            continue; // no point checking further
        }

        // Rule 2: ConflictingWedges — Up and Down bonds from same center.
        let mut has_up = false;
        let mut has_down = false;
        for (_, bid) in mol.neighbors(idx) {
            let bond = mol.bond(bid);
            if bond.atom1 == idx {
                match bond.order {
                    BondOrder::Up   => has_up   = true,
                    BondOrder::Down => has_down = true,
                    _ => {}
                }
            }
        }
        if has_up && has_down {
            errors.push(StereoError { atom_idx: i, kind: StereoErrorKind::ConflictingWedges });
        }

        // Rule 3: RedundantStereo — all heavy neighbours have the same rank.
        if !heavy_neighbors.is_empty() {
            let first_rank = ranks[heavy_neighbors[0].0 as usize];
            let all_same = heavy_neighbors.iter().all(|nb| ranks[nb.0 as usize] == first_rank);
            // Also check the center itself doesn't break ties via implicit H.
            if all_same && implicit_h == 0 {
                errors.push(StereoError { atom_idx: i, kind: StereoErrorKind::RedundantStereo });
            }
        }
    }

    errors
}

/// Summarise how many stereocenters in `mol` have been specified vs left open.
///
/// A potential stereocenter is an sp3 atom with 4 distinct heavy-atom
/// neighbours (counting one implicit H as a distinct group when present).
pub fn stereo_completeness(mol: &Molecule) -> StereoCompleteness {
    let ranks = simple_morgan_ranks(mol);
    let mut specified = 0usize;
    let mut unspecified = 0usize;

    for (idx, atom) in mol.atoms() {
        // Skip aromatics and obvious non-centers.
        if atom.aromatic { continue; }

        let heavy_nbs: Vec<AtomIdx> = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .map(|(nb, _)| nb)
            .collect();
        let implicit_h = chematic_core::implicit_hcount(mol, idx) as usize;
        let groups = heavy_nbs.len() + implicit_h;

        if groups != 4 { continue; } // only tetrahedral candidates

        // Check all neighbours have distinct ranks (including implicit H as rank 0).
        let mut nb_ranks: Vec<u64> = heavy_nbs.iter().map(|nb| ranks[nb.0 as usize]).collect();
        if implicit_h > 0 { nb_ranks.push(0); }

        let mut sorted = nb_ranks.clone();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() < 4 { continue; } // symmetric neighbours — not a stereocenter

        if atom.chirality != Chirality::None {
            specified += 1;
        } else {
            unspecified += 1;
        }
    }

    StereoCompleteness {
        specified,
        unspecified,
        total_centers: specified + unspecified,
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
        assert!(errors.is_empty(), "L-alanine should have no stereo errors: {errors:?}");
    }

    #[test]
    fn test_impossible_center_explicit_h_zero() {
        // A carbon with chirality annotation, 1 heavy bond, and explicit H=0
        // gives total_groups = 1 → ImpossibleCenter.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder, Chirality};
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
            errors.iter().any(|e| e.atom_idx == 0 && e.kind == StereoErrorKind::ImpossibleCenter),
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
}
