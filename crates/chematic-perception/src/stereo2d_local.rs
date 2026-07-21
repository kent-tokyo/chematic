//! Local (CIP-independent) tetrahedral parity from wedge bonds and 2D coordinates.
//!
//! [`local_parity_from_wedges`] computes `Atom.chirality` plus a matching
//! `stereo_neighbor_order` directly from wedge/hash bonds and a 2D layout,
//! using only the neighbor order already present in `mol` (the order the
//! reader built the adjacency list in). It never consults CIP priority --
//! unlike [`crate::stereo2d::assign_stereo_from_2d`], which ranks substituents
//! by CIP to produce an R/S label. A CIP tie must not prevent a molecule from
//! having *some* known chirality: the wedge drawing itself is a complete,
//! CIP-independent statement of local parity, and CIP ranking is a separate,
//! strictly later, optional stage (mirroring how RDKit's own pipeline keeps
//! geometry-to-parity and parity-to-label as two independent stages).
//!
//! **Convention** (measured against RDKit's `CHI_TETRAHEDRAL_CW`/`CCW` tag on
//! frame-aligned fixtures -- i.e. fixtures where chematic's `mol.neighbors()`
//! order and RDKit's `GetNeighbors()` order are confirmed identical, not
//! derived by analogy to the CIP R/S convention):
//! - 4 explicit neighbors: viewing *from* the first neighbor (in `mol`'s
//!   adjacency order) toward the center, the remaining three go
//!   counterclockwise for [`Chirality::CounterClockwise`] (`@`), clockwise for
//!   [`Chirality::Clockwise`] (`@@`) -- `chematic_core::Chirality`'s own
//!   documented convention.
//! - 3 heavy neighbors + 1 implicit H: no synthetic 3D position is invented
//!   for the H, matching RDKit's own `atomChiralTypeFromBondDirPseudo3D`
//!   nNbrs==3 path, which computes purely from the 3 real bond vectors *from
//!   the center atom* (mathematically the same triple product as the 4-real
//!   case with the center itself, rather than a neighbor, as the pivot). The
//!   resulting `stereo_neighbor_order` places [`STEREO_H_SENTINEL`] *last*;
//!   the sign flips relative to the 4-explicit case because the pivot moved
//!   from "first neighbor" to "center" (an independently-calibrated, not
//!   assumed, relationship -- confirmed by an odd-permutation cross-check
//!   against RDKit's own root-atom SMILES output for the same fixture).
//!
//! This module never calls into `cip_priority` and never touches
//! `Atom.cip_code`. Nothing in the reader crates calls this yet -- integration
//! is a separate, later step.

use chematic_core::{AtomIdx, Chirality, Molecule, STEREO_H_SENTINEL};

use crate::stereo2d::{P3, signed_volume, wedge_z};

/// Tolerance below which a signed volume is treated as coplanar/degenerate.
const VOLUME_EPS: f64 = 1e-6;

/// Compute local tetrahedral parity for `center` from wedge bonds and 2D
/// coordinates, without any CIP ranking.
///
/// Returns `(chirality, stereo_neighbor_order)` on success, using the
/// neighbor order already recorded in `mol` (see module docs for the exact
/// sign convention). Returns `None` when:
/// - `center` doesn't have exactly 4 explicit neighbors, or exactly 3
///   explicit neighbors with `implicit_hcount(mol, center) == 1` (any other
///   neighbor count, or a 3-neighbor atom whose "missing" valence isn't a
///   single implicit H, isn't a tetrahedral stereocenter this function
///   handles);
/// - any required neighbor's or the center's 2D coordinate is missing;
/// - more than one neighbor bond carries a wedge/hash (contradictory
///   notation -- a valid wedge drawing has at most one non-planar bond per
///   stereocenter);
/// - the resulting signed volume is (near-)zero, i.e. the drawing is
///   coplanar/degenerate -- including the case where no bond is wedged at all.
pub fn local_parity_from_wedges(
    mol: &Molecule,
    coords: &[(f64, f64)],
    center: AtomIdx,
) -> Option<(Chirality, Vec<u32>)> {
    let nbs: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();

    match nbs.len() {
        4 => tetrahedral_4(mol, coords, center, &nbs),
        3 if chematic_core::implicit_hcount(mol, center) == 1 => {
            tetrahedral_3_implicit_h(mol, coords, center, &nbs)
        }
        _ => None,
    }
}

fn point_for(coords: &[(f64, f64)], mol: &Molecule, center: AtomIdx, nb: AtomIdx) -> Option<P3> {
    let (x, y) = coords.get(nb.0 as usize).copied()?;
    Some(P3 {
        x,
        y,
        z: wedge_z(mol, center, nb),
    })
}

/// At most one wedge/hash bond may originate from `center`; more than one is
/// contradictory notation.
fn has_at_most_one_wedge(mol: &Molecule, center: AtomIdx, nbs: &[AtomIdx]) -> bool {
    nbs.iter()
        .filter(|&&nb| wedge_z(mol, center, nb) != 0.0)
        .count()
        <= 1
}

fn tetrahedral_4(
    mol: &Molecule,
    coords: &[(f64, f64)],
    center: AtomIdx,
    nbs: &[AtomIdx],
) -> Option<(Chirality, Vec<u32>)> {
    if !has_at_most_one_wedge(mol, center, nbs) {
        return None;
    }
    let pts: Vec<P3> = nbs
        .iter()
        .map(|&nb| point_for(coords, mol, center, nb))
        .collect::<Option<_>>()?;

    // Apex = first-listed neighbor; viewed = the other three, in order.
    let vol = signed_volume(pts[1], pts[2], pts[3], pts[0]);
    if vol.abs() < VOLUME_EPS {
        return None;
    }
    let chirality = if vol < 0.0 {
        Chirality::CounterClockwise
    } else {
        Chirality::Clockwise
    };
    let order = nbs.iter().map(|a| a.0).collect();
    Some((chirality, order))
}

fn tetrahedral_3_implicit_h(
    mol: &Molecule,
    coords: &[(f64, f64)],
    center: AtomIdx,
    nbs: &[AtomIdx],
) -> Option<(Chirality, Vec<u32>)> {
    if !has_at_most_one_wedge(mol, center, nbs) {
        return None;
    }
    let pts: Vec<P3> = nbs
        .iter()
        .map(|&nb| point_for(coords, mol, center, nb))
        .collect::<Option<_>>()?;
    let (cx, cy) = coords.get(center.0 as usize).copied()?;
    let center_pt = P3 {
        x: cx,
        y: cy,
        z: 0.0,
    };

    // No synthetic position for the implicit H: the triple product of the
    // three real bond vectors from `center` already carries full parity.
    let vol = signed_volume(pts[0], pts[1], pts[2], center_pt);
    if vol.abs() < VOLUME_EPS {
        return None;
    }
    let chirality = if vol < 0.0 {
        Chirality::Clockwise
    } else {
        Chirality::CounterClockwise
    };
    let mut order: Vec<u32> = nbs.iter().map(|a| a.0).collect();
    order.push(STEREO_H_SENTINEL);
    Some((chirality, order))
}

/// Apply [`local_parity_from_wedges`] to every eligible atom in `mol`,
/// writing `Atom.chirality` and `stereo_neighbor_order` in-place.
///
/// Does not touch `Atom.cip_code`. Not called by any reader yet -- callers
/// opt in explicitly.
pub fn apply_local_parity_from_wedges(mol: &mut Molecule, coords: &[(f64, f64)]) {
    let atom_indices: Vec<AtomIdx> = mol.atoms().map(|(idx, _)| idx).collect();
    for idx in atom_indices {
        if let Some((chirality, order)) = local_parity_from_wedges(mol, coords, idx) {
            mol.set_chirality(idx, chirality);
            mol.set_stereo_neighbor_order(idx, order);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

    /// Asymmetric, non-degenerate 4-position layout (same shape used by the
    /// P1-A0 diagnosis fixtures) so no accidental coplanarity sneaks in.
    fn quad_positions() -> [(f64, f64); 4] {
        [(-1.0, 0.4), (0.9, 0.7), (-0.5, -1.1), (0.8, -0.6)]
    }

    fn chfclbr(wedge_on_first: bool) -> (Molecule, Vec<(f64, f64)>, AtomIdx) {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        let order = if wedge_on_first {
            BondOrder::Up
        } else {
            BondOrder::Single
        };
        b.add_bond(c, f, order).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        (b.build(), coords, c)
    }

    #[test]
    fn tetrahedral_4heavy_wedge_gives_counterclockwise() {
        // Matches the calibration fixture measured against RDKit directly:
        // wedge on the first-declared bond -> negative volume -> CCW,
        // cross-checked against RDKit's CHI_TETRAHEDRAL_CCW on the same
        // frame-aligned neighbor order.
        let (mol, coords, c) = chfclbr(true);
        let (chirality, order) = local_parity_from_wedges(&mol, &coords, c).unwrap();
        assert_eq!(chirality, Chirality::CounterClockwise);
        assert_eq!(order, vec![1, 2, 3, 4]);
    }

    #[test]
    fn tetrahedral_4heavy_no_h_all_explicit() {
        // C(F)(Cl)(Br)(I): zero H anywhere, still just the 4-explicit path.
        let (mol, coords, c) = chfclbr(true);
        assert!(local_parity_from_wedges(&mol, &coords, c).is_some());
    }

    #[test]
    fn tetrahedral_4neighbors_explicit_h() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let h = b.add_atom(Atom::new(Element::H));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, h, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        let mol = b.build();
        let (chirality, order) = local_parity_from_wedges(&mol, &coords, c).unwrap();
        assert_eq!(chirality, Chirality::CounterClockwise);
        assert_eq!(order, vec![1, 2, 3, 4]);
    }

    #[test]
    fn wedge_hash_inversion_flips_chirality() {
        let (mol_wedge, coords, c) = chfclbr(true);
        let (wedge_chirality, _) = local_parity_from_wedges(&mol_wedge, &coords, c).unwrap();

        let mut b = MoleculeBuilder::new();
        let c2 = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c2, f, BondOrder::Down).unwrap();
        b.add_bond(c2, cl, BondOrder::Single).unwrap();
        b.add_bond(c2, br, BondOrder::Single).unwrap();
        b.add_bond(c2, i, BondOrder::Single).unwrap();
        let mol_hash = b.build();
        let (hash_chirality, _) = local_parity_from_wedges(&mol_hash, &coords, c2).unwrap();

        assert_ne!(wedge_chirality, hash_chirality);
    }

    #[test]
    fn bond_atom_order_inversion_flips_chirality() {
        // Same physical molecule, same wedge (still on C-F, positions
        // unchanged) -- only the BOND BLOCK lists the four bonds in reverse
        // (I, Br, Cl, F instead of F, Cl, Br, I). Reversing a 4-element list
        // is an even permutation, so per the calibration this must give the
        // SAME sign, matching RDKit's own behavior on the equivalent
        // reordering (confirmed against real RDKit output, not assumed).
        let (mol1, coords, c1) = chfclbr(true);
        let (chirality1, _) = local_parity_from_wedges(&mol1, &coords, c1).unwrap();

        let mut b = MoleculeBuilder::new();
        let c2 = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        b.add_bond(c2, i, BondOrder::Single).unwrap();
        b.add_bond(c2, br, BondOrder::Single).unwrap();
        b.add_bond(c2, cl, BondOrder::Single).unwrap();
        b.add_bond(c2, f, BondOrder::Up).unwrap();
        let mol2 = b.build();
        let (chirality2, order2) = local_parity_from_wedges(&mol2, &coords, c2).unwrap();

        assert_eq!(chirality1, chirality2);
        assert_eq!(order2, vec![4, 3, 2, 1]);
    }

    #[test]
    fn multiple_stereocenters_both_assigned() {
        // Two independent CHFClBr-like centers joined by a bond; each should
        // get its own chirality without interference.
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let f1 = b.add_atom(Atom::new(Element::F));
        let cl1 = b.add_atom(Atom::new(Element::CL));
        let br1 = b.add_atom(Atom::new(Element::BR));
        let c2 = b.add_atom(Atom::new(Element::C));
        let f2 = b.add_atom(Atom::new(Element::F));
        let cl2 = b.add_atom(Atom::new(Element::CL));
        let br2 = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c1, f1, BondOrder::Up).unwrap();
        b.add_bond(c1, cl1, BondOrder::Single).unwrap();
        b.add_bond(c1, br1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, f2, BondOrder::Down).unwrap();
        b.add_bond(c2, cl2, BondOrder::Single).unwrap();
        b.add_bond(c2, br2, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![
            (0.0, 0.0),
            quad[0],
            quad[1],
            quad[2],
            (3.0, 0.0),
            (3.0 + quad[0].0, quad[0].1),
            (3.0 + quad[1].0, quad[1].1),
            (3.0 + quad[2].0, quad[2].1),
        ];
        let mut mol = b.build();
        apply_local_parity_from_wedges(&mut mol, &coords);
        assert_eq!(mol.atom(c1).chirality, Chirality::CounterClockwise);
        assert_eq!(mol.atom(c2).chirality, Chirality::Clockwise);
        // cip_code must stay untouched -- this module never assigns it.
        assert_eq!(mol.atom(c1).cip_code, None);
        assert_eq!(mol.atom(c2).cip_code, None);
    }

    #[test]
    fn cip_priority_tie_still_gets_chirality() {
        // Two neighbors (F and Cl) are placed so the local geometry is
        // unambiguous, but with a tied CIP-relevant substituent pair the
        // CIP-based assign_stereo_from_2d would refuse (see rank4's
        // Ordering::Equal -> None branch in stereo2d.rs). This module never
        // calls rank4/cip_priority, so a CIP tie must not matter at all.
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        // Two identical -CH2-CH3 branches: tied under any CIP ranking, but
        // geometrically perfectly resolvable from wedge + coordinates.
        let et1 = b.add_atom(Atom::new(Element::C));
        let et1b = b.add_atom(Atom::new(Element::C));
        let et2 = b.add_atom(Atom::new(Element::C));
        let et2b = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        b.add_bond(c, et1, BondOrder::Single).unwrap();
        b.add_bond(et1, et1b, BondOrder::Single).unwrap();
        b.add_bond(c, et2, BondOrder::Single).unwrap();
        b.add_bond(et2, et2b, BondOrder::Single).unwrap();
        b.add_bond(c, f, BondOrder::Up).unwrap();
        let h = b.add_atom(Atom::new(Element::H));
        b.add_bond(c, h, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![
            (0.0, 0.0),
            quad[0],
            (quad[0].0 + 0.3, quad[0].1 + 1.0),
            quad[1],
            (quad[1].0 + 0.3, quad[1].1 - 1.0),
            quad[2],
            quad[3],
        ];
        let mol = b.build();

        // Sanity check: the CIP-dependent path really does refuse this atom
        // (tied ethyl branches), so the two functions are genuinely being
        // exercised on the same tie condition.
        let cip_result = crate::stereo2d::assign_stereo_from_2d(&mol, &coords);
        assert!(cip_result.get(c).is_none(), "CIP-based path should tie");

        let (chirality, _) = local_parity_from_wedges(&mol, &coords, c).unwrap();
        assert_ne!(chirality, Chirality::None);
    }

    #[test]
    fn missing_coordinates_no_assignment() {
        let (mol, mut coords, c) = chfclbr(true);
        coords.truncate(3); // drop the last neighbor's coordinate entirely
        assert!(local_parity_from_wedges(&mol, &coords, c).is_none());
    }

    #[test]
    fn degenerate_coplanar_no_assignment() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        // No wedge bonds at all -- every z is 0, so the volume is exactly 0.
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        let mol = b.build();
        assert!(local_parity_from_wedges(&mol, &coords, c).is_none());
    }

    #[test]
    fn contradictory_wedges_no_assignment() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let i = b.add_atom(Atom::new(Element::I));
        // Two wedges from the same center: contradictory notation.
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Up).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, i, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2], quad[3]];
        let mol = b.build();
        assert!(local_parity_from_wedges(&mol, &coords, c).is_none());
    }

    #[test]
    fn tetrahedral_3heavy_implicit_h_wedge() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Up).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2]];
        let mol = b.build();
        // Matches the calibration fixture measured against RDKit: 3 heavy
        // neighbors + implicit H, wedge on the first bond -> RDKit's raw
        // CHI_TETRAHEDRAL_CW (and independently confirmed via RDKit's own
        // root-atom SMILES [C@H](F)(Cl)Br for the wedge case, translated
        // through the H-last vs H-first permutation parity).
        let (chirality, order) = local_parity_from_wedges(&mol, &coords, c).unwrap();
        assert_eq!(chirality, Chirality::Clockwise);
        assert_eq!(order, vec![1, 2, 3, STEREO_H_SENTINEL]);
    }

    #[test]
    fn tetrahedral_3heavy_implicit_h_hash_inverts() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Down).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2]];
        let mol = b.build();
        let (chirality, _) = local_parity_from_wedges(&mol, &coords, c).unwrap();
        assert_eq!(chirality, Chirality::CounterClockwise);
    }

    #[test]
    fn tetrahedral_3heavy_bond_order_reversed_inverts() {
        // Same physical molecule and wedge as tetrahedral_3heavy_implicit_h_wedge
        // (positions unchanged) -- only the bond block lists the three bonds
        // in reverse (Br, Cl, F instead of F, Cl, Br). Reversing a 3-element
        // list is an odd permutation (unlike the 4-neighbor case), so per
        // calibration the sign must flip.
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, f, BondOrder::Up).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2]];
        let mol = b.build();
        let (chirality, order) = local_parity_from_wedges(&mol, &coords, c).unwrap();
        assert_eq!(order, vec![3, 2, 1, STEREO_H_SENTINEL]);
        assert_eq!(chirality, Chirality::CounterClockwise);
    }

    #[test]
    fn only_three_heavy_no_implicit_h_no_assignment() {
        // 3 explicit neighbors but valence implies 0 implicit H (e.g. a
        // charged center) -- not the 3-heavy-plus-H shape this function
        // handles, so it must refuse rather than guess.
        let mut b = MoleculeBuilder::new();
        let mut n_atom = Atom::new(Element::N);
        n_atom.charge = 1;
        let n = b.add_atom(n_atom);
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(n, c1, BondOrder::Double).unwrap();
        b.add_bond(n, c2, BondOrder::Single).unwrap();
        b.add_bond(n, c3, BondOrder::Single).unwrap();
        let quad = quad_positions();
        let coords = vec![(0.0, 0.0), quad[0], quad[1], quad[2]];
        let mol = b.build();
        assert!(local_parity_from_wedges(&mol, &coords, n).is_none());
    }
}
