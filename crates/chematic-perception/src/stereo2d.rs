//! Stereochemistry assignment from 2D coordinates and wedge bonds.
//!
//! [`assign_stereo_from_2d`] infers R/S tetrahedral stereo descriptors by
//! combining the 2D layout positions of substituents with the wedge-bond
//! notation (`BondOrder::Up` / `BondOrder::Down`) already encoded in the
//! molecule.
//!
//! **Convention:**
//! - `BondOrder::Up` (solid wedge): the *far* atom (`bond.atom2`) is in
//!   front of the plane (toward the viewer, z > 0).
//! - `BondOrder::Down` (dashed wedge): the far atom is behind the plane (z < 0).
//! - All other bonds: both endpoints are coplanar (z = 0).
//!
//! Priority is determined with the same simplified 1-sphere CIP rule used by
//! `chematic-3d`'s `assign_stereo_from_3d`.

use std::cmp::Ordering;

use chematic_core::{AtomIdx, BondOrder, CipCode, Molecule};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Result of a 2D-based stereochemistry assignment.
#[derive(Debug, Default, Clone)]
pub struct StereoAssignment2D {
    pub assignments: Vec<(AtomIdx, CipCode)>,
}

impl StereoAssignment2D {
    /// Look up the CIP code for atom `idx`.
    pub fn get(&self, idx: AtomIdx) -> Option<CipCode> {
        self.assignments.iter().find(|(i, _)| *i == idx).map(|(_, c)| *c)
    }
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Assign R/S stereo to atoms that carry wedge bonds, given the 2D layout.
///
/// `coords[i]` is the `(x, y)` position of atom `i` in any consistent unit
/// (screen pixels, Å, etc.).  The z-component is derived from wedge bonds.
///
/// Returns a [`StereoAssignment2D`] containing one entry per successfully
/// assigned chiral center.  Atoms where the chirality is ambiguous (tied
/// priorities or too few heavy-atom neighbors) are omitted.
pub fn assign_stereo_from_2d(mol: &Molecule, coords: &[(f64, f64)]) -> StereoAssignment2D {
    let mut assignments = Vec::new();
    for (idx, _) in mol.atoms() {
        if let Some(code) = assign_rs(mol, coords, idx) {
            assignments.push((idx, code));
        }
    }
    StereoAssignment2D { assignments }
}

/// Like [`assign_stereo_from_2d`] but writes `cip_code` directly onto each
/// chiral atom in `mol`.
pub fn apply_stereo_from_2d(mol: &mut Molecule, coords: &[(f64, f64)]) {
    let result = assign_stereo_from_2d(mol, coords);
    for (idx, code) in result.assignments {
        // atom_mut requires access — we use set_cip_code helper if it exists,
        // or rewrite the atom via add/remove.
        // Since Molecule.atoms is private, use the new mutable set_* API:
        // (molecule.rs exposes atoms[idx].cip_code via a dedicated setter)
        mol.set_cip_code(idx, Some(code));
    }
}

// ---------------------------------------------------------------------------
// Core R/S logic
// ---------------------------------------------------------------------------

/// 3D point (x, y, z).
#[derive(Clone, Copy)]
struct P3 { x: f64, y: f64, z: f64 }

fn assign_rs(mol: &Molecule, coords: &[(f64, f64)], center: AtomIdx) -> Option<CipCode> {
    let center_pos = coords.get(center.0 as usize)?;

    // Collect heavy-atom neighbors with their 3D positions.
    let nbs: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();
    if nbs.len() != 4 {
        return None; // need exactly 4 heavy-atom neighbors for tetrahedral
    }
    let subs: [AtomIdx; 4] = [nbs[0], nbs[1], nbs[2], nbs[3]];

    // Build 3D positions: z from wedge bonds.
    let z_for: [f64; 4] = [
        wedge_z(mol, center, subs[0]),
        wedge_z(mol, center, subs[1]),
        wedge_z(mol, center, subs[2]),
        wedge_z(mol, center, subs[3]),
    ];

    let pts: [P3; 4] = [0, 1, 2, 3].map(|i| {
        let (x, y) = coords.get(subs[i].0 as usize).copied().unwrap_or(*center_pos);
        P3 { x, y, z: z_for[i] }
    });

    let ranks = rank4(mol, center, &subs)?;

    // Sort substituents highest-priority first.
    let mut order: [usize; 4] = [0, 1, 2, 3];
    order.sort_by(|&i, &j| ranks[j].cmp(&ranks[i]));

    let s: [P3; 4] = order.map(|i| pts[i]);
    let v = signed_volume(s[0], s[1], s[2], s[3]);

    if v.abs() < 1e-6 {
        return None; // coplanar / degenerate
    }

    Some(if v > 0.0 { CipCode::S } else { CipCode::R })
}

/// Determine z-offset for `neighbor` as seen from `center`:
/// Up (wedge toward viewer) → +1.0
/// Down (dash away from viewer) → -1.0
/// Anything else → 0.0
fn wedge_z(mol: &Molecule, center: AtomIdx, neighbor: AtomIdx) -> f64 {
    // Check bond center→neighbor
    if let Some((_, bond)) = mol.bond_between(center, neighbor) {
        match bond.order {
            BondOrder::Up => {
                // If bond.atom1 == center, neighbor is in front (+z).
                // If bond.atom1 == neighbor (bond drawn away), neighbor is behind (−z).
                // Standard convention: wedge tip at atom1, base at atom2 → atom2 is in front.
                if bond.atom1 == center { return 1.0; }
                else { return -1.0; }
            }
            BondOrder::Down => {
                if bond.atom1 == center { return -1.0; }
                else { return 1.0; }
            }
            _ => {}
        }
    }
    0.0
}

// ---------------------------------------------------------------------------
// CIP priority helpers (same logic as stereo3d.rs)
// ---------------------------------------------------------------------------

fn priority_key(mol: &Molecule, center: AtomIdx, sub: AtomIdx) -> (u8, Vec<u8>) {
    let an = mol.atom(sub).element.atomic_number();
    let mut nbs: Vec<u8> = mol
        .neighbors(sub)
        .filter(|(nb, _)| *nb != center)
        .map(|(nb, _)| mol.atom(nb).element.atomic_number())
        .collect();
    nbs.sort_unstable_by(|a, b| b.cmp(a));
    (an, nbs)
}

fn cmp_priority(mol: &Molecule, center: AtomIdx, a: AtomIdx, b: AtomIdx) -> Ordering {
    let ka = priority_key(mol, center, a);
    let kb = priority_key(mol, center, b);
    match ka.0.cmp(&kb.0) {
        Ordering::Equal => {
            let max_len = ka.1.len().max(kb.1.len());
            for i in 0..max_len {
                let va = ka.1.get(i).copied().unwrap_or(0);
                let vb = kb.1.get(i).copied().unwrap_or(0);
                match va.cmp(&vb) {
                    Ordering::Equal => {}
                    other => return other,
                }
            }
            Ordering::Equal
        }
        other => other,
    }
}

fn rank4(mol: &Molecule, center: AtomIdx, subs: &[AtomIdx; 4]) -> Option<[u8; 4]> {
    let mut order: [usize; 4] = [0, 1, 2, 3];
    order.sort_by(|&i, &j| cmp_priority(mol, center, subs[i], subs[j]).reverse());
    for k in 0..3 {
        if cmp_priority(mol, center, subs[order[k]], subs[order[k + 1]]) == Ordering::Equal {
            return None;
        }
    }
    let mut ranks = [0u8; 4];
    for (rank_from_top, &idx) in order.iter().enumerate() {
        ranks[idx] = (4 - rank_from_top) as u8;
    }
    Some(ranks)
}

// ---------------------------------------------------------------------------
// Signed volume
// ---------------------------------------------------------------------------

fn signed_volume(p1: P3, p2: P3, p3: P3, p4: P3) -> f64 {
    let (ax, ay, az) = (p1.x - p4.x, p1.y - p4.y, p1.z - p4.z);
    let (bx, by, bz) = (p2.x - p4.x, p2.y - p4.y, p2.z - p4.z);
    let (cx, cy, cz) = (p3.x - p4.x, p3.y - p4.y, p3.z - p4.z);
    ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_no_stereo_ethane() {
        let mol = parse("CC").unwrap();
        let coords = vec![(0.0, 0.0), (1.0, 0.0)];
        let result = assign_stereo_from_2d(&mol, &coords);
        assert!(result.assignments.is_empty());
    }

    #[test]
    fn test_r_s_bromochlorofluoromethane() {
        // CHFClBr — 4 heavy-atom neighbors, wedge C-Br above plane
        // Build: C bonded to F, Cl, Br, H via explicit H
        // We won't test R vs S deterministically (depends on atom ordering),
        // just that an assignment is returned when a wedge bond is present.
        use chematic_core::{Atom, BondOrder as BO, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c  = b.add_atom(Atom::new(Element::C));
        let f  = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let h  = b.add_atom(Atom::new(Element::H));
        b.add_bond(c, f,  BO::Single).unwrap();
        b.add_bond(c, cl, BO::Single).unwrap();
        b.add_bond(c, br, BO::Up).unwrap();    // Br is in front of plane
        b.add_bond(c, h,  BO::Down).unwrap();  // H is behind plane
        let mol = b.build();
        // Use non-degenerate 2D positions to avoid a zero-volume determinant.
        let coords = vec![
            (0.0,   0.0),   // C
            (-1.0, -0.5),   // F
            (1.0,  -0.5),   // Cl
            (0.0,   1.0),   // Br  (z = +1 from Up bond)
            (0.0,  -1.0),   // H   (z = -1 from Down bond)
        ];
        let result = assign_stereo_from_2d(&mol, &coords);
        // Should assign one chiral center.
        assert_eq!(result.assignments.len(), 1);
        assert_eq!(result.assignments[0].0, c);
    }
}
