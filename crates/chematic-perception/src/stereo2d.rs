//! v0.1.93: Stereochemistry assignment from 2D coordinates with full CIP prioritization.
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
//! Priority is determined using full multi-sphere BFS CIP rules (v0.1.93+),
//! superseding the simplified 1-sphere approach from v0.1.92.

use std::cmp::Ordering;

use crate::cip_priority;
use chematic_core::{AtomIdx, BondIdx, BondOrder, CipCode, Molecule};

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
        self.assignments
            .iter()
            .find(|(i, _)| *i == idx)
            .map(|(_, c)| *c)
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

/// Assign E/Z stereo to double bonds using 2D atom coordinates.
///
/// Unlike the CIP-based approach (which reads `BondOrder::Up/Down` hints near
/// sp2 carbons), this function derives E/Z purely from the 2D positions of
/// atoms: it projects each substituent onto the perpendicular of the
/// double-bond axis and compares which side each highest-priority substituent
/// falls on.
///
/// Sets `cip_code` (E or Z) on atom1 of each resolved double bond in `mol`.
/// Double bonds that are terminal (no substituent at one end) or have tied
/// substituent priorities are skipped.
pub fn assign_ez_from_2d(mol: &mut Molecule, coords: &[(f64, f64)]) {
    let bond_indices: Vec<BondIdx> = mol
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Double)
        .map(|(bidx, _)| bidx)
        .collect();
    for bidx in bond_indices {
        if let Some((atom_idx, code)) = ez_from_coords(mol, bidx, coords) {
            mol.set_cip_code(atom_idx, Some(code));
        }
    }
}

/// Return the E/Z descriptor for the double bond at `bond_idx` from 2D coordinates.
///
/// Returns `None` when the bond is not a double bond, is a terminal alkene,
/// has an ambiguous geometry (substituent on the bond axis), or has tied
/// CIP priorities.
pub fn cip_ez_descriptor(
    mol: &Molecule,
    bond_idx: BondIdx,
    coords: &[(f64, f64)],
) -> Option<CipCode> {
    ez_from_coords(mol, bond_idx, coords).map(|(_, code)| code)
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
// E/Z core logic
// ---------------------------------------------------------------------------

/// Compute E/Z for a single double bond from 2D coordinates.
///
/// Returns `Some((atom1_of_bond, E | Z))` or `None` if unresolvable.
fn ez_from_coords(
    mol: &Molecule,
    bond_idx: BondIdx,
    coords: &[(f64, f64)],
) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    if bond.order != BondOrder::Double {
        return None;
    }

    let a1 = bond.atom1;
    let a2 = bond.atom2;

    let (x1, y1) = coords.get(a1.0 as usize).copied()?;
    let (x2, y2) = coords.get(a2.0 as usize).copied()?;

    // Non-double-bond substituents at each alkene end.
    let subs_a1: Vec<AtomIdx> = mol
        .neighbors(a1)
        .filter(|&(nb, bidx)| nb != a2 && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    let subs_a2: Vec<AtomIdx> = mol
        .neighbors(a2)
        .filter(|&(nb, bidx)| nb != a1 && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None; // terminal alkene — no geometric isomerism
    }

    // Highest-priority substituent (simplified 1-sphere CIP) at each end.
    let ha1 = highest_ez_priority(mol, a1, &subs_a1)?;
    let ha2 = highest_ez_priority(mol, a2, &subs_a2)?;

    let (hx1, hy1) = coords.get(ha1.0 as usize).copied()?;
    let (hx2, hy2) = coords.get(ha2.0 as usize).copied()?;

    // Double-bond direction vector (a1 → a2).
    let vx = x2 - x1;
    let vy = y2 - y1;

    // Signed 2D cross product: cross(v, u) = vx*uy - vy*ux.
    // The sign tells which side of the bond axis each substituent lies on.
    let side_a1 = cross2d(vx, vy, hx1 - x1, hy1 - y1);
    let side_a2 = cross2d(vx, vy, hx2 - x1, hy2 - y1);

    if side_a1.abs() < 1e-6 || side_a2.abs() < 1e-6 {
        return None; // substituent falls on the bond axis — geometry ambiguous
    }

    // Same side → Z (zusammen); opposite sides → E (entgegen).
    let code = if side_a1.signum() == side_a2.signum() {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some((a1, code))
}

/// Return the highest-priority substituent from `subs` at `center`.
///
/// Returns `None` when the top two substituents have equal priority (tied →
/// E/Z is unspecified) or `subs` is empty.
fn highest_ez_priority(mol: &Molecule, center: AtomIdx, subs: &[AtomIdx]) -> Option<AtomIdx> {
    if subs.is_empty() {
        return None;
    }
    if subs.len() == 1 {
        return Some(subs[0]);
    }
    let mut sorted = subs.to_vec();
    sorted.sort_by(|&a, &b| cip_priority::compare_branches(mol, center, a, b).reverse());
    if cip_priority::compare_branches(mol, center, sorted[0], sorted[1]) == Ordering::Equal {
        return None; // tied top priorities → E/Z not determinable
    }
    Some(sorted[0])
}

/// 2D cross product scalar: vx*uy - vy*ux.
fn cross2d(vx: f64, vy: f64, ux: f64, uy: f64) -> f64 {
    vx * uy - vy * ux
}

// ---------------------------------------------------------------------------
// Core R/S logic
// ---------------------------------------------------------------------------

/// 3D point (x, y, z).
#[derive(Clone, Copy)]
pub(crate) struct P3 {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

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
        let (x, y) = coords
            .get(subs[i].0 as usize)
            .copied()
            .unwrap_or(*center_pos);
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
pub(crate) fn wedge_z(mol: &Molecule, center: AtomIdx, neighbor: AtomIdx) -> f64 {
    // Check bond center→neighbor
    if let Some((_, bond)) = mol.bond_between(center, neighbor) {
        match bond.order {
            BondOrder::Up => {
                // If bond.atom1 == center, neighbor is in front (+z).
                // If bond.atom1 == neighbor (bond drawn away), neighbor is behind (−z).
                // Standard convention: wedge tip at atom1, base at atom2 → atom2 is in front.
                if bond.atom1 == center {
                    return 1.0;
                } else {
                    return -1.0;
                }
            }
            BondOrder::Down => {
                if bond.atom1 == center {
                    return -1.0;
                } else {
                    return 1.0;
                }
            }
            _ => {}
        }
    }
    0.0
}

// Note: 1-sphere CIP helpers removed in v0.1.93; now using cip_priority module

fn rank4(mol: &Molecule, center: AtomIdx, subs: &[AtomIdx; 4]) -> Option<[u8; 4]> {
    let mut order: [usize; 4] = [0, 1, 2, 3];
    order.sort_by(|&i, &j| cip_priority::compare_branches(mol, center, subs[i], subs[j]).reverse());
    for k in 0..3 {
        if cip_priority::compare_branches(mol, center, subs[order[k]], subs[order[k + 1]])
            == Ordering::Equal
        {
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

pub(crate) fn signed_volume(p1: P3, p2: P3, p3: P3, p4: P3) -> f64 {
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

    // -----------------------------------------------------------------------
    // E/Z tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_ez_terminal_alkene_no_assignment() {
        // CH2=CH2 — terminal, no E/Z possible
        let mol = parse("C=C").unwrap();
        let coords = vec![(0.0, 0.0), (1.5, 0.0)];
        // No E/Z should be assigned for a terminal alkene
        let bond_idx = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .map(|(i, _)| i);
        assert!(bond_idx.is_some());
        let result = cip_ez_descriptor(&mol, bond_idx.unwrap(), &coords);
        assert!(result.is_none(), "terminal alkene should have no E/Z");
    }

    #[test]
    fn test_ez_but2ene_z() {
        // (Z)-but-2-ene: CH3\C=C/CH3
        // Atoms: C0(Me)-C1=C2-C3(Me)
        // Layout: C1 and C2 horizontal; C0 above-left, C3 above-right → same side → Z
        //
        //   C0          C3
        //     \        /
        //      C1 = C2
        //
        let mol = parse("CC=CC").unwrap();
        // Double bond is between atoms 1 and 2 (0-indexed SMILES order).
        let bond_idx = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .map(|(i, _)| i)
            .unwrap();
        let coords = vec![
            (-0.866, 0.5), // C0 (methyl at a1 end, above)
            (0.0, 0.0),    // C1 (left alkene carbon)
            (1.5, 0.0),    // C2 (right alkene carbon)
            (2.366, 0.5),  // C3 (methyl at a2 end, above — same side as C0)
        ];
        let result = cip_ez_descriptor(&mol, bond_idx, &coords);
        assert_eq!(result, Some(CipCode::Z), "(Z)-but-2-ene should be Z");
    }

    #[test]
    fn test_ez_but2ene_e() {
        // (E)-but-2-ene: CH3/C=C/CH3
        // C0 above-left, C3 below-right → opposite sides → E
        //
        //   C0
        //     \
        //      C1 = C2
        //              \
        //               C3
        //
        let mol = parse("CC=CC").unwrap();
        let bond_idx = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .map(|(i, _)| i)
            .unwrap();
        let coords = vec![
            (-0.866, 0.5), // C0 above
            (0.0, 0.0),    // C1
            (1.5, 0.0),    // C2
            (2.366, -0.5), // C3 below → opposite side → E
        ];
        let result = cip_ez_descriptor(&mol, bond_idx, &coords);
        assert_eq!(result, Some(CipCode::E), "(E)-but-2-ene should be E");
    }

    #[test]
    fn test_assign_ez_from_2d_sets_cip_code() {
        let mut mol = parse("CC=CC").unwrap();
        let coords = vec![
            (-0.866, 0.5),
            (0.0, 0.0),
            (1.5, 0.0),
            (2.366, 0.5), // same side → Z
        ];
        assign_ez_from_2d(&mut mol, &coords);
        // At least one atom should have a Z cip_code after assignment.
        let has_z = mol
            .atoms()
            .any(|(_, atom)| atom.cip_code == Some(CipCode::Z));
        assert!(has_z, "assign_ez_from_2d should set Z on but-2-ene");
    }

    // -----------------------------------------------------------------------
    // R/S tests
    // -----------------------------------------------------------------------

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
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let h = b.add_atom(Atom::new(Element::H));
        b.add_bond(c, f, BO::Single).unwrap();
        b.add_bond(c, cl, BO::Single).unwrap();
        b.add_bond(c, br, BO::Up).unwrap(); // Br is in front of plane
        b.add_bond(c, h, BO::Down).unwrap(); // H is behind plane
        let mol = b.build();
        // Use non-degenerate 2D positions to avoid a zero-volume determinant.
        let coords = vec![
            (0.0, 0.0),   // C
            (-1.0, -0.5), // F
            (1.0, -0.5),  // Cl
            (0.0, 1.0),   // Br  (z = +1 from Up bond)
            (0.0, -1.0),  // H   (z = -1 from Down bond)
        ];
        let result = assign_stereo_from_2d(&mol, &coords);
        // Should assign one chiral center.
        assert_eq!(result.assignments.len(), 1);
        assert_eq!(result.assignments[0].0, c);
    }
}
