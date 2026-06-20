//! v0.1.93: Stereochemistry assignment from 3D coordinates with full CIP prioritization.
//!
//! [`assign_stereo_from_3d`] determines R/S (tetrahedral) and E/Z (alkene)
//! stereo descriptors from 3D atom positions, without relying on SMILES
//! wedge/dash bond annotations.
//!
//! **Priority rule**: v0.1.93+ uses full multi-sphere BFS CIP priority rules,
//! superseding the simplified 1-sphere approach from prior versions.

use std::cmp::Ordering;

use chematic_core::{AtomIdx, BondIdx, BondOrder, CipCode, Molecule};
use chematic_perception::cip_priority;

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Public result type
// ---------------------------------------------------------------------------

/// Result of a 3D-based stereochemistry assignment.
#[derive(Debug, Default)]
pub struct StereoAssignment3D {
    pub assignments: Vec<(AtomIdx, CipCode)>,
}

impl StereoAssignment3D {
    /// Look up the CIP code for a given atom index.
    pub fn get(&self, idx: AtomIdx) -> Option<CipCode> {
        self.assignments
            .iter()
            .find(|(i, _)| *i == idx)
            .map(|(_, c)| *c)
    }
}

// Note: priority helpers moved to chematic_perception::cip_priority (v0.1.93+)

/// Rank 4 substituents by priority (rank 1 = lowest, rank 4 = highest).
/// Returns `None` if any two substituents are tied (ambiguous).
/// Uses full multi-sphere BFS CIP from chematic_perception.
fn rank4(mol: &Molecule, center: AtomIdx, subs: &[AtomIdx; 4]) -> Option<[u8; 4]> {
    let mut order: [usize; 4] = [0, 1, 2, 3];
    order.sort_by(|&i, &j| cip_priority::compare_branches(mol, center, subs[i], subs[j]).reverse());

    // Check for ties.
    for k in 0..3 {
        if cip_priority::compare_branches(mol, center, subs[order[k]], subs[order[k + 1]])
            == Ordering::Equal
        {
            return None;
        }
    }

    // ranks[i] = rank of subs[i] (1 = lowest, 4 = highest).
    let mut ranks = [0u8; 4];
    for (rank_from_top, &idx) in order.iter().enumerate() {
        ranks[idx] = (4 - rank_from_top) as u8;
    }
    Some(ranks)
}

// ---------------------------------------------------------------------------
// Signed volume
// ---------------------------------------------------------------------------

/// Signed volume of the tetrahedron with vertices p1, p2, p3, p4.
///
/// V = det([p1-p4, p2-p4, p3-p4])
/// V > 0 → p1→p2→p3 appears CCW when viewed from p4 → S
/// V < 0 → p1→p2→p3 appears CW when viewed from p4 → R
fn signed_volume(p1: Point3, p2: Point3, p3: Point3, p4: Point3) -> f64 {
    let ax = p1.x - p4.x;
    let ay = p1.y - p4.y;
    let az = p1.z - p4.z;
    let bx = p2.x - p4.x;
    let by = p2.y - p4.y;
    let bz = p2.z - p4.z;
    let cx = p3.x - p4.x;
    let cy = p3.y - p4.y;
    let cz = p3.z - p4.z;
    // det([a, b, c])
    ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
}

// ---------------------------------------------------------------------------
// R/S assignment
// ---------------------------------------------------------------------------

fn assign_rs(mol: &Molecule, coords: &Coords3D, idx: AtomIdx) -> Option<CipCode> {
    // Collect heavy-atom neighbors.
    let nbs: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();
    if nbs.len() != 4 {
        return None;
    }
    let subs: [AtomIdx; 4] = [nbs[0], nbs[1], nbs[2], nbs[3]];

    let ranks = rank4(mol, idx, &subs)?;

    // Sort subs so subs_sorted[0]=highest priority, subs_sorted[3]=lowest.
    let mut order: [usize; 4] = [0, 1, 2, 3];
    order.sort_by(|&i, &j| ranks[j].cmp(&ranks[i])); // descending rank
    let s: [Point3; 4] = [
        coords.get(subs[order[0]]),
        coords.get(subs[order[1]]),
        coords.get(subs[order[2]]),
        coords.get(subs[order[3]]),
    ];

    // V = signed_volume(highest, 2nd, 3rd; viewed from lowest)
    let v = signed_volume(s[0], s[1], s[2], s[3]);

    if v.abs() < 1e-6 {
        return None; // degenerate / coplanar
    }

    // V > 0 → CCW sequence of 1→2→3 when viewed from 4 → S
    // V < 0 → CW → R
    Some(if v > 0.0 { CipCode::S } else { CipCode::R })
}

// ---------------------------------------------------------------------------
// E/Z assignment
// ---------------------------------------------------------------------------

/// Dihedral angle (radians) for the sequence h1–a1=a2–h2.
/// Returns `None` if vectors are degenerate.
pub(crate) fn dihedral(pa1: Point3, pa2: Point3, ph1: Point3, ph2: Point3) -> Option<f64> {
    // b1 = a1→h1, b2 = a1→a2 (bond axis), b3 = a2→h2
    let b1 = ph1.sub(&pa1);
    let b2 = pa2.sub(&pa1);
    let b3 = ph2.sub(&pa2);

    let n1 = b1.cross(&b2);
    let n2 = b3.cross(&b2);

    let d1 = n1.norm();
    let d2 = n2.norm();
    if d1 < 1e-10 || d2 < 1e-10 {
        return None;
    }

    let cos_a = n1.dot(&n2) / (d1 * d2);
    let angle = cos_a.clamp(-1.0, 1.0).acos();

    // Sign: (n1 × n2) · b2 > 0 → positive (same sense as b2)
    let sign = n1.cross(&n2).dot(&b2);
    Some(if sign < 0.0 { -angle } else { angle })
}

fn assign_ez(mol: &Molecule, coords: &Coords3D, bond_idx: BondIdx) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    if bond.order != BondOrder::Double {
        return None;
    }

    let a1 = bond.atom1;
    let a2 = bond.atom2;

    // Substituents at each end (exclude the other alkene carbon).
    let subs_a1: Vec<AtomIdx> = mol
        .neighbors(a1)
        .filter(|(nb, _)| *nb != a2)
        .map(|(nb, _)| nb)
        .collect();
    let subs_a2: Vec<AtomIdx> = mol
        .neighbors(a2)
        .filter(|(nb, _)| *nb != a1)
        .map(|(nb, _)| nb)
        .collect();

    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None; // terminal alkene
    }

    // Highest-priority substituent at each end (by full multi-sphere CIP).
    let h1 = *subs_a1
        .iter()
        .max_by(|&&a, &&b| cip_priority::compare_branches(mol, a1, a, b))?;
    let h2 = *subs_a2
        .iter()
        .max_by(|&&a, &&b| cip_priority::compare_branches(mol, a2, a, b))?;

    // If either end has two equal-priority substituents, skip.
    if subs_a1.len() == 2
        && cip_priority::compare_branches(mol, a1, subs_a1[0], subs_a1[1]) == Ordering::Equal
    {
        return None;
    }
    if subs_a2.len() == 2
        && cip_priority::compare_branches(mol, a2, subs_a2[0], subs_a2[1]) == Ordering::Equal
    {
        return None;
    }

    let pa1 = coords.get(a1);
    let pa2 = coords.get(a2);
    let ph1 = coords.get(h1);
    let ph2 = coords.get(h2);

    let angle = dihedral(pa1, pa2, ph1, ph2)?;

    // |angle| < π/2 → same side → Z (cis); |angle| ≥ π/2 → E (trans).
    let code = if angle.abs() < std::f64::consts::FRAC_PI_2 {
        CipCode::Z
    } else {
        CipCode::E
    };

    Some((a1, code))
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Assign R/S and E/Z stereochemistry from 3D coordinates.
///
/// Uses full multi-sphere BFS CIP priority rules (v0.1.93+).
/// Centers where priorities cannot be resolved are omitted from the result.
///
/// # Example
/// ```rust,ignore
/// let mol = parse("N[C@@H](C)C(=O)O").unwrap();
/// let coords = assign_stereo_from_3d(&mol, &coords_3d);
/// // returns S for the chiral center if the 3D coords encode L-alanine
/// ```
pub fn assign_stereo_from_3d(mol: &Molecule, coords: &Coords3D) -> StereoAssignment3D {
    let mut assignments = Vec::new();

    // R/S for sp3 centers (degree 4, all bonds single/up/down).
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        // Only carbon and heteroatoms with exactly 4 heavy-atom neighbors.
        if mol.atom(idx).element.atomic_number() == 1 {
            continue;
        }
        let nb_count = mol.neighbors(idx).count();
        if nb_count != 4 {
            continue;
        }
        // Verify all bonds from this atom are single (sp3-like).
        let all_single = mol.neighbors(idx).all(|(_, bidx)| {
            matches!(
                mol.bond(bidx).order,
                BondOrder::Single | BondOrder::Up | BondOrder::Down
            )
        });
        if !all_single {
            continue;
        }
        if let Some(code) = assign_rs(mol, coords, idx) {
            assignments.push((idx, code));
        }
    }

    // E/Z for double bonds.
    for j in 0..mol.bond_count() {
        if let Some((atom_idx, code)) = assign_ez(mol, coords, BondIdx(j as u32)) {
            assignments.push((atom_idx, code));
        }
    }

    StereoAssignment3D { assignments }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    use crate::coords::Coords3D;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    /// Build a known tetrahedral center manually for testing.
    /// Center at origin; substituents placed at given positions.
    #[allow(dead_code)]
    fn make_coords4(
        center: Point3,
        s1: Point3,
        s2: Point3,
        s3: Point3,
        s4: Point3,
        n: usize,
    ) -> Coords3D {
        let mut c = Coords3D::new_zeroed(n);
        c.set(AtomIdx(0), center);
        if n > 1 {
            c.set(AtomIdx(1), s1);
        }
        if n > 2 {
            c.set(AtomIdx(2), s2);
        }
        if n > 3 {
            c.set(AtomIdx(3), s3);
        }
        if n > 4 {
            c.set(AtomIdx(4), s4);
        }
        c
    }

    // -------------------------------------------------------------------------
    // signed_volume
    // -------------------------------------------------------------------------

    #[test]
    fn signed_volume_positive_ccw() {
        // p1=(1,0,0), p2=(0,1,0), p3=(0,0,1), p4=(-1,-1,-1)
        // V = det([p1-p4, p2-p4, p3-p4]) = det([(2,1,1),(1,2,1),(1,1,2)]) = 4 > 0
        let v = signed_volume(
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(-1.0, -1.0, -1.0),
        );
        assert!(v > 0.0, "expected positive signed volume, got {v}");
    }

    #[test]
    fn signed_volume_negative_cw() {
        // Swap p1 and p2 → negative.
        let v = signed_volume(
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(-1.0, -1.0, -1.0),
        );
        assert!(v < 0.0, "expected negative signed volume, got {v}");
    }

    // -------------------------------------------------------------------------
    // R/S assignment with manual Coords3D
    // -------------------------------------------------------------------------
    //
    // Molecule: BrC(F)(Cl)I — 4 different halogens + C center.
    // C=atom0, Br=atom1(an=35), F=atom2(an=9), Cl=atom3(an=17), I=atom4(an=53).
    // CIP priority: I(53) > Br(35) > Cl(17) > F(9).
    //
    // Place them so that I→Br→Cl appears CCW viewed from F → should be S.

    #[test]
    fn rs_manual_s_center() {
        let m = mol("[C]([Br])([F])([Cl])[I]");
        assert_eq!(m.atom_count(), 5, "C + 4 halogens");
        // center=atom0 at origin.
        // I(atom4)  = (0, 0,  1)   → highest priority (53)
        // Br(atom1) = (1, 0, -0.3) → 2nd (35)
        // Cl(atom3) = (-0.5,  0.87, -0.3) → 3rd (17)
        // F(atom2)  = (-0.5, -0.87, -0.3) → lowest (9)
        // Using signed_volume(I, Br, Cl; viewed from F) > 0 → S
        let mut coords = Coords3D::new_zeroed(5);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0)); // C
        coords.set(AtomIdx(1), Point3::new(1.0, 0.0, -0.3)); // Br
        coords.set(AtomIdx(2), Point3::new(-0.5, -0.87, -0.3)); // F
        coords.set(AtomIdx(3), Point3::new(-0.5, 0.87, -0.3)); // Cl
        coords.set(AtomIdx(4), Point3::new(0.0, 0.0, 1.0)); // I

        let result = assign_stereo_from_3d(&m, &coords);
        let code = result.get(AtomIdx(0));
        assert_eq!(code, Some(CipCode::S), "expected S, got {code:?}");
    }

    #[test]
    fn rs_manual_r_center() {
        let m = mol("[C]([Br])([F])([Cl])[I]");
        // Swap I and Br positions → inverts stereo → R
        let mut coords = Coords3D::new_zeroed(5);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0)); // C
        coords.set(AtomIdx(1), Point3::new(0.0, 0.0, 1.0)); // Br (now at top)
        coords.set(AtomIdx(2), Point3::new(-0.5, -0.87, -0.3)); // F
        coords.set(AtomIdx(3), Point3::new(-0.5, 0.87, -0.3)); // Cl
        coords.set(AtomIdx(4), Point3::new(1.0, 0.0, -0.3)); // I (now at side)
        // Now: I(53)=side, Br(35)=top, Cl(17)=side2, F(9)=bottom
        // signed_volume(I, Br, Cl; viewed from F) should be < 0 → R
        let result = assign_stereo_from_3d(&m, &coords);
        let code = result.get(AtomIdx(0));
        assert_eq!(code, Some(CipCode::R), "expected R, got {code:?}");
    }

    // -------------------------------------------------------------------------
    // E/Z assignment
    // -------------------------------------------------------------------------

    #[test]
    fn ez_trans_but2ene() {
        // (E)-but-2-ene: Cl-CH=CH-Br type for easy priority.
        // Use ClCH=CHBr (trans = E).
        let m = mol("ClC=CBr");
        // Place trans: Cl and Br on opposite sides.
        // C0=Cl(an=17), C1=C(alkene,atom1), C2=C(alkene,atom2), C3=Br(an=35)
        // C1=C2 double bond; subs: Cl at C1, Br at C2.
        // Trans arrangement: Cl-C1=C2-Br dihedral ≈ 180°
        let mut coords = Coords3D::new_zeroed(4);
        coords.set(AtomIdx(0), Point3::new(-1.5, 0.0, 0.0)); // Cl
        coords.set(AtomIdx(1), Point3::new(-0.67, 0.0, 0.0)); // C1 (alkene)
        coords.set(AtomIdx(2), Point3::new(0.67, 0.0, 0.0)); // C2 (alkene)
        coords.set(AtomIdx(3), Point3::new(2.0, 0.5, 0.0)); // Br (same side as Cl → trans? no...)

        // For E (trans), Cl and Br should be on opposite sides of the C=C axis.
        // Let's define: Cl at (-1.5, 0.5, 0), C1 at (-0.67, 0, 0), C2 at (0.67, 0, 0), Br at (1.5, -0.5, 0)
        // Dihedral Cl-C1=C2-Br:
        //   b1 = C1→Cl = (-0.83, 0.5, 0)
        //   b2 = C1→C2 = (1.34, 0, 0)
        //   b3 = C2→Br = (0.83, -0.5, 0)
        //   n1 = b1 × b2 = (0, 0, 0.5*0 - 0*0) ... let's just trust the code.
        coords.set(AtomIdx(0), Point3::new(-1.5, 0.5, 0.0)); // Cl (above)
        coords.set(AtomIdx(1), Point3::new(-0.67, 0.0, 0.0)); // C1
        coords.set(AtomIdx(2), Point3::new(0.67, 0.0, 0.0)); // C2
        coords.set(AtomIdx(3), Point3::new(1.5, -0.5, 0.0)); // Br (below) → opposite = E

        let result = assign_stereo_from_3d(&m, &coords);
        // The double bond C1=C2 should be labeled E.
        let found_e = result
            .assignments
            .iter()
            .any(|(_, code)| *code == CipCode::E);
        assert!(
            found_e,
            "expected E assignment for trans ClCH=CHBr, got {:?}",
            result.assignments
        );
    }

    #[test]
    fn ez_cis_arrangement() {
        // Z: Cl and Br on the same side.
        let m = mol("ClC=CBr");
        let mut coords = Coords3D::new_zeroed(4);
        coords.set(AtomIdx(0), Point3::new(-1.5, 0.5, 0.0)); // Cl (above)
        coords.set(AtomIdx(1), Point3::new(-0.67, 0.0, 0.0)); // C1
        coords.set(AtomIdx(2), Point3::new(0.67, 0.0, 0.0)); // C2
        coords.set(AtomIdx(3), Point3::new(1.5, 0.5, 0.0)); // Br (also above) → same side = Z

        let result = assign_stereo_from_3d(&m, &coords);
        let found_z = result
            .assignments
            .iter()
            .any(|(_, code)| *code == CipCode::Z);
        assert!(
            found_z,
            "expected Z assignment for cis ClCH=CHBr, got {:?}",
            result.assignments
        );
    }

    // -------------------------------------------------------------------------
    // Empty / degenerate inputs
    // -------------------------------------------------------------------------

    #[test]
    fn no_stereo_for_benzene() {
        let m = mol("c1ccccc1");
        let mut c = Coords3D::new_zeroed(m.atom_count());
        for i in 0..m.atom_count() {
            c.set(AtomIdx(i as u32), Point3::zero());
        }
        let result = assign_stereo_from_3d(&m, &c);
        assert!(
            result.assignments.is_empty(),
            "benzene should have no stereo assignments"
        );
    }
}
