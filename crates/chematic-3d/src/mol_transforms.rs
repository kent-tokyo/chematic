use crate::coords::{Coords3D, Point3};
use chematic_core::{AtomIdx, Molecule};
use std::collections::HashSet;
use std::fmt;

/// Failure returned by the checked conformer-editing operations.
///
/// The variants mirror the preconditions of RDKit's `rdMolTransforms`
/// setters while keeping invalid indices and degenerate coordinates out of
/// panic paths in Rust, Python, and WebAssembly callers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MolTransformError {
    CoordinateCountMismatch { expected: usize, actual: usize },
    AtomIndexOutOfRange { atom: u32, atom_count: usize },
    NonFiniteCoordinates,
    NonFiniteTarget,
    AtomsNotBonded { atom1: u32, atom2: u32 },
    RingBond { atom1: u32, atom2: u32 },
    BothAngleBondsInRing,
    CoincidentAtoms { atom1: u32, atom2: u32 },
    CollinearAngle,
    UndefinedDihedral,
}

impl fmt::Display for MolTransformError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CoordinateCountMismatch { expected, actual } => {
                write!(f, "expected {expected} coordinates, got {actual}")
            }
            Self::AtomIndexOutOfRange { atom, atom_count } => {
                write!(
                    f,
                    "atom index {atom} is out of range for {atom_count} atoms"
                )
            }
            Self::NonFiniteCoordinates => write!(f, "coordinates must all be finite"),
            Self::NonFiniteTarget => write!(f, "target value must be finite"),
            Self::AtomsNotBonded { atom1, atom2 } => {
                write!(f, "atoms {atom1} and {atom2} must be bonded")
            }
            Self::RingBond { atom1, atom2 } => {
                write!(f, "bond ({atom1},{atom2}) must not belong to a ring")
            }
            Self::BothAngleBondsInRing => {
                write!(f, "both bonds defining the angle must not belong to a ring")
            }
            Self::CoincidentAtoms { atom1, atom2 } => {
                write!(f, "atoms {atom1} and {atom2} have identical 3D coordinates")
            }
            Self::CollinearAngle => {
                write!(f, "angle rotation axis is undefined for collinear atoms")
            }
            Self::UndefinedDihedral => write!(f, "dihedral is undefined for collinear atoms"),
        }
    }
}

impl std::error::Error for MolTransformError {}

/// Get bond length in Ångströms between atoms a and b.
pub fn get_bond_length(coords: &Coords3D, a: AtomIdx, b: AtomIdx) -> f64 {
    coords.get(a).distance(&coords.get(b))
}

/// Checked bond-length measurement for external bindings.
pub fn try_get_bond_length(
    coords: &Coords3D,
    mol: &Molecule,
    a: AtomIdx,
    b: AtomIdx,
) -> Result<f64, MolTransformError> {
    validate_coords_and_atoms(coords, mol, &[a, b])?;
    Ok(get_bond_length(coords, a, b))
}

/// Get bond angle A—Center—B in radians.
pub fn get_bond_angle(coords: &Coords3D, a: AtomIdx, center: AtomIdx, b: AtomIdx) -> f64 {
    super::constraints::compute_angle(coords, a, center, b)
}

/// Get bond angle A—Center—B in degrees.
pub fn get_bond_angle_deg(coords: &Coords3D, a: AtomIdx, center: AtomIdx, b: AtomIdx) -> f64 {
    get_bond_angle(coords, a, center, b).to_degrees()
}

/// Checked angle measurement in degrees for external bindings.
pub fn try_get_bond_angle_deg(
    coords: &Coords3D,
    mol: &Molecule,
    a: AtomIdx,
    center: AtomIdx,
    b: AtomIdx,
) -> Result<f64, MolTransformError> {
    validate_coords_and_atoms(coords, mol, &[a, center, b])?;
    if coords.get(a).distance(&coords.get(center)) <= 1e-8 {
        return Err(MolTransformError::CoincidentAtoms {
            atom1: a.0,
            atom2: center.0,
        });
    }
    if coords.get(center).distance(&coords.get(b)) <= 1e-8 {
        return Err(MolTransformError::CoincidentAtoms {
            atom1: center.0,
            atom2: b.0,
        });
    }
    Ok(get_bond_angle_deg(coords, a, center, b))
}

/// Set the distance between bonded atoms `fixed` and `moved`.
///
/// Atom `fixed` stays in place. Atom `moved` and every atom reachable on its
/// side of the bond are translated together, matching RDKit's
/// `SetBondLength` contract. Ring bonds are rejected because they do not have
/// a unique movable side.
pub fn set_bond_length(
    coords: &Coords3D,
    mol: &Molecule,
    fixed: AtomIdx,
    moved: AtomIdx,
    target_length: f64,
) -> Result<Coords3D, MolTransformError> {
    validate_transform_input(coords, mol, &[fixed, moved], target_length)?;
    let (bond_idx, _) = require_bond(mol, fixed, moved)?;
    if chematic_perception::ring_bond_flags_shared(mol)[bond_idx.0 as usize] {
        return Err(MolTransformError::RingBond {
            atom1: fixed.0,
            atom2: moved.0,
        });
    }

    let fixed_pos = coords.get(fixed);
    let moved_pos = coords.get(moved);
    let fixed_to_moved = moved_pos.sub(&fixed_pos);
    let current_length = fixed_to_moved.norm();
    if current_length <= 1e-8 {
        return Err(MolTransformError::CoincidentAtoms {
            atom1: fixed.0,
            atom2: moved.0,
        });
    }

    let translation = fixed_to_moved.scale(target_length / current_length - 1.0);
    let mut result = coords.clone();
    for atom in find_rotated_atoms(mol, fixed, moved) {
        result.set(atom, coords.get(atom).add(&translation));
    }
    Ok(result)
}

/// Set angle `a-center-moved` in radians.
///
/// `a` and `center` remain fixed. `moved` and every atom on its side of the
/// `center-moved` bond are rotated together, matching RDKit's `SetAngleRad`
/// contract. If both angle bonds are ring bonds there is no independent side
/// to rotate and the operation is rejected.
pub fn set_bond_angle(
    coords: &Coords3D,
    mol: &Molecule,
    a: AtomIdx,
    center: AtomIdx,
    moved: AtomIdx,
    target_angle: f64,
) -> Result<Coords3D, MolTransformError> {
    validate_transform_input(coords, mol, &[a, center, moved], target_angle)?;
    let (first_bond, _) = require_bond(mol, a, center)?;
    let (second_bond, _) = require_bond(mol, center, moved)?;
    let ring_bonds = chematic_perception::ring_bond_flags_shared(mol);
    if ring_bonds[first_bond.0 as usize] && ring_bonds[second_bond.0 as usize] {
        return Err(MolTransformError::BothAngleBondsInRing);
    }

    let center_pos = coords.get(center);
    let center_to_a = coords.get(a).sub(&center_pos);
    let center_to_moved = coords.get(moved).sub(&center_pos);
    if center_to_a.norm() <= 1e-8 {
        return Err(MolTransformError::CoincidentAtoms {
            atom1: a.0,
            atom2: center.0,
        });
    }
    if center_to_moved.norm() <= 1e-8 {
        return Err(MolTransformError::CoincidentAtoms {
            atom1: center.0,
            atom2: moved.0,
        });
    }

    let current_angle = get_bond_angle(coords, a, center, moved);
    let delta = target_angle - current_angle;
    if delta.abs() <= f64::EPSILON {
        return Ok(coords.clone());
    }
    let axis = center_to_a
        .cross(&center_to_moved)
        .try_normalize()
        .ok_or(MolTransformError::CollinearAngle)?;

    let mut result = coords.clone();
    for atom in find_rotated_atoms(mol, center, moved) {
        let relative = coords.get(atom).sub(&center_pos);
        let rotated = super::constraints::rotate_around_axis(relative, axis, delta);
        result.set(atom, center_pos.add(&rotated));
    }
    Ok(result)
}

/// Get dihedral angle A—B—C—D in radians.
/// Returns None if any three atoms are collinear.
pub fn get_dihedral(
    coords: &Coords3D,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    d: AtomIdx,
) -> Option<f64> {
    let pa = coords.get(a);
    let pb = coords.get(b);
    let pc = coords.get(c);
    let pd = coords.get(d);
    super::stereo3d::dihedral(pb, pc, pa, pd)
}

/// Get dihedral angle A—B—C—D in degrees.
pub fn get_dihedral_deg(
    coords: &Coords3D,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    d: AtomIdx,
) -> Option<f64> {
    get_dihedral(coords, a, b, c, d).map(|rad| rad.to_degrees())
}

/// Checked dihedral measurement in degrees for external bindings.
pub fn try_get_dihedral_deg(
    coords: &Coords3D,
    mol: &Molecule,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    d: AtomIdx,
) -> Result<f64, MolTransformError> {
    validate_coords_and_atoms(coords, mol, &[a, b, c, d])?;
    for (first, second) in [(a, b), (b, c), (c, d)] {
        if coords.get(first).distance(&coords.get(second)) <= 1e-8 {
            return Err(MolTransformError::CoincidentAtoms {
                atom1: first.0,
                atom2: second.0,
            });
        }
    }
    get_dihedral_deg(coords, a, b, c, d).ok_or(MolTransformError::UndefinedDihedral)
}

/// Compute centroid (center of mass) of all atoms.
pub fn compute_centroid(coords: &Coords3D) -> [f64; 3] {
    let all_coords: Vec<[f64; 3]> = (0..coords.atom_count())
        .map(|i| {
            let p = coords.get(AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();
    super::usr::centroid(&all_coords)
}

/// Set dihedral angle A—B—C—D to target_angle (radians) by rotating the B-C bond.
/// Returns new Coords3D with atoms on the C side of the bond rotated.
pub fn set_dihedral(
    coords: &Coords3D,
    mol: &Molecule,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    d: AtomIdx,
    target_angle: f64,
) -> Coords3D {
    try_set_dihedral(coords, mol, a, b, c, d, target_angle).unwrap_or_else(|_| coords.clone())
}

/// Checked counterpart to [`set_dihedral`], matching RDKit's requirement
/// that the central bond exists and is not part of a ring.
pub fn try_set_dihedral(
    coords: &Coords3D,
    mol: &Molecule,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    d: AtomIdx,
    target_angle: f64,
) -> Result<Coords3D, MolTransformError> {
    validate_transform_input(coords, mol, &[a, b, c, d], target_angle)?;
    let (central_bond, _) = require_bond(mol, b, c)?;
    if chematic_perception::ring_bond_flags_shared(mol)[central_bond.0 as usize] {
        return Err(MolTransformError::RingBond {
            atom1: b.0,
            atom2: c.0,
        });
    }
    let current = get_dihedral(coords, a, b, c, d).ok_or(MolTransformError::UndefinedDihedral)?;

    let delta = target_angle - current;
    let pb = coords.get(b);
    let pc = coords.get(c);

    // Guard against degenerate geometry: B and C define the rotation axis.
    if pc.distance(&pb) < 1e-10 {
        return Err(MolTransformError::CoincidentAtoms {
            atom1: b.0,
            atom2: c.0,
        });
    }

    let axis = pc.sub(&pb).normalize();

    let mut new_coords = coords.clone();

    // Find atoms on the C side of the B-C bond.
    let rotated = find_rotated_atoms(mol, b, c);

    // Rotate C-side atoms around the B-C axis. Atom C lies on the axis and is unchanged.
    for idx in rotated {
        let p = coords.get(idx);
        let relative = p.sub(&pc);
        let rotated_relative = super::constraints::rotate_around_axis(relative, axis, delta);
        let new_p = pc.add(&rotated_relative);
        new_coords.set(idx, new_p);
    }

    Ok(new_coords)
}

fn validate_transform_input(
    coords: &Coords3D,
    mol: &Molecule,
    atoms: &[AtomIdx],
    target: f64,
) -> Result<(), MolTransformError> {
    validate_coords_and_atoms(coords, mol, atoms)?;
    if !target.is_finite() {
        return Err(MolTransformError::NonFiniteTarget);
    }
    Ok(())
}

fn validate_coords_and_atoms(
    coords: &Coords3D,
    mol: &Molecule,
    atoms: &[AtomIdx],
) -> Result<(), MolTransformError> {
    if coords.atom_count() != mol.atom_count() {
        return Err(MolTransformError::CoordinateCountMismatch {
            expected: mol.atom_count(),
            actual: coords.atom_count(),
        });
    }
    for &atom in atoms {
        if atom.0 as usize >= mol.atom_count() {
            return Err(MolTransformError::AtomIndexOutOfRange {
                atom: atom.0,
                atom_count: mol.atom_count(),
            });
        }
    }
    if !coords.is_finite() {
        return Err(MolTransformError::NonFiniteCoordinates);
    }
    Ok(())
}

fn require_bond(
    mol: &Molecule,
    a: AtomIdx,
    b: AtomIdx,
) -> Result<(chematic_core::BondIdx, &chematic_core::BondEntry), MolTransformError> {
    mol.bond_between(a, b)
        .ok_or(MolTransformError::AtomsNotBonded {
            atom1: a.0,
            atom2: b.0,
        })
}

/// Center conformer on origin by translating centroid to (0, 0, 0).
pub fn center_on_origin(coords: &Coords3D) -> Coords3D {
    let [cx, cy, cz] = compute_centroid(coords);
    let mut result = coords.clone();
    for i in 0..coords.atom_count() {
        let idx = AtomIdx(i as u32);
        let p = coords.get(idx);
        result.set(
            idx,
            Point3 {
                x: p.x - cx,
                y: p.y - cy,
                z: p.z - cz,
            },
        );
    }
    result
}

/// Apply 4×4 homogeneous transformation matrix to all atom positions.
pub fn transform_conformer(coords: &Coords3D, matrix: &[[f64; 4]; 4]) -> Coords3D {
    let mut result = coords.clone();
    for i in 0..coords.atom_count() {
        let idx = AtomIdx(i as u32);
        let p = coords.get(idx);

        let x = matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2] * p.z + matrix[0][3];
        let y = matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2] * p.z + matrix[1][3];
        let z = matrix[2][0] * p.x + matrix[2][1] * p.y + matrix[2][2] * p.z + matrix[2][3];

        result.set(idx, Point3 { x, y, z });
    }
    result
}

/// Find all atoms in the subtree rooted at `target` when coming from `parent`.
fn find_rotated_atoms(mol: &Molecule, parent: AtomIdx, target: AtomIdx) -> Vec<AtomIdx> {
    let mut result = vec![target];
    let mut visited = HashSet::from([target]);
    let mut stack = vec![target];

    while let Some(current) = stack.pop() {
        for (neighbor, _) in mol.neighbors(current) {
            if neighbor != parent && visited.insert(neighbor) {
                result.push(neighbor);
                stack.push(neighbor);
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_get_bond_length() {
        let mol = parse("CC").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        let len = get_bond_length(&coords, AtomIdx(0), AtomIdx(1));
        // C-C single bond ≈ 1.54 Å
        assert!((len - 1.54).abs() < 0.1);
    }

    #[test]
    fn test_get_bond_angle() {
        let mol = parse("CCC").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        let angle = get_bond_angle(&coords, AtomIdx(0), AtomIdx(1), AtomIdx(2));
        // C-C-C sp3 tetrahedral ≈ 109.5° ≈ 1.91 rad
        assert!((angle - std::f64::consts::PI / 2.0).abs() < 0.5);
    }

    #[test]
    fn test_compute_centroid() {
        let mol = parse("C").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        let [cx, cy, cz] = compute_centroid(&coords);
        // Single carbon — centroid should be close to its position
        assert!(cx.is_finite() && cy.is_finite() && cz.is_finite());
    }

    #[test]
    fn test_center_on_origin() {
        let mol = parse("CC").unwrap();
        let mut coords = crate::dg::generate_coords(&mol);
        // Translate first
        let p0 = coords.get(AtomIdx(0));
        coords.set(
            AtomIdx(0),
            Point3 {
                x: p0.x + 10.0,
                y: p0.y + 20.0,
                z: p0.z + 30.0,
            },
        );

        let centered = center_on_origin(&coords);
        let [cx, cy, cz] = compute_centroid(&centered);
        assert!(cx.abs() < 1e-10 && cy.abs() < 1e-10 && cz.abs() < 1e-10);
    }

    #[test]
    fn test_transform_matrix() {
        let mol = parse("C").unwrap();
        let coords = crate::dg::generate_coords(&mol);

        // Identity matrix
        let identity = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let result = transform_conformer(&coords, &identity);
        let p_orig = coords.get(AtomIdx(0));
        let p_result = result.get(AtomIdx(0));
        assert!((p_orig.x - p_result.x).abs() < 1e-10);

        // Translation by (5, 10, 15)
        let translation = [
            [1.0, 0.0, 0.0, 5.0],
            [0.0, 1.0, 0.0, 10.0],
            [0.0, 0.0, 1.0, 15.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let result = transform_conformer(&coords, &translation);
        let p_result = result.get(AtomIdx(0));
        assert!((p_result.x - p_orig.x - 5.0).abs() < 1e-10);
        assert!((p_result.y - p_orig.y - 10.0).abs() < 1e-10);
        assert!((p_result.z - p_orig.z - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_rotated_atoms_linear() {
        // Linear chain: C-C-C-C
        let mol = parse("CCCC").unwrap();
        let rotated = find_rotated_atoms(&mol, AtomIdx(1), AtomIdx(2));
        // Rotating around bond 1-2, atom 1 stays, atom 2+ rotate
        assert!(rotated.contains(&AtomIdx(2)));
        assert!(rotated.contains(&AtomIdx(3)));
        assert!(!rotated.contains(&AtomIdx(0)));
        assert!(!rotated.contains(&AtomIdx(1)));
    }

    #[test]
    fn test_find_rotated_atoms_branched() {
        // Branched: C(C)C with bond 0-1
        let mol = parse("CC(C)C").unwrap();
        let rotated = find_rotated_atoms(&mol, AtomIdx(0), AtomIdx(1));
        // Rotating around bond 0-1, atom 0 stays, atoms 1+ rotate
        assert!(rotated.contains(&AtomIdx(1)));
        assert!(rotated.contains(&AtomIdx(2)));
        assert!(rotated.contains(&AtomIdx(3)));
        assert!(!rotated.contains(&AtomIdx(0)));
    }

    #[test]
    fn test_set_dihedral_linear() {
        let mol = parse("CCCC").unwrap();
        let mut coords = Coords3D::new_zeroed(4);
        coords.set(AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(1.0, 0.0, 0.0));
        coords.set(AtomIdx(2), Point3::new(1.0, 1.0, 0.0));
        coords.set(AtomIdx(3), Point3::new(2.0, 1.0, 1.0));
        let original_p3 = coords.get(AtomIdx(3));

        let target = 60.0_f64.to_radians();
        let result = set_dihedral(
            &coords,
            &mol,
            AtomIdx(0),
            AtomIdx(1),
            AtomIdx(2),
            AtomIdx(3),
            target,
        );
        let new_p3 = result.get(AtomIdx(3));

        let dist_change = original_p3.distance(&new_p3);
        assert!(
            dist_change > 0.1,
            "Atom position should change for different dihedral"
        );

        let dihedral =
            get_dihedral(&result, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)).unwrap();
        assert!((dihedral - target).abs() < 1e-8);
    }

    #[test]
    fn test_set_dihedral_degenerate() {
        let mol = parse("CCC").unwrap();
        let coords = crate::dg::generate_coords(&mol);

        // Degenerate case: C and D are too close (not a valid 4-atom chain)
        // Should return coords unchanged or early
        let result = set_dihedral(
            &coords,
            &mol,
            AtomIdx(0),
            AtomIdx(1),
            AtomIdx(2),
            AtomIdx(2),
            180.0,
        );
        let p = result.get(AtomIdx(0));
        let p_orig = coords.get(AtomIdx(0));

        // Degenerate input should not crash; coords may be unchanged
        assert!((p.x - p_orig.x).abs() < 1e-8);
    }

    #[test]
    fn set_bond_length_moves_the_selected_side_as_a_rigid_fragment() {
        let mol = parse("CCCC").unwrap();
        let mut coords = Coords3D::new_zeroed(4);
        for i in 0..4 {
            coords.set(AtomIdx(i), Point3::new(i as f64, 0.0, 0.0));
        }

        let result = set_bond_length(&coords, &mol, AtomIdx(0), AtomIdx(1), 2.0).unwrap();
        assert!((get_bond_length(&result, AtomIdx(0), AtomIdx(1)) - 2.0).abs() < 1e-12);
        assert_eq!(result.get(AtomIdx(0)), Point3::new(0.0, 0.0, 0.0));
        assert_eq!(result.get(AtomIdx(1)), Point3::new(2.0, 0.0, 0.0));
        assert_eq!(result.get(AtomIdx(2)), Point3::new(3.0, 0.0, 0.0));
        assert_eq!(result.get(AtomIdx(3)), Point3::new(4.0, 0.0, 0.0));
    }

    #[test]
    fn set_bond_length_rejects_ring_bonds() {
        let mol = parse("C1CCCCC1").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        assert!(matches!(
            set_bond_length(&coords, &mol, AtomIdx(0), AtomIdx(1), 2.0),
            Err(MolTransformError::RingBond { .. })
        ));
    }

    #[test]
    fn set_bond_angle_moves_the_selected_side_as_a_rigid_fragment() {
        let mol = parse("CCCC").unwrap();
        let mut coords = Coords3D::new_zeroed(4);
        coords.set(AtomIdx(0), Point3::new(1.0, 0.0, 0.0));
        coords.set(AtomIdx(1), Point3::new(0.0, 0.0, 0.0));
        coords.set(AtomIdx(2), Point3::new(0.0, 1.0, 0.0));
        coords.set(AtomIdx(3), Point3::new(0.0, 2.0, 0.0));

        let result = set_bond_angle(
            &coords,
            &mol,
            AtomIdx(0),
            AtomIdx(1),
            AtomIdx(2),
            60.0_f64.to_radians(),
        )
        .unwrap();
        assert!(
            (get_bond_angle_deg(&result, AtomIdx(0), AtomIdx(1), AtomIdx(2)) - 60.0).abs() < 1e-10
        );
        assert!((get_bond_length(&result, AtomIdx(2), AtomIdx(3)) - 1.0).abs() < 1e-12);
        assert_eq!(result.get(AtomIdx(0)), coords.get(AtomIdx(0)));
        assert_eq!(result.get(AtomIdx(1)), coords.get(AtomIdx(1)));
    }

    #[test]
    fn set_bond_angle_rejects_two_ring_bonds() {
        let mol = parse("C1CCCCC1").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        assert_eq!(
            set_bond_angle(
                &coords,
                &mol,
                AtomIdx(0),
                AtomIdx(1),
                AtomIdx(2),
                90.0_f64.to_radians(),
            ),
            Err(MolTransformError::BothAngleBondsInRing)
        );
    }

    #[test]
    fn checked_dihedral_rejects_a_ring_central_bond() {
        let mol = parse("C1CCCCC1").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        assert!(matches!(
            try_set_dihedral(
                &coords,
                &mol,
                AtomIdx(0),
                AtomIdx(1),
                AtomIdx(2),
                AtomIdx(3),
                60.0_f64.to_radians(),
            ),
            Err(MolTransformError::RingBond { .. })
        ));
    }

    #[test]
    fn checked_transforms_reject_coordinate_count_and_non_finite_targets() {
        let mol = parse("CC").unwrap();
        let coords = Coords3D::new_zeroed(1);
        assert!(matches!(
            set_bond_length(&coords, &mol, AtomIdx(0), AtomIdx(1), 1.5),
            Err(MolTransformError::CoordinateCountMismatch { .. })
        ));

        let coords = crate::dg::generate_coords(&mol);
        assert_eq!(
            set_bond_length(&coords, &mol, AtomIdx(0), AtomIdx(1), f64::NAN),
            Err(MolTransformError::NonFiniteTarget)
        );
    }

    #[test]
    fn checked_measurements_reject_degenerate_coordinates() {
        let mol = parse("CCC").unwrap();
        let mut coords = Coords3D::new_zeroed(3);
        coords.set(AtomIdx(2), Point3::new(1.0, 0.0, 0.0));
        assert!(matches!(
            try_get_bond_angle_deg(&coords, &mol, AtomIdx(0), AtomIdx(1), AtomIdx(2)),
            Err(MolTransformError::CoincidentAtoms { .. })
        ));
    }
}
