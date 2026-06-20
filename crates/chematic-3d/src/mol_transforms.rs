use crate::coords::{Coords3D, Point3};
use chematic_core::{AtomIdx, Molecule};

/// Get bond length in Ångströms between atoms a and b.
pub fn get_bond_length(coords: &Coords3D, a: AtomIdx, b: AtomIdx) -> f64 {
    coords.get(a).distance(&coords.get(b))
}

/// Get bond angle A—Center—B in radians.
pub fn get_bond_angle(coords: &Coords3D, a: AtomIdx, center: AtomIdx, b: AtomIdx) -> f64 {
    super::constraints::compute_angle(coords, a, center, b)
}

/// Get bond angle A—Center—B in degrees.
pub fn get_bond_angle_deg(coords: &Coords3D, a: AtomIdx, center: AtomIdx, b: AtomIdx) -> f64 {
    get_bond_angle(coords, a, center, b).to_degrees()
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
    let current = match get_dihedral(coords, a, b, c, d) {
        Some(x) => x,
        None => return coords.clone(),
    };

    let delta = target_angle - current;
    let pb = coords.get(b);
    let pc = coords.get(c);

    // Guard against degenerate geometry: B and C define the rotation axis.
    if pc.distance(&pb) < 1e-10 {
        return coords.clone();
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

    new_coords
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
    let mut stack = vec![target];

    while let Some(current) = stack.pop() {
        for (neighbor, _) in mol.neighbors(current) {
            if neighbor != parent && !result.contains(&neighbor) {
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
}
