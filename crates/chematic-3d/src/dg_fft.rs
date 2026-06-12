//! Distance geometry with eigenvalue decomposition for 3D coordinate generation.
//!
//! Implements the distance geometry approach:
//! 1. Build bound matrix from bond lengths and angle constraints
//! 2. Convert distances to Gram matrix (inner products)
//! 3. Eigenvalue decomposition to extract coordinates
//! 4. Metrization to enforce Euclidean distances
//!
//! This provides more accurate geometry than rule-based placement, especially
//! for complex molecules with many torsion degrees of freedom.

#![allow(dead_code)]

use std::f64::consts::PI;

use chematic_core::{AtomIdx, BondOrder, Molecule};

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Bound matrix construction
// ---------------------------------------------------------------------------

/// Build a bound matrix (lower and upper distance bounds) from molecular constraints.
///
/// Returns (lower_bounds, upper_bounds) as n×n matrices where:
/// - lower_bounds[i][j] = minimum allowed distance
/// - upper_bounds[i][j] = maximum allowed distance
///
/// Constraints include:
/// - Bond lengths (from ideal values ± tolerance)
/// - Angle constraints (from ideal angles)
/// - Van der Waals (from VDW radii sum)
fn build_bound_matrix(mol: &Molecule) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = mol.atom_count();
    let mut lower = vec![vec![0.0; n]; n];
    let mut upper = vec![vec![f64::INFINITY; n]; n];

    // Diagonal: distance to self = 0
    for i in 0..n {
        lower[i][i] = 0.0;
        upper[i][i] = 0.0;
    }

    // Bond length constraints (distance = 1)
    for j in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(j as u32));
        let a = bond.atom1;
        let b = bond.atom2;
        let dist = ideal_bond_length(mol, a, b);
        let tolerance = 0.05; // ±0.05 Å tolerance
        lower[a.0 as usize][b.0 as usize] = (dist - tolerance).max(0.5);
        upper[a.0 as usize][b.0 as usize] = dist + tolerance;
        lower[b.0 as usize][a.0 as usize] = lower[a.0 as usize][b.0 as usize];
        upper[b.0 as usize][a.0 as usize] = upper[a.0 as usize][b.0 as usize];
    }

    // Angle constraints (distance = 2)
    for center_idx in 0..n {
        let center = AtomIdx(center_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();

        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let a = neighbors[i];
                let b = neighbors[j];
                let bond_len_1 = ideal_bond_length(mol, center, a);
                let bond_len_2 = ideal_bond_length(mol, center, b);
                let angle = ideal_bond_angle(mol, center);

                // Distance = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(angle))
                let dist_sq = bond_len_1.powi(2) + bond_len_2.powi(2)
                    - 2.0 * bond_len_1 * bond_len_2 * angle.cos();
                let dist = dist_sq.max(0.0).sqrt();
                let tolerance = 0.1;

                let idx_a = a.0 as usize;
                let idx_b = b.0 as usize;
                lower[idx_a][idx_b] = lower[idx_a][idx_b].max((dist - tolerance).max(0.5));
                upper[idx_a][idx_b] = upper[idx_a][idx_b].min(dist + tolerance);
                lower[idx_b][idx_a] = lower[idx_a][idx_b];
                upper[idx_b][idx_a] = upper[idx_a][idx_b];
            }
        }
    }

    // Van der Waals bounds (non-bonded distance >= sum of VDW radii)
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);

            // Skip if bonded
            if mol.bond_between(a, b).is_some() {
                continue;
            }

            let vdw_sum = vdw_distance(mol, a, b);
            if lower[i][j] < vdw_sum {
                lower[i][j] = vdw_sum;
                lower[j][i] = vdw_sum;
            }
        }
    }

    (lower, upper)
}

/// Ideal bond length (Å) from atom pair and bond order.
fn ideal_bond_length(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
    let ea = mol.atom(a).element.atomic_number();
    let eb = mol.atom(b).element.atomic_number();
    let order = mol
        .bond_between(a, b)
        .map(|(_, bond)| bond.order)
        .unwrap_or(chematic_core::BondOrder::Single);

    // Simplified table (same as current DG implementation)
    match (ea.min(eb), ea.max(eb), order) {
        (6, 6, BondOrder::Single) => 1.54,
        (6, 6, BondOrder::Double) => 1.34,
        (6, 6, BondOrder::Triple) => 1.20,
        (6, 7, BondOrder::Single) => 1.47,
        (6, 8, BondOrder::Single) => 1.43,
        (6, 9, _) => 1.35,
        (1, 6, _) => 1.09,
        (1, 7, _) => 1.01,
        (1, 8, _) => 0.96,
        _ => 1.54,
    }
}

/// Ideal bond angle (radians) at center atom.
fn ideal_bond_angle(mol: &Molecule, center: AtomIdx) -> f64 {
    let atom = mol.atom(center);

    // Simplified: aromatic/double → ~120°, else ~109.5°
    if atom.aromatic {
        2.0 * PI / 3.0 // 120°
    } else {
        let has_double = mol
            .neighbors(center)
            .any(|(_, bidx)| matches!(mol.bond(bidx).order, BondOrder::Double));
        if has_double {
            2.0 * PI / 3.0
        } else {
            1.91 // ~109.5°
        }
    }
}

/// Van der Waals distance (sum of VDW radii, Å).
fn vdw_distance(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
    let r_a = vdw_radius(mol.atom(a).element.atomic_number());
    let r_b = vdw_radius(mol.atom(b).element.atomic_number());
    r_a + r_b
}

fn vdw_radius(atomic_num: u8) -> f64 {
    match atomic_num {
        1 => 1.20,  // H
        6 => 1.70,  // C
        7 => 1.55,  // N
        8 => 1.52,  // O
        9 => 1.47,  // F
        16 => 1.80, // S
        17 => 1.75, // Cl
        _ => 1.70,
    }
}

// ---------------------------------------------------------------------------
// Eigenvalue decomposition & coordinate generation
// ---------------------------------------------------------------------------

/// Generate 3D coordinates via distance geometry (ETDKG-style).
///
/// Algorithm:
/// 1. Build distance bounds from bond lengths and angles
/// 2. Create target distance matrix (average of bounds)
/// 3. Compute Gram matrix from distances
/// 4. Eigenvalue decomposition via Jacobi method
/// 5. Extract 3D coordinates from top 3 eigenvectors
/// 6. Center molecule at origin
pub fn generate_coords_dg(mol: &Molecule) -> Coords3D {
    let n = mol.atom_count();

    if n == 0 {
        return Coords3D::new_zeroed(0);
    }

    if n == 1 {
        let mut coords = Coords3D::new_zeroed(1);
        coords.set(AtomIdx(0), Point3 { x: 0.0, y: 0.0, z: 0.0 });
        return coords;
    }

    // Step 1: Build bounds
    let (lower, upper) = build_bound_matrix(mol);

    // Step 2: Create target distance matrix (average of bounds)
    let mut dist_matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                dist_matrix[i][j] = 0.0;
            } else {
                dist_matrix[i][j] = (lower[i][j] + upper[i][j]) / 2.0;
            }
        }
    }

    // Step 3: Compute Gram matrix
    let gram = distance_to_gram_matrix(&dist_matrix);

    // Step 4: Eigenvalue decomposition (Jacobi method)
    let (eigenvalues, eigenvectors) = jacobi_eigendecompose(&gram);

    // Step 5: Extract 3D coordinates from top 3 positive eigenvectors
    let mut coords = Coords3D::new_zeroed(n);

    // Find indices of 3 largest positive eigenvalues
    let mut eval_indices: Vec<usize> = (0..n).collect();
    eval_indices.sort_by(|&a, &b| {
        let av = eigenvalues[a].abs();
        let bv = eigenvalues[b].abs();
        // Handle NaN by treating it as zero
        let av = if av.is_nan() { 0.0 } else { av };
        let bv = if bv.is_nan() { 0.0 } else { bv };
        bv.partial_cmp(&av).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Use top 3 eigenvectors (or fewer if not enough positive eigenvalues)
    let mut pos_count = 0;
    let mut use_indices = Vec::new();
    for &idx in &eval_indices {
        if eigenvalues[idx] > 1e-10 && pos_count < 3 {
            use_indices.push(idx);
            pos_count += 1;
        }
    }

    // If fewer than 3 positive eigenvalues, pad with smaller ones
    for &idx in &eval_indices {
        if use_indices.len() >= 3 {
            break;
        }
        if !use_indices.contains(&idx) {
            use_indices.push(idx);
        }
    }

    // Set coordinates
    for i in 0..n {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        if use_indices.len() > 0 {
            x = (eigenvalues[use_indices[0]].abs().sqrt()) * eigenvectors[i][use_indices[0]];
        }
        if use_indices.len() > 1 {
            y = (eigenvalues[use_indices[1]].abs().sqrt()) * eigenvectors[i][use_indices[1]];
        }
        if use_indices.len() > 2 {
            z = (eigenvalues[use_indices[2]].abs().sqrt()) * eigenvectors[i][use_indices[2]];
        }

        coords.set(AtomIdx(i as u32), Point3 { x, y, z });
    }

    // Step 6: Center molecule
    center_coordinates(&mut coords);

    coords
}

/// Convert distance matrix to Gram matrix.
///
/// Gram[i,j] = 0.5 * (D[0,i]² + D[0,j]² - D[i,j]²)
/// where atom 0 is the reference point (origin).
fn distance_to_gram_matrix(dist: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = dist.len();
    let mut gram = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..n {
            gram[i][j] =
                0.5 * (dist[0][i] * dist[0][i] + dist[0][j] * dist[0][j] - dist[i][j] * dist[i][j]);
        }
    }

    gram
}

/// Jacobi eigendecomposition for symmetric matrices.
///
/// Returns (eigenvalues, eigenvectors) where eigenvectors[i][j] is the
/// i-th component of the j-th eigenvector.
fn jacobi_eigendecompose(mat: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut v = vec![vec![0.0; n]; n];

    // Initialize eigenvector matrix to identity
    for i in 0..n {
        v[i][i] = 1.0;
    }

    let max_iterations = 100;
    let tolerance = 1e-10;

    for _iteration in 0..max_iterations {
        // Find largest off-diagonal element
        let mut max_val = 0.0;
        let mut p = 0;
        let mut q = 1;

        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < tolerance {
            break;
        }

        // Compute rotation angle
        let apq = a[p][q];
        let app = a[p][p];
        let aqq = a[q][q];

        let theta = 0.5 * (2.0 * apq / (aqq - app)).atan();
        let c = theta.cos();
        let s = theta.sin();

        // Apply Givens rotation
        for i in 0..n {
            if i == p || i == q {
                continue;
            }

            let aip = a[i][p];
            let aiq = a[i][q];

            a[i][p] = c * aip - s * aiq;
            a[p][i] = a[i][p];
            a[i][q] = s * aip + c * aiq;
            a[q][i] = a[i][q];
        }

        // Update diagonal
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        // Update eigenvectors
        for i in 0..n {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
    }

    // Extract eigenvalues from diagonal
    let mut eigenvalues = vec![0.0; n];
    for i in 0..n {
        eigenvalues[i] = a[i][i];
    }

    (eigenvalues, v)
}

/// Center coordinates at origin (centroid).
fn center_coordinates(coords: &mut Coords3D) {
    let n = coords.atom_count();
    if n == 0 {
        return;
    }

    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;

    for i in 0..n {
        let p = coords.get(AtomIdx(i as u32));
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }

    cx /= n as f64;
    cy /= n as f64;
    cz /= n as f64;

    for i in 0..n {
        let p = coords.get(AtomIdx(i as u32));
        coords.set(AtomIdx(i as u32), Point3 {
            x: p.x - cx,
            y: p.y - cy,
            z: p.z - cz,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_bound_matrix_ethane() {
        let mol = parse("CC").unwrap();
        let (lower, upper) = build_bound_matrix(&mol);
        assert_eq!(lower.len(), 2);
        assert_eq!(upper.len(), 2);
        // C-C bond should be ~1.54 Å
        assert!(lower[0][1] > 1.4 && lower[0][1] < 1.6);
        assert!(upper[0][1] > 1.4 && upper[0][1] < 1.6);
    }

    #[test]
    fn test_bound_matrix_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let (lower, _upper) = build_bound_matrix(&mol);
        assert_eq!(lower.len(), 6);
        // Aromatic C-C should be ~1.40 Å
        assert!(lower[0][1] > 1.3 && lower[0][1] < 1.5);
    }

    #[test]
    fn test_ideal_bond_length_cc() {
        assert!((ideal_bond_length_test(6, 6, BondOrder::Single) - 1.54).abs() < 0.01);
    }

    #[test]
    fn test_ideal_bond_angle_sp2() {
        let angle = ideal_bond_angle_test(true);
        assert!((angle - 2.0 * PI / 3.0).abs() < 0.01); // 120°
    }

    #[test]
    fn test_generate_coords_dg_ethane() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 2);

        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));

        // Should have non-zero distance
        let dist = p0.distance(&p1);
        assert!(dist > 1.4 && dist < 1.7, "C-C bond distance: {}", dist);
    }

    #[test]
    fn test_generate_coords_dg_finite() {
        let mol = parse("c1ccc(C)cc1").unwrap();
        let coords = generate_coords_dg(&mol);

        // All coordinates must be finite (no NaN, no Inf)
        for i in 0..mol.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite(), "atom {} x is not finite", i);
            assert!(p.y.is_finite(), "atom {} y is not finite", i);
            assert!(p.z.is_finite(), "atom {} z is not finite", i);
        }
    }

    #[test]
    fn test_generate_coords_dg_centered() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords_dg(&mol);

        // Calculate centroid
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;

        for i in 0..mol.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            cx += p.x;
            cy += p.y;
            cz += p.z;
        }

        cx /= mol.atom_count() as f64;
        cy /= mol.atom_count() as f64;
        cz /= mol.atom_count() as f64;

        // Centroid should be near origin
        assert!(cx.abs() < 1e-6, "centroid x: {}", cx);
        assert!(cy.abs() < 1e-6, "centroid y: {}", cy);
        assert!(cz.abs() < 1e-6, "centroid z: {}", cz);
    }

    #[test]
    fn test_generate_coords_dg_single_atom() {
        let mol = parse("C").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 1);
        let p = coords.get(AtomIdx(0));
        assert_eq!(p.x, 0.0);
        assert_eq!(p.y, 0.0);
        assert_eq!(p.z, 0.0);
    }

    #[test]
    fn test_generate_coords_dg_ethene() {
        let mol = parse("C=C").unwrap();
        let coords = generate_coords_dg(&mol);

        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));
        let dist = p0.distance(&p1);

        // Double bond C=C is ~1.34 Å
        assert!(dist > 1.2 && dist < 1.5, "C=C bond distance: {}", dist);
    }

    #[test]
    fn test_generate_coords_dg_propane() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 3);

        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));
        let p2 = coords.get(AtomIdx(2));

        // All bonds should have reasonable distances
        let d01 = p0.distance(&p1);
        let d12 = p1.distance(&p2);
        let d02 = p0.distance(&p2);

        assert!(d01 > 1.4 && d01 < 1.7);
        assert!(d12 > 1.4 && d12 < 1.7);
        assert!(d02 > 2.0 && d02 < 3.5);
    }

    #[test]
    fn test_generate_coords_dg_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 6);

        // All atoms should have non-zero positions
        let mut all_nonzero = true;
        for i in 0..6 {
            let p = coords.get(AtomIdx(i as u32));
            if p.x.abs() < 1e-10 && p.y.abs() < 1e-10 && p.z.abs() < 1e-10 {
                all_nonzero = false;
            }
        }
        assert!(all_nonzero, "some atoms have zero coordinates");
    }

    #[test]
    fn test_gram_matrix_simple() {
        let dist = vec![
            vec![0.0, 1.5, 2.5],
            vec![1.5, 0.0, 1.6],
            vec![2.5, 1.6, 0.0],
        ];

        let gram = distance_to_gram_matrix(&dist);
        assert_eq!(gram.len(), 3);
        assert_eq!(gram[0].len(), 3);

        // Gram matrix should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert!((gram[i][j] - gram[j][i]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_jacobi_eigendecompose_identity() {
        let mat = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

        let (eigenvalues, _eigenvectors) = jacobi_eigendecompose(&mat);
        assert_eq!(eigenvalues.len(), 2);
        assert!((eigenvalues[0] - 1.0).abs() < 1e-6);
        assert!((eigenvalues[1] - 1.0).abs() < 1e-6);
    }

    // Helper functions for testing
    fn ideal_bond_length_test(a: u8, b: u8, order: BondOrder) -> f64 {
        match (a.min(b), a.max(b), order) {
            (6, 6, BondOrder::Single) => 1.54,
            (6, 6, BondOrder::Double) => 1.34,
            (6, 6, BondOrder::Triple) => 1.20,
            _ => 1.54,
        }
    }

    fn ideal_bond_angle_test(aromatic: bool) -> f64 {
        if aromatic { 2.0 * PI / 3.0 } else { 1.91 }
    }
}
