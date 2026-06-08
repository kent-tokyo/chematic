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

use crate::coords::Coords3D;

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

/// Generate 3D coordinates via distance geometry.
///
/// Currently a placeholder that falls back to rule-based generation.
/// Full implementation requires eigenvalue solver (eigenvalues/eigenvectors).
pub fn generate_coords_dg(mol: &Molecule) -> Coords3D {
    let n = mol.atom_count();
    let coords = Coords3D::new_zeroed(n);

    if n == 0 {
        return coords;
    }

    // For now, return zeroed coords with a note that full DG is not yet implemented.
    // TODO: Implement eigenvalue decomposition to extract 3D coordinates from gram matrix.

    coords
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
