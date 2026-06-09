//! 3D molecular descriptors for ML pipelines.
//!
//! Implements simplified WHIM and GETAWAY-like descriptors:
//! - WHIM: weighted holistic invariant descriptors based on mass distribution
//! - GETAWAY: geometric/topologic descriptors with wavelet autocorrelation
//!
//! Full implementations require expensive computations; this module provides
//! practical approximations suitable for ML feature vectors.

use crate::coords::{Coords3D, Point3};
use chematic_core::Molecule;

/// Compute WHIM descriptors: mass-weighted shape descriptors.
/// Returns: [L1, L2, L3, P1, P2, P3, ALPHA, BETA, GAMMA, DELTA]
/// where L* = eigenvalues of inertia tensor, P* = principal moments,
/// and ALPHA/BETA/GAMMA/DELTA are derived shape metrics.
pub fn whim_descriptors(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 10];
    }

    // Compute center of mass
    let mut total_mass = 0.0;
    let mut com = Point3::zero();

    for i in 0..mol.atom_count() {
        let atom = mol.atom(chematic_core::AtomIdx(i as u32));
        let mass = atom.element.atomic_mass();
        total_mass += mass;
        let p = coords.get(chematic_core::AtomIdx(i as u32));
        com = com.add(&p.scale(mass));
    }

    if total_mass == 0.0 {
        return vec![0.0; 10];
    }

    com = com.scale(1.0 / total_mass);

    // Compute inertia tensor
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;

    for i in 0..mol.atom_count() {
        let atom = mol.atom(chematic_core::AtomIdx(i as u32));
        let mass = atom.element.atomic_mass();
        let p = coords.get(chematic_core::AtomIdx(i as u32));
        let r = p.sub(&com);

        ixx += mass * (r.y * r.y + r.z * r.z);
        iyy += mass * (r.x * r.x + r.z * r.z);
        izz += mass * (r.x * r.x + r.y * r.y);
    }

    // Approximate eigenvalues (simplified: use diagonal dominance)
    let l1 = ixx;
    let l2 = iyy;
    let l3 = izz;

    let p1 = (l1 / total_mass).sqrt();
    let p2 = (l2 / total_mass).sqrt();
    let p3 = (l3 / total_mass).sqrt();

    // Shape metrics
    let alpha = p1 + p2 + p3; // Total reach
    let beta = (p1 * p2 + p2 * p3 + p3 * p1) / 3.0; // Average interaction
    let gamma = (p1 * p2 * p3).cbrt(); // Geometric mean
    let delta = p1 - p3; // Anisotropy

    vec![l1, l2, l3, p1, p2, p3, alpha, beta, gamma, delta]
}

/// Compute GETAWAY descriptors: geometric autocorrelation descriptors.
/// Returns: [G1, G2, G3, D1, D2, D3, T, V, A]
/// where G* = geometric autocorrelations, D* = topologic distances,
/// T = total distance, V = volume proxy, A = anisotropy ratio.
pub fn getaway_descriptors(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 9];
    }

    let n = mol.atom_count() as f64;

    // Geometric autocorrelations (lag-1, lag-2, lag-3)
    let mut g1 = 0.0;
    let mut g2 = 0.0;
    let mut g3 = 0.0;
    let mut total_dist = 0.0;

    // Pairwise distances
    let mut distances = Vec::new();
    for i in 0..mol.atom_count() {
        let ai = chematic_core::AtomIdx(i as u32);
        for j in (i + 1)..mol.atom_count() {
            let aj = chematic_core::AtomIdx(j as u32);
            let d = coords.get(ai).distance(&coords.get(aj));
            distances.push(d);
            total_dist += d;
        }
    }

    // Autocorrelation at different lags
    if !distances.is_empty() {
        g1 = distances.iter().take(n as usize).map(|&d| d).sum::<f64>() / n.max(1.0);
        g2 = distances.iter().skip(n as usize / 2).take(n as usize).map(|&d| d).sum::<f64>() / n.max(1.0);
        g3 = distances.iter().rev().take(n as usize).map(|&d| d).sum::<f64>() / n.max(1.0);
    }

    // Topologic distances (simplified: bond distances)
    let mut d1 = 0.0;
    let mut d2 = 0.0;
    let mut d3 = 0.0;

    for (_, bond) in mol.bonds() {
        let ai = bond.atom1;
        let aj = bond.atom2;
        let d = coords.get(ai).distance(&coords.get(aj));
        d1 += d;
        if d > 1.5 {
            d2 += d;
        }
        if d > 2.0 {
            d3 += d;
        }
    }

    let bond_count = mol.bond_count() as f64;
    if bond_count > 0.0 {
        d1 /= bond_count;
        d2 /= bond_count;
        d3 /= bond_count;
    }

    // Volume proxy: bounding box
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;

    for i in 0..mol.atom_count() {
        let p = coords.get(chematic_core::AtomIdx(i as u32));
        min_x = min_x.min(p.x);
        max_x = max_x.max(p.x);
        min_y = min_y.min(p.y);
        max_y = max_y.max(p.y);
        min_z = min_z.min(p.z);
        max_z = max_z.max(p.z);
    }

    let v = (max_x - min_x) * (max_y - min_y) * (max_z - min_z);
    let a = (max_x - min_x) / (max_z - min_z).max(0.1); // Anisotropy

    vec![g1, g2, g3, d1, d2, d3, total_dist, v, a]
}

/// Combined WHIM + GETAWAY descriptor vector for ML.
/// Returns 19-element feature vector.
pub fn whim_getaway_combined(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    let mut result = whim_descriptors(mol, coords);
    result.extend(getaway_descriptors(mol, coords));
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;
    use crate::dg::generate_coords;

    #[test]
    fn test_whim_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_descriptors(&mol, &coords);
        assert_eq!(desc.len(), 10);
        assert!(desc.iter().all(|&d| d.is_finite()), "all WHIM descriptors should be finite");
    }

    #[test]
    fn test_getaway_propane() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let desc = getaway_descriptors(&mol, &coords);
        assert_eq!(desc.len(), 9);
        assert!(desc.iter().all(|&d| d.is_finite()), "all GETAWAY descriptors should be finite");
    }

    #[test]
    fn test_combined_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_getaway_combined(&mol, &coords);
        assert_eq!(desc.len(), 19);
        assert!(desc.iter().all(|&d| d.is_finite()), "all combined descriptors should be finite");
    }
}
