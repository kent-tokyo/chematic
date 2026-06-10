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

/// AutoCorr3D: Moreau-Broto Self-Correlation (Euclidean Distance).
///
/// Compute self-correlation for 3D coordinates using binned Euclidean distances.
/// Each lag corresponds to a distance bin (k * 1Å):
/// - lag 1: 0-1 Å
/// - lag 2: 1-2 Å
/// - lag 3: 2-3 Å
/// - ... lag 8: 7-8 Å
///
/// For each lag, sum over all atom pairs (i,j) with distance in bin of v(i) * v(j),
/// where v(i) is the atomic mass (simplified feature for 3D self-correlation).
pub fn autocorr_3d(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 8];
    }

    let n = mol.atom_count();
    let mut result = vec![0.0; 8];

    for lag in 1..=8 {
        let lower = (lag - 1) as f64;
        let upper = lag as f64;
        let mut sum = 0.0;

        for i in 0..n {
            let idx_i = chematic_core::AtomIdx(i as u32);
            let atom_i = mol.atom(idx_i);
            let mass_i = atom_i.element.atomic_mass();
            let p_i = coords.get(idx_i);

            for j in (i + 1)..n {
                let idx_j = chematic_core::AtomIdx(j as u32);
                let atom_j = mol.atom(idx_j);
                let mass_j = atom_j.element.atomic_mass();
                let p_j = coords.get(idx_j);

                let dist = p_i.distance(&p_j);
                if dist >= lower && dist < upper {
                    sum += mass_i * mass_j;
                }
            }
        }
        result[lag - 1] = sum;
    }

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

    #[test]
    fn test_autocorr_3d_single_atom() {
        let mol = parse("C").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Single atom: no pairs → all zeros
        for val in ac {
            assert!((val - 0.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_autocorr_3d_ethane() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Ethane: C-C distance ≈ 1.54 Å → lag 2 (1-2 Å)
        // Mass of C ≈ 12.0, so product ≈ 144
        assert!(ac[0] < 1.0, "lag 1 (0-1Å) should be minimal: {}", ac[0]);
        assert!(ac[1] > 100.0, "lag 2 (1-2Å) should be ~144: {}", ac[1]);
    }

    #[test]
    fn test_autocorr_3d_propane() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Should have non-zero values in appropriate distance bins
        assert!(ac.iter().any(|&x| x > 0.0), "should have non-zero autocorr values");
        assert!(ac.iter().all(|&x| x.is_finite()), "all values should be finite");
    }

    #[test]
    fn test_autocorr_3d_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Benzene: ring structure with various distances
        assert!(ac.iter().any(|&x| x > 0.0), "benzene should have non-zero autocorr");
        assert!(ac.iter().all(|&x| x.is_finite()), "all values should be finite");
    }
}
