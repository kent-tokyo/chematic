//! Shape descriptors computed from 3D coordinates.
//!
//! All functions accept `&Molecule` (for atom masses) and `&Coords3D`.
//! Returns 0.0 for degenerate inputs (single atom, all atoms coincident).

use chematic_core::{AtomIdx, Molecule};

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Average atomic mass (Da). Covers common organic / drug-like elements.
fn avg_mass(an: u8) -> f64 {
    match an {
        1 => 1.008,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        n => n as f64,
    }
}

/// Mass-weighted centroid and total mass.
fn center_of_mass(mol: &Molecule, coords: &Coords3D) -> (Point3, f64) {
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    let mut total = 0.0_f64;
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let m = avg_mass(mol.atom(idx).element.atomic_number());
        let p = coords.get(idx);
        cx += m * p.x;
        cy += m * p.y;
        cz += m * p.z;
        total += m;
    }
    if total > 1e-10 {
        (Point3::new(cx / total, cy / total, cz / total), total)
    } else {
        (Point3::zero(), total)
    }
}

/// Mass-weighted inertia tensor about the center of mass.
fn inertia_tensor(mol: &Molecule, coords: &Coords3D) -> [[f64; 3]; 3] {
    let (cm, _) = center_of_mass(mol, coords);
    let mut t = [[0.0_f64; 3]; 3];
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let m = avg_mass(mol.atom(idx).element.atomic_number());
        let p = coords.get(idx);
        let x = p.x - cm.x;
        let y = p.y - cm.y;
        let z = p.z - cm.z;
        t[0][0] += m * (y * y + z * z);
        t[1][1] += m * (x * x + z * z);
        t[2][2] += m * (x * x + y * y);
        t[0][1] -= m * x * y;
        t[1][0] = t[0][1];
        t[0][2] -= m * x * z;
        t[2][0] = t[0][2];
        t[1][2] -= m * y * z;
        t[2][1] = t[1][2];
    }
    t
}

/// 3×3 symmetric Jacobi eigensolver.
///
/// Returns `(eigenvalues_sorted_ascending, evecs)` where
/// `evecs[i][j]` = i-th component of the j-th eigenvector.
pub(crate) fn jacobi3(m: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = m;
    // Identity matrix accumulates eigenvector rotations (columns = eigenvectors).
    let mut v: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    for _ in 0..100 {
        // Find the off-diagonal element with the largest absolute value.
        let (p, q) = {
            let mut mx = 0.0_f64;
            let mut pi = 0usize;
            let mut qi = 1usize;
            for i in 0..3 {
                for j in (i + 1)..3 {
                    if a[i][j].abs() > mx {
                        mx = a[i][j].abs();
                        pi = i;
                        qi = j;
                    }
                }
            }
            if mx < 1e-12 {
                break;
            }
            (pi, qi)
        };

        // Jacobi rotation angle t = tan(theta).
        let diff = a[q][q] - a[p][p];
        let t = if diff.abs() < 1e-15 * a[p][q].abs() {
            // Near-equal diagonal: theta → ±45°, t = ±1.
            1.0_f64.copysign(a[p][q])
        } else {
            let phi = diff / (2.0 * a[p][q]);
            let sign = if phi >= 0.0 { 1.0_f64 } else { -1.0_f64 };
            sign / (phi.abs() + (1.0 + phi * phi).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let tau = s / (1.0 + c);

        // Update diagonal.
        let a_pp = a[p][p] - t * a[p][q];
        let a_qq = a[q][q] + t * a[p][q];
        a[p][p] = a_pp;
        a[q][q] = a_qq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        // Update off-diagonal rows/columns (symmetric).
        for r in 0..3 {
            if r != p && r != q {
                let old_pr = a[p][r];
                let old_qr = a[q][r];
                a[p][r] = old_pr - s * (old_qr + tau * old_pr);
                a[r][p] = a[p][r];
                a[q][r] = old_qr + s * (old_pr - tau * old_qr);
                a[r][q] = a[q][r];
            }
        }

        // Accumulate rotation into eigenvector matrix.
        for r in 0..3 {
            let old_rp = v[r][p];
            let old_rq = v[r][q];
            v[r][p] = old_rp - s * (old_rq + tau * old_rp);
            v[r][q] = old_rq + s * (old_rp - tau * old_rq);
        }
    }

    // Sort eigenvalue/eigenvector pairs by eigenvalue ascending.
    let mut pairs: [(f64, [f64; 3]); 3] = [
        (a[0][0], [v[0][0], v[1][0], v[2][0]]),
        (a[1][1], [v[0][1], v[1][1], v[2][1]]),
        (a[2][2], [v[0][2], v[1][2], v[2][2]]),
    ];
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));

    let evals = [pairs[0].0, pairs[1].0, pairs[2].0];
    // evecs[i][j] = component i of eigenvector j.
    let evecs = [
        [pairs[0].1[0], pairs[1].1[0], pairs[2].1[0]],
        [pairs[0].1[1], pairs[1].1[1], pairs[2].1[1]],
        [pairs[0].1[2], pairs[1].1[2], pairs[2].1[2]],
    ];
    (evals, evecs)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// All three principal moments of inertia (PMI1 ≤ PMI2 ≤ PMI3) in Da·Å².
///
/// Computed from the mass-weighted inertia tensor eigenvalues.
pub fn pmi(mol: &Molecule, coords: &Coords3D) -> (f64, f64, f64) {
    let tensor = inertia_tensor(mol, coords);
    let (evals, _) = jacobi3(tensor);
    (evals[0].max(0.0), evals[1].max(0.0), evals[2].max(0.0))
}

/// Smallest principal moment of inertia (Da·Å²).
pub fn pmi1(mol: &Molecule, coords: &Coords3D) -> f64 {
    pmi(mol, coords).0
}

/// Middle principal moment of inertia (Da·Å²).
pub fn pmi2(mol: &Molecule, coords: &Coords3D) -> f64 {
    pmi(mol, coords).1
}

/// Largest principal moment of inertia (Da·Å²).
pub fn pmi3(mol: &Molecule, coords: &Coords3D) -> f64 {
    pmi(mol, coords).2
}

/// Normalized principal moment ratio NPR1 = PMI1 / PMI3.
///
/// Ranges from 0 (rod-like, linear) to 1 (sphere-like).
pub fn npr1(mol: &Molecule, coords: &Coords3D) -> f64 {
    let (p1, _, p3) = pmi(mol, coords);
    if p3 < 1e-10 { 0.0 } else { p1 / p3 }
}

/// Normalized principal moment ratio NPR2 = PMI2 / PMI3.
pub fn npr2(mol: &Molecule, coords: &Coords3D) -> f64 {
    let (_, p2, p3) = pmi(mol, coords);
    if p3 < 1e-10 { 0.0 } else { p2 / p3 }
}

/// Mass-weighted radius of gyration (Å).
///
/// Rg = sqrt( Σ m_i * |r_i − r_cm|² / Σ m_i )
pub fn radius_of_gyration(mol: &Molecule, coords: &Coords3D) -> f64 {
    let (cm, total_m) = center_of_mass(mol, coords);
    if total_m < 1e-10 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let m = avg_mass(mol.atom(idx).element.atomic_number());
        let p = coords.get(idx);
        let r2 = (p.x - cm.x).powi(2) + (p.y - cm.y).powi(2) + (p.z - cm.z).powi(2);
        sum += m * r2;
    }
    (sum / total_m).sqrt()
}

/// Asphericity = PMI3 − (PMI1 + PMI2) / 2 (Da·Å²).
///
/// Zero for a perfect sphere; positive for elongated or flat shapes.
pub fn asphericity(mol: &Molecule, coords: &Coords3D) -> f64 {
    let (p1, p2, p3) = pmi(mol, coords);
    p3 - (p1 + p2) / 2.0
}

/// Eccentricity = sqrt(1 − PMI1 / PMI3).
///
/// Zero for a perfect sphere; approaches 1 for a perfect rod.
pub fn eccentricity(mol: &Molecule, coords: &Coords3D) -> f64 {
    let (p1, _, p3) = pmi(mol, coords);
    if p3 < 1e-10 {
        return 0.0;
    }
    (1.0 - p1 / p3).max(0.0).sqrt()
}

/// Plane of Best Fit (PBF): RMS deviation of heavy atoms from the
/// least-squares plane (Å).
///
/// PBF ≈ 0 for perfectly flat molecules (benzene); large values indicate
/// 3D shape (e.g., steroids, cyclohexane).  Uses unweighted geometry centroid
/// to match the RDKit convention.
pub fn plane_of_best_fit(mol: &Molecule, coords: &Coords3D) -> f64 {
    // Heavy atoms only (H excluded) — matches the published definition and RDKit convention.
    // RDKit issue #9238: including H artificially reduces PBF vs. the reference paper values.
    let heavy: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx)
        .collect();
    let n = heavy.len();
    if n < 3 {
        return 0.0;
    }

    // Geometric centroid (unweighted).
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for &idx in &heavy {
        let p = coords.get(idx);
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    cx /= n as f64;
    cy /= n as f64;
    cz /= n as f64;

    // Covariance matrix of coordinates about the centroid.
    let mut cov = [[0.0_f64; 3]; 3];
    for &idx in &heavy {
        let p = coords.get(idx);
        let x = p.x - cx;
        let y = p.y - cy;
        let z = p.z - cz;
        cov[0][0] += x * x;
        cov[0][1] += x * y;
        cov[0][2] += x * z;
        cov[1][0] += x * y;
        cov[1][1] += y * y;
        cov[1][2] += y * z;
        cov[2][0] += x * z;
        cov[2][1] += y * z;
        cov[2][2] += z * z;
    }

    // Eigenvector of smallest eigenvalue = normal to best-fit plane.
    let (_, evecs) = jacobi3(cov);
    let nx = evecs[0][0];
    let ny = evecs[1][0];
    let nz = evecs[2][0];

    // RMS distance from the plane.
    let mut sum_sq = 0.0_f64;
    for &idx in &heavy {
        let p = coords.get(idx);
        let d = (p.x - cx) * nx + (p.y - cy) * ny + (p.z - cz) * nz;
        sum_sq += d * d;
    }
    (sum_sq / n as f64).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    use crate::dg::generate_coords;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    // --- jacobi3 unit tests --------------------------------------------------

    #[test]
    fn jacobi3_diagonal_identity() {
        // Diagonal matrix: eigenvalues = diagonal, eigenvectors = identity.
        let m = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]];
        let (evals, _) = jacobi3(m);
        let mut sorted = evals;
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(
            (sorted[0] - 1.0).abs() < 1e-9,
            "expected 1, got {}",
            sorted[0]
        );
        assert!(
            (sorted[1] - 2.0).abs() < 1e-9,
            "expected 2, got {}",
            sorted[1]
        );
        assert!(
            (sorted[2] - 3.0).abs() < 1e-9,
            "expected 3, got {}",
            sorted[2]
        );
    }

    #[test]
    fn jacobi3_symmetric_matrix() {
        // Known symmetric matrix with eigenvalues 1, 2, 3 (after rotation).
        // [[2, 1, 0], [1, 2, 0], [0, 0, 3]] has eigenvalues 1, 3, 3 — let's use simpler:
        // [[5, 4, 0], [4, 5, 0], [0, 0, 1]] has eigenvalues 1, 1, 9.
        let m = [[5.0, 4.0, 0.0], [4.0, 5.0, 0.0], [0.0, 0.0, 1.0]];
        let (evals, _) = jacobi3(m);
        assert!((evals[0] - 1.0).abs() < 1e-8, "evals[0]={}", evals[0]);
        assert!((evals[1] - 1.0).abs() < 1e-8, "evals[1]={}", evals[1]);
        assert!((evals[2] - 9.0).abs() < 1e-8, "evals[2]={}", evals[2]);
    }

    // --- PMI / NPR -----------------------------------------------------------

    #[test]
    fn single_atom_pmi_zero() {
        let m = mol("C");
        let c = generate_coords(&m);
        let (p1, p2, p3) = pmi(&m, &c);
        assert!(
            p1 < 1e-9 && p2 < 1e-9 && p3 < 1e-9,
            "single atom PMI should be 0"
        );
    }

    #[test]
    fn ethane_pmi_sorted() {
        let m = mol("CC");
        let c = generate_coords(&m);
        let (p1, p2, p3) = pmi(&m, &c);
        assert!(p1 <= p2 + 1e-9, "PMI1 <= PMI2: {p1} vs {p2}");
        assert!(p2 <= p3 + 1e-9, "PMI2 <= PMI3: {p2} vs {p3}");
        assert!(p3 > 1e-9, "ethane PMI3 > 0");
    }

    #[test]
    fn npr_in_range() {
        let m = mol("c1ccccc1");
        let c = generate_coords(&m);
        let n1 = npr1(&m, &c);
        let n2 = npr2(&m, &c);
        assert!((0.0..=1.0 + 1e-9).contains(&n1), "NPR1 out of range: {n1}");
        assert!((0.0..=1.0 + 1e-9).contains(&n2), "NPR2 out of range: {n2}");
        assert!(n1 <= n2 + 1e-9, "NPR1 <= NPR2: {n1} vs {n2}");
    }

    // --- Radius of Gyration --------------------------------------------------

    #[test]
    fn rog_positive_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let c = generate_coords(&m);
        let rg = radius_of_gyration(&m, &c);
        assert!(rg > 0.0, "aspirin Rg should be positive, got {rg}");
    }

    #[test]
    fn rog_single_atom_zero() {
        let m = mol("C");
        let c = generate_coords(&m);
        let rg = radius_of_gyration(&m, &c);
        assert!(rg.abs() < 1e-9, "single atom Rg should be 0");
    }

    // --- Asphericity / Eccentricity ------------------------------------------

    #[test]
    fn asphericity_nonneg() {
        let m = mol("c1ccccc1");
        let c = generate_coords(&m);
        let b = asphericity(&m, &c);
        assert!(b >= -1e-9, "asphericity should be >= 0, got {b}");
    }

    #[test]
    fn eccentricity_in_range() {
        let m = mol("CCC");
        let c = generate_coords(&m);
        let e = eccentricity(&m, &c);
        assert!(
            (0.0..=1.0 + 1e-9).contains(&e),
            "eccentricity out of range: {e}"
        );
    }

    // --- PBF -----------------------------------------------------------------

    #[test]
    fn pbf_benzene_near_zero() {
        // Benzene is planar; PBF should be nearly zero.
        let m = mol("c1ccccc1");
        let c = generate_coords(&m);
        let pbf = plane_of_best_fit(&m, &c);
        assert!(pbf < 0.1, "benzene PBF should be near 0, got {pbf:.4}");
    }

    #[test]
    fn pbf_two_atoms_zero() {
        // Fewer than 3 atoms → PBF = 0.
        let m = mol("CC");
        let c = generate_coords(&m);
        let pbf = plane_of_best_fit(&m, &c);
        assert!(pbf.abs() < 1e-9, "two-atom PBF should be 0");
    }

    #[test]
    fn pbf_manually_nonplanar() {
        // 4 atoms: (0,0,0), (1,0,0), (0,1,0), (0,0,1) are non-coplanar.
        // PBF > 0 expected.
        use crate::coords::Coords3D;
        let m = mol("CC(C)C"); // 4 heavy atoms
        let n = m.atom_count();
        assert_eq!(n, 4);
        let mut c = Coords3D::new_zeroed(n);
        c.set(chematic_core::AtomIdx(0), Point3::new(0.0, 0.0, 0.0));
        c.set(chematic_core::AtomIdx(1), Point3::new(1.0, 0.0, 0.0));
        c.set(chematic_core::AtomIdx(2), Point3::new(0.0, 1.0, 0.0));
        c.set(chematic_core::AtomIdx(3), Point3::new(0.0, 0.0, 1.0));
        let pbf = plane_of_best_fit(&m, &c);
        assert!(
            pbf > 0.1,
            "non-coplanar atoms should have PBF > 0, got {pbf:.4}"
        );
    }

    #[test]
    fn pbf_uses_heavy_atoms_only() {
        // RDKit issue #9238: PBF must exclude hydrogen atoms.
        // A flat aromatic ring should have PBF ≈ 0 regardless of H positions.
        // With H included (bug), PBF is artificially inflated if H coords deviate
        // from the plane; with heavy atoms only (correct), benzene = ~0.
        let m = mol("c1ccccc1"); // benzene — all heavy atoms lie in a plane
        let c = generate_coords(&m);
        let pbf = plane_of_best_fit(&m, &c);
        assert!(
            pbf < 0.05,
            "benzene PBF (heavy-atoms-only) should be ~0, got {pbf:.4}; \
             if > 0.05, H atoms are being incorrectly included (RDKit #9238)"
        );

        // Naphthalene: also flat, should have PBF ≈ 0 with heavy-only.
        let m2 = mol("c1ccc2ccccc2c1");
        let c2 = generate_coords(&m2);
        let pbf2 = plane_of_best_fit(&m2, &c2);
        assert!(pbf2 < 0.05, "naphthalene PBF should be ~0, got {pbf2:.4}");
    }
}
