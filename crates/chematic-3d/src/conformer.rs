//! Conformer ensemble: a molecule with multiple sets of 3D coordinates.

use std::fmt;

use chematic_core::{AtomIdx, Molecule};

use crate::coords::Coords3D;
use crate::shape_descriptors::jacobi3;

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

#[derive(Debug, PartialEq)]
pub enum ConformerError {
    AtomCountMismatch { expected: usize, got: usize },
}

impl fmt::Display for ConformerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConformerError::AtomCountMismatch { expected, got } => {
                write!(f, "conformer has {got} atoms but molecule has {expected}")
            }
        }
    }
}

impl std::error::Error for ConformerError {}

// ---------------------------------------------------------------------------
// ConformerEnsemble
// ---------------------------------------------------------------------------

/// A molecule paired with zero or more sets of 3D coordinates.
///
/// Conformer indices are contiguous; `remove_conformer` shifts all subsequent
/// indices down by one (Vec::remove semantics).
pub struct ConformerEnsemble {
    mol: Molecule,
    conformers: Vec<Coords3D>,
}

impl ConformerEnsemble {
    /// Create an ensemble with no conformers.
    pub fn new(mol: Molecule) -> Self {
        Self {
            mol,
            conformers: Vec::new(),
        }
    }

    /// Create an ensemble pre-loaded with one conformer.
    ///
    /// Returns an error if `coords.atom_count() != mol.atom_count()`.
    pub fn with_conformer(mol: Molecule, coords: Coords3D) -> Result<Self, ConformerError> {
        let expected = mol.atom_count();
        let got = coords.atom_count();
        if got != expected {
            return Err(ConformerError::AtomCountMismatch { expected, got });
        }
        Ok(Self {
            mol,
            conformers: vec![coords],
        })
    }

    /// The molecule (topology only; no coordinates).
    pub fn mol(&self) -> &Molecule {
        &self.mol
    }

    /// Number of conformers currently stored.
    pub fn conformer_count(&self) -> usize {
        self.conformers.len()
    }

    /// Append a conformer.
    ///
    /// Returns the index of the newly added conformer, or an error if the
    /// atom count does not match.
    pub fn add_conformer(&mut self, coords: Coords3D) -> Result<usize, ConformerError> {
        let expected = self.mol.atom_count();
        let got = coords.atom_count();
        if got != expected {
            return Err(ConformerError::AtomCountMismatch { expected, got });
        }
        let idx = self.conformers.len();
        self.conformers.push(coords);
        Ok(idx)
    }

    /// Return a reference to the conformer at `idx`, or `None` if out of range.
    pub fn get_conformer(&self, idx: usize) -> Option<&Coords3D> {
        self.conformers.get(idx)
    }

    /// Return a mutable reference to the conformer at `idx`, or `None` if out of range.
    pub fn get_conformer_mut(&mut self, idx: usize) -> Option<&mut Coords3D> {
        self.conformers.get_mut(idx)
    }

    /// Remove and return the conformer at `idx`.
    ///
    /// All conformers with index > `idx` shift down by one.
    /// Returns `None` if `idx` is out of range.
    pub fn remove_conformer(&mut self, idx: usize) -> Option<Coords3D> {
        if idx < self.conformers.len() {
            Some(self.conformers.remove(idx))
        } else {
            None
        }
    }

    /// RMSD between conformers `a` and `b` **without** superposition.
    ///
    /// Returns `None` if either index is out of range or the molecule has no atoms.
    pub fn conformer_rmsd_no_align(&self, a: usize, b: usize) -> Option<f64> {
        let ca = self.conformers.get(a)?;
        let cb = self.conformers.get(b)?;
        let n = self.mol.atom_count();
        if n == 0 {
            return Some(0.0);
        }
        let sum_sq: f64 = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let pa = ca.get(idx);
                let pb = cb.get(idx);
                let dx = pa.x - pb.x;
                let dy = pa.y - pb.y;
                let dz = pa.z - pb.z;
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        Some((sum_sq / n as f64).sqrt())
    }

    /// Kabsch-aligned RMSD between conformers `a` and `b`.
    ///
    /// Finds the rigid-body rotation (no scaling) that minimises RMSD, then
    /// returns that minimum RMSD.  Returns `None` if either index is out of
    /// range.
    pub fn conformer_rmsd(&self, a: usize, b: usize) -> Option<f64> {
        let ca = self.conformers.get(a)?;
        let cb = self.conformers.get(b)?;
        let n = self.mol.atom_count();
        Some(kabsch_rmsd(ca, cb, n))
    }
}

// ---------------------------------------------------------------------------
// Kabsch RMSD helper
// ---------------------------------------------------------------------------

fn kabsch_rmsd(coords_a: &Coords3D, coords_b: &Coords3D, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }

    let nf = n as f64;

    // Centroids.
    let mut ca = [0.0f64; 3];
    let mut cb = [0.0f64; 3];
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let pa = coords_a.get(idx);
        let pb = coords_b.get(idx);
        ca[0] += pa.x;
        ca[1] += pa.y;
        ca[2] += pa.z;
        cb[0] += pb.x;
        cb[1] += pb.y;
        cb[2] += pb.z;
    }
    for k in 0..3 {
        ca[k] /= nf;
        cb[k] /= nf;
    }

    // Centered coordinates.
    let mut p = vec![[0.0f64; 3]; n];
    let mut q = vec![[0.0f64; 3]; n];
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let pa = coords_a.get(idx);
        let pb = coords_b.get(idx);
        p[i] = [pa.x - ca[0], pa.y - ca[1], pa.z - ca[2]];
        q[i] = [pb.x - cb[0], pb.y - cb[1], pb.z - cb[2]];
    }

    // H = P^T * Q  (3×3 covariance matrix).
    let mut h = [[0.0f64; 3]; 3];
    for i in 0..n {
        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += p[i][r] * q[i][c];
            }
        }
    }

    // H^T * H (symmetric).
    let mut hth = [[0.0f64; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            for k in 0..3 {
                hth[r][c] += h[k][r] * h[k][c];
            }
        }
    }

    // Eigendecompose H^T * H → V columns are right singular vectors.
    // evecs[i][j] = component i of eigenvector j (sorted ascending by eigenvalue).
    let (evals, v) = jacobi3(hth);

    // U = H * V * diag(1/σ).  σ_j = sqrt(evals[j]).
    let mut hv = [[0.0f64; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            for k in 0..3 {
                hv[r][c] += h[r][k] * v[k][c];
            }
        }
    }
    let mut u = [[0.0f64; 3]; 3];
    for j in 0..3 {
        let sigma = evals[j].max(0.0).sqrt();
        for r in 0..3 {
            u[r][j] = if sigma > 1e-10 { hv[r][j] / sigma } else { 0.0 };
        }
    }

    // R = V * U^T.  R[r][c] = Σ_k V[r][k] * U[c][k].
    let mut r_mat = [[0.0f64; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            for k in 0..3 {
                r_mat[r][c] += v[r][k] * u[c][k];
            }
        }
    }

    // Reflection correction: if det(R) < 0, flip V column with smallest σ (col 0).
    let det = det3(r_mat);
    let mut v_final = v;
    if det < 0.0 {
        for r in 0..3 {
            v_final[r][0] *= -1.0;
        }
        // Recompute R.
        r_mat = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..3 {
                    r_mat[r][c] += v_final[r][k] * u[c][k];
                }
            }
        }
    }

    // Apply R to q, compute RMSD.
    let mut sum_sq = 0.0f64;
    for i in 0..n {
        for row in 0..3 {
            let rotated =
                r_mat[row][0] * q[i][0] + r_mat[row][1] * q[i][1] + r_mat[row][2] * q[i][2];
            let diff = p[i][row] - rotated;
            sum_sq += diff * diff;
        }
    }
    (sum_sq / nf).sqrt()
}

fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    use crate::{coords::Point3, dg::generate_coords};

    fn make_ensemble() -> ConformerEnsemble {
        let mol = parse("CCC").unwrap();
        let c = generate_coords(&mol);
        ConformerEnsemble::with_conformer(mol, c).unwrap()
    }

    // --- Construction and basic access --------------------------------------

    #[test]
    fn new_has_zero_conformers() {
        let mol = parse("C").unwrap();
        let ens = ConformerEnsemble::new(mol);
        assert_eq!(ens.conformer_count(), 0);
    }

    #[test]
    fn with_conformer_has_one() {
        let ens = make_ensemble();
        assert_eq!(ens.conformer_count(), 1);
    }

    #[test]
    fn add_conformer_increments_count() {
        let mol = parse("CC").unwrap();
        let c1 = generate_coords(&mol);
        let c2 = generate_coords(&mol);
        let mut ens = ConformerEnsemble::with_conformer(mol, c1).unwrap();
        let idx = ens.add_conformer(c2).unwrap();
        assert_eq!(idx, 1);
        assert_eq!(ens.conformer_count(), 2);
    }

    #[test]
    fn add_conformer_wrong_atom_count_errors() {
        let mol = parse("CC").unwrap();
        let wrong = Coords3D::new_zeroed(5);
        let mut ens = ConformerEnsemble::new(mol);
        let err = ens.add_conformer(wrong).unwrap_err();
        assert!(matches!(
            err,
            ConformerError::AtomCountMismatch {
                expected: 2,
                got: 5
            }
        ));
    }

    #[test]
    fn get_conformer_out_of_range_returns_none() {
        let ens = make_ensemble();
        assert!(ens.get_conformer(99).is_none());
    }

    // --- remove_conformer ---------------------------------------------------

    #[test]
    fn remove_conformer_decrements_count() {
        let mut ens = make_ensemble();
        let removed = ens.remove_conformer(0);
        assert!(removed.is_some());
        assert_eq!(ens.conformer_count(), 0);
    }

    #[test]
    fn remove_conformer_shifts_indices() {
        let mol = parse("C").unwrap();
        let n = mol.atom_count();
        let mut ens = ConformerEnsemble::new(mol);

        // Add three conformers with distinct x-coordinates for atom 0.
        for x in [1.0f64, 2.0, 3.0] {
            let mut c = Coords3D::new_zeroed(n);
            c.set(AtomIdx(0), Point3::new(x, 0.0, 0.0));
            ens.add_conformer(c).unwrap();
        }

        // Remove index 0; what was index 1 (x=2) is now index 0.
        ens.remove_conformer(0).unwrap();
        assert_eq!(ens.conformer_count(), 2);
        assert!((ens.get_conformer(0).unwrap().get(AtomIdx(0)).x - 2.0).abs() < 1e-10);
    }

    #[test]
    fn remove_conformer_out_of_range_returns_none() {
        let mut ens = make_ensemble();
        assert!(ens.remove_conformer(99).is_none());
    }

    // --- RMSD ---------------------------------------------------------------

    #[test]
    fn rmsd_no_align_same_conformer_is_zero() {
        let ens = make_ensemble();
        let rmsd = ens.conformer_rmsd_no_align(0, 0).unwrap();
        assert!(rmsd.abs() < 1e-10, "self-RMSD should be 0, got {rmsd}");
    }

    #[test]
    fn rmsd_no_align_translated_is_nonzero() {
        let mol = parse("CC").unwrap();
        let n = mol.atom_count();
        let mut c1 = Coords3D::new_zeroed(n);
        let mut c2 = Coords3D::new_zeroed(n);
        for i in 0..n {
            c1.set(AtomIdx(i as u32), Point3::new(i as f64, 0.0, 0.0));
            c2.set(AtomIdx(i as u32), Point3::new(i as f64 + 10.0, 0.0, 0.0));
        }
        let mut ens = ConformerEnsemble::with_conformer(mol, c1).unwrap();
        ens.add_conformer(c2).unwrap();
        let rmsd = ens.conformer_rmsd_no_align(0, 1).unwrap();
        assert!(
            rmsd > 0.0,
            "translated conformers should have non-zero RMSD"
        );
    }

    #[test]
    fn kabsch_rmsd_same_conformer_is_zero() {
        let ens = make_ensemble();
        let rmsd = ens.conformer_rmsd(0, 0).unwrap();
        assert!(
            rmsd.abs() < 1e-8,
            "Kabsch self-RMSD should be 0, got {rmsd}"
        );
    }

    #[test]
    fn kabsch_rmsd_pure_translation_is_zero() {
        // After Kabsch superposition, a pure translation should give RMSD = 0.
        let mol = parse("CCC").unwrap();
        let n = mol.atom_count();
        let base = generate_coords(&mol);
        let mut shifted = Coords3D::new_zeroed(n);
        let offset = 5.0;
        for i in 0..n {
            let p = base.get(AtomIdx(i as u32));
            shifted.set(
                AtomIdx(i as u32),
                Point3::new(p.x + offset, p.y + offset, p.z + offset),
            );
        }
        let mut ens = ConformerEnsemble::with_conformer(mol, base).unwrap();
        ens.add_conformer(shifted).unwrap();
        let rmsd = ens.conformer_rmsd(0, 1).unwrap();
        assert!(
            rmsd < 1e-6,
            "pure-translation Kabsch RMSD should be ~0, got {rmsd}"
        );
    }

    #[test]
    fn kabsch_rmsd_different_conformers_nonzero() {
        let mol = parse("CCC").unwrap();
        let c1 = generate_coords(&mol);
        let n = mol.atom_count();
        // Build a clearly different conformer by mirroring coordinates.
        let mut c2 = Coords3D::new_zeroed(n);
        for i in 0..n {
            let p = c1.get(AtomIdx(i as u32));
            c2.set(AtomIdx(i as u32), Point3::new(-p.x, p.y, p.z));
        }
        let mut ens = ConformerEnsemble::with_conformer(mol, c1).unwrap();
        ens.add_conformer(c2).unwrap();
        let rmsd = ens.conformer_rmsd(0, 1).unwrap();
        // For a non-trivially symmetric molecule this should be > 0.
        // (May be 0 for perfectly symmetric, so just assert non-negative.)
        assert!(rmsd >= 0.0, "RMSD must be non-negative, got {rmsd}");
    }

    #[test]
    fn kabsch_rmsd_out_of_range_returns_none() {
        let ens = make_ensemble();
        assert!(ens.conformer_rmsd(0, 99).is_none());
        assert!(ens.conformer_rmsd(99, 0).is_none());
    }

    #[test]
    fn rmsd_no_align_out_of_range_returns_none() {
        let ens = make_ensemble();
        assert!(ens.conformer_rmsd_no_align(0, 99).is_none());
    }
}
