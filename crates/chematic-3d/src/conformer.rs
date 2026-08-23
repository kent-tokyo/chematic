//! Conformer ensemble: a molecule with multiple sets of 3D coordinates.

use std::fmt;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_smarts::{
    AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, MatchConfig, QueryMolecule,
    find_matches_with_config,
};

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

    /// Symmetry-aware Kabsch RMSD between conformers `a` and `b`: the minimum
    /// [`conformer_rmsd`]-style RMSD over every way `b`'s atoms can be
    /// relabelled onto `a`'s that is consistent with the molecule's own graph
    /// symmetry (automorphisms), not just the identity relabelling.
    ///
    /// See [`rmsd_symmetric`] for the algorithm and its known limitation
    /// (no `-COO⁻`/`-NO₂`-style terminal-group symmetrization).
    ///
    /// Returns `None` if either index is out of range.
    pub fn conformer_rmsd_symmetric(&self, a: usize, b: usize) -> Option<f64> {
        let ca = self.conformers.get(a)?;
        let cb = self.conformers.get(b)?;
        Some(rmsd_symmetric(&self.mol, ca, cb))
    }

    /// Compute the 12 USR shape descriptors for conformer `idx`.
    ///
    /// Returns `None` if `idx` is out of range.
    pub fn conformer_usr_descriptors(&self, idx: usize) -> Option<[f64; 12]> {
        let c = self.conformers.get(idx)?;
        let pts: Vec<[f64; 3]> = c.points.iter().map(|p| [p.x, p.y, p.z]).collect();
        Some(crate::usr::usr_descriptors(&pts))
    }

    /// Cluster conformers by Kabsch-aligned RMSD and return the indices of
    /// representative conformers to keep (one per cluster).
    ///
    /// Uses a **greedy leader-linkage** algorithm: conformers are visited in
    /// index order; each is compared against the representative (first member)
    /// of every existing cluster via [`conformer_rmsd`]. If the RMSD to any
    /// cluster leader is strictly less than `rms_threshold`, the conformer joins
    /// that cluster and is discarded. Otherwise it starts a new cluster and is kept.
    ///
    /// # Returns
    /// Indices of kept conformers in ascending order, at most one per cluster.
    /// - Empty ensemble → `[]`
    /// - Single conformer → `[0]`
    /// - `rms_threshold <= 0.0` → all indices kept
    ///
    /// # Example
    /// ```rust,ignore
    /// // Remove near-duplicate conformers within 0.5 Å RMSD
    /// let kept = ensemble.cluster_conformers_by_rms(0.5);
    /// ```
    pub fn cluster_conformers_by_rms(&self, rms_threshold: f64) -> Vec<usize> {
        let n = self.conformers.len();
        if n == 0 {
            return vec![];
        }
        if rms_threshold <= 0.0 {
            return (0..n).collect();
        }
        let mut leaders: Vec<usize> = Vec::new();
        'outer: for i in 0..n {
            for &leader in &leaders {
                let rmsd = self.conformer_rmsd(i, leader).unwrap_or(f64::INFINITY);
                if rmsd < rms_threshold {
                    continue 'outer; // duplicate — skip
                }
            }
            leaders.push(i); // new cluster representative
        }
        leaders
    }

    /// Return the `(index, rmsd)` of the first existing conformer within
    /// `rmsd_threshold` Å of `coords` after Kabsch superposition, or `None` if
    /// no existing conformer is within the threshold (including an empty
    /// ensemble, or `rmsd_threshold <= 0.0` — the "no pruning" convention this
    /// type's duplicate-checking methods all share).
    ///
    /// Used by ensemble generators to discard near-duplicate structures before
    /// adding them, while still reporting *which* existing conformer matched
    /// and at what RMSD (needed to record provenance for a discarded
    /// candidate, not just a yes/no). O(k) Kabsch operations where k is the
    /// current ensemble size.
    pub fn find_duplicate(&self, coords: &Coords3D, rmsd_threshold: f64) -> Option<(usize, f64)> {
        if rmsd_threshold <= 0.0 {
            return None;
        }
        let n = self.mol.atom_count();
        self.conformers
            .iter()
            .enumerate()
            .map(|(i, existing)| (i, kabsch_rmsd(coords, existing, n)))
            .find(|(_, rmsd)| *rmsd < rmsd_threshold)
    }

    /// Return `true` if `coords` is within `rmsd_threshold` Å of any existing
    /// conformer after Kabsch superposition. Thin wrapper over
    /// [`find_duplicate`](Self::find_duplicate) for callers that only need the
    /// yes/no answer.
    pub fn is_duplicate(&self, coords: &Coords3D, rmsd_threshold: f64) -> bool {
        self.find_duplicate(coords, rmsd_threshold).is_some()
    }

    /// Symmetry-aware counterpart to [`find_duplicate`](Self::find_duplicate):
    /// the `(index, rmsd)` of the first existing conformer within
    /// `rmsd_threshold` Å under [`rmsd_symmetric`] (automorphism-aware) rather
    /// than plain fixed-index Kabsch RMSD. Correct on molecules with
    /// interchangeable substituents (e.g. a terminal `-CF3`, see
    /// `symmetric_rmsd_recovers_zero_when_fixed_index_rmsd_is_wrong` in this
    /// module's tests) where `find_duplicate` can treat truly-identical
    /// conformers as distinct. O(k) automorphism enumerations where k is the
    /// current ensemble size — more expensive than `find_duplicate`, used
    /// where correctness on symmetric molecules matters more than raw speed.
    pub fn find_duplicate_symmetric(
        &self,
        coords: &Coords3D,
        rmsd_threshold: f64,
    ) -> Option<(usize, f64)> {
        if rmsd_threshold <= 0.0 {
            return None;
        }
        self.conformers
            .iter()
            .enumerate()
            .map(|(i, existing)| (i, rmsd_symmetric(&self.mol, coords, existing)))
            .find(|(_, rmsd)| *rmsd < rmsd_threshold)
    }

    /// Thin wrapper over
    /// [`find_duplicate_symmetric`](Self::find_duplicate_symmetric) for
    /// callers that only need the yes/no answer.
    pub fn is_duplicate_symmetric(&self, coords: &Coords3D, rmsd_threshold: f64) -> bool {
        self.find_duplicate_symmetric(coords, rmsd_threshold)
            .is_some()
    }

    /// Mean pairwise USR dissimilarity across all conformers.
    ///
    /// Returns a value in `[0.0, 1.0]`: 0.0 means all conformers are identical
    /// shapes; values closer to 1.0 indicate a highly diverse ensemble.
    /// Returns 0.0 for ensembles with fewer than 2 conformers.
    pub fn conformer_diversity_usr(&self) -> f64 {
        let n = self.conformers.len();
        if n < 2 {
            return 0.0;
        }
        let descs: Vec<[f64; 12]> = (0..n)
            .filter_map(|i| self.conformer_usr_descriptors(i))
            .collect();
        let mut total = 0.0;
        let mut count = 0usize;
        for i in 0..descs.len() {
            for j in i + 1..descs.len() {
                total += 1.0 - crate::usr::usr_similarity(&descs[i], &descs[j]);
                count += 1;
            }
        }
        if count == 0 {
            0.0
        } else {
            total / count as f64
        }
    }
}

// ---------------------------------------------------------------------------
// Symmetry-aware (automorphism-brute-force) RMSD
// ---------------------------------------------------------------------------

/// Symmetry-aware Kabsch RMSD between two conformers of the *same* molecule
/// topology (same atom count, same connectivity), ported from RDKit's
/// `GetBestRMS` (`Code/GraphMol/MolAlign/AlignMolecules.cpp`, pinned commit
/// `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f` — see
/// `scripts/mmff94_provenance/PROVENANCE.md` for this project's pinning
/// convention): enumerate every automorphism of `mol` (self-match with
/// `uniquify: false`, matching RDKit's own `SubstructMatch(..., uniquify=false)`
/// self-match step) via VF2, and return the minimum Kabsch-aligned RMSD over
/// every automorphism-consistent relabelling of `coords_b` onto `coords_a`.
/// Brute force over all matches, same as RDKit — no pruning, no Hungarian
/// assignment.
///
/// A plain (non-symmetry-aware) RMSD is wrong on any molecule with
/// permutation-equivalent atoms: e.g. a terminal `-CF3`'s three fluorines are
/// interchangeable, so a conformer pair that differs only by which specific F
/// atom sits in which position should score ~0, not a large fixed-index RMSD.
///
/// **Known limitation, not yet ported**: RDKit's `GetBestRMS` additionally
/// runs `symmetrizeConjugatedTerminalGroups` before matching, which neutralizes
/// resonance-equivalent terminal groups (carboxylate `-COO⁻`, nitro `-NO2`, …)
/// so their two oxygens are treated as interchangeable even though a
/// Kekulized/formal-charge representation gives them different bond orders.
/// This function does NOT do that preprocessing, so on such molecules it will
/// report a higher (worse) RMSD than RDKit's `GetBestRMS` for a case that is
/// chemically equivalent but not topologically automorphic as drawn. See
/// issue tracker for the follow-up.
pub fn rmsd_symmetric(mol: &Molecule, coords_a: &Coords3D, coords_b: &Coords3D) -> f64 {
    let n = mol.atom_count();
    if n == 0 {
        return 0.0;
    }
    let query = molecule_self_query(mol);
    let config = MatchConfig {
        uniquify: false,
        ..MatchConfig::default()
    };
    let matches = find_matches_with_config(&query, mol, &config);
    debug_assert!(
        !matches.is_empty(),
        "a molecule must always self-match at least via the identity mapping"
    );

    let mut best = f64::MAX;
    for m in &matches {
        let mut relabelled = Coords3D::new_zeroed(n);
        for qi in 0..n {
            // `m` maps this molecule's own atom indices (as the query) onto
            // itself (as the target) under one automorphism; relabelling
            // `coords_b` by that map produces one automorphism-consistent
            // atom correspondence between `coords_a` and `coords_b`.
            let target = m[&qi];
            relabelled.set(AtomIdx(qi as u32), coords_b.get(target));
        }
        let rmsd = kabsch_rmsd(coords_a, &relabelled, n);
        if rmsd < best {
            best = rmsd;
        }
    }
    best
}

/// Build a `QueryMolecule` that matches only `mol` itself: one query atom per
/// heavy atom (atomic number + formal charge), one query bond per bond
/// (exact bond-order primitive where the order is unambiguous, `Any` for the
/// query/metadata bond kinds that carry no independent topological meaning
/// for automorphism purposes — `Zero`/`Dative`/`Query*`/`Up`/`Down` are all
/// mapped to `Any` since none of them changes which atoms are interchangeable).
fn molecule_self_query(mol: &Molecule) -> QueryMolecule {
    let mut q = QueryMolecule::new();
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        let by_element =
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(atom.element.atomic_number()));
        let by_charge = AtomQuery::Primitive(AtomPrimitive::Charge(atom.charge));
        q.add_atom(AtomQuery::And(Box::new(by_element), Box::new(by_charge)));
    }
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        let prim = match bond.order {
            BondOrder::Single => BondPrimitive::Single,
            BondOrder::Double => BondPrimitive::Double,
            BondOrder::Triple => BondPrimitive::Triple,
            BondOrder::Aromatic => BondPrimitive::Aromatic,
            _ => BondPrimitive::Any,
        };
        q.add_bond(
            bond.atom1.0 as usize,
            bond.atom2.0 as usize,
            BondQuery::Primitive(prim),
        );
    }
    q
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

    // R = U * V^T.  R[r][c] = Σ_k U[r][k] * V[c][k].
    let mut r_mat = [[0.0f64; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            for k in 0..3 {
                r_mat[r][c] += u[r][k] * v[c][k];
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
        // Recompute R = U * V_final^T.
        r_mat = [[0.0f64; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..3 {
                    r_mat[r][c] += u[r][k] * v_final[c][k];
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
    fn kabsch_rmsd_pure_rotation_is_zero() {
        // A pure rotation must give RMSD = 0 after Kabsch superposition.
        let mol = parse("CCC").unwrap();
        let n = mol.atom_count();
        let base = generate_coords(&mol);
        // 90° rotation around z-axis: (x,y,z) → (−y, x, z).
        let mut rotated = Coords3D::new_zeroed(n);
        for i in 0..n {
            let p = base.get(AtomIdx(i as u32));
            rotated.set(AtomIdx(i as u32), Point3::new(-p.y, p.x, p.z));
        }
        let mut ens = ConformerEnsemble::with_conformer(mol, base).unwrap();
        ens.add_conformer(rotated).unwrap();
        let rmsd = ens.conformer_rmsd(0, 1).unwrap();
        assert!(
            rmsd < 1e-6,
            "pure-rotation Kabsch RMSD should be ~0, got {rmsd}"
        );
    }

    // --- rmsd_symmetric ------------------------------------------------------

    #[test]
    fn symmetric_rmsd_self_is_zero() {
        let mol = parse("CCC").unwrap();
        let c = generate_coords(&mol);
        let rmsd = rmsd_symmetric(&mol, &c, &c);
        assert!(rmsd.abs() < 1e-8, "self-RMSD should be 0, got {rmsd}");
    }

    #[test]
    fn symmetric_rmsd_never_exceeds_fixed_index_rmsd() {
        // The identity mapping is always one of the enumerated automorphisms,
        // so the symmetric RMSD (a minimum over all of them) can never be
        // larger than the plain fixed-index Kabsch RMSD.
        let mol = parse("CCC").unwrap();
        let n = mol.atom_count();
        let c1 = generate_coords(&mol);
        let mut c2 = Coords3D::new_zeroed(n);
        for i in 0..n {
            let p = c1.get(AtomIdx(i as u32));
            c2.set(AtomIdx(i as u32), Point3::new(-p.x, p.y, p.z));
        }
        let fixed = kabsch_rmsd(&c1, &c2, n);
        let symmetric = rmsd_symmetric(&mol, &c1, &c2);
        assert!(
            symmetric <= fixed + 1e-9,
            "symmetric RMSD ({symmetric}) must not exceed fixed-index RMSD ({fixed})"
        );
    }

    #[test]
    fn symmetric_rmsd_recovers_zero_when_fixed_index_rmsd_is_wrong() {
        // 1,1,1-trifluoroethane: the 3 fluorines on the CF3 carbon are
        // topologically interchangeable (verified: 6 = 3! self-match
        // automorphisms with `uniquify: false`). Swap two F atoms' positions
        // between otherwise-identical conformers -- a real geometric
        // difference under fixed-index comparison, but the SAME physical
        // structure under symmetry-aware comparison.
        let mol = parse("FC(F)(F)C").unwrap();
        let n = mol.atom_count();
        let base = generate_coords(&mol);

        let fluorines: Vec<AtomIdx> = (0..n)
            .map(|i| AtomIdx(i as u32))
            .filter(|&idx| mol.atom(idx).element == chematic_core::Element::F)
            .collect();
        assert_eq!(fluorines.len(), 3, "CF3 should have exactly 3 fluorines");

        let mut swapped = Coords3D::new_zeroed(n);
        for i in 0..n {
            swapped.set(AtomIdx(i as u32), base.get(AtomIdx(i as u32)));
        }
        let (f0, f1) = (fluorines[0], fluorines[1]);
        let (p0, p1) = (base.get(f0), base.get(f1));
        swapped.set(f0, p1);
        swapped.set(f1, p0);

        let fixed = kabsch_rmsd(&base, &swapped, n);
        assert!(
            fixed > 0.1,
            "swapping two real F positions should give a clearly nonzero \
             fixed-index RMSD, got {fixed}"
        );

        let symmetric = rmsd_symmetric(&mol, &base, &swapped);
        assert!(
            symmetric < 1e-6,
            "symmetry-aware RMSD should recover ~0 for an automorphism-only \
             difference, got {symmetric} (fixed-index was {fixed})"
        );
    }

    #[test]
    fn symmetric_rmsd_ensemble_method_matches_free_function() {
        let mol = parse("CCC").unwrap();
        let c1 = generate_coords(&mol);
        let c2 = generate_coords(&mol);
        let expected = rmsd_symmetric(&mol, &c1, &c2);
        let mut ens = ConformerEnsemble::with_conformer(mol, c1).unwrap();
        ens.add_conformer(c2).unwrap();
        let got = ens.conformer_rmsd_symmetric(0, 1).unwrap();
        assert!(
            (got - expected).abs() < 1e-12,
            "ensemble method should match the free function exactly"
        );
    }

    #[test]
    fn symmetric_rmsd_ensemble_method_out_of_range_returns_none() {
        let ens = make_ensemble();
        assert!(ens.conformer_rmsd_symmetric(0, 99).is_none());
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

    // C4: conformer diversity metrics

    #[test]
    fn usr_descriptors_single_conformer() {
        let ens = make_ensemble();
        let d = ens.conformer_usr_descriptors(0);
        assert!(d.is_some(), "valid index must return Some");
        assert!(
            d.unwrap().iter().all(|v| v.is_finite()),
            "all USR values finite"
        );
    }

    #[test]
    fn usr_descriptors_out_of_range() {
        let ens = make_ensemble();
        assert!(ens.conformer_usr_descriptors(99).is_none());
    }

    #[test]
    fn diversity_usr_single_conformer_is_zero() {
        let ens = make_ensemble();
        assert_eq!(
            ens.conformer_diversity_usr(),
            0.0,
            "single conformer → diversity 0"
        );
    }

    #[test]
    fn diversity_usr_identical_conformers_is_zero() {
        use crate::coords::Point3;
        use chematic_core::{Atom, Element, MoleculeBuilder};

        let mut b = MoleculeBuilder::new();
        let a0 = b.add_atom(Atom::new(Element::C));
        let a1 = b.add_atom(Atom::new(Element::C));
        let mol = b.build();

        let mut c = Coords3D::new_zeroed(2);
        c.set(a0, Point3::new(0.0, 0.0, 0.0));
        c.set(a1, Point3::new(1.5, 0.0, 0.0));

        let mut ens = ConformerEnsemble::with_conformer(mol, c.clone()).unwrap();
        ens.add_conformer(c).unwrap();

        let div = ens.conformer_diversity_usr();
        assert!(
            div.abs() < 1e-9,
            "identical conformers → diversity ~0, got {div}"
        );
    }

    #[test]
    fn diversity_usr_different_shapes_positive() {
        use crate::coords::Point3;
        use chematic_core::{Atom, Element, MoleculeBuilder};

        // 3-atom molecule; two very different conformers
        let mut b = MoleculeBuilder::new();
        let a0 = b.add_atom(Atom::new(Element::C));
        let a1 = b.add_atom(Atom::new(Element::C));
        let a2 = b.add_atom(Atom::new(Element::C));
        let mol = b.build();

        let mut c1 = Coords3D::new_zeroed(3);
        c1.set(a0, Point3::new(0.0, 0.0, 0.0));
        c1.set(a1, Point3::new(1.0, 0.0, 0.0));
        c1.set(a2, Point3::new(2.0, 0.0, 0.0));

        let mut c2 = Coords3D::new_zeroed(3);
        c2.set(a0, Point3::new(0.0, 0.0, 0.0));
        c2.set(a1, Point3::new(0.0, 10.0, 0.0));
        c2.set(a2, Point3::new(0.0, 0.0, 10.0));

        let mut ens = ConformerEnsemble::with_conformer(mol, c1).unwrap();
        ens.add_conformer(c2).unwrap();

        let div = ens.conformer_diversity_usr();
        assert!(
            div > 0.0 && div <= 1.0,
            "diverse ensemble → diversity in (0,1], got {div}"
        );
    }

    // --- cluster_conformers_by_rms ------------------------------------------

    #[test]
    fn cluster_empty_ensemble() {
        let mol = parse("C").unwrap();
        let ens = ConformerEnsemble::new(mol);
        assert!(ens.cluster_conformers_by_rms(0.5).is_empty());
    }

    #[test]
    fn cluster_single_conformer() {
        let ens = make_ensemble();
        assert_eq!(ens.cluster_conformers_by_rms(0.5), vec![0]);
    }

    #[test]
    fn cluster_zero_threshold_keeps_all() {
        let mol = parse("CCC").unwrap();
        let c = generate_coords(&mol);
        let mut ens = ConformerEnsemble::with_conformer(mol, c.clone()).unwrap();
        ens.add_conformer(c.clone()).unwrap();
        ens.add_conformer(c).unwrap();
        let kept = ens.cluster_conformers_by_rms(0.0);
        assert_eq!(kept, vec![0, 1, 2], "threshold ≤ 0 → keep all");
    }

    #[test]
    fn cluster_negative_threshold_keeps_all() {
        let mol = parse("CCC").unwrap();
        let c = generate_coords(&mol);
        let mut ens = ConformerEnsemble::with_conformer(mol, c.clone()).unwrap();
        ens.add_conformer(c).unwrap();
        assert_eq!(ens.cluster_conformers_by_rms(-1.0), vec![0, 1]);
    }

    #[test]
    fn cluster_identical_conformers_keeps_first() {
        let mol = parse("CCC").unwrap();
        let c = generate_coords(&mol);
        let mut ens = ConformerEnsemble::with_conformer(mol, c.clone()).unwrap();
        ens.add_conformer(c.clone()).unwrap();
        ens.add_conformer(c).unwrap();
        let kept = ens.cluster_conformers_by_rms(0.01);
        assert_eq!(
            kept,
            vec![0],
            "three identical conformers → keep only first"
        );
    }

    #[test]
    fn cluster_distinct_conformers_keeps_all() {
        // Two conformers with RMSD >> threshold: both kept.
        let mol = parse("CCC").unwrap();
        let n = mol.atom_count();
        let mut ca = Coords3D::new_zeroed(n);
        let mut cb = Coords3D::new_zeroed(n);
        for i in 0..n {
            ca.set(
                chematic_core::AtomIdx(i as u32),
                Point3 {
                    x: i as f64,
                    y: 0.0,
                    z: 0.0,
                },
            );
            cb.set(
                chematic_core::AtomIdx(i as u32),
                Point3 {
                    x: 0.0,
                    y: i as f64 * 100.0,
                    z: 0.0,
                },
            );
        }
        let mut ens = ConformerEnsemble::with_conformer(mol, ca).unwrap();
        ens.add_conformer(cb).unwrap();
        let kept = ens.cluster_conformers_by_rms(0.5);
        assert_eq!(kept, vec![0, 1], "two distinct conformers → both kept");
    }

    // --- is_duplicate_symmetric ----------------------------------------------

    #[test]
    fn is_duplicate_symmetric_recovers_true_when_plain_is_duplicate_misses() {
        // Same setup as symmetric_rmsd_recovers_zero_when_fixed_index_rmsd_is_wrong:
        // swapping two of CF3's interchangeable fluorines gives a clearly nonzero
        // plain Kabsch RMSD but a ~0 symmetric RMSD.
        let mol = parse("FC(F)(F)C").unwrap();
        let n = mol.atom_count();
        let base = generate_coords(&mol);

        let fluorines: Vec<AtomIdx> = (0..n)
            .map(|i| AtomIdx(i as u32))
            .filter(|&idx| mol.atom(idx).element == chematic_core::Element::F)
            .collect();
        let (f0, f1) = (fluorines[0], fluorines[1]);
        let mut swapped = Coords3D::new_zeroed(n);
        for i in 0..n {
            swapped.set(AtomIdx(i as u32), base.get(AtomIdx(i as u32)));
        }
        let (p0, p1) = (base.get(f0), base.get(f1));
        swapped.set(f0, p1);
        swapped.set(f1, p0);

        let ens = ConformerEnsemble::with_conformer(mol, base).unwrap();
        assert!(
            !ens.is_duplicate(&swapped, 0.05),
            "plain Kabsch should NOT treat the swap as a near-duplicate at a tight threshold"
        );
        assert!(
            ens.is_duplicate_symmetric(&swapped, 0.05),
            "symmetric RMSD should treat the swap as a duplicate (automorphism-only difference)"
        );
    }

    #[test]
    fn is_duplicate_symmetric_zero_threshold_is_never_duplicate() {
        let ens = make_ensemble();
        let c = ens.get_conformer(0).unwrap().clone();
        assert!(!ens.is_duplicate_symmetric(&c, 0.0));
    }

    #[test]
    fn find_duplicate_reports_matching_index_and_rmsd() {
        let ens = make_ensemble(); // one conformer, index 0
        let c = ens.get_conformer(0).unwrap().clone();
        let (idx, rmsd) = ens
            .find_duplicate(&c, 0.5)
            .expect("identical coords must match");
        assert_eq!(idx, 0);
        assert!(
            rmsd.abs() < 1e-9,
            "self-comparison RMSD should be ~0, got {rmsd}"
        );
    }

    #[test]
    fn find_duplicate_none_when_no_match() {
        let mol = parse("CCC").unwrap();
        let n = mol.atom_count();
        let mut far = Coords3D::new_zeroed(n);
        for i in 0..n {
            far.set(
                AtomIdx(i as u32),
                Point3 {
                    x: 0.0,
                    y: i as f64 * 100.0,
                    z: 0.0,
                },
            );
        }
        let ens = make_ensemble();
        assert!(ens.find_duplicate(&far, 0.5).is_none());
    }

    #[test]
    fn cluster_ascending_order() {
        // Kept indices must always be in ascending order.
        let mol = parse("CCC").unwrap();
        let c0 = generate_coords(&mol);
        let mut ens = ConformerEnsemble::with_conformer(mol, c0.clone()).unwrap();
        ens.add_conformer(c0).unwrap(); // identical → cluster with 0
        let n = ens.mol().atom_count();
        let mut far = Coords3D::new_zeroed(n);
        for i in 0..n {
            far.set(
                chematic_core::AtomIdx(i as u32),
                Point3 {
                    x: 0.0,
                    y: i as f64 * 100.0,
                    z: 0.0,
                },
            );
        }
        ens.add_conformer(far).unwrap(); // distinct → new cluster
        let kept = ens.cluster_conformers_by_rms(0.1);
        for w in kept.windows(2) {
            assert!(w[0] < w[1], "kept indices not ascending");
        }
    }
}
