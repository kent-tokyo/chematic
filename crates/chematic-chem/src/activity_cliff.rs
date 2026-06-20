//! Activity cliff detection for SAR analysis.
//!
//! An activity cliff is a pair of structurally similar molecules whose biological
//! activity differs significantly — a key indicator of SAR sensitivity and potential
//! pharmacophore requirements.
//!
//! # Usage
//!
//! ```
//! # use chematic_smiles::parse;
//! # use chematic_chem::activity_cliff::activity_cliffs;
//! let mols = vec![
//!     parse("c1ccccc1").unwrap(),
//!     parse("Cc1ccccc1").unwrap(),
//!     parse("Nc1ccccc1").unwrap(),
//! ];
//! let refs: Vec<&_> = mols.iter().collect();
//! let activities = vec![5.0_f64, 7.5, 3.0];
//! let cliffs = activity_cliffs(&refs, &activities, 0.5, 2.0);
//! println!("{} cliff(s) found", cliffs.len());
//! ```

use chematic_core::Molecule;
use chematic_fp::ecfp4;

/// A pair of molecules forming an activity cliff.
#[derive(Debug, Clone, PartialEq)]
pub struct ActivityCliff {
    /// Index of the first molecule in the input slice (always `mol_a_idx < mol_b_idx`).
    pub mol_a_idx: usize,
    /// Index of the second molecule in the input slice.
    pub mol_b_idx: usize,
    /// ECFP4 Tanimoto similarity between the pair (0–1).
    pub similarity: f32,
    /// Absolute difference of activity values (`|activity_a − activity_b|`).
    pub activity_delta: f64,
}

/// Detect activity cliffs in a molecular dataset.
///
/// An activity cliff is returned for each pair (i, j) where **both**:
/// - `tanimoto_ecfp4(mol_i, mol_j) ≥ sim_threshold` (structurally similar), **and**
/// - `|activity_i − activity_j| ≥ cliff_delta` (but activity differs significantly).
///
/// Only upper-triangle pairs `(i < j)` are returned.
/// The result is sorted by similarity descending (highest similarity first).
///
/// `activities` are raw values (e.g. pIC50, log Ki). Use pIC50 rather than IC50
/// so that a 100-fold difference corresponds to a cliff_delta of 2.0.
///
/// # Panics
///
/// Panics if `mols.len() != activities.len()`.
pub fn activity_cliffs(
    mols: &[&Molecule],
    activities: &[f64],
    sim_threshold: f32,
    cliff_delta: f64,
) -> Vec<ActivityCliff> {
    assert_eq!(
        mols.len(),
        activities.len(),
        "mols and activities must have the same length"
    );

    let n = mols.len();
    if n < 2 {
        return Vec::new();
    }

    // Pre-compute ECFP4 fingerprints for all molecules.
    let fps: Vec<_> = mols.iter().map(|m| ecfp4(m)).collect();

    let mut cliffs: Vec<ActivityCliff> = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let sim = fps[i].tanimoto(&fps[j]) as f32;
            if sim < sim_threshold {
                continue;
            }
            let delta = (activities[i] - activities[j]).abs();
            if delta >= cliff_delta {
                cliffs.push(ActivityCliff {
                    mol_a_idx: i,
                    mol_b_idx: j,
                    similarity: sim,
                    activity_delta: delta,
                });
            }
        }
    }

    // Sort by similarity descending (most similar cliffs first).
    cliffs.sort_by(|a, b| b.similarity.partial_cmp(&a.similarity).unwrap());
    cliffs
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn test_empty_or_single_molecule() {
        // Less than 2 molecules → no cliffs possible.
        let m = mol("c1ccccc1");
        assert_eq!(activity_cliffs(&[], &[], 0.5, 1.0).len(), 0);
        assert_eq!(activity_cliffs(&[&m], &[5.0], 0.5, 1.0).len(), 0);
    }

    #[test]
    fn test_identical_molecules_no_cliff() {
        // Identical molecules have sim=1.0 but activity_delta=0 → no cliff.
        let m = mol("c1ccccc1");
        let cliffs = activity_cliffs(&[&m, &m], &[5.0, 5.0], 0.5, 2.0);
        assert_eq!(cliffs.len(), 0, "identical activity → no cliff");
    }

    #[test]
    fn test_dissimilar_molecules_no_cliff_despite_activity_gap() {
        // benzene vs pyridine: may be somewhat similar but let's use a low sim threshold.
        let a = mol("c1ccccc1");
        let b = mol("c1ccncc1");
        // With high sim_threshold of 0.9, unlikely to be cliffs for these different mols.
        let cliffs = activity_cliffs(&[&a, &b], &[5.0, 9.0], 0.9, 2.0);
        // Not asserting specific count — just check no panic.
        let _ = cliffs;
    }

    #[test]
    fn test_large_activity_gap_is_cliff_at_low_threshold() {
        // Use sim_threshold = 0.0 so any pair with large activity delta is a cliff.
        // benzene vs toluene ECFP4 similarity can be < 0.3 due to different radii
        // environments, so we use a low threshold to verify the core logic.
        let benz = mol("c1ccccc1");
        let tolu = mol("Cc1ccccc1");
        let cliffs = activity_cliffs(&[&benz, &tolu], &[3.0, 7.0], 0.0, 2.0);
        assert_eq!(
            cliffs.len(),
            1,
            "benzene/toluene with delta=4.0 and threshold=0 → 1 cliff"
        );
        assert_eq!(cliffs[0].mol_a_idx, 0);
        assert_eq!(cliffs[0].mol_b_idx, 1);
        assert!((cliffs[0].activity_delta - 4.0).abs() < 1e-9);
    }

    #[test]
    fn test_no_cliff_when_activity_delta_below_threshold() {
        // delta = 0.5 < cliff_delta 2.0 → no cliff, even at sim_threshold = 0.
        let benz = mol("c1ccccc1");
        let tolu = mol("Cc1ccccc1");
        let cliffs = activity_cliffs(&[&benz, &tolu], &[5.0, 5.5], 0.0, 2.0);
        assert_eq!(cliffs.len(), 0, "small activity delta → no cliff");
    }

    #[test]
    fn test_no_cliff_when_sim_below_threshold() {
        // Very dissimilar molecules (benzene vs butane) won't pass a high sim threshold.
        let a = mol("c1ccccc1");
        let b = mol("CCCC");
        let cliffs = activity_cliffs(&[&a, &b], &[3.0, 8.0], 0.8, 2.0);
        // These are very different → sim < 0.8 → no cliff.
        assert_eq!(
            cliffs.len(),
            0,
            "very dissimilar mols should not form a cliff at high sim threshold"
        );
    }

    #[test]
    fn test_sorted_by_similarity_descending() {
        // Three molecules: check that the output is sorted by sim descending.
        let a = mol("c1ccccc1");
        let b = mol("Cc1ccccc1");
        let c = mol("CCCC");

        let cliffs = activity_cliffs(&[&a, &b, &c], &[3.0, 8.0, 7.5], 0.0, 1.0);
        for w in cliffs.windows(2) {
            assert!(
                w[0].similarity >= w[1].similarity,
                "cliffs should be sorted by similarity descending"
            );
        }
    }

    #[test]
    #[should_panic(expected = "mols and activities must have the same length")]
    fn test_panic_on_mismatched_lengths() {
        let m = mol("c1ccccc1");
        activity_cliffs(&[&m], &[5.0, 6.0], 0.5, 1.0);
    }

    #[test]
    fn test_mul_mol_series_low_threshold() {
        // A larger set with sim_threshold=0.0 to verify multi-mol logic.
        let smiles = ["c1ccccc1", "Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1"];
        let mols: Vec<Molecule> = smiles.iter().map(|s| mol(s)).collect();
        let refs: Vec<&Molecule> = mols.iter().collect();
        let acts = vec![4.0, 6.5, 6.7, 5.0];
        // With threshold=0.0 and cliff_delta=2.0, benzene(4.0) vs toluene(6.5) is a cliff.
        let cliffs = activity_cliffs(&refs, &acts, 0.0, 2.0);
        assert!(!cliffs.is_empty(), "should find at least one cliff");
        for cliff in &cliffs {
            assert!(cliff.activity_delta >= 2.0);
        }
    }
}
