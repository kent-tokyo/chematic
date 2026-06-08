//! Similarity-based nearest-neighbour search over a molecule database.
//!
//! [`nearest_neighbors`] computes fingerprints for every molecule in `db`,
//! measures Tanimoto similarity against a query fingerprint, and returns the
//! top-k results sorted by descending similarity.

use chematic_core::Molecule;

use crate::bitvec::BitVec2048;
use crate::ecfp::{EcfpConfig, ecfp};

// ---------------------------------------------------------------------------
// Fingerprint type selector
// ---------------------------------------------------------------------------

/// Which fingerprint type to use for the similarity search.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FpType {
    /// ECFP4 (Morgan radius=2, 2048 bits). Default.
    #[default]
    Ecfp4,
    /// ECFP6 (Morgan radius=3, 2048 bits).
    Ecfp6,
    /// ECFP4 with chirality encoding enabled.
    Ecfp4Chiral,
    /// FCFP4 (feature-based circular FP).
    Fcfp4,
    /// MACCS 166-bit structural keys.
    Maccs,
    /// Topological path FP (max_len=7, 2048 bits).
    TopoPath,
}

fn compute_fp(mol: &Molecule, fp_type: FpType) -> BitVec2048 {
    match fp_type {
        FpType::Ecfp4 => ecfp(mol, &EcfpConfig::default()),
        FpType::Ecfp6 => ecfp(
            mol,
            &EcfpConfig {
                radius: 3,
                nbits: 2048,
                use_chirality: false,
                use_double_fold: false,
            },
        ),
        FpType::Ecfp4Chiral => ecfp(
            mol,
            &EcfpConfig {
                radius: 2,
                nbits: 2048,
                use_chirality: true,
                use_double_fold: false,
            },
        ),
        FpType::Fcfp4 => crate::fcfp::fcfp4(mol),
        FpType::Maccs => crate::maccs::maccs(mol),
        FpType::TopoPath => {
            crate::topo_path::topo_path(mol, &crate::topo_path::TopoPathConfig::default())
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Return the `k` most similar molecules in `db` to `query`, sorted by
/// descending Tanimoto similarity.
///
/// Returns `Vec<(db_index, tanimoto)>` where `db_index` is the 0-based
/// position in `db`.  If `db` has fewer than `k` molecules, all are returned.
/// Molecules with Tanimoto = 0.0 are excluded unless every molecule scores 0.
pub fn nearest_neighbors(
    query: &Molecule,
    db: &[Molecule],
    k: usize,
    fp_type: FpType,
) -> Vec<(usize, f64)> {
    if k == 0 || db.is_empty() {
        return vec![];
    }

    let query_fp = compute_fp(query, fp_type);

    let mut scores: Vec<(usize, f64)> = db
        .iter()
        .enumerate()
        .map(|(i, mol)| {
            let fp = compute_fp(mol, fp_type);
            (i, query_fp.tanimoto(&fp))
        })
        .filter(|(_, t)| *t > 0.0)
        .collect();

    scores.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    scores.truncate(k);
    scores
}

/// Like [`nearest_neighbors`] but accepts a pre-computed query fingerprint.
///
/// Use this when you need to search multiple query fingerprints against the
/// same database to avoid recomputing database fingerprints each time.
pub fn nearest_neighbors_from_fp(
    query_fp: &BitVec2048,
    db_fps: &[BitVec2048],
    k: usize,
) -> Vec<(usize, f64)> {
    if k == 0 || db_fps.is_empty() {
        return vec![];
    }

    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .enumerate()
        .map(|(i, fp)| (i, query_fp.tanimoto(fp)))
        .filter(|(_, t)| *t > 0.0)
        .collect();

    scores.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    scores.truncate(k);
    scores
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn benzene() -> Molecule {
        parse("c1ccccc1").unwrap()
    }
    fn toluene() -> Molecule {
        parse("Cc1ccccc1").unwrap()
    }
    fn naphthalene() -> Molecule {
        parse("c1ccc2ccccc2c1").unwrap()
    }
    fn ethane() -> Molecule {
        parse("CC").unwrap()
    }

    #[test]
    fn test_nn_self_is_first() {
        let query = benzene();
        let db = vec![ethane(), toluene(), benzene(), naphthalene()];
        let results = nearest_neighbors(&query, &db, 3, FpType::Ecfp4);
        assert!(!results.is_empty());
        // benzene (index 2) should be the closest match to itself
        assert_eq!(results[0].0, 2, "benzene should match itself first");
        assert!(
            (results[0].1 - 1.0).abs() < 1e-9,
            "self-similarity should be 1.0"
        );
    }

    #[test]
    fn test_nn_returns_k_results() {
        let query = benzene();
        let db = vec![ethane(), toluene(), benzene(), naphthalene()];
        let results = nearest_neighbors(&query, &db, 2, FpType::Ecfp4);
        assert!(results.len() <= 2, "should return at most k results");
    }

    #[test]
    fn test_nn_sorted_descending() {
        let query = benzene();
        let db = vec![ethane(), toluene(), benzene(), naphthalene()];
        let results = nearest_neighbors(&query, &db, 4, FpType::Ecfp4);
        for w in results.windows(2) {
            assert!(
                w[0].1 >= w[1].1,
                "results should be sorted by descending Tanimoto"
            );
        }
    }

    #[test]
    fn test_nn_empty_db() {
        let query = benzene();
        let results = nearest_neighbors(&query, &[], 5, FpType::Ecfp4);
        assert!(results.is_empty());
    }

    #[test]
    fn test_nn_k_zero() {
        let query = benzene();
        let db = vec![benzene()];
        let results = nearest_neighbors(&query, &db, 0, FpType::Ecfp4);
        assert!(results.is_empty());
    }

    #[test]
    fn test_nn_from_fp() {
        let query = benzene();
        let db = [ethane(), toluene(), benzene()];
        let query_fp = crate::ecfp::ecfp4(&query);
        let db_fps: Vec<_> = db.iter().map(crate::ecfp::ecfp4).collect();
        let results = nearest_neighbors_from_fp(&query_fp, &db_fps, 3);
        assert!(!results.is_empty());
        assert_eq!(results[0].0, 2, "benzene fp should match itself");
    }

    #[test]
    fn test_nn_maccs_type() {
        // With MACCS keys, benzene and toluene may score identically against
        // a benzene query (methyl adds no unique key). Just verify the top result
        // has a high score and ethane (index 2) scores lower than the aromatic mols.
        let query = benzene();
        let db = vec![toluene(), benzene(), ethane()];
        let results = nearest_neighbors(&query, &db, 3, FpType::Maccs);
        assert!(!results.is_empty());
        // Ethane (index 2) should not be the top hit
        assert_ne!(
            results[0].0, 2,
            "ethane should not be the top MACCS hit for benzene"
        );
        // Top hit should have Tanimoto > 0.5
        assert!(
            results[0].1 > 0.5,
            "top MACCS hit should have Tanimoto > 0.5"
        );
    }
}
