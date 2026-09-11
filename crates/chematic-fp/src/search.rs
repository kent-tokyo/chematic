//! Similarity-based nearest-neighbour search over a molecule database.
//!
//! [`nearest_neighbors`] computes fingerprints for every molecule in `db`,
//! measures Tanimoto similarity against a query fingerprint, and returns the
//! top-k results sorted by descending similarity.
//!
//! [`PreparedFingerprintIndex`] is the reusable path for repeated queries: it
//! computes the database fingerprints once and reuses them for every search.

use chematic_core::Molecule;

use crate::bitvec::BitVec2048;
use crate::ecfp::{EcfpConfig, ecfp};
use crate::rdkit_morgan_ecfp4::{RdkitMorganError, rdkit_morgan_ecfp4_experimental};

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
    /// RDKit Morgan radius-2, 2048-bit compatible ECFP4. This profile is
    /// fallible because RDKit-parity aromaticity preprocessing can fail.
    RdkitEcfp4,
}

fn compute_fp(mol: &Molecule, fp_type: FpType) -> Result<BitVec2048, RdkitMorganError> {
    let fingerprint = match fp_type {
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
        FpType::RdkitEcfp4 => return Ok(rdkit_morgan_ecfp4_experimental(mol)?.fingerprint),
    };
    Ok(fingerprint)
}

/// A fingerprint preparation or query failure. The index never falls back to
/// a different fingerprint profile after this error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreparedFingerprintError {
    Database {
        index: usize,
        source: RdkitMorganError,
    },
    Query {
        source: RdkitMorganError,
    },
}

impl std::fmt::Display for PreparedFingerprintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Database { index, source } => {
                write!(
                    f,
                    "fingerprint preparation failed at database index {index}: {source}"
                )
            }
            Self::Query { source } => write!(f, "fingerprint query failed: {source}"),
        }
    }
}

impl std::error::Error for PreparedFingerprintError {}

/// A reusable, in-memory fingerprint index for repeated nearest-neighbour
/// queries over the same molecule database.
///
/// The index deliberately stores only fingerprints, not molecules. This keeps
/// the hot query path independent of parsing and molecular graph traversal.
/// Returned indices always refer to the order of the slice passed to [`new`].
#[derive(Debug, Clone)]
pub struct PreparedFingerprintIndex {
    fp_type: FpType,
    fingerprints: Vec<BitVec2048>,
    popcounts: Vec<u32>,
}

impl PreparedFingerprintIndex {
    /// Build an index by computing one fingerprint for each database molecule.
    pub fn new(db: &[Molecule], fp_type: FpType) -> Self {
        Self::try_new(db, fp_type).expect("non-fallible fingerprint profile failed")
    }

    /// Build an index with explicit failure reporting for fallible profiles.
    pub fn try_new(db: &[Molecule], fp_type: FpType) -> Result<Self, PreparedFingerprintError> {
        let fingerprints: Vec<_> = db
            .iter()
            .enumerate()
            .map(|(index, mol)| {
                compute_fp(mol, fp_type)
                    .map_err(|source| PreparedFingerprintError::Database { index, source })
            })
            .collect::<Result<_, _>>()?;
        let popcounts = fingerprints.iter().map(BitVec2048::popcount).collect();
        Ok(Self {
            fp_type,
            fingerprints,
            popcounts,
        })
    }

    /// Fingerprint family used by this index.
    pub fn fp_type(&self) -> FpType {
        self.fp_type
    }

    /// Number of indexed molecules.
    pub fn len(&self) -> usize {
        self.fingerprints.len()
    }

    /// Whether the index contains no molecules.
    pub fn is_empty(&self) -> bool {
        self.fingerprints.is_empty()
    }

    /// Search the prepared database with a molecule query.
    pub fn search(&self, query: &Molecule, k: usize) -> Vec<(usize, f64)> {
        let query_fp =
            compute_fp(query, self.fp_type).expect("non-fallible fingerprint profile failed");
        self.search_fp(&query_fp, k)
    }

    /// Search using a fallible fingerprint profile without fallback.
    pub fn try_search(
        &self,
        query: &Molecule,
        k: usize,
    ) -> Result<Vec<(usize, f64)>, PreparedFingerprintError> {
        let query_fp = compute_fp(query, self.fp_type)
            .map_err(|source| PreparedFingerprintError::Query { source })?;
        Ok(self.search_fp(&query_fp, k))
    }

    /// Search the prepared database with an already computed fingerprint.
    pub fn search_fp(&self, query_fp: &BitVec2048, k: usize) -> Vec<(usize, f64)> {
        nearest_neighbors_from_prepared_fp(query_fp, &self.fingerprints, &self.popcounts, k)
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

    PreparedFingerprintIndex::new(db, fp_type).search(query, k)
}

/// Fallible nearest-neighbour search for profiles such as RDKit-compatible
/// ECFP4. No alternate fingerprint profile is used after a failure.
pub fn try_nearest_neighbors(
    query: &Molecule,
    db: &[Molecule],
    k: usize,
    fp_type: FpType,
) -> Result<Vec<(usize, f64)>, PreparedFingerprintError> {
    if k == 0 || db.is_empty() {
        return Ok(vec![]);
    }
    PreparedFingerprintIndex::try_new(db, fp_type)?.try_search(query, k)
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

    let query_popcount = query_fp.popcount();
    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .enumerate()
        .map(|(i, fp)| {
            (
                i,
                query_fp.tanimoto_with_counts_f64(fp, query_popcount, fp.popcount()),
            )
        })
        .filter(|(_, t)| *t > 0.0)
        .collect();

    rank_top_k(&mut scores, k);
    scores
}

fn nearest_neighbors_from_prepared_fp(
    query_fp: &BitVec2048,
    db_fps: &[BitVec2048],
    db_popcounts: &[u32],
    k: usize,
) -> Vec<(usize, f64)> {
    if k == 0 || db_fps.is_empty() {
        return vec![];
    }

    let query_popcount = query_fp.popcount();
    let mut scores: Vec<(usize, f64)> = db_fps
        .iter()
        .zip(db_popcounts.iter().copied())
        .enumerate()
        .map(|(i, (fp, popcount))| {
            (
                i,
                query_fp.tanimoto_with_counts_f64(fp, query_popcount, popcount),
            )
        })
        .filter(|(_, t)| *t > 0.0)
        .collect();
    rank_top_k(&mut scores, k);
    scores
}

/// Keep the exact top-k set while avoiding a full sort of the candidate list.
/// Both selection and final ordering use the public deterministic contract:
/// similarity descending, then database index ascending.
fn rank_top_k(scores: &mut Vec<(usize, f64)>, k: usize) {
    if scores.len() > k {
        scores.select_nth_unstable_by(k - 1, compare_rank);
        scores.truncate(k);
    }
    scores.sort_unstable_by(compare_rank);
}

fn compare_rank(a: &(usize, f64), b: &(usize, f64)) -> std::cmp::Ordering {
    b.1.partial_cmp(&a.1)
        .unwrap_or(std::cmp::Ordering::Equal)
        .then_with(|| a.0.cmp(&b.0))
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
    fn prepared_index_matches_one_shot_search() {
        let query = benzene();
        let db = vec![ethane(), toluene(), benzene(), naphthalene()];
        let expected = nearest_neighbors(&query, &db, 3, FpType::Ecfp4);
        let index = PreparedFingerprintIndex::new(&db, FpType::Ecfp4);
        assert_eq!(index.len(), db.len());
        assert!(!index.is_empty());
        assert_eq!(index.fp_type(), FpType::Ecfp4);
        assert_eq!(index.search(&query, 3), expected);
    }

    #[test]
    fn prepared_index_handles_empty_and_zero_k() {
        let query = benzene();
        let empty = PreparedFingerprintIndex::new(&[], FpType::Ecfp4);
        assert!(empty.is_empty());
        assert!(empty.search(&query, 5).is_empty());

        let db = vec![benzene()];
        let index = PreparedFingerprintIndex::new(&db, FpType::Ecfp4);
        assert!(index.search(&query, 0).is_empty());
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

    #[test]
    fn top_k_ties_are_ordered_by_database_index() {
        let query = benzene();
        let db = vec![ethane(), benzene(), benzene(), toluene()];
        let results = nearest_neighbors(&query, &db, 3, FpType::Ecfp4);
        assert_eq!(results[0].0, 1);
        assert_eq!(results[1].0, 2);
    }

    #[test]
    fn prepared_index_supports_fallible_rdkit_ecfp4_profile() {
        let query = benzene();
        let db = vec![ethane(), benzene(), toluene()];
        let index = PreparedFingerprintIndex::try_new(&db, FpType::RdkitEcfp4).unwrap();
        assert_eq!(index.fp_type(), FpType::RdkitEcfp4);
        let results = index.try_search(&query, 3).unwrap();
        assert_eq!(results[0].0, 1);
        assert!((results[0].1 - 1.0).abs() < 1e-9);
    }
}
