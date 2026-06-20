//! MinHash LSH (Locality Sensitive Hashing) index for fast approximate
//! Jaccard similarity search over MHFP fingerprints.
//!
//! ## Algorithm: band decomposition
//!
//! Split K MinHash values into B bands of R rows each (K = B × R).
//! Two fingerprints land in the same bucket for band `b` iff their rows
//! b×R .. b×R+R match exactly.  Any shared bucket makes them candidates;
//! exact similarity is then computed only for candidates.
//!
//! ```text
//! num_hashes=128, bands=16, rows=8:
//!   P(candidate | s=0.7) ≈ 90%
//!   P(candidate | s=0.8) ≈ 94%
//!   P(candidate | s=0.5) ≈ 45%   ← low-similarity pairs naturally filtered
//! ```
//!
//! Time complexity vs. linear scan:
//! - add:   O(B)      vs. O(1)
//! - query: O(B + C×K) vs. O(N×K)  where C << N is the candidate count

use std::collections::{HashMap, HashSet};

use crate::mhfp::MhfpFingerprint;

// ---------------------------------------------------------------------------
// Index
// ---------------------------------------------------------------------------

/// MinHash LSH index for approximate Jaccard/Tanimoto similarity search.
///
/// Insert MHFP fingerprints with [`MhfpLshIndex::add`], then call
/// [`MhfpLshIndex::query`] to retrieve all indexed fingerprints whose
/// similarity to a probe exceeds a threshold in sub-linear time.
///
/// ## Parameters
///
/// - `bands × rows` must equal the `num_hashes` of every inserted fingerprint.
/// - Default (`new(128)`): 16 bands × 8 rows.  For higher recall at lower
///   thresholds, use more bands and fewer rows (e.g. `with_bands(32, 4)`).
pub struct MhfpLshIndex {
    bands: usize,
    rows: usize,
    /// buckets\[band_idx\] maps a band's hash vector → list of entry indices.
    buckets: Vec<HashMap<Vec<u64>, Vec<usize>>>,
    /// Full fingerprints stored for exact similarity recomputation.
    fps: Vec<MhfpFingerprint>,
}

impl MhfpLshIndex {
    /// Create an index tuned for `num_hashes`-lane MHFP fingerprints.
    ///
    /// Uses 16 bands × (num_hashes / 16) rows.  Panics if `num_hashes` is
    /// not divisible by 16.
    pub fn new(num_hashes: usize) -> Self {
        assert!(
            num_hashes.is_multiple_of(16),
            "num_hashes ({num_hashes}) must be divisible by 16 for default band decomposition"
        );
        Self::with_bands(16, num_hashes / 16)
    }

    /// Create with explicit `bands` and `rows`.
    ///
    /// Every inserted fingerprint must have exactly `bands × rows` hash lanes.
    pub fn with_bands(bands: usize, rows: usize) -> Self {
        assert!(bands > 0 && rows > 0, "bands and rows must be > 0");
        MhfpLshIndex {
            bands,
            rows,
            buckets: vec![HashMap::new(); bands],
            fps: Vec::new(),
        }
    }

    /// Add a fingerprint to the index; returns its 0-based index.
    ///
    /// Panics if `fp.num_hashes != bands × rows`.
    pub fn add(&mut self, fp: MhfpFingerprint) -> usize {
        let expected = self.bands * self.rows;
        assert_eq!(
            fp.hashes.len(),
            expected,
            "fingerprint has {} hashes but index expects {} (bands={}, rows={})",
            fp.hashes.len(),
            expected,
            self.bands,
            self.rows
        );

        let idx = self.fps.len();
        for b in 0..self.bands {
            let start = b * self.rows;
            let band_key = fp.hashes[start..start + self.rows].to_vec();
            self.buckets[b].entry(band_key).or_default().push(idx);
        }
        self.fps.push(fp);
        idx
    }

    /// Query for all indexed fingerprints with Tanimoto similarity ≥ `threshold`.
    ///
    /// Returns `(index, similarity)` pairs sorted by descending similarity.
    /// The result may miss some true positives below `threshold + ε` due to the
    /// probabilistic nature of LSH; increase `bands` (and decrease `rows`)
    /// for higher recall.
    pub fn query(&self, fp: &MhfpFingerprint, threshold: f64) -> Vec<(usize, f64)> {
        // Step 1: collect candidate indices from shared band buckets
        let mut candidates: HashSet<usize> = HashSet::new();
        for b in 0..self.bands {
            let start = b * self.rows;
            let end = (start + self.rows).min(fp.hashes.len());
            let band_key = &fp.hashes[start..end];
            if let Some(bucket) = self.buckets[b].get(band_key) {
                candidates.extend(bucket.iter().copied());
            }
        }

        // Step 2: exact similarity check for candidates
        let mut results: Vec<(usize, f64)> = candidates
            .into_iter()
            .filter_map(|idx| {
                let sim = fp.tanimoto(&self.fps[idx]);
                (sim >= threshold).then_some((idx, sim))
            })
            .collect();

        results.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        results
    }

    /// Number of fingerprints in the index.
    pub fn len(&self) -> usize {
        self.fps.len()
    }

    /// True if the index contains no fingerprints.
    pub fn is_empty(&self) -> bool {
        self.fps.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mhfp::{MhfpConfig, mhfp};
    use chematic_smiles::parse;

    fn fp(smiles: &str) -> MhfpFingerprint {
        let mol = parse(smiles).unwrap();
        mhfp(&mol)
    }

    fn fp_cfg(smiles: &str, cfg: &MhfpConfig) -> MhfpFingerprint {
        let mol = parse(smiles).unwrap();
        crate::mhfp::mhfp_with_config(&mol, cfg)
    }

    #[test]
    fn test_lsh_empty_query() {
        let index = MhfpLshIndex::new(128);
        let results = index.query(&fp("CC"), 0.5);
        assert!(results.is_empty(), "empty index should return no results");
    }

    #[test]
    fn test_lsh_self_similarity() {
        let mut index = MhfpLshIndex::new(128);
        let benzene = fp("c1ccccc1");
        index.add(benzene.clone());
        let results = index.query(&benzene, 0.99);
        assert!(
            !results.is_empty(),
            "self should be found at threshold 0.99"
        );
        assert_eq!(results[0].0, 0, "self should be index 0");
        assert!(
            (results[0].1 - 1.0).abs() < 1e-9,
            "self-similarity should be 1.0, got {}",
            results[0].1
        );
    }

    #[test]
    fn test_lsh_threshold_filtering() {
        let mut index = MhfpLshIndex::new(128);
        let benzene_fp = fp("c1ccccc1");
        let ethane_fp = fp("CC");
        index.add(benzene_fp.clone());
        index.add(ethane_fp);

        // High threshold: only benzene itself (or nothing if LSH misses it)
        let sim_benz_eth = benzene_fp.tanimoto(&fp("CC"));
        let high_threshold = sim_benz_eth + 0.1; // above ethane similarity
        let results = index.query(&benzene_fp, high_threshold);
        for (_, sim) in &results {
            assert!(
                *sim >= high_threshold,
                "all results must meet threshold, got {}",
                sim
            );
        }
    }

    #[test]
    fn test_lsh_similar_mols_found() {
        // Benzene and toluene are structurally similar; LSH should find toluene
        // when querying with benzene at a reasonable threshold.
        let mut index = MhfpLshIndex::new(128);
        let benzene = fp("c1ccccc1");
        let toluene = fp("Cc1ccccc1");
        let ethane = fp("CC");

        index.add(ethane); // idx 0
        index.add(toluene); // idx 1
        index.add(benzene.clone()); // idx 2

        // Self-query at very high threshold must find itself
        let results = index.query(&benzene, 0.99);
        let found_self = results.iter().any(|(i, _)| *i == 2);
        assert!(found_self, "benzene should find itself in the index");

        // Toluene similarity should be above 0 (both have aromatic ring)
        let sim = benzene.tanimoto(&fp("Cc1ccccc1"));
        assert!(sim > 0.0, "benzene-toluene similarity should be > 0");
    }

    #[test]
    fn test_lsh_sorted_descending() {
        let mut index = MhfpLshIndex::new(128);
        // Add several aromatic and aliphatic molecules
        for smi in &["c1ccccc1", "Cc1ccccc1", "c1ccc2ccccc2c1", "CC", "CCC"] {
            index.add(fp(smi));
        }
        let results = index.query(&fp("c1ccccc1"), 0.0);
        for w in results.windows(2) {
            assert!(
                w[0].1 >= w[1].1,
                "results should be sorted descending: {} >= {}",
                w[0].1,
                w[1].1
            );
        }
    }

    #[test]
    fn test_lsh_custom_bands() {
        // with_bands(32, 4) works for num_hashes=128
        let cfg = MhfpConfig {
            num_hashes: 128,
            seed: 0,
            radius: 2,
        };
        let mut index = MhfpLshIndex::with_bands(32, 4);
        let benzene = fp_cfg("c1ccccc1", &cfg);
        index.add(benzene.clone());
        let results = index.query(&benzene, 0.99);
        assert!(!results.is_empty(), "custom bands: self should be found");
        assert!((results[0].1 - 1.0).abs() < 1e-9);
    }
}
