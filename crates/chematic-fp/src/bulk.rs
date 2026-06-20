//! Bulk Tanimoto computation over BitVec2048 fingerprint slices.
//!
//! All functions return `f32` scores (vs the `f64` used by `search`)
//! to halve matrix memory footprint while retaining full precision for
//! popcounts ≤ 2048.
//!
//! # Performance
//! - Query popcount is hoisted outside the inner loop (computed once per query).
//! - No heap allocations inside the hot path (`intersection_popcount` operates
//!   in-place over the word arrays).
//! - The tight `zip + AND + count_ones` loop autovectorizes to POPCNT on x86
//!   and to `vcnt` on ARM without unsafe code.

use crate::bitvec::BitVec2048;

/// Compute Tanimoto similarity between `query` and every entry in `db`.
///
/// Returns a **dense** `Vec<f32>` of length `db.len()` with no zero-filtering.
/// Position `i` in the output corresponds directly to `db[i]`.
/// Returns `1.0` for any pair where both fingerprints are all-zero.
///
/// # Complexity
/// O(N × 32) word operations; autovectorizable.
pub fn tanimoto_slice(query: &BitVec2048, db: &[BitVec2048]) -> Vec<f32> {
    if db.is_empty() {
        return vec![];
    }
    let qa = query.popcount();
    db.iter()
        .map(|fp| query.tanimoto_with_counts(fp, qa, fp.popcount()))
        .collect()
}

/// Compute an M×N Tanimoto similarity matrix, returned as a flat row-major `Vec<f32>`.
///
/// `result[i * db.len() + j]` is the Tanimoto of `queries[i]` against `db[j]`.
///
/// All query and db popcounts are precomputed before the inner loop,
/// so each popcount is computed exactly once.
///
/// Returns an empty vec if either slice is empty.
///
/// # Complexity
/// O(M×N×32) word operations (M=queries, N=db).
pub fn tanimoto_matrix(queries: &[BitVec2048], db: &[BitVec2048]) -> Vec<f32> {
    let m = queries.len();
    let n = db.len();
    if m == 0 || n == 0 {
        return vec![];
    }
    let q_counts: Vec<u32> = queries.iter().map(|q| q.popcount()).collect();
    let d_counts: Vec<u32> = db.iter().map(|d| d.popcount()).collect();

    let mut out = Vec::with_capacity(m * n);
    for (i, q) in queries.iter().enumerate() {
        for (j, d) in db.iter().enumerate() {
            out.push(q.tanimoto_with_counts(d, q_counts[i], d_counts[j]));
        }
    }
    out
}

/// Return the `k` most similar entries in `db` to `query`, sorted by descending Tanimoto.
///
/// Returns `Vec<(usize, f32)>` of `(db_index, tanimoto)` pairs.
/// Returns `min(k, db.len())` results — **no zero-filtering** (unlike `nearest_neighbors_from_fp`).
/// Returns `[]` if `k == 0` or `db` is empty.
///
/// # Complexity
/// O(N) average for partial selection, then O(k log k) to sort the front.
pub fn top_k_similar(query: &BitVec2048, db: &[BitVec2048], k: usize) -> Vec<(usize, f32)> {
    if k == 0 || db.is_empty() {
        return vec![];
    }
    let k = k.min(db.len());
    let scores = tanimoto_slice(query, db);
    let mut indexed: Vec<(usize, f32)> = scores.into_iter().enumerate().collect();

    // Partial selection: move top-k to front in O(N) average
    indexed.select_nth_unstable_by(k - 1, |a, b| {
        b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    indexed.truncate(k);
    indexed.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    indexed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecfp::ecfp4;
    use chematic_smiles::parse;

    fn fp(smi: &str) -> BitVec2048 {
        ecfp4(&parse(smi).unwrap())
    }

    #[test]
    fn test_tanimoto_slice_parity_with_single() {
        let q = fp("c1ccccc1");
        let db = vec![fp("CC"), fp("c1ccccc1N"), fp("CCO"), fp("c1ccccc1")];
        let bulk = tanimoto_slice(&q, &db);
        for (i, bv) in db.iter().enumerate() {
            let single = q.tanimoto(bv) as f32;
            assert!(
                (bulk[i] - single).abs() < 1e-5,
                "slice[{i}] = {:.6} ≠ single {:.6}",
                bulk[i],
                single
            );
        }
    }

    #[test]
    fn test_tanimoto_slice_dense_no_zero_filter() {
        let q = fp("c1ccccc1");
        // FNV hash-mapped ECFP4 for "CC" and "c1ccccc1" are disjoint → Tanimoto = 0
        // The result must still be present (no filtering).
        let db = vec![fp("CC")];
        let result = tanimoto_slice(&q, &db);
        assert_eq!(result.len(), 1);
        // Score may be 0 or small — just ensure it's in [0, 1]
        assert!((0.0..=1.0).contains(&result[0]));
    }

    #[test]
    fn test_tanimoto_slice_all_zero_convention() {
        let zero = BitVec2048::new();
        let result = tanimoto_slice(&zero, std::slice::from_ref(&zero));
        assert_eq!(result.len(), 1);
        assert!((result[0] - 1.0).abs() < 1e-6, "both zero → 1.0");
    }

    #[test]
    fn test_tanimoto_slice_empty_db() {
        let q = fp("CC");
        let result = tanimoto_slice(&q, &[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_tanimoto_matrix_row_major() {
        let queries = vec![fp("CC"), fp("c1ccccc1")];
        let db = vec![fp("CCO"), fp("c1ccccc1N"), fp("CCCC")];
        let matrix = tanimoto_matrix(&queries, &db);
        assert_eq!(matrix.len(), queries.len() * db.len());

        for (i, q) in queries.iter().enumerate() {
            let row = tanimoto_slice(q, &db);
            for (j, &expected) in row.iter().enumerate() {
                assert!(
                    (matrix[i * db.len() + j] - expected).abs() < 1e-6,
                    "matrix[{i}][{j}] mismatch"
                );
            }
        }
    }

    #[test]
    fn test_tanimoto_matrix_empty() {
        let q = vec![fp("CC")];
        assert!(tanimoto_matrix(&q, &[]).is_empty());
        assert!(tanimoto_matrix(&[], &q).is_empty());
    }

    #[test]
    fn test_top_k_similar_sorted_descending() {
        let query = fp("c1ccccc1");
        let db = vec![fp("CC"), fp("c1ccccc1"), fp("c1ccccc1N"), fp("CCCC")];
        let k = 3;
        let hits = top_k_similar(&query, &db, k);
        assert_eq!(hits.len(), k);
        for w in hits.windows(2) {
            assert!(
                w[0].1 >= w[1].1,
                "not sorted: {:.4} < {:.4}",
                w[0].1,
                w[1].1
            );
        }
    }

    #[test]
    fn test_top_k_similar_matches_sorted_slice() {
        let query = fp("c1ccccc1");
        let db = vec![
            fp("CC"),
            fp("c1ccccc1"),
            fp("CCO"),
            fp("c1ccccc1N"),
            fp("CCCC"),
        ];
        let k = 3;
        let hits = top_k_similar(&query, &db, k);

        // Independent reference: sort the full slice
        let mut all: Vec<(usize, f32)> = tanimoto_slice(&query, &db)
            .into_iter()
            .enumerate()
            .collect();
        all.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        all.truncate(k);

        for (i, (hit_idx, hit_score)) in hits.iter().enumerate() {
            assert_eq!(*hit_idx, all[i].0, "index mismatch at rank {i}");
            assert!(
                (hit_score - all[i].1).abs() < 1e-6,
                "score mismatch at rank {i}"
            );
        }
    }

    #[test]
    fn test_top_k_similar_k_zero() {
        let query = fp("CC");
        let db = vec![fp("CCC")];
        assert!(top_k_similar(&query, &db, 0).is_empty());
    }

    #[test]
    fn test_top_k_similar_k_exceeds_db() {
        let query = fp("CC");
        let db = vec![fp("CCC"), fp("CCCC")];
        let hits = top_k_similar(&query, &db, 100);
        assert_eq!(hits.len(), db.len(), "returns all when k > db.len()");
    }

    #[test]
    fn test_top_k_similar_empty_db() {
        let query = fp("CC");
        assert!(top_k_similar(&query, &[], 5).is_empty());
    }

    #[test]
    fn test_self_similarity_is_one() {
        let query = fp("c1ccccc1N");
        let db = vec![fp("c1ccccc1N")];
        let hits = top_k_similar(&query, &db, 1);
        assert_eq!(hits.len(), 1);
        assert!(
            (hits[0].1 - 1.0).abs() < 1e-6,
            "self-similarity must be 1.0"
        );
    }
}
