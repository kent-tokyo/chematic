//! Ultrafast Shape Recognition (USR) 3D similarity descriptors.
//!
//! Implements Ballester & Richards (2007): each molecule is summarised by 12
//! moments (mean, variance, skewness) of four distance distributions measured
//! from four reference points:
//!
//! 1. **ctd** — centroid of all atoms
//! 2. **cst** — atom closest to ctd
//! 3. **fct** — atom farthest from ctd
//! 4. **ftf** — atom farthest from fct
//!
//! Similarity is computed as a Soergel distance over the 12 moments:
//! `S(a, b) = 1 / (1 + (1/12) Σ|a_i − b_i|)`.
//!
//! Reference: Ballester, P.J.; Richards, W.G. (2007). *Ultrafast shape
//! recognition to search compound databases for similar molecular shapes*.
//! J. Comput. Chem. 28, 1711–1723.

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the 12 USR shape descriptors for a set of 3D coordinates.
///
/// `coords` is a slice of `[x, y, z]` positions in any consistent unit (Å).
/// Returns `[f64; 12]` — three statistical moments for each of four reference
/// points (ctd, cst, fct, ftf).  Returns all-zeros for empty or
/// single-atom inputs.
pub fn usr_descriptors(coords: &[[f64; 3]]) -> [f64; 12] {
    let n = coords.len();
    if n <= 1 {
        return [0.0; 12];
    }

    // 1. Centroid (ctd).
    let ctd = centroid(coords);

    // 2. Distances from ctd.
    let d_ctd: Vec<f64> = coords.iter().map(|p| dist(p, &ctd)).collect();

    // 3. Closest atom to ctd (cst).
    let cst_idx = extreme_idx(coords, &d_ctd, Extreme::Min);
    let cst = coords[cst_idx];

    // 4. Farthest atom from ctd (fct).
    //
    // The original `min_by`/`max_by` calls had NO tie-break at all
    // (`max_by` returns the LAST element on an exact tie per its own docs,
    // `min_by` the FIRST) -- an inconsistent default that made the reference
    // point (and thus all 12 descriptors) depend on atom input order for
    // two chemically-different atoms tying on distance (confirmed
    // empirically: two points placed symmetrically so they tie exactly on
    // distance from a third, unambiguous reference point -- swapping their
    // ARRAY POSITIONS changed which one became the reference point).
    //
    // An index-based tie-break ("smallest index wins") is NOT sufficient to
    // fix this: swapping two points' array positions is exactly the
    // operation that changes which physical point occupies the "smallest
    // index" -- so an index-based rule stays internally consistent but
    // still varies with array order. `extreme_idx` tie-breaks on the
    // coordinate VALUES themselves (content, not position) instead, which
    // is invariant to which order the SAME set of points was supplied in.
    let fct_idx = extreme_idx(coords, &d_ctd, Extreme::Max);
    let fct = coords[fct_idx];

    // 5. Farthest atom from fct (ftf). Same content-based tie-break.
    let d_fct: Vec<f64> = coords.iter().map(|p| dist(p, &fct)).collect();
    let ftf_idx = extreme_idx(coords, &d_fct, Extreme::Max);
    let ftf = coords[ftf_idx];

    // 6. Compute moment triplets for each reference point.
    let d_cst: Vec<f64> = coords.iter().map(|p| dist(p, &cst)).collect();
    let d_ftf: Vec<f64> = coords.iter().map(|p| dist(p, &ftf)).collect();

    let mut out = [0.0f64; 12];
    let (m0, m1, m2) = moments(&d_ctd);
    out[0] = m0;
    out[1] = m1;
    out[2] = m2;
    let (m0, m1, m2) = moments(&d_cst);
    out[3] = m0;
    out[4] = m1;
    out[5] = m2;
    let (m0, m1, m2) = moments(&d_fct);
    out[6] = m0;
    out[7] = m1;
    out[8] = m2;
    let (m0, m1, m2) = moments(&d_ftf);
    out[9] = m0;
    out[10] = m1;
    out[11] = m2;
    out
}

/// Compute the USR similarity between two descriptor vectors.
///
/// `S = 1 / (1 + (1/12) Σ|a_i − b_i|)`
///
/// Returns 1.0 for identical descriptors and approaches 0.0 for very
/// different shapes.
pub fn usr_similarity(a: &[f64; 12], b: &[f64; 12]) -> f64 {
    let sum: f64 = a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).sum();
    1.0 / (1.0 + sum / 12.0)
}

// ---------------------------------------------------------------------------
// Convenience wrappers using distance-geometry 3D coordinates
// ---------------------------------------------------------------------------

/// Compute USR descriptors for a molecule using auto-generated 3D coordinates.
///
/// Generates rule-based 3D coordinates via distance geometry, then computes
/// 12 USR moments. Heavy atoms only (H excluded).
pub fn usr_from_dg(mol: &chematic_core::Molecule) -> [f64; 12] {
    use chematic_core::AtomIdx;
    let coords = crate::dg::generate_coords(mol);
    let pts: Vec<[f64; 3]> = (0..mol.atom_count())
        .filter(|&i| mol.atom(AtomIdx(i as u32)).element.atomic_number() != 1)
        .map(|i| {
            let p = coords.get(AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();
    usr_descriptors(&pts)
}

/// Screen a library of molecules by 3D shape similarity to a query molecule.
///
/// Returns `(index, similarity)` pairs sorted by decreasing similarity.
/// Both query and library molecules have 3D coordinates auto-generated.
pub fn shape_screen(
    query: &chematic_core::Molecule,
    library: &[&chematic_core::Molecule],
) -> Vec<(usize, f64)> {
    let q_desc = usr_from_dg(query);
    let mut results: Vec<(usize, f64)> = library
        .iter()
        .enumerate()
        .map(|(i, mol)| {
            let desc = usr_from_dg(mol);
            (i, usr_similarity(&q_desc, &desc))
        })
        .collect();
    results.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    results
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

pub(crate) fn centroid(coords: &[[f64; 3]]) -> [f64; 3] {
    let n = coords.len() as f64;
    let mut c = [0.0f64; 3];
    for p in coords {
        c[0] += p[0];
        c[1] += p[1];
        c[2] += p[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

enum Extreme {
    Min,
    Max,
}

/// Index of the min/max `values[i]`, tie-breaking on the lexicographic
/// (x, then y, then z) order of `coords[i]` itself -- content, not array
/// position, so the result is the same regardless of which order the same
/// set of points was supplied in. An index-based tie-break (e.g. "smallest
/// index wins") is NOT sufficient here: swapping two tied points' array
/// positions is exactly the operation that changes which point occupies
/// the smallest index, so it would still vary with input order.
fn extreme_idx(coords: &[[f64; 3]], values: &[f64], which: Extreme) -> usize {
    let mut best = 0usize;
    for i in 1..values.len() {
        let better = match which {
            Extreme::Min => values[i] < values[best],
            Extreme::Max => values[i] > values[best],
        };
        let tied = values[i] == values[best];
        if better || (tied && coords[i] < coords[best]) {
            best = i;
        }
    }
    best
}

/// First three standardised moments (mean, variance, skewness) of a distance distribution.
fn moments(ds: &[f64]) -> (f64, f64, f64) {
    let n = ds.len() as f64;
    let mean = ds.iter().sum::<f64>() / n;
    let var = ds.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / n;
    let skew = ds.iter().map(|d| (d - mean).powi(3)).sum::<f64>() / n;
    (mean, var, skew)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn square() -> Vec<[f64; 3]> {
        // Unit square in the XY plane.
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    }

    #[test]
    fn test_reference_point_ties_are_order_independent() {
        // An outlier far enough to be the unambiguous fct (farthest from
        // centroid, no tie there), plus two points P and Q placed
        // symmetrically in y so they tie EXACTLY on distance from the
        // outlier (the ftf pick) -- but a filler F off that mirror plane
        // (nonzero z) breaks any symmetry of the point set AS A WHOLE, so
        // swapping P and Q is a genuine non-automorphic tie, not a
        // harmless one. Before the fix, swapping which of P/Q came first
        // in the input flipped which one became ftf, changing d_ftf (and
        // thus 3 of the 12 descriptor values) since P != Q in position.
        let outlier = [100.0, 0.0, 0.0];
        let p = [0.0, 5.0, 0.0];
        let q = [0.0, -5.0, 0.0];
        let filler = [0.0, 2.0, 3.0];
        let a = vec![outlier, p, q, filler];
        let mut b = a.clone();
        b.swap(1, 2); // same point SET, P and Q swapped in input order
        let da = usr_descriptors(&a);
        let db = usr_descriptors(&b);
        for (x, y) in da.iter().zip(db.iter()) {
            // Tolerance, not exact equality: summing the same set of floats
            // in a different order is expected to differ in the last few
            // bits (IEEE 754 non-associativity) even when the SAME reference
            // atoms are picked -- this test's job is to catch the previous
            // bug's ~6-unit divergence from picking a DIFFERENT atom, not to
            // demand bit-identical floating-point summation order.
            assert!(
                (x - y).abs() < 1e-9,
                "usr_descriptors depends on the input order of tied reference-point candidates: {da:?} vs {db:?}"
            );
        }
    }

    #[test]
    fn test_self_similarity_is_one() {
        let coords = square();
        let d = usr_descriptors(&coords);
        assert!((usr_similarity(&d, &d) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_different_shapes_less_than_one() {
        let square = square();
        // Elongated shape (line).
        let line = vec![
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [15.0, 0.0, 0.0],
        ];
        let d_sq = usr_descriptors(&square);
        let d_ln = usr_descriptors(&line);
        let sim = usr_similarity(&d_sq, &d_ln);
        assert!((0.0..1.0).contains(&sim), "sim={sim:.4} should be in [0,1)");
    }

    #[test]
    fn test_empty_input() {
        let d = usr_descriptors(&[]);
        assert_eq!(d, [0.0f64; 12]);
    }

    #[test]
    fn test_single_atom() {
        let d = usr_descriptors(&[[1.0, 2.0, 3.0]]);
        assert_eq!(d, [0.0f64; 12]);
    }

    #[test]
    fn test_usr_from_dg_benzene() {
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let desc = usr_from_dg(&mol);
        // Self-similarity should be 1.0
        assert!((usr_similarity(&desc, &desc) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_usr_from_dg_different_molecules() {
        let benzene = chematic_smiles::parse("c1ccccc1").unwrap();
        let ethanol = chematic_smiles::parse("CCO").unwrap();
        let d_benz = usr_from_dg(&benzene);
        let d_eth = usr_from_dg(&ethanol);
        let sim = usr_similarity(&d_benz, &d_eth);
        assert!(
            sim < 1.0,
            "benzene and ethanol should not have identical shape"
        );
        assert!(sim >= 0.0, "similarity must be non-negative");
    }

    #[test]
    fn test_shape_screen_returns_sorted() {
        let query = chematic_smiles::parse("c1ccccc1").unwrap();
        let toluene = chematic_smiles::parse("Cc1ccccc1").unwrap();
        let ethanol = chematic_smiles::parse("CCO").unwrap();
        let library = vec![&ethanol as &chematic_core::Molecule, &toluene];
        let results = shape_screen(&query, &library);
        assert_eq!(results.len(), 2);
        // Results should be sorted descending
        assert!(
            results[0].1 >= results[1].1,
            "results should be sorted by decreasing similarity"
        );
        // Toluene (index 1) should score higher than ethanol (index 0) vs benzene
        let toluene_hit = results.iter().find(|(i, _)| *i == 1);
        assert!(toluene_hit.is_some());
    }

    #[test]
    fn test_translation_invariance() {
        let c1 = square();
        let offset = 100.0;
        let c2: Vec<[f64; 3]> = c1
            .iter()
            .map(|p| [p[0] + offset, p[1] + offset, p[2] + offset])
            .collect();
        // USR is NOT translation-invariant by design (it uses absolute distances
        // from reference points, not pairwise distances). However, centroid-based
        // reference points make it centroid-relative, so pure translation DOES
        // preserve the descriptor.
        let d1 = usr_descriptors(&c1);
        let d2 = usr_descriptors(&c2);
        let sim = usr_similarity(&d1, &d2);
        assert!(
            (sim - 1.0).abs() < 1e-6,
            "USR should be translation-invariant, sim={sim}"
        );
    }
}
