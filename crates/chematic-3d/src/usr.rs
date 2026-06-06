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
    let cst_idx = d_ctd
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let cst = coords[cst_idx];

    // 4. Farthest atom from ctd (fct).
    let fct_idx = d_ctd
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let fct = coords[fct_idx];

    // 5. Farthest atom from fct (ftf).
    let ftf_idx = coords
        .iter()
        .enumerate()
        .max_by(|a, b| dist(a.1, &fct).partial_cmp(&dist(b.1, &fct)).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    let ftf = coords[ftf_idx];

    // 6. Compute moment triplets for each reference point.
    let d_cst: Vec<f64> = coords.iter().map(|p| dist(p, &cst)).collect();
    let d_fct: Vec<f64> = coords.iter().map(|p| dist(p, &fct)).collect();
    let d_ftf: Vec<f64> = coords.iter().map(|p| dist(p, &ftf)).collect();

    let mut out = [0.0f64; 12];
    let (m0, m1, m2) = moments(&d_ctd);
    out[0] = m0; out[1] = m1; out[2] = m2;
    let (m0, m1, m2) = moments(&d_cst);
    out[3] = m0; out[4] = m1; out[5] = m2;
    let (m0, m1, m2) = moments(&d_fct);
    out[6] = m0; out[7] = m1; out[8] = m2;
    let (m0, m1, m2) = moments(&d_ftf);
    out[9] = m0; out[10] = m1; out[11] = m2;
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
// Helpers
// ---------------------------------------------------------------------------

fn centroid(coords: &[[f64; 3]]) -> [f64; 3] {
    let n = coords.len() as f64;
    let mut c = [0.0f64; 3];
    for p in coords {
        c[0] += p[0]; c[1] += p[1]; c[2] += p[2];
    }
    [c[0] / n, c[1] / n, c[2] / n]
}

fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0]-b[0]).powi(2) + (a[1]-b[1]).powi(2) + (a[2]-b[2]).powi(2)).sqrt()
}

/// First three standardised moments (mean, variance, skewness) of a distance distribution.
fn moments(ds: &[f64]) -> (f64, f64, f64) {
    let n = ds.len() as f64;
    let mean = ds.iter().sum::<f64>() / n;
    let var  = ds.iter().map(|d| (d - mean).powi(2)).sum::<f64>() / n;
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
        vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0]]
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
        let line = vec![[0.0,0.0,0.0],[5.0,0.0,0.0],[10.0,0.0,0.0],[15.0,0.0,0.0]];
        let d_sq = usr_descriptors(&square);
        let d_ln = usr_descriptors(&line);
        let sim = usr_similarity(&d_sq, &d_ln);
        assert!(sim < 1.0 && sim >= 0.0, "sim={sim:.4} should be in [0,1)");
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
    fn test_translation_invariance() {
        let c1 = square();
        let offset = 100.0;
        let c2: Vec<[f64;3]> = c1.iter().map(|p| [p[0]+offset, p[1]+offset, p[2]+offset]).collect();
        // USR is NOT translation-invariant by design (it uses absolute distances
        // from reference points, not pairwise distances). However, centroid-based
        // reference points make it centroid-relative, so pure translation DOES
        // preserve the descriptor.
        let d1 = usr_descriptors(&c1);
        let d2 = usr_descriptors(&c2);
        let sim = usr_similarity(&d1, &d2);
        assert!((sim - 1.0).abs() < 1e-6, "USR should be translation-invariant, sim={sim}");
    }
}
