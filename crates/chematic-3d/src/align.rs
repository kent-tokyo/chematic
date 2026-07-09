//! Standalone 3D coordinate alignment using the Kabsch algorithm.
//!
//! [`align_coords`] finds the rotation and translation that minimises the RMSD
//! between two sets of paired 3D coordinates.  The implementation reuses the
//! Jacobi eigensolver from [`crate::shape_descriptors`].
//!
//! [`rmsd_no_align`] computes the raw RMSD without any superposition.

use crate::shape_descriptors::jacobi3;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Result of a Kabsch alignment.
#[derive(Debug, Clone)]
pub struct AlignResult {
    /// RMSD after optimal superposition (Å or whatever unit `reference` uses).
    pub rmsd: f64,
    /// 3×3 rotation matrix R such that `mobile_centred @ R ≈ reference_centred`.
    pub rotation: [[f64; 3]; 3],
    /// Translation to apply to mobile *after* rotation to match reference centroid.
    pub translation: [f64; 3],
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute the RMSD between two sets of paired coordinates **without** any
/// rotation or translation.
///
/// `a` and `b` must have the same length.  Returns `0.0` for empty slices.
pub fn rmsd_no_align(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = a
        .iter()
        .zip(b.iter())
        .map(|(pa, pb)| (0..3).map(|i| (pa[i] - pb[i]).powi(2)).sum::<f64>())
        .sum();
    (sum_sq / n as f64).sqrt()
}

/// Find the optimal superposition of `mobile` onto `reference` using the
/// Kabsch algorithm and return the RMSD plus the rotation matrix and
/// translation vector.
///
/// The two slices must have the same length (atom-correspondence is assumed).
/// Returns an [`AlignResult`] with `rmsd = 0.0` for empty or single-atom inputs.
///
/// To obtain the aligned coordinates, apply:
/// ```text
/// aligned[i][j] = Σ_k rotation[j][k] * (mobile[i][k] - mobile_centroid[k])
///                + reference_centroid[j]
/// ```
/// Or use [`apply_alignment`].
pub fn align_coords(reference: &[[f64; 3]], mobile: &[[f64; 3]]) -> AlignResult {
    let n = reference.len().min(mobile.len());
    if n == 0 {
        return AlignResult {
            rmsd: 0.0,
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0, 0.0, 0.0],
        };
    }

    let nf = n as f64;

    // Centroids.
    let mut cr = [0.0f64; 3];
    let mut cm = [0.0f64; 3];
    for i in 0..n {
        for k in 0..3 {
            cr[k] += reference[i][k];
            cm[k] += mobile[i][k];
        }
    }
    for k in 0..3 {
        cr[k] /= nf;
        cm[k] /= nf;
    }

    // Centered coordinates.
    let p: Vec<[f64; 3]> = reference
        .iter()
        .map(|v| [v[0] - cr[0], v[1] - cr[1], v[2] - cr[2]])
        .collect();
    let q: Vec<[f64; 3]> = mobile
        .iter()
        .map(|v| [v[0] - cm[0], v[1] - cm[1], v[2] - cm[2]])
        .collect();

    // H = P^T * Q  (3×3 cross-covariance).
    let mut h = [[0.0f64; 3]; 3];
    for i in 0..n {
        for r in 0..3 {
            for c in 0..3 {
                h[r][c] += p[i][r] * q[i][c];
            }
        }
    }

    // SVD via eigendecomposition of H^T * H.
    let mut hth = [[0.0f64; 3]; 3];
    for r in 0..3 {
        for c in 0..3 {
            for k in 0..3 {
                hth[r][c] += h[k][r] * h[k][c];
            }
        }
    }
    let (evals, v) = jacobi3(hth);

    // U = H * V * diag(1/σ).
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

    // R = U * V^T (maximises trace(H^T R), the Kabsch objective).
    let mut rot = [[0.0f64; 3]; 3];
    let mut v_final = v;
    for r in 0..3 {
        for c in 0..3 {
            for k in 0..3 {
                rot[r][c] += u[r][k] * v_final[c][k];
            }
        }
    }

    // Reflection correction: flip the column of V paired with the smallest
    // singular value (index 0, ascending order) when det(R) < 0.
    let det = det3(rot);
    if det < 0.0 {
        for r in 0..3 {
            v_final[r][0] *= -1.0;
        }
        rot = [[0.0; 3]; 3];
        for r in 0..3 {
            for c in 0..3 {
                for k in 0..3 {
                    rot[r][c] += u[r][k] * v_final[c][k];
                }
            }
        }
    }

    // RMSD after optimal rotation.
    let mut sum_sq = 0.0f64;
    for i in 0..n {
        for row in 0..3 {
            let rotated = (0..3).map(|k| rot[row][k] * q[i][k]).sum::<f64>();
            let diff = p[i][row] - rotated;
            sum_sq += diff * diff;
        }
    }
    let rmsd = (sum_sq / nf).sqrt();

    // Translation: after rotating mobile around its centroid, shift to reference centroid.
    let translation = [cr[0] - cm[0], cr[1] - cm[1], cr[2] - cm[2]];

    AlignResult {
        rmsd,
        rotation: rot,
        translation,
    }
}

/// Apply an [`AlignResult`] to transform `mobile` coordinates.
///
/// Returns a new `Vec` of aligned coordinates.
pub fn apply_alignment(mobile: &[[f64; 3]], result: &AlignResult) -> Vec<[f64; 3]> {
    // Compute mobile centroid.
    let n = mobile.len();
    if n == 0 {
        return Vec::new();
    }
    let mut cm = [0.0f64; 3];
    for v in mobile {
        for k in 0..3 {
            cm[k] += v[k];
        }
    }
    for k in 0..3 {
        cm[k] /= n as f64;
    }

    // Compute reference centroid from translation.
    // translation = cr - cm_orig, so cr = cm_orig + translation.
    // We need the reference centroid to shift into, but it's embedded in translation.
    // Simpler: rotate (mobile - cm) then add (cm + translation).
    let cr = [
        cm[0] + result.translation[0],
        cm[1] + result.translation[1],
        cm[2] + result.translation[2],
    ];

    mobile
        .iter()
        .map(|v| {
            let centered = [v[0] - cm[0], v[1] - cm[1], v[2] - cm[2]];
            let mut out = [0.0f64; 3];
            for row in 0..3 {
                out[row] = (0..3)
                    .map(|k| result.rotation[row][k] * centered[k])
                    .sum::<f64>()
                    + cr[row];
            }
            out
        })
        .collect()
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

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn test_rmsd_no_align_identical() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        assert!(approx_eq(rmsd_no_align(&coords, &coords), 0.0, 1e-10));
    }

    #[test]
    fn test_rmsd_no_align_empty() {
        assert_eq!(rmsd_no_align(&[], &[]), 0.0);
    }

    #[test]
    fn test_rmsd_no_align_translated() {
        let a = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let b = vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]; // translated by 1 Å along x
        // RMSD without align = 1.0
        assert!(approx_eq(rmsd_no_align(&a, &b), 1.0, 1e-9));
    }

    #[test]
    fn test_align_identical() {
        let coords = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.75, 1.3, 0.0]];
        let result = align_coords(&coords, &coords);
        assert!(
            approx_eq(result.rmsd, 0.0, 1e-9),
            "identical coords → RMSD 0"
        );
    }

    #[test]
    fn test_align_pure_translation() {
        // Two identical shapes offset by a constant translation → RMSD after align = 0.
        let reference = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let mobile: Vec<[f64; 3]> = reference
            .iter()
            .map(|v| [v[0] + 3.0, v[1] - 2.0, v[2] + 1.0])
            .collect();
        let result = align_coords(&reference, &mobile);
        assert!(
            approx_eq(result.rmsd, 0.0, 1e-6),
            "pure translation → RMSD 0 after Kabsch"
        );
    }

    #[test]
    fn test_align_different_shapes_nonzero_rmsd() {
        let reference = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let mobile = vec![[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [0.5, 1.1, 0.0]]; // perturbed
        let result = align_coords(&reference, &mobile);
        assert!(result.rmsd > 0.0, "different shapes → RMSD > 0");
    }

    #[test]
    fn test_apply_alignment_reduces_rmsd() {
        let reference = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let mobile: Vec<[f64; 3]> = reference
            .iter()
            .map(|v| [v[0] + 2.0, v[1] + 2.0, v[2]])
            .collect();
        let result = align_coords(&reference, &mobile);
        let aligned = apply_alignment(&mobile, &result);
        let rmsd_after = rmsd_no_align(&reference, &aligned);
        assert!(
            approx_eq(rmsd_after, result.rmsd, 1e-6),
            "apply_alignment should match reported RMSD"
        );
    }

    /// Regression test for a rotation-direction bug: `align_coords` once
    /// returned the reference→mobile rotation instead of mobile→reference,
    /// which is invisible under pure translation (H is symmetric there, so
    /// U == V and the two directions coincide) but produces a large bogus
    /// RMSD once an actual rotation is involved.
    #[test]
    fn test_align_pure_rotation_recovers_zero_rmsd() {
        let reference = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [1.0, 1.0, 4.0],
        ];
        // 90° rotation about z (x,y,z) -> (-y,x,z), plus a translation.
        let mobile: Vec<[f64; 3]> = reference
            .iter()
            .map(|p| [-p[1] + 5.0, p[0] - 3.0, p[2] + 2.0])
            .collect();
        let result = align_coords(&reference, &mobile);
        assert!(
            approx_eq(result.rmsd, 0.0, 1e-6),
            "pure rotation + translation → RMSD ~0 after Kabsch, got {}",
            result.rmsd
        );
    }

    #[test]
    fn test_apply_alignment_under_rotation_matches_reference_pointwise() {
        let reference = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [1.0, 1.0, 4.0],
        ];
        let mobile: Vec<[f64; 3]> = reference
            .iter()
            .map(|p| [-p[1] + 5.0, p[0] - 3.0, p[2] + 2.0])
            .collect();
        let result = align_coords(&reference, &mobile);
        let aligned = apply_alignment(&mobile, &result);
        for (r, a) in reference.iter().zip(aligned.iter()) {
            for k in 0..3 {
                assert!(
                    approx_eq(r[k], a[k], 1e-6),
                    "aligned point {a:?} should match reference point {r:?}"
                );
            }
        }
    }

    /// Forces the det<0 reflection-correction branch: a mirror image of a
    /// genuinely non-planar (chiral) point set cannot be superposed by any
    /// proper rotation, so the corrected result must still have det(R) = +1
    /// and a nonzero residual RMSD (never a silent improper "rotation").
    #[test]
    fn test_align_mirror_image_forces_reflection_correction() {
        let reference = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [1.0, 1.0, 4.0],
        ];
        // Reflect through the xy-plane (z -> -z): improper for a non-planar set.
        let mobile: Vec<[f64; 3]> = reference.iter().map(|p| [p[0], p[1], -p[2]]).collect();
        let result = align_coords(&reference, &mobile);
        assert!(
            approx_eq(det3(result.rotation), 1.0, 1e-6),
            "corrected rotation must be proper (det=+1), got det={}",
            det3(result.rotation)
        );
        assert!(
            result.rmsd > 0.5,
            "a chiral point set can't be superposed onto its mirror image by any rotation, rmsd={}",
            result.rmsd
        );
    }
}
