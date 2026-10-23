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
    /// 3×3 rotation matrix R such that `R * mobile_centred ≈ reference_centred`.
    pub rotation: [[f64; 3]; 3],
    /// Translation in `aligned = rotation * mobile + translation`.
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
/// aligned[i][j] = Σ_k rotation[j][k] * mobile[i][k] + translation[j]
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

    let reference_centroid = centroid(&reference[..n]);
    let mobile_centroid = centroid(&mobile[..n]);
    let reference_centered = center(&reference[..n], reference_centroid);
    let mobile_centered = center(&mobile[..n], mobile_centroid);
    let rotation = kabsch_rotation(&reference_centered, &mobile_centered);
    let rmsd = aligned_rmsd(&reference_centered, &mobile_centered, rotation);
    let translation = std::array::from_fn(|axis| {
        reference_centroid[axis]
            - (0..3)
                .map(|column| rotation[axis][column] * mobile_centroid[column])
                .sum::<f64>()
    });

    AlignResult {
        rmsd,
        rotation,
        translation,
    }
}

fn centroid(coords: &[[f64; 3]]) -> [f64; 3] {
    let mut centroid = [0.0; 3];
    for point in coords {
        for axis in 0..3 {
            centroid[axis] += point[axis];
        }
    }
    for value in &mut centroid {
        *value /= coords.len() as f64;
    }
    centroid
}

fn center(coords: &[[f64; 3]], centroid: [f64; 3]) -> Vec<[f64; 3]> {
    coords
        .iter()
        .map(|point| {
            [
                point[0] - centroid[0],
                point[1] - centroid[1],
                point[2] - centroid[2],
            ]
        })
        .collect()
}

fn multiply_transpose(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut product = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            for inner in 0..3 {
                product[row][column] += left[row][inner] * right[column][inner];
            }
        }
    }
    product
}

fn kabsch_rotation(reference: &[[f64; 3]], mobile: &[[f64; 3]]) -> [[f64; 3]; 3] {
    let mut covariance = [[0.0; 3]; 3];
    for (reference, mobile) in reference.iter().zip(mobile) {
        for row in 0..3 {
            for column in 0..3 {
                covariance[row][column] += reference[row] * mobile[column];
            }
        }
    }
    let mut square = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            for inner in 0..3 {
                square[row][column] += covariance[inner][row] * covariance[inner][column];
            }
        }
    }
    let (values, right) = jacobi3(square);
    let mut left = [[0.0; 3]; 3];
    for column in 0..3 {
        let sigma = values[column].max(0.0).sqrt();
        for row in 0..3 {
            let projected = (0..3)
                .map(|inner| covariance[row][inner] * right[inner][column])
                .sum::<f64>();
            left[row][column] = if sigma > 1e-10 {
                projected / sigma
            } else {
                0.0
            };
        }
    }
    // An exactly planar or linear cloud has one or two zero singular values.
    // The corresponding projected left columns above are zero: using them
    // directly produces a projection matrix, not a rigid-body rotation.
    let rank = values
        .iter()
        .filter(|&&value| value.max(0.0).sqrt() > 1e-10)
        .count();
    if rank == 0 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    if rank < 3 {
        let dominant = [left[0][2], left[1][2], left[2][2]];
        let dominant = normalize3(dominant);
        for axis in 0..3 {
            left[axis][2] = dominant[axis];
        }
        let second = if rank == 2 {
            normalize3([left[0][1], left[1][1], left[2][1]])
        } else {
            // Choose a deterministic axis furthest from the dominant one.
            let axis = dominant
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.abs().total_cmp(&b.abs()))
                .map(|(index, _)| index)
                .unwrap();
            let mut basis = [0.0; 3];
            basis[axis] = 1.0;
            normalize3(cross3(dominant, basis))
        };
        for axis in 0..3 {
            left[axis][1] = second[axis];
        }
        let first = cross3(second, dominant);
        for axis in 0..3 {
            left[axis][0] = first[axis];
        }
    }
    let mut rotation = multiply_transpose(left, right);
    if det3(rotation) < 0.0 {
        for row in &mut left {
            row[0] *= -1.0;
        }
        rotation = multiply_transpose(left, right);
    }
    rotation
}

fn aligned_rmsd(reference: &[[f64; 3]], mobile: &[[f64; 3]], rotation: [[f64; 3]; 3]) -> f64 {
    let sum = reference
        .iter()
        .zip(mobile)
        .map(|(reference, mobile)| {
            (0..3)
                .map(|row| {
                    let rotated = (0..3)
                        .map(|column| rotation[row][column] * mobile[column])
                        .sum::<f64>();
                    (reference[row] - rotated).powi(2)
                })
                .sum::<f64>()
        })
        .sum::<f64>();
    (sum / reference.len() as f64).sqrt()
}

/// Apply an [`AlignResult`] to transform `mobile` coordinates.
///
/// Returns a new `Vec` of aligned coordinates.
pub fn apply_alignment(mobile: &[[f64; 3]], result: &AlignResult) -> Vec<[f64; 3]> {
    mobile
        .iter()
        .map(|v| {
            let mut out = [0.0f64; 3];
            for row in 0..3 {
                out[row] = (0..3).map(|k| result.rotation[row][k] * v[k]).sum::<f64>()
                    + result.translation[row];
            }
            out
        })
        .collect()
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    v.map(|x| x / norm)
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

    #[test]
    fn planar_fit_is_a_proper_rotation_and_transforms_extra_atoms() {
        let reference = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let mobile = [[5.0, -3.0, 2.0], [5.0, -1.0, 2.0], [2.0, -3.0, 2.0]];
        let fit = align_coords(&reference, &mobile);
        assert!(fit.rmsd < 1e-8);
        assert!((det3(fit.rotation) - 1.0).abs() < 1e-8);
        let extra = [4.0, -2.0, 4.0];
        let transformed = apply_alignment(&[extra], &fit)[0];
        let expected = [1.0, 1.0, 2.0];
        for axis in 0..3 {
            assert!((transformed[axis] - expected[axis]).abs() < 1e-8);
        }
        for (found, expected) in apply_alignment(&mobile, &fit).iter().zip(reference) {
            for axis in 0..3 {
                assert!((found[axis] - expected[axis]).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn linear_fit_preserves_distances_off_axis() {
        let reference = [[3.0, 4.0, 5.0], [3.0, 6.0, 5.0]];
        let mobile = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let fit = align_coords(&reference, &mobile);
        assert!(fit.rmsd < 1e-8);
        assert!((det3(fit.rotation) - 1.0).abs() < 1e-8);
        let transformed = apply_alignment(&[[0.0, 1.0, 0.0]], &fit)[0];
        let dx = transformed[0] - reference[0][0];
        let dy = transformed[1] - reference[0][1];
        let dz = transformed[2] - reference[0][2];
        assert!((dx * dx + dy * dy + dz * dz - 1.0).abs() < 1e-8);
    }
}
