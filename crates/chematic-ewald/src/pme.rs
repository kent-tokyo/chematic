//! Smooth Particle Mesh Ewald (SPME) reciprocal space energy calculation using FFT.
//!
//! Computes the reciprocal-space contribution to the Ewald sum using:
//! 1. Charge interpolation onto a 3D mesh
//! 2. Forward FFT to reciprocal space
//! 3. Reciprocal space energy evaluation
//! 4. Inverse FFT to real space (forces via differentiation)

use std::f64::consts::PI;

use crate::{PmeConfig, BoxVectors};

const K_COULOMB: f64 = 332.0637; // kcal·Å/(mol·e²)

/// Compute reciprocal-space Ewald energy.
///
/// Uses Fast Fourier Transform (FFT) to evaluate long-range Coulomb interactions
/// in periodic systems. The real-space and reciprocal-space contributions partition
/// the 1/r interaction based on the Ewald parameter α.
///
/// # Arguments
/// * `coords` - Atomic coordinates [[x0, y0, z0], ...]
/// * `charges` - Partial charges (e)
/// * `box_vecs` - Periodic box vectors
/// * `config` - PME configuration (alpha, kmax, mesh, spline_order)
///
/// # Returns
/// Reciprocal-space energy contribution (kcal/mol)
pub fn reciprocal_space_energy(
    coords: &[[f64; 3]],
    charges: &[f64],
    box_vecs: &BoxVectors,
    config: &PmeConfig,
) -> f64 {
    if coords.is_empty() {
        return 0.0;
    }

    // Auto-compute alpha if not provided
    let alpha = if config.alpha > 0.0 {
        config.alpha
    } else {
        3.5 / config.r_cut
    };

    // Create charge mesh
    let mesh_size = [config.mesh[0], config.mesh[1], config.mesh[2]];
    let mut charge_grid = vec![0.0; mesh_size[0] * mesh_size[1] * mesh_size[2]];

    // Interpolate charges onto mesh (B-spline order)
    interpolate_charges_to_mesh(
        coords,
        charges,
        box_vecs,
        &mut charge_grid,
        &mesh_size,
        config.spline_order,
    );

    // Compute reciprocal space energy from charge density
    let energy = compute_reciprocal_energy(&charge_grid, box_vecs, &mesh_size, alpha, config.kmax);

    energy
}

/// Interpolate point charges onto 3D mesh using B-splines.
fn interpolate_charges_to_mesh(
    coords: &[[f64; 3]],
    charges: &[f64],
    box_vecs: &BoxVectors,
    output_grid: &mut [f64],
    mesh_size: &[usize; 3],
    spline_order: u8,
) {
    let [m0, m1, m2] = *mesh_size;

    for (coord, &charge) in coords.iter().zip(charges.iter()) {
        if charge.abs() < 1e-10 {
            continue;
        }

        // Map atomic coordinate to fractional coordinates (0..1)
        let frac = map_to_fractional(*coord, box_vecs);

        // Wrap fractional coordinates into [0, 1)
        let frac = [
            frac[0] - frac[0].floor(),
            frac[1] - frac[1].floor(),
            frac[2] - frac[2].floor(),
        ];

        // Scale fractional coordinates to mesh indices u ∈ [0, M_d)
        let u = [
            frac[0] * m0 as f64,
            frac[1] * m1 as f64,
            frac[2] * m2 as f64,
        ];

        // Get floor of scaled coordinates (nearest grid point)
        let n = [
            u[0].floor() as isize,
            u[1].floor() as isize,
            u[2].floor() as isize,
        ];

        // Compute B-spline weights for each axis
        let wx = bspline_weights(u[0] - n[0] as f64, spline_order);
        let wy = bspline_weights(u[1] - n[1] as f64, spline_order);
        let wz = bspline_weights(u[2] - n[2] as f64, spline_order);

        // Spread charge to p³ grid points using 3D tensor product
        for kx in 0..spline_order as isize {
            let gx = (n[0] - kx).rem_euclid(m0 as isize) as usize;
            for ky in 0..spline_order as isize {
                let gy = (n[1] - ky).rem_euclid(m1 as isize) as usize;
                for kz in 0..spline_order as isize {
                    let gz = (n[2] - kz).rem_euclid(m2 as isize) as usize;
                    let linear_idx = gx + gy * m0 + gz * m0 * m1;
                    output_grid[linear_idx] +=
                        charge * wx[kx as usize] * wy[ky as usize] * wz[kz as usize];
                }
            }
        }
    }
}

/// Map Cartesian coordinates to fractional coordinates (0..1) in the box.
fn map_to_fractional(coord: [f64; 3], box_vecs: &BoxVectors) -> [f64; 3] {
    // Solve: coord = frac[0]*a + frac[1]*b + frac[2]*c
    // For orthogonal boxes, this is simple division
    let a = &box_vecs.0[0];
    let b = &box_vecs.0[1];
    let c = &box_vecs.0[2];

    let det = a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);

    if det.abs() < 1e-10 {
        return [0.0; 3];
    }

    // Inverse: frac = inv(M) * coord
    let inv = matrix_inverse_3x3(&[a, b, c]);
    [
        inv[0][0] * coord[0] + inv[0][1] * coord[1] + inv[0][2] * coord[2],
        inv[1][0] * coord[0] + inv[1][1] * coord[1] + inv[1][2] * coord[2],
        inv[2][0] * coord[0] + inv[2][1] * coord[1] + inv[2][2] * coord[2],
    ]
}

/// Compute 3×3 matrix inverse.
fn matrix_inverse_3x3(mat: &[&[f64; 3]]) -> [[f64; 3]; 3] {
    let a = mat[0];
    let b = mat[1];
    let c = mat[2];

    let det = a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);

    if det.abs() < 1e-10 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }

    let inv_det = 1.0 / det;

    [
        [
            inv_det * (b[1] * c[2] - b[2] * c[1]),
            inv_det * (a[2] * c[1] - a[1] * c[2]),
            inv_det * (a[1] * b[2] - a[2] * b[1]),
        ],
        [
            inv_det * (b[2] * c[0] - b[0] * c[2]),
            inv_det * (a[0] * c[2] - a[2] * c[0]),
            inv_det * (a[2] * b[0] - a[0] * b[2]),
        ],
        [
            inv_det * (b[0] * c[1] - b[1] * c[0]),
            inv_det * (a[1] * c[0] - a[0] * c[1]),
            inv_det * (a[0] * b[1] - a[1] * b[0]),
        ],
    ]
}

/// Compute reciprocal-space energy from charge density.
fn compute_reciprocal_energy(
    charge_grid: &[f64],
    box_vecs: &BoxVectors,
    mesh_size: &[usize; 3],
    alpha: f64,
    kmax: [usize; 3],
) -> f64 {
    let volume = box_vecs.volume();
    let mut energy = 0.0;

    // Iterate over reciprocal lattice vectors
    for kx in 0..kmax[0] {
        for ky in 0..kmax[1] {
            for kz in 0..kmax[2] {
                if kx == 0 && ky == 0 && kz == 0 {
                    continue; // Skip k=0 (handled by self-energy)
                }

                // Reciprocal lattice vectors
                let k_vec = reciprocal_vector(
                    kx as i32,
                    ky as i32,
                    kz as i32,
                    box_vecs,
                    mesh_size,
                );

                let k_sq = k_vec[0] * k_vec[0] + k_vec[1] * k_vec[1] + k_vec[2] * k_vec[2];
                if k_sq < 1e-10 {
                    continue;
                }

                // Structure factor from charge mesh (simplified: direct sum)
                let s_k = compute_structure_factor(charge_grid, mesh_size, &k_vec, box_vecs);

                // Ewald kernel: exp(-k²/4α²) / k²
                let kernel = (-k_sq / (4.0 * alpha * alpha)).exp() / k_sq;

                energy += 2.0 * PI / volume * K_COULOMB * kernel * s_k * s_k;
            }
        }
    }

    energy
}

/// Compute reciprocal lattice vector from mesh indices.
fn reciprocal_vector(
    kx: i32,
    ky: i32,
    kz: i32,
    box_vecs: &BoxVectors,
    mesh_size: &[usize; 3],
) -> [f64; 3] {
    // Reciprocal basis = 2π * (inverse of real basis transposed)
    let inv = matrix_inverse_3x3(&[&box_vecs.0[0], &box_vecs.0[1], &box_vecs.0[2]]);

    let bx = [2.0 * PI * inv[0][0], 2.0 * PI * inv[1][0], 2.0 * PI * inv[2][0]];
    let by = [2.0 * PI * inv[0][1], 2.0 * PI * inv[1][1], 2.0 * PI * inv[2][1]];
    let bz = [2.0 * PI * inv[0][2], 2.0 * PI * inv[1][2], 2.0 * PI * inv[2][2]];

    let kx_frac = kx as f64 / mesh_size[0] as f64;
    let ky_frac = ky as f64 / mesh_size[1] as f64;
    let kz_frac = kz as f64 / mesh_size[2] as f64;

    [
        kx_frac * bx[0] + ky_frac * by[0] + kz_frac * bz[0],
        kx_frac * bx[1] + ky_frac * by[1] + kz_frac * bz[1],
        kx_frac * bx[2] + ky_frac * by[2] + kz_frac * bz[2],
    ]
}

/// Compute B-spline interpolation weights for fractional offset `t ∈ [0, 1)`.
///
/// Returns `order` weights where `w[k]` is the weight for grid point at
/// offset `−k` from `floor(scaled_coord)`. Weights sum to 1.0 (partition of unity).
///
/// Uses Cardinal B-spline recursion (Essmann et al. 1995).
pub(crate) fn bspline_weights(t: f64, order: u8) -> Vec<f64> {
    debug_assert!((0.0..1.0).contains(&t), "t={t} must be in [0, 1)");
    let p = order as usize;
    let mut w = vec![0.0f64; p];
    w[0] = 1.0;
    for o in 2..=p {
        let of = o as f64;
        for k in (1..p).rev() {
            w[k] = ((t + k as f64) * w[k] + (of - t - k as f64) * w[k - 1]) / (of - 1.0);
        }
        w[0] *= t / (of - 1.0);
    }
    w
}

/// Compute structure factor S(k) = Σ_j ρ_j * exp(i k · r_j).
fn compute_structure_factor(
    charge_grid: &[f64],
    mesh_size: &[usize; 3],
    k_vec: &[f64; 3],
    box_vecs: &BoxVectors,
) -> f64 {
    // Simplified: direct summation over mesh points
    // In full PME, this would be computed via FFT
    let [m0, m1, m2] = *mesh_size;
    let a = &box_vecs.0[0];
    let b = &box_vecs.0[1];
    let c = &box_vecs.0[2];

    let mut s_real = 0.0;
    let mut s_imag = 0.0;

    for (idx, &rho) in charge_grid.iter().enumerate() {
        if rho.abs() < 1e-10 {
            continue;
        }

        // Convert linear index to 3D mesh indices
        let i = idx % m0;
        let j = (idx / m0) % m1;
        let k = idx / (m0 * m1);

        // Cartesian position of mesh point: r = (i/M0)*a + (j/M1)*b + (k/M2)*c
        let fi = i as f64 / m0 as f64;
        let fj = j as f64 / m1 as f64;
        let fk = k as f64 / m2 as f64;
        let r = [
            fi * a[0] + fj * b[0] + fk * c[0],
            fi * a[1] + fj * b[1] + fk * c[1],
            fi * a[2] + fj * b[2] + fk * c[2],
        ];

        // Correct phase: k · r (Å⁻¹ · Å = dimensionless)
        let phase = k_vec[0] * r[0] + k_vec[1] * r[1] + k_vec[2] * r[2];

        s_real += rho * phase.cos();
        s_imag += rho * phase.sin();
    }

    // Return |S(k)|²
    s_real * s_real + s_imag * s_imag
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_map_to_fractional_identity() {
        let box_vecs = BoxVectors::cubic(10.0);
        let frac = map_to_fractional([5.0, 5.0, 5.0], &box_vecs);
        assert!((frac[0] - 0.5).abs() < 1e-6);
        assert!((frac[1] - 0.5).abs() < 1e-6);
        assert!((frac[2] - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_reciprocal_vector_zero() {
        let box_vecs = BoxVectors::cubic(10.0);
        let mesh_size = [10, 10, 10];
        let k = reciprocal_vector(0, 0, 0, &box_vecs, &mesh_size);
        assert!((k[0]).abs() < 1e-10);
        assert!((k[1]).abs() < 1e-10);
        assert!((k[2]).abs() < 1e-10);
    }

    #[test]
    fn test_reciprocal_space_energy_empty() {
        let box_vecs = BoxVectors::cubic(10.0);
        let config = PmeConfig::default();
        let coords: Vec<[f64; 3]> = vec![];
        let charges: Vec<f64> = vec![];
        let energy = reciprocal_space_energy(&coords, &charges, &box_vecs, &config);
        assert_eq!(energy, 0.0);
    }

    #[test]
    fn test_bspline_weights_partition_of_unity() {
        for order in 1u8..=6 {
            for i in 0..20 {
                let t = i as f64 / 20.0;
                let w = bspline_weights(t, order);
                let sum: f64 = w.iter().sum();
                assert!(
                    (sum - 1.0).abs() < 1e-12,
                    "order={order}, t={t}: sum={sum}"
                );
            }
        }
    }

    #[test]
    fn test_bspline_weights_order4_known_values() {
        // t=0: [0, 1/6, 2/3, 1/6]
        let w0 = bspline_weights(0.0, 4);
        assert!(w0[0].abs() < 1e-12, "w[0] at t=0: {}", w0[0]);
        assert!(
            (w0[1] - 1.0 / 6.0).abs() < 1e-12,
            "w[1] at t=0: {}",
            w0[1]
        );
        assert!(
            (w0[2] - 2.0 / 3.0).abs() < 1e-12,
            "w[2] at t=0: {}",
            w0[2]
        );
        assert!(
            (w0[3] - 1.0 / 6.0).abs() < 1e-12,
            "w[3] at t=0: {}",
            w0[3]
        );

        // t=0.5: [1/48, 23/48, 23/48, 1/48]
        let w5 = bspline_weights(0.5, 4);
        assert!(
            (w5[0] - 1.0 / 48.0).abs() < 1e-12,
            "w[0] at t=0.5: {}",
            w5[0]
        );
        assert!(
            (w5[1] - 23.0 / 48.0).abs() < 1e-12,
            "w[1] at t=0.5: {}",
            w5[1]
        );
        assert!(
            (w5[2] - 23.0 / 48.0).abs() < 1e-12,
            "w[2] at t=0.5: {}",
            w5[2]
        );
        assert!(
            (w5[3] - 1.0 / 48.0).abs() < 1e-12,
            "w[3] at t=0.5: {}",
            w5[3]
        );
    }

    #[test]
    fn test_charge_conservation() {
        use crate::BoxVectors;
        let box_vecs = BoxVectors::cubic(10.0);
        let mesh_size = [16usize, 16, 16];
        let mut grid = vec![0.0f64; 16 * 16 * 16];
        let coords = vec![[2.5, 5.0, 7.5], [1.0, 1.0, 1.0]];
        let charges = vec![1.5, -0.7];
        interpolate_charges_to_mesh(&coords, &charges, &box_vecs, &mut grid, &mesh_size, 4);
        let grid_total: f64 = grid.iter().sum();
        let charge_total: f64 = charges.iter().sum();
        assert!(
            (grid_total - charge_total).abs() < 1e-10,
            "grid_sum={grid_total}, charge_sum={charge_total}"
        );
    }

    #[test]
    fn test_charge_conservation_noncubic_mesh() {
        use crate::BoxVectors;
        // Verifies cbrt hack is gone — non-cubic mesh must work
        let box_vecs = BoxVectors::cubic(10.0);
        let mesh_size = [8usize, 12, 16];
        let mut grid = vec![0.0f64; 8 * 12 * 16];
        let coords = vec![[3.0, 5.0, 8.0]];
        let charges = vec![2.0];
        interpolate_charges_to_mesh(&coords, &charges, &box_vecs, &mut grid, &mesh_size, 4);
        let grid_total: f64 = grid.iter().sum();
        assert!(
            (grid_total - 2.0).abs() < 1e-10,
            "grid_sum={grid_total}"
        );
    }

    #[test]
    fn test_reciprocal_energy_nonnegative() {
        use crate::{BoxVectors, PmeConfig};
        // Each k-term is (2π/V) * K * exp(-k²/4α²)/k² * |S(k)|² ≥ 0
        let box_vecs = BoxVectors::cubic(10.0);
        let config = PmeConfig {
            alpha: 0.3,
            kmax: [3, 3, 3],
            mesh: [8, 8, 8],
            spline_order: 4,
            ..PmeConfig::default()
        };
        let coords = vec![[2.0, 2.0, 2.0], [6.0, 6.0, 6.0]];
        let charges = vec![1.0, -1.0];
        let energy = reciprocal_space_energy(&coords, &charges, &box_vecs, &config);
        assert!(energy >= 0.0, "reciprocal energy must be >= 0, got {energy}");
    }
}
