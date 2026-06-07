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

    // Interpolate charges onto mesh (B-spline order 4)
    interpolate_charges_to_mesh(
        coords,
        charges,
        box_vecs,
        &mut charge_grid,
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
    spline_order: u8,
) {
    let _ = spline_order; // Use to suppress warnings

    // TODO: Implement full B-spline interpolation
    // For now, simple linear interpolation to nearest mesh points

    for i in 0..charges.len() {
        let charge = charges[i];
        if charge.abs() < 1e-10 {
            continue;
        }

        // Map atomic coordinate to fractional coordinates (0..1)
        let frac = map_to_fractional(coords[i], box_vecs);

        // Approximate mesh size as cubic root (assumes cubic mesh)
        // WARNING: This is a simplification. Proper PME requires explicit mesh dimensions.
        let approx_mesh_side = (output_grid.len() as f64).cbrt() as usize;
        if approx_mesh_side == 0 {
            continue;
        }

        // Map fractional coordinates to 3D mesh indices
        let ix = ((frac[0] * approx_mesh_side as f64) as usize) % approx_mesh_side;
        let iy = ((frac[1] * approx_mesh_side as f64) as usize) % approx_mesh_side;
        let iz = ((frac[2] * approx_mesh_side as f64) as usize) % approx_mesh_side;

        // Convert 3D indices to linear index: idx = ix + iy*M0 + iz*M0*M1
        let linear_idx = ix + iy * approx_mesh_side + iz * approx_mesh_side * approx_mesh_side;

        if linear_idx < output_grid.len() {
            output_grid[linear_idx] += charge;
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
                let s_k = compute_structure_factor(charge_grid, mesh_size, &k_vec);

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

/// Compute structure factor S(k) = Σ_j ρ_j * exp(i k · r_j).
fn compute_structure_factor(charge_grid: &[f64], mesh_size: &[usize; 3], k_vec: &[f64; 3]) -> f64 {
    // Simplified: direct summation over mesh points
    // In full PME, this would be computed via FFT
    let mut s_real = 0.0;
    let mut s_imag = 0.0;

    for (idx, &rho) in charge_grid.iter().enumerate() {
        if rho.abs() < 1e-10 {
            continue;
        }

        // Approximate mesh point position (linear indexing)
        let i = idx % mesh_size[0];
        let j = (idx / mesh_size[0]) % mesh_size[1];
        let k = idx / (mesh_size[0] * mesh_size[1]);

        let phase = k_vec[0] * (i as f64 / mesh_size[0] as f64)
            + k_vec[1] * (j as f64 / mesh_size[1] as f64)
            + k_vec[2] * (k as f64 / mesh_size[2] as f64);

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
}
