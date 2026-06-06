//! Smooth Particle Mesh Ewald (SPME) for long-range electrostatic interactions.
//!
//! Provides efficient computation of long-range Coulomb forces in periodic and non-periodic systems.
//! Uses FFT-based reciprocal space summation combined with real-space cutoff.

#![forbid(unsafe_code)]

pub mod real;
pub mod pme;

pub use real::{direct_coulomb, direct_coulomb_cutoff, direct_coulomb_damped, K_COULOMB};
pub use pme::reciprocal_space_energy;

/// Configuration for PME/SPME calculations.
#[derive(Clone, Debug)]
pub struct PmeConfig {
    /// Ewald splitting parameter α (default: auto-computed from r_cut).
    pub alpha: f64,
    /// Reciprocal space cutoff: kmax per axis (default: [5, 5, 5]).
    pub kmax: [usize; 3],
    /// Real-space cutoff in Ångströms (default: 8.0).
    pub r_cut: f64,
    /// PME grid dimensions (default: [16, 16, 16]).
    pub mesh: [usize; 3],
    /// B-spline interpolation order (default: 4).
    pub spline_order: u8,
}

impl Default for PmeConfig {
    fn default() -> Self {
        Self {
            alpha: 0.0, // Will be auto-computed
            kmax: [5, 5, 5],
            r_cut: 8.0,
            mesh: [16, 16, 16],
            spline_order: 4,
        }
    }
}

/// Periodic box vectors for MD simulation (3x3 matrix, rows are box vectors).
#[derive(Clone, Debug)]
pub struct BoxVectors(pub [[f64; 3]; 3]);

impl BoxVectors {
    /// Create a cubic box with side length `a` Ångströms.
    pub fn cubic(a: f64) -> Self {
        BoxVectors([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    }

    /// Compute box volume (Ų).
    pub fn volume(&self) -> f64 {
        let a = self.0[0];
        let b = self.0[1];
        let c = self.0[2];
        (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            .abs()
    }
}

/// Result of Ewald energy calculation.
#[derive(Clone, Debug)]
pub struct EwaldResult {
    /// Real-space Coulomb energy (kcal/mol).
    pub real_energy: f64,
    /// Reciprocal-space (PME) contribution (kcal/mol).
    pub reciprocal_energy: f64,
    /// Self-energy correction for Gaussian clouds.
    pub self_energy: f64,
    /// Total Ewald energy = real + reciprocal + self.
    pub total_energy: f64,
}


/// SPME Ewald energy for periodic system.
///
/// Computes long-range electrostatic energy using Smooth Particle Mesh Ewald (SPME).
/// Splits calculation into real-space (short-range, direct sum with damping) and
/// reciprocal-space (long-range, FFT-based).
///
/// # Arguments
/// * `coords` - Atomic coordinates [[x0, y0, z0], [x1, y1, z1], ...]
/// * `charges` - Partial charges (e)
/// * `box_vecs` - Periodic box vectors
/// * `config` - PME configuration (alpha, kmax, mesh, spline_order)
///
/// # Returns
/// `EwaldResult` with breakdown:
/// - real_energy: Short-range damped Coulomb (r < r_cut)
/// - reciprocal_energy: Long-range FFT-based contribution
/// - self_energy: Self-interaction correction
/// - total_energy: Sum of all contributions
///
/// # Notes
/// - Requires coordinates to be within the periodic box
/// - Auto-computes alpha if not provided (α = 3.5 / r_cut)
/// - Accurate for N up to several thousand atoms with standard mesh
pub fn spme_energy(
    coords: &[[f64; 3]],
    charges: &[f64],
    box_vecs: &BoxVectors,
    config: &PmeConfig,
) -> EwaldResult {
    let alpha = if config.alpha > 0.0 {
        config.alpha
    } else {
        3.5 / config.r_cut
    };

    // Real-space: short-range damped Coulomb
    let real = real::direct_coulomb_damped(coords, charges, alpha);

    // Reciprocal-space: FFT-based long-range
    let reciprocal = pme::reciprocal_space_energy(coords, charges, box_vecs, config);

    // Self-energy: correction for Gaussian charge distribution
    let self_corr = compute_self_energy(charges, alpha);

    EwaldResult {
        real_energy: real,
        reciprocal_energy: reciprocal,
        self_energy: self_corr,
        total_energy: real + reciprocal + self_corr,
    }
}

/// Helper: Direct Coulomb with cutoff using flat coordinate array.
fn direct_coulomb_cutoff_coords(coords: &[[f64; 3]], charges: &[f64], r_cut: f64) -> f64 {
    let mut energy = 0.0;
    let n = coords.len();

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r > 1e-6 && r < r_cut {
                energy += K_COULOMB * charges[i] * charges[j] / r;
            }
        }
    }

    energy
}

fn compute_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let alpha = if alpha > 0.0 { alpha } else { 3.5 / 8.0 };
    let sqrt_pi = std::f64::consts::PI.sqrt();
    let mut self_e = 0.0;
    for &q in charges {
        self_e -= K_COULOMB * alpha / sqrt_pi * q * q;
    }
    self_e
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_box_volume() {
        let box_vecs = BoxVectors::cubic(10.0);
        assert!((box_vecs.volume() - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_direct_coulomb_flat_coords() {
        // Two unit charges 1 Å apart
        let coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let energy = direct_coulomb(&coords, &charges);
        // Expected: K * 1.0 * (-1.0) / 1.0 = -332.0637 kcal/mol (approx)
        let expected = -K_COULOMB;
        assert!((energy - expected).abs() < 1.0, "energy={}, expected={}", energy, expected);
    }
}
