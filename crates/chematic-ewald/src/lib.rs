//! Smooth Particle Mesh Ewald (SPME) for long-range electrostatic interactions.
//!
//! Provides efficient computation of long-range Coulomb forces in periodic and non-periodic systems.
//! Uses FFT-based reciprocal space summation combined with real-space cutoff.

#![forbid(unsafe_code)]

use std::fmt;

pub mod pme;
pub mod real;

pub use pme::reciprocal_space_energy;
pub use real::{
    K_COULOMB, direct_coulomb, direct_coulomb_cutoff, direct_coulomb_damped,
    direct_coulomb_damped_cutoff,
};

/// Error type for SPME calculations.
#[derive(Debug, Clone, PartialEq)]
pub enum EwaldError {
    /// The periodic box or PME configuration is invalid or exceeds the
    /// implementation's bounded work limits.
    ///
    /// The original public API exposed only this variant. It remains the
    /// single variant so rejecting newly validated inputs does not break
    /// exhaustive matches in downstream Rust code.
    SingularBoxMatrix,
}

impl fmt::Display for EwaldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EwaldError::SingularBoxMatrix => {
                write!(
                    f,
                    "invalid PME input, configuration, or singular box matrix"
                )
            }
        }
    }
}

impl std::error::Error for EwaldError {}

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

pub(crate) const MAX_PME_GRID_POINTS: usize = 16_777_216;
pub(crate) const MAX_PME_RECIPROCAL_POINTS: usize = 2_000_000;
pub(crate) const MAX_PME_SPLINE_ORDER: u8 = 16;

pub(crate) fn validate_pme_inputs(
    coords: &[[f64; 3]],
    charges: &[f64],
    box_vecs: &BoxVectors,
    config: &PmeConfig,
) -> Result<usize, EwaldError> {
    if coords.len() != charges.len() {
        return Err(EwaldError::SingularBoxMatrix);
    }
    if coords
        .iter()
        .flatten()
        .chain(charges.iter())
        .chain(box_vecs.0.iter().flatten())
        .any(|value| !value.is_finite())
    {
        return Err(EwaldError::SingularBoxMatrix);
    }
    let volume = box_vecs.volume();
    if !volume.is_finite() || volume < 1e-10 {
        return Err(EwaldError::SingularBoxMatrix);
    }
    if !config.r_cut.is_finite() || config.r_cut <= 0.0 {
        return Err(EwaldError::SingularBoxMatrix);
    }
    if !config.alpha.is_finite() || config.alpha < 0.0 {
        return Err(EwaldError::SingularBoxMatrix);
    }
    if config.mesh.contains(&0) {
        return Err(EwaldError::SingularBoxMatrix);
    }
    if config.spline_order == 0 || config.spline_order > MAX_PME_SPLINE_ORDER {
        return Err(EwaldError::SingularBoxMatrix);
    }
    let grid_points = config
        .mesh
        .iter()
        .try_fold(1usize, |product, &size| product.checked_mul(size))
        .ok_or(EwaldError::SingularBoxMatrix)?;
    if grid_points > MAX_PME_GRID_POINTS {
        return Err(EwaldError::SingularBoxMatrix);
    }
    let reciprocal_points = config
        .kmax
        .iter()
        .enumerate()
        .try_fold(1usize, |product, (axis, &maximum)| {
            let width = if axis == 2 {
                maximum.checked_add(1)
            } else {
                maximum.checked_mul(2)?.checked_add(1)
            }?;
            product.checked_mul(width)
        })
        .ok_or(EwaldError::SingularBoxMatrix)?;
    if reciprocal_points > MAX_PME_RECIPROCAL_POINTS
        || config
            .kmax
            .iter()
            .any(|&maximum| maximum > i32::MAX as usize)
    {
        return Err(EwaldError::SingularBoxMatrix);
    }
    Ok(grid_points)
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
        (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            .abs()
    }
}

/// Result of Ewald energy calculation.
#[derive(Clone, Debug, PartialEq)]
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
/// `Result<EwaldResult, EwaldError>` with breakdown:
/// - real_energy: Short-range damped Coulomb (r < r_cut)
/// - reciprocal_energy: Long-range FFT-based contribution
/// - self_energy: Self-interaction correction
/// - total_energy: Sum of all contributions
///
/// # Errors
/// Returns `EwaldError::SingularBoxMatrix` if the box matrix is singular.
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
) -> Result<EwaldResult, EwaldError> {
    validate_pme_inputs(coords, charges, box_vecs, config)?;

    let alpha = if config.alpha > 0.0 {
        config.alpha
    } else {
        3.5 / config.r_cut
    };

    // Real-space: short-range damped Coulomb
    let real = real::direct_coulomb_damped_cutoff(coords, charges, alpha, config.r_cut);

    // Reciprocal-space: FFT-based long-range
    let reciprocal = pme::reciprocal_space_energy(coords, charges, box_vecs, config)?;

    // Self-energy: correction for Gaussian charge distribution
    let self_corr = compute_self_energy(charges, alpha);

    Ok(EwaldResult {
        real_energy: real,
        reciprocal_energy: reciprocal,
        self_energy: self_corr,
        total_energy: real + reciprocal + self_corr,
    })
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
        assert!(
            (energy - expected).abs() < 1.0,
            "energy={}, expected={}",
            energy,
            expected
        );
    }

    #[test]
    fn test_spme_singular_box_returns_error() {
        let degenerate_box = BoxVectors([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let coords = [[0.0, 0.0, 0.0]];
        let charges = [1.0];
        let result = spme_energy(&coords, &charges, &degenerate_box, &PmeConfig::default());
        assert_eq!(result, Err(EwaldError::SingularBoxMatrix));
    }

    #[test]
    fn test_spme_valid_box() {
        let box_vecs = BoxVectors::cubic(10.0);
        let coords = [[5.0, 5.0, 5.0]];
        let charges = [1.0];
        let result = spme_energy(&coords, &charges, &box_vecs, &PmeConfig::default());
        assert!(result.is_ok());
    }

    #[test]
    fn spme_rejects_coordinate_charge_length_mismatch() {
        let result = spme_energy(
            &[[0.0, 0.0, 0.0]],
            &[],
            &BoxVectors::cubic(10.0),
            &PmeConfig::default(),
        );
        assert_eq!(result, Err(EwaldError::SingularBoxMatrix));
    }

    #[test]
    fn spme_rejects_zero_mesh_without_panicking() {
        let config = PmeConfig {
            mesh: [0, 16, 16],
            ..PmeConfig::default()
        };
        let result = spme_energy(
            &[[0.0, 0.0, 0.0]],
            &[1.0],
            &BoxVectors::cubic(10.0),
            &config,
        );
        assert_eq!(result, Err(EwaldError::SingularBoxMatrix));
    }

    #[test]
    fn spme_rejects_excessive_grid_before_allocation() {
        let config = PmeConfig {
            mesh: [usize::MAX, 2, 2],
            ..PmeConfig::default()
        };
        let result = spme_energy(
            &[[0.0, 0.0, 0.0]],
            &[1.0],
            &BoxVectors::cubic(10.0),
            &config,
        );
        assert_eq!(result, Err(EwaldError::SingularBoxMatrix));
    }

    #[test]
    fn spme_rejects_zero_spline_order_without_panicking() {
        let config = PmeConfig {
            spline_order: 0,
            ..PmeConfig::default()
        };
        let result = spme_energy(
            &[[0.0, 0.0, 0.0]],
            &[1.0],
            &BoxVectors::cubic(10.0),
            &config,
        );
        assert_eq!(result, Err(EwaldError::SingularBoxMatrix));
    }

    #[test]
    fn spme_rejects_non_finite_and_invalid_numeric_inputs() {
        let box_vecs = BoxVectors::cubic(10.0);
        for (coords, charges) in [
            (vec![[f64::NAN, 0.0, 0.0]], vec![1.0]),
            (vec![[0.0, 0.0, 0.0]], vec![f64::INFINITY]),
        ] {
            assert_eq!(
                spme_energy(&coords, &charges, &box_vecs, &PmeConfig::default()),
                Err(EwaldError::SingularBoxMatrix)
            );
        }

        let invalid_box = BoxVectors([[f64::NAN, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]);
        assert_eq!(
            spme_energy(&[], &[], &invalid_box, &PmeConfig::default()),
            Err(EwaldError::SingularBoxMatrix)
        );

        for config in [
            PmeConfig {
                r_cut: 0.0,
                ..PmeConfig::default()
            },
            PmeConfig {
                r_cut: f64::NAN,
                ..PmeConfig::default()
            },
            PmeConfig {
                alpha: -1.0,
                ..PmeConfig::default()
            },
            PmeConfig {
                alpha: f64::INFINITY,
                ..PmeConfig::default()
            },
            PmeConfig {
                spline_order: MAX_PME_SPLINE_ORDER + 1,
                ..PmeConfig::default()
            },
        ] {
            assert_eq!(
                spme_energy(&[], &[], &box_vecs, &config),
                Err(EwaldError::SingularBoxMatrix)
            );
        }
    }

    #[test]
    fn spme_rejects_excessive_reciprocal_work_before_iteration() {
        for kmax in [
            [usize::MAX, 1, 1],
            [1, 1, usize::MAX],
            [1_000, 1_000, 1_000],
            [i32::MAX as usize + 1, 0, 0],
        ] {
            let config = PmeConfig {
                kmax,
                ..PmeConfig::default()
            };
            assert_eq!(
                spme_energy(&[], &[], &BoxVectors::cubic(10.0), &config),
                Err(EwaldError::SingularBoxMatrix)
            );
        }
    }
}
