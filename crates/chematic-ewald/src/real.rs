//! Real-space Coulomb energy calculations with optional cutoff.

pub const K_COULOMB: f64 = 332.0637; // kcal·Å/(mol·e²)

/// Direct Coulomb sum without cutoff (non-periodic).
///
/// Computes the sum of all pairwise Coulomb interactions:
/// E = K_e * Σ_{i<j} q_i * q_j / r_ij
///
/// # Arguments
/// * `coords` - Atomic coordinates as flat array [x0, y0, z0, x1, y1, z1, ...]
/// * `charges` - Partial charges (in units of electron charge)
///
/// # Returns
/// Total Coulomb energy in kcal/mol
pub fn direct_coulomb(coords: &[[f64; 3]], charges: &[f64]) -> f64 {
    let mut energy = 0.0;
    let n = coords.len();

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r > 1e-6 {
                energy += K_COULOMB * charges[i] * charges[j] / r;
            }
        }
    }

    energy
}

/// Direct Coulomb sum with real-space cutoff.
///
/// Only includes pairwise interactions where r_ij < r_cut.
///
/// # Arguments
/// * `coords` - Atomic coordinates [[x0, y0, z0], [x1, y1, z1], ...]
/// * `charges` - Partial charges
/// * `r_cut` - Real-space cutoff in Ångströms
///
/// # Returns
/// Coulomb energy for pairs within cutoff (kcal/mol)
pub fn direct_coulomb_cutoff(coords: &[[f64; 3]], charges: &[f64], r_cut: f64) -> f64 {
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

/// Coulomb energy with damped potential for short-range overlap.
///
/// Uses complementary error function erfc for smooth damping:
/// E = K_e * Σ_{i<j} q_i * q_j * erfc(α * r_ij) / r_ij
///
/// where α is the Ewald splitting parameter.
///
/// # Arguments
/// * `coords` - Atomic coordinates [[x0, y0, z0], [x1, y1, z1], ...]
/// * `charges` - Partial charges
/// * `alpha` - Ewald splitting parameter (typical: 3.5 / r_cut)
///
/// # Returns
/// Damped Coulomb energy (kcal/mol)
pub fn direct_coulomb_damped(coords: &[[f64; 3]], charges: &[f64], alpha: f64) -> f64 {
    let mut energy = 0.0;
    let n = coords.len();

    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();

            if r > 1e-6 {
                let ar = alpha * r;
                let erfc_ar = complementary_error_function(ar);
                energy += K_COULOMB * charges[i] * charges[j] * erfc_ar / r;
            }
        }
    }

    energy
}

/// Complementary error function erfc(x) = 1 - erf(x).
///
/// Approximation using Abramowitz and Stegun (1964) formula 7.1.26 for erf(x),
/// then erfc(x) = 1 - erf(x).
fn complementary_error_function(x: f64) -> f64 {
    const A1: f64 = 0.254829592;
    const A2: f64 = -0.284496736;
    const A3: f64 = 1.421413741;
    const A4: f64 = -1.453152027;
    const A5: f64 = 1.061405429;
    const P: f64 = 0.3275911;

    let x_abs = x.abs();
    let t = 1.0 / (1.0 + P * x_abs);

    // erf(x) approximation
    let erf = 1.0 - (((((A5 * t + A4) * t + A3) * t + A2) * t + A1) * t * (-x_abs * x_abs).exp());

    // erfc(x) = 1 - erf(x), but erf(-x) = -erf(x)
    if x >= 0.0 {
        1.0 - erf
    } else {
        1.0 + erf  // erf(x) = -erf(-x), so erfc(x) = 1 - erf(x) = 1 + erf(-x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erfc_zero() {
        let val = complementary_error_function(0.0);
        assert!((val - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_erfc_large() {
        let val = complementary_error_function(5.0);
        assert!(val < 1e-6); // erfc(5) ≈ 0
    }

    #[test]
    fn test_erfc_symmetry() {
        let p = complementary_error_function(2.0);
        let n = complementary_error_function(-2.0);
        assert!((p + n - 2.0).abs() < 1e-6); // erfc(x) + erfc(-x) = 2
    }
}
