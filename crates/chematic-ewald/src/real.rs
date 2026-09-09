//! Real-space Coulomb energy calculations with optional cutoff.

use std::collections::HashMap;

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
    // A positive finite cutoff permits a conservative uniform-cell traversal.
    // Keep the original all-pair path for invalid coordinates/cutoffs so NaN
    // propagation and unusual caller inputs are unchanged.
    if r_cut.is_finite()
        && r_cut > 0.0
        && coords
            .iter()
            .all(|coord| coord.iter().all(|value| value.is_finite()))
    {
        let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
        for (index, coord) in coords.iter().enumerate() {
            let key = [
                (coord[0] / r_cut).floor() as i64,
                (coord[1] / r_cut).floor() as i64,
                (coord[2] / r_cut).floor() as i64,
            ];
            cells.entry(key).or_default().push(index);
        }

        let mut energy = 0.0;
        for (i, coord_i) in coords.iter().enumerate() {
            let key = [
                (coord_i[0] / r_cut).floor() as i64,
                (coord_i[1] / r_cut).floor() as i64,
                (coord_i[2] / r_cut).floor() as i64,
            ];
            for dx_cell in -1..=1 {
                for dy_cell in -1..=1 {
                    for dz_cell in -1..=1 {
                        let neighbor = [key[0] + dx_cell, key[1] + dy_cell, key[2] + dz_cell];
                        let Some(indices) = cells.get(&neighbor) else {
                            continue;
                        };
                        for &j in indices {
                            if j <= i {
                                continue;
                            }
                            let dx = coord_i[0] - coords[j][0];
                            let dy = coord_i[1] - coords[j][1];
                            let dz = coord_i[2] - coords[j][2];
                            let r = (dx * dx + dy * dy + dz * dz).sqrt();
                            // Preserve the public function's strict cutoff.
                            if r > 1e-6 && r < r_cut {
                                energy += K_COULOMB * charges[i] * charges[j] / r;
                            }
                        }
                    }
                }
            }
        }
        return energy;
    }

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

/// Damped Coulomb energy using a validated uniform cell list.
///
/// This is the real-space path used by PME.  A cell has an edge length equal
/// to `r_cut`, so only the 27 neighboring cells can contain an interacting
/// pair.  The strict `r < r_cut` boundary and the overlap guard intentionally
/// match [`direct_coulomb_cutoff`].  Non-finite inputs fall back to the direct
/// implementation instead of being used as hash-map coordinates.
pub fn direct_coulomb_damped_cutoff(
    coords: &[[f64; 3]],
    charges: &[f64],
    alpha: f64,
    r_cut: f64,
) -> f64 {
    if !r_cut.is_finite() || r_cut <= 0.0 || coords.iter().any(|c| c.iter().any(|v| !v.is_finite()))
    {
        return direct_coulomb_damped(coords, charges, alpha);
    }

    let mut cells: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (index, coord) in coords.iter().enumerate() {
        let key = [
            (coord[0] / r_cut).floor() as i64,
            (coord[1] / r_cut).floor() as i64,
            (coord[2] / r_cut).floor() as i64,
        ];
        cells.entry(key).or_default().push(index);
    }

    let mut energy = 0.0;
    for (i, coord_i) in coords.iter().enumerate() {
        let key = [
            (coord_i[0] / r_cut).floor() as i64,
            (coord_i[1] / r_cut).floor() as i64,
            (coord_i[2] / r_cut).floor() as i64,
        ];
        for dx_cell in -1..=1 {
            for dy_cell in -1..=1 {
                for dz_cell in -1..=1 {
                    let neighbor = [key[0] + dx_cell, key[1] + dy_cell, key[2] + dz_cell];
                    let Some(indices) = cells.get(&neighbor) else {
                        continue;
                    };
                    for &j in indices {
                        if j <= i {
                            continue;
                        }
                        let dx = coord_i[0] - coords[j][0];
                        let dy = coord_i[1] - coords[j][1];
                        let dz = coord_i[2] - coords[j][2];
                        let r = (dx * dx + dy * dy + dz * dz).sqrt();
                        if r > 1e-6 && r < r_cut {
                            let ar = alpha * r;
                            energy += K_COULOMB
                                * charges[i]
                                * charges[j]
                                * complementary_error_function(ar)
                                / r;
                        }
                    }
                }
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
        1.0 + erf // erf(x) = -erf(-x), so erfc(x) = 1 - erf(x) = 1 + erf(-x)
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

    #[test]
    fn cell_list_matches_all_pair_cutoff_including_boundaries() {
        let coords = [
            [-4.0_f64, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, -4.0],
            [7.2, 0.0, 0.0],
        ];
        let charges = [1.0, -0.5, 0.75, -1.0, 0.25, 0.5];
        let r_cut = 4.0;
        let mut expected = 0.0;
        for i in 0..coords.len() {
            for j in (i + 1)..coords.len() {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r > 1e-6 && r < r_cut {
                    expected += K_COULOMB * charges[i] * charges[j] / r;
                }
            }
        }
        let actual = direct_coulomb_cutoff(&coords, &charges, r_cut);
        assert!((actual - expected).abs() < 1e-10);
    }

    #[test]
    fn cell_list_matches_direct_damped_cutoff() {
        let coords = [
            [-2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
            [9.0, 0.0, 0.0],
            [0.0, 2.9, 0.0],
            [0.0, 0.0, -3.1],
        ];
        let charges = [1.0, -0.5, 0.75, -1.0, 0.25, 0.5];
        let alpha = 3.5 / 3.0;
        let expected = direct_coulomb_cutoff(&coords, &charges, 3.0);
        let actual = direct_coulomb_damped_cutoff(&coords, &charges, alpha, 3.0);
        let direct = direct_coulomb_damped(&coords, &charges, alpha);
        assert!(actual.is_finite());
        assert!(direct.is_finite());
        assert!((actual - damped_reference(&coords, &charges, alpha, 3.0)).abs() < 1e-10);
        assert!(expected.is_finite());
    }

    #[test]
    fn cell_list_handles_negative_coordinates_and_boundary() {
        let coords = [[-3.0, 0.0, 0.0], [0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0, 1.0];
        let alpha = 1.0;
        let actual = direct_coulomb_damped_cutoff(&coords, &charges, alpha, 3.0);
        let expected = damped_reference(&coords, &charges, alpha, 3.0);
        assert!((actual - expected).abs() < 1e-10);
    }

    fn damped_reference(coords: &[[f64; 3]], charges: &[f64], alpha: f64, r_cut: f64) -> f64 {
        let mut energy = 0.0;
        for i in 0..coords.len() {
            for j in i + 1..coords.len() {
                let dx = coords[i][0] - coords[j][0];
                let dy = coords[i][1] - coords[j][1];
                let dz = coords[i][2] - coords[j][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r > 1e-6 && r < r_cut {
                    energy += K_COULOMB
                        * charges[i]
                        * charges[j]
                        * complementary_error_function(alpha * r)
                        / r;
                }
            }
        }
        energy
    }
}
