//! ETKDG (Experimental Torsion Knowledge Distance Geometry) v3 conformer generation.
//!
//! This module augments rule-based 3D coordinate generation with experimental
//! torsion angle preferences derived from crystal structure databases.
//!
//! Algorithm:
//! 1. Generate initial 3D coordinates via rule-based DG
//! 2. Apply preferred torsion angles for common patterns (carbonyl, aromatic, etc.)
//! 3. Minimize with constraints to preserve ring geometry

use crate::coords::Coords3D;
use crate::etkdg_knowledge::{default_torsion_preference, get_torsion_preference};
use chematic_core::{AtomIdx, Molecule};
use fastrand;

/// Generate 3D coordinates using ETKDG with torsion angle preferences.
///
/// This function improves upon basic rule-based 3D generation by applying
/// experimental torsion angle preferences to common structural patterns.
pub fn generate_coords_etkdg(mol: &Molecule) -> Coords3D {
    generate_coords_etkdg_with_noise(mol, 0.0)
}

/// Generate 3D coordinates using ETKDG with optional torsion noise for
/// conformer ensemble diversity.
///
/// `noise_sigma_deg`: if > 0, adds uniform random noise in
/// `[-noise_sigma_deg, +noise_sigma_deg]` to each preferred torsion angle
/// before applying it.  Use 0.0 for deterministic single-conformer generation;
/// use ~30° for ensemble sampling to break out of local minima.
pub fn generate_coords_etkdg_with_noise(mol: &Molecule, noise_sigma_deg: f64) -> Coords3D {
    let mut coords = super::dg::generate_coords(mol);

    if mol.atom_count() < 4 {
        return coords;
    }

    apply_torsion_preferences_with_noise(mol, &mut coords, noise_sigma_deg);

    let constraints = super::constraints::build_constraints(mol);
    coords = super::constraints::satisfy_constraints(&coords, mol, &constraints, 3);

    coords
}

/// Apply torsion angle preferences to coordinates using the knowledge base.
///
/// When `noise_sigma_deg > 0`, uniform random noise in
/// `[-noise_sigma_deg, +noise_sigma_deg]` is added to each preferred angle to
/// generate diverse conformers for ensemble sampling.
fn apply_torsion_preferences_with_noise(
    mol: &Molecule,
    coords: &mut Coords3D,
    noise_sigma_deg: f64,
) {
    let n = mol.atom_count();
    let mut applied = std::collections::HashSet::new();

    // Scan for 4-atom chains (A-B-C-D torsions)
    for b in 0..n {
        for c in 0..n {
            if b == c {
                continue;
            }

            let b_idx = AtomIdx(b as u32);
            let c_idx = AtomIdx(c as u32);
            if mol.bond_between(b_idx, c_idx).is_none() {
                continue;
            }

            let b_neighbors: Vec<usize> = mol
                .neighbors(b_idx)
                .filter(|(nb, _)| nb.0 as usize != c)
                .map(|(nb, _)| nb.0 as usize)
                .collect();

            let c_neighbors: Vec<usize> = mol
                .neighbors(c_idx)
                .filter(|(nb, _)| nb.0 as usize != b)
                .map(|(nb, _)| nb.0 as usize)
                .collect();

            for &a in &b_neighbors {
                for &d in &c_neighbors {
                    let key = (a.min(d), a.max(d), b.min(c), b.max(c));
                    if applied.contains(&key) {
                        continue;
                    }

                    let a_idx = AtomIdx(a as u32);
                    let d_idx = AtomIdx(d as u32);

                    let current =
                        super::mol_transforms::get_dihedral(coords, a_idx, b_idx, c_idx, d_idx);
                    if current.is_none() {
                        continue;
                    }

                    let current_deg = current.unwrap() * 180.0 / std::f64::consts::PI;

                    let preference = get_torsion_preference(mol, a_idx, b_idx, c_idx, d_idx)
                        .unwrap_or_else(default_torsion_preference);

                    // Add uniform noise when generating ensemble conformers.
                    let noise = if noise_sigma_deg > 0.0 {
                        (fastrand::f64() * 2.0 - 1.0) * noise_sigma_deg
                    } else {
                        0.0
                    };
                    // Normalize target to [-180, 180] so that preferences near ±180°
                    // plus noise (e.g. 180° + 30° = 210°) stay in the canonical range.
                    let raw = preference.angle_deg + noise;
                    let target_deg = ((raw + 180.0).rem_euclid(360.0)) - 180.0;

                    // Use circular (shortest-path) distance for the threshold check.
                    // Without wrapping, a current of −150° vs target of 210° would
                    // give |diff| = 360° and trigger a spurious full-circle rotation.
                    let diff = {
                        let d = (target_deg - current_deg).rem_euclid(360.0);
                        if d > 180.0 { d - 360.0 } else { d }
                    };

                    // Only apply if difference is significant (> 20°)
                    if diff.abs() > 20.0 {
                        let target_rad = target_deg * std::f64::consts::PI / 180.0;
                        *coords = super::mol_transforms::set_dihedral(
                            coords, mol, a_idx, b_idx, c_idx, d_idx, target_rad,
                        );
                        applied.insert(key);
                    }
                }
            }
        }
    }
}
