//! ETKDG (Experimental Torsion Knowledge Distance Geometry) v3 conformer generation.
//!
//! This module augments rule-based 3D coordinate generation with experimental
//! torsion angle preferences derived from crystal structure databases.
//!
//! Algorithm:
//! 1. Generate initial 3D coordinates via rule-based DG
//! 2. Apply preferred torsion angles for common patterns (carbonyl, aromatic, etc.)
//! 3. Minimize with constraints to preserve ring geometry

use chematic_core::{Molecule, AtomIdx};
use crate::coords::Coords3D;
use crate::etkdg_knowledge::{get_torsion_preference, default_torsion_preference};

/// Generate 3D coordinates using ETKDG with torsion angle preferences.
///
/// This function improves upon basic rule-based 3D generation by applying
/// experimental torsion angle preferences to common structural patterns.
pub fn generate_coords_etkdg(mol: &Molecule) -> Coords3D {
    // Step 1: Generate initial coordinates via rule-based DG
    let mut coords = super::dg::generate_coords(mol);

    if mol.atom_count() < 4 {
        // Need at least 4 atoms for meaningful torsions
        return coords;
    }

    // Step 2: Apply torsion preferences to aliphatic C-C-C-C rotations
    // Heuristic: prefer 180° (anti/staggered) for most aliphatic chains
    apply_torsion_preferences(mol, &mut coords);

    // Step 3: Satisfy distance constraints to preserve geometry
    let constraints = super::constraints::build_constraints(mol);
    coords = super::constraints::satisfy_constraints(&coords, mol, &constraints, 3);

    coords
}

/// Apply torsion angle preferences to coordinates using knowledge base.
///
/// Uses experimental torsion angle preferences based on atom types and
/// chemical patterns to refine the initial coordinates.
fn apply_torsion_preferences(mol: &Molecule, coords: &mut Coords3D) {
    let n = mol.atom_count();
    let mut applied = std::collections::HashSet::new();

    // Scan for 4-atom chains (A-B-C-D torsions)
    for b in 0..n {
        for c in 0..n {
            if b == c {
                continue;
            }

            // Check if B-C bond exists
            let b_idx = AtomIdx(b as u32);
            let c_idx = AtomIdx(c as u32);
            if mol.bond_between(b_idx, c_idx).is_none() {
                continue;
            }

            // Find neighbors of B (excluding C) and C (excluding B)
            let b_neighbors: Vec<usize> = (0..n)
                .filter(|&a| a != c && mol.bond_between(b_idx, AtomIdx(a as u32)).is_some())
                .collect();

            let c_neighbors: Vec<usize> = (0..n)
                .filter(|&a| a != b && mol.bond_between(c_idx, AtomIdx(a as u32)).is_some())
                .collect();

            // Try all A-B-C-D combinations
            for &a in &b_neighbors {
                for &d in &c_neighbors {
                    let key = (a.min(d), a.max(d), b.min(c), b.max(c));
                    if applied.contains(&key) {
                        continue;
                    }

                    let a_idx = AtomIdx(a as u32);
                    let d_idx = AtomIdx(d as u32);

                    // Get current dihedral
                    let current = super::mol_transforms::get_dihedral(
                        coords,
                        a_idx,
                        b_idx,
                        c_idx,
                        d_idx,
                    );

                    if current.is_none() {
                        continue;
                    }

                    let current_deg = current.unwrap() * 180.0 / std::f64::consts::PI;

                    // Get torsion preference from knowledge base
                    let preference = get_torsion_preference(mol, a_idx, b_idx, c_idx, d_idx)
                        .unwrap_or_else(|| default_torsion_preference());

                    let target_deg = preference.angle_deg;

                    // Only apply if difference is significant (> 20°)
                    if (current_deg - target_deg).abs() > 20.0 {
                        let target_rad = target_deg * std::f64::consts::PI / 180.0;
                        *coords = super::mol_transforms::set_dihedral(
                            &coords,
                            mol,
                            a_idx,
                            b_idx,
                            c_idx,
                            d_idx,
                            target_rad,
                        );
                        applied.insert(key);
                    }
                }
            }
        }
    }
}
