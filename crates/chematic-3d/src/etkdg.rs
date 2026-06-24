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
use crate::etkdg_knowledge::{
    AtomType, build_smarts_torsion_map, classify_atom_type, default_torsion_preference,
    get_torsion_preference,
};
use crate::prng::Prng;
use chematic_core::{AtomIdx, Molecule};

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

    let mut prng = Prng::new();
    apply_torsion_preferences_with_noise(mol, &mut coords, noise_sigma_deg, &mut prng);

    let constraints = super::constraints::build_constraints(mol);
    coords = super::constraints::satisfy_constraints(&coords, mol, &constraints, 3);

    // Hard-snap amide omega angles to 0° or 180° after constraint satisfaction.
    // Addresses RDKit issue #9266: ETKDGv3 generates twisted tertiary amides.
    snap_amide_torsions(mol, &mut coords);

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
    prng: &mut Prng,
) {
    let n = mol.atom_count();
    let mut applied = std::collections::HashSet::new();

    // Build a set of ring bonds so SMARTS rules can skip them.
    let ring_set = chematic_perception::find_sssr(mol);
    let ring_bond_set: std::collections::HashSet<(u32, u32)> = ring_set
        .rings()
        .iter()
        .flat_map(|ring| {
            let len = ring.len();
            (0..len).flat_map(move |i| {
                let a = ring[i].0;
                let b = ring[(i + 1) % len].0;
                [(a, b), (b, a)]
            })
        })
        .collect();

    // SMARTS-derived bond preferences (more precise than atom-type rules).
    // Computed once per molecule and consulted before the atom-type fallback.
    let smarts_map = build_smarts_torsion_map(mol, &ring_bond_set);

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

                    // SMARTS map is consulted first (B-C bond context);
                    // atom-type rules are the fallback.
                    let preference = smarts_map
                        .get(&(b as u32, c as u32))
                        .cloned()
                        .or_else(|| get_torsion_preference(mol, a_idx, b_idx, c_idx, d_idx))
                        .unwrap_or_else(default_torsion_preference);

                    // Add Gaussian noise scaled by bond flexibility.
                    // Rigid bonds (amide, biaryl) receive less noise than freely
                    // rotatable single bonds, preserving chemically important
                    // constraints while still exploring diverse conformers.
                    let scale = bond_flexibility_scale(mol, b_idx, c_idx);
                    let noise = if noise_sigma_deg > 0.0 && scale > 0.0 {
                        prng.gaussian_f64() * noise_sigma_deg * scale
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

/// Snap all amide omega angles (A–N–C(=O)–D) to the nearest 0° or 180°.
///
/// Addresses RDKit issue #9266: ETKDGv3 can leave tertiary amides in a twisted
/// gauche (~60° / ~120°) geometry because soft torsion penalties are insufficient
/// when constraint satisfaction later distorts the amide bond. This post-processing
/// step enforces planarity (|ω| ≤ 30° or |ω − 180°| ≤ 30°) for all amide bonds.
fn snap_amide_torsions(mol: &Molecule, coords: &mut super::coords::Coords3D) {
    use super::mol_transforms::{get_dihedral, set_dihedral};
    let n = mol.atom_count();

    for b in 0..n {
        let b_idx = AtomIdx(b as u32);
        let b_atom = mol.atom(b_idx);
        // Match any non-aromatic nitrogen — catches both secondary (NSp2) and
        // tertiary (NSp3) amide N. classify_atom_type returns NSp3 for degree-3 N,
        // which would silently skip all tertiary amides if checked against NSp2 only.
        if b_atom.element.atomic_number() != 7 || b_atom.aromatic {
            continue;
        }
        // Require a C(=O) neighbor — the defining feature of amide N.
        let c_idx = mol
            .neighbors(b_idx)
            .find(|(nb, _)| classify_atom_type(mol, *nb) == AtomType::CCarbonyl)
            .map(|(nb, _)| nb);
        let Some(c_idx) = c_idx else { continue };

        let b_neighbors: Vec<AtomIdx> = mol
            .neighbors(b_idx)
            .filter(|(nb, _)| *nb != c_idx)
            .map(|(nb, _)| nb)
            .collect();
        let c_neighbors: Vec<AtomIdx> = mol
            .neighbors(c_idx)
            .filter(|(nb, _)| *nb != b_idx)
            .map(|(nb, _)| nb)
            .collect();

        // Snap using ONE (a, d) pair per N-C(=O) bond. All torsions around the
        // same bond are rotated identically by set_dihedral, so snapping additional
        // pairs after the first would read stale coords and may overcorrect.
        'snap: for &a_idx in &b_neighbors {
            for &d_idx in &c_neighbors {
                let Some(omega_rad) = get_dihedral(coords, a_idx, b_idx, c_idx, d_idx) else {
                    continue;
                };
                let omega_deg = omega_rad.to_degrees();
                // Distance to nearest planar angle (0° or ±180°).
                let to_0 = omega_deg.abs();
                let to_180 = (omega_deg.abs() - 180.0).abs();
                let min_dist = to_0.min(to_180);
                if min_dist > 30.0 {
                    let target_deg = if to_0 < to_180 { 0.0_f64 } else { 180.0_f64 };
                    *coords = set_dihedral(
                        coords,
                        mol,
                        a_idx,
                        b_idx,
                        c_idx,
                        d_idx,
                        target_deg.to_radians(),
                    );
                }
                // One snap per bond is sufficient — subsequent (a, d) pairs would
                // read already-modified coords and risk double-correction.
                break 'snap;
            }
        }
    }
}

/// Return a noise scale factor [0, 1] for the B–C central bond of a torsion.
///
/// Freely rotatable single bonds receive 1.0 (full sigma); bonds with partial
/// double-bond character or aromatic conjugation receive a reduced scale to
/// preserve chemically important geometry.
fn bond_flexibility_scale(mol: &Molecule, b_idx: AtomIdx, c_idx: AtomIdx) -> f64 {
    let b = classify_atom_type(mol, b_idx);
    let c = classify_atom_type(mol, c_idx);

    // Check if the bond itself is a double or triple bond — no torsion noise.
    if let Some((_, bond)) = mol.bond_between(b_idx, c_idx) {
        match bond.order {
            chematic_core::BondOrder::Double
            | chematic_core::BondOrder::Triple
            | chematic_core::BondOrder::Aromatic => return 0.0,
            _ => {}
        }
    }

    match (b, c) {
        // Amide/urea N–C(=O): partial double bond, highly restricted.
        (AtomType::NSp2, AtomType::CCarbonyl) | (AtomType::CCarbonyl, AtomType::NSp2) => 0.20,
        // Biaryl / hetaryl–aryl: moderate restriction (biphenyl-like twist).
        (AtomType::CAromatic, AtomType::CAromatic)
        | (AtomType::NAromatic, AtomType::CAromatic)
        | (AtomType::CAromatic, AtomType::NAromatic)
        | (AtomType::OAromatic, AtomType::CAromatic)
        | (AtomType::CAromatic, AtomType::OAromatic)
        | (AtomType::SAromatic, AtomType::CAromatic)
        | (AtomType::CAromatic, AtomType::SAromatic) => 0.50,
        // Vinyl/alkene adjacent: semi-rigid conjugation.
        (AtomType::CSp2Alkene, AtomType::CCarbonyl)
        | (AtomType::CCarbonyl, AtomType::CSp2Alkene)
        | (AtomType::CSp2Alkene, AtomType::NSp2)
        | (AtomType::NSp2, AtomType::CSp2Alkene)
        | (AtomType::CSp2Alkene, AtomType::CAromatic)
        | (AtomType::CAromatic, AtomType::CSp2Alkene) => 0.30,
        // Aryl–sp3: slightly restricted.
        (AtomType::CAromatic, AtomType::CSp3) | (AtomType::CSp3, AtomType::CAromatic) => 0.70,
        // All other single bonds: fully flexible.
        _ => 1.0,
    }
}
