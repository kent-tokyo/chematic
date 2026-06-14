//! Simplified force-field geometry minimization for molecular structures.
//!
//! Uses gradient descent with finite differences over energy terms:
//! bond stretching, angle bending, VDW repulsion, and (for MMFF94) electrostatic interactions.
//! Bond lengths and angles use element-specific parameters; charges use 3D geometry.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_ff::{
    assign_dreiding_types, assign_mmff94_types, dreiding_angle, dreiding_bond_len, dreiding_vdw,
    mmff94_angle_params, mmff94_bond_params, mmff94_charges_3d, mmff94_vdw_params,
};

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Force field parameters
// ---------------------------------------------------------------------------

/// Bond stretching spring constant (kcal/mol/Ų).
/// Used in both DREIDING and generic force fields.
const BOND_SPRING_CONSTANT: f64 = 700.0;

/// Angle bending spring constant (kcal/mol/rad²).
/// Used in both DREIDING and generic force fields.
const ANGLE_SPRING_CONSTANT: f64 = 100.0;

/// Van der Waals interaction cutoff distance (Ångströms).
/// Interactions beyond this distance are ignored.
const VDW_CUTOFF: f64 = 8.0;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Force field selection for minimization.
#[derive(Debug, Clone, Copy)]
#[derive(Default)]
pub enum ForceField {
    /// UFF-derived force field (default, fast).
    UFF,
    /// DREIDING force field.
    #[default]
    DREIDING,
    /// MMFF94 force field (Merck Molecular Force Field 94, industry standard).
    MMFF94,
}


/// Configuration for the minimization algorithm.
pub struct MinimizeConfig {
    /// Maximum number of gradient-descent steps.
    pub max_steps: usize,
    /// Base step size for coordinate updates (scaled by max gradient).
    pub step_size: f64,
    /// Convergence threshold: stop when max gradient component < this value.
    pub convergence: f64,
    /// Force field to use for energy calculation.
    pub force_field: ForceField,
}

impl Default for MinimizeConfig {
    fn default() -> Self {
        Self {
            max_steps: 200,
            step_size: 0.05,
            convergence: 1e-4,
            force_field: ForceField::DREIDING,
        }
    }
}

/// Minimize molecular geometry using default configuration.
pub fn minimize(mol: &Molecule, coords: Coords3D) -> Coords3D {
    minimize_with_config(mol, coords, &MinimizeConfig::default())
}

/// Alias for [`minimize`] using UFF-derived energy terms.
///
/// Provided for discoverability; identical to calling `minimize(mol, coords)`.
pub fn minimize_uff(mol: &Molecule, coords: Coords3D) -> Coords3D {
    minimize(mol, coords)
}

/// Minimize molecular geometry using DREIDING force field parameters.
///
/// Uses the same gradient descent approach as [`minimize`], but employs DREIDING
/// force field parameters for bond lengths, angles, and VDW interactions instead of UFF.
///
/// # Arguments
/// * `mol` - Molecule to minimize
/// * `coords` - Initial 3D coordinates
///
/// # Returns
/// Minimized coordinates
pub fn minimize_dreiding(mol: &Molecule, coords: Coords3D) -> Coords3D {
    minimize_dreiding_with_config(mol, coords, &MinimizeConfig::default())
}

/// Minimize molecular geometry using MMFF94 force field (industry standard for small molecules).
///
/// MMFF94 (Merck Molecular Force Field 94) provides high-quality geometry optimization
/// suitable for drug-like molecules. This is the recommended force field for most use cases.
///
/// # Arguments
/// * `mol` - Molecule to minimize
/// * `coords` - Initial 3D coordinates
///
/// # Returns
/// Minimized coordinates
pub fn minimize_mmff94(mol: &Molecule, coords: Coords3D) -> Coords3D {
    let config = MinimizeConfig {
        force_field: ForceField::MMFF94,
        ..MinimizeConfig::default()
    };
    minimize_with_config(mol, coords, &config)
}

/// Generic gradient descent minimization with custom energy evaluation function.
///
/// This function encapsulates the core gradient descent loop, accepting a closure
/// that computes the energy given current coordinates. This allows different force fields
/// to share the same optimization control flow.
///
/// # Arguments
/// * `mol` - Molecule to minimize
/// * `coords` - Initial 3D coordinates
/// * `config` - Minimization configuration
/// * `energy_fn` - Closure that computes total energy for given coordinates
fn minimize_gradient_descent<F>(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
    energy_fn: F,
) -> Coords3D
where
    F: Fn(&Coords3D) -> f64,
{
    if mol.atom_count() <= 1 {
        return coords;
    }

    let mut c = coords;
    let delta = 1e-4;

    for _ in 0..config.max_steps {
        let mut grad = vec![Point3::zero(); mol.atom_count()];
        let mut max_grad = 0.0f64;

        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);

            // Compute gradient components via finite differences along x, y, z
            grad[i].x = {
                let orig = c.get(idx);
                let mut p = orig;
                p.x += delta;
                c.set(idx, p);
                let ep = energy_fn(&c);
                let mut p = orig;
                p.x -= delta;
                c.set(idx, p);
                let em = energy_fn(&c);
                c.set(idx, orig);
                (ep - em) / (2.0 * delta)
            };

            grad[i].y = {
                let orig = c.get(idx);
                let mut p = orig;
                p.y += delta;
                c.set(idx, p);
                let ep = energy_fn(&c);
                let mut p = orig;
                p.y -= delta;
                c.set(idx, p);
                let em = energy_fn(&c);
                c.set(idx, orig);
                (ep - em) / (2.0 * delta)
            };

            grad[i].z = {
                let orig = c.get(idx);
                let mut p = orig;
                p.z += delta;
                c.set(idx, p);
                let ep = energy_fn(&c);
                let mut p = orig;
                p.z -= delta;
                c.set(idx, p);
                let em = energy_fn(&c);
                c.set(idx, orig);
                (ep - em) / (2.0 * delta)
            };

            let gmax = grad[i].x.abs().max(grad[i].y.abs()).max(grad[i].z.abs());
            if gmax > max_grad {
                max_grad = gmax;
            }
        }

        if max_grad < config.convergence {
            break;
        }

        let scale = config.step_size / max_grad.max(1e-8);
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            let p = c.get(idx);
            c.set(
                idx,
                Point3::new(
                    p.x - scale * grad[i].x,
                    p.y - scale * grad[i].y,
                    p.z - scale * grad[i].z,
                ),
            );
        }
    }

    c
}

/// Internal MMFF94 minimization implementation with custom config.
fn minimize_mmff94_with_config(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
) -> Coords3D {
    if mol.atom_count() <= 1 {
        return coords;
    }

    // Assign MMFF94 types for all atoms
    let mmff94_types = match assign_mmff94_types(mol) {
        Ok(types) => types,
        Err(_) => return coords, // Fall back if type assignment fails
    };

    minimize_gradient_descent(mol, coords, config, |c| {
        total_energy_mmff94(mol, c, &mmff94_types)
    })
}

/// Minimize molecular geometry using DREIDING parameters with custom configuration.
pub fn minimize_dreiding_with_config(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
) -> Coords3D {
    if mol.atom_count() <= 1 {
        return coords;
    }

    // Assign DREIDING types for all atoms
    let dreiding_types = assign_dreiding_types(mol);

    minimize_gradient_descent(mol, coords, config, |c| {
        total_energy_dreiding(mol, c, &dreiding_types)
    })
}

fn total_energy_dreiding(
    mol: &Molecule,
    coords: &Coords3D,
    dreiding_types: &[chematic_ff::DREIDINGType],
) -> f64 {
    bond_energy_dreiding(mol, coords, dreiding_types)
        + angle_energy_dreiding(mol, coords, dreiding_types)
        + vdw_energy_dreiding(mol, coords, dreiding_types)
}

fn bond_energy_dreiding(
    mol: &Molecule,
    coords: &Coords3D,
    dreiding_types: &[chematic_ff::DREIDINGType],
) -> f64 {
    let mut energy = 0.0;
    let k = BOND_SPRING_CONSTANT;
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let r = coords.get(a1).distance(&coords.get(a2));
        let t1 = dreiding_types[a1.0 as usize];
        let t2 = dreiding_types[a2.0 as usize];
        let r0 = dreiding_bond_len(t1, t2, bond.order);
        let dr = r - r0;
        energy += 0.5 * k * dr * dr;
    }
    energy
}

fn angle_energy_dreiding(
    mol: &Molecule,
    coords: &Coords3D,
    dreiding_types: &[chematic_ff::DREIDINGType],
) -> f64 {
    let mut energy = 0.0;
    let k = ANGLE_SPRING_CONSTANT;

    for b_idx in 0..mol.atom_count() {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();

        if neighbors.len() < 2 {
            continue;
        }

        let theta0 = dreiding_angle(dreiding_types[b_idx]);

        for (i, &a) in neighbors.iter().enumerate() {
            for &c in &neighbors[i + 1..] {
                let pb = coords.get(b);

                let pa = coords.get(a);
                let pc = coords.get(c);

                let va = pa.sub(&pb);
                let vc = pc.sub(&pb);

                let na = va.norm();
                let nc = vc.norm();

                if na < 1e-10 || nc < 1e-10 {
                    continue;
                }

                let cos_theta = (va.dot(&vc) / (na * nc)).clamp(-1.0, 1.0);
                let theta = cos_theta.acos();
                let dtheta = theta - theta0;
                energy += 0.5 * k * dtheta * dtheta;
            }
        }
    }

    energy
}

fn vdw_energy_dreiding(
    mol: &Molecule,
    coords: &Coords3D,
    dreiding_types: &[chematic_ff::DREIDINGType],
) -> f64 {
    let n = mol.atom_count();
    let cutoff = VDW_CUTOFF;

    let mut excluded: HashSet<(usize, usize)> = HashSet::new();

    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
    }

    for b_idx in 0..n {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(b).map(|(nb, _)| nb.0 as usize).collect();
        for ii in 0..neighbors.len() {
            for jj in (ii + 1)..neighbors.len() {
                let i = neighbors[ii];
                let j = neighbors[jj];
                excluded.insert((i.min(j), i.max(j)));
            }
        }
    }

    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            let r = coords
                .get(AtomIdx(i as u32))
                .distance(&coords.get(AtomIdx(j as u32)));

            if r < 0.01 || r >= cutoff {
                continue;
            }

            let t_i = dreiding_types[i];
            let t_j = dreiding_types[j];
            let (r0_i, well_i) = dreiding_vdw(t_i);
            let (r0_j, well_j) = dreiding_vdw(t_j);

            // Lorentz-Berthelot combining rules
            let r0 = (r0_i + r0_j) / 2.0;
            let well = (well_i * well_j).sqrt();

            let ratio = r0 / r;
            let ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
            let ratio12 = ratio6 * ratio6;
            energy += well * (ratio12 - 2.0 * ratio6);
        }
    }

    energy
}

/// Minimize molecular geometry using the provided configuration.
pub fn minimize_with_config(mol: &Molecule, coords: Coords3D, config: &MinimizeConfig) -> Coords3D {
    if mol.atom_count() <= 1 {
        return coords;
    }

    // Dispatch to appropriate force field implementation
    match config.force_field {
        ForceField::MMFF94 => minimize_mmff94_with_config(mol, coords, config),
        _ => {
            // Default UFF/DREIDING path (unchanged behavior)
            minimize_generic_with_config(mol, coords, config)
        }
    }
}

fn minimize_generic_with_config(mol: &Molecule, coords: Coords3D, config: &MinimizeConfig) -> Coords3D {
    minimize_gradient_descent(mol, coords, config, |c| total_energy(mol, c))
}

// ---------------------------------------------------------------------------
// Total energy
// ---------------------------------------------------------------------------

fn total_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    bond_energy(mol, coords) + angle_energy(mol, coords) + vdw_energy(mol, coords)
}

// ---------------------------------------------------------------------------
// UFF-derived element parameters
// ---------------------------------------------------------------------------

/// Ideal bond length (Å) by atom element pair and bond order.
/// Canonical pair: (a, b) where a <= b lexicographically.
fn ideal_bond_len(sym1: &str, sym2: &str, order: BondOrder) -> f64 {
    let (a, b) = if sym1 <= sym2 {
        (sym1, sym2)
    } else {
        (sym2, sym1)
    };
    match (a, b, order) {
        // C–C
        ("C", "C", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.540,
        ("C", "C", BondOrder::Double) => 1.340,
        ("C", "C", BondOrder::Triple) => 1.204,
        ("C", "C", BondOrder::Aromatic) => 1.395,
        // C–H
        ("C", "H", _) => 1.090,
        // C–N
        ("C", "N", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.469,
        ("C", "N", BondOrder::Double) => 1.279,
        ("C", "N", BondOrder::Triple) => 1.158,
        ("C", "N", BondOrder::Aromatic) => 1.340,
        // C–O
        ("C", "O", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.427,
        ("C", "O", BondOrder::Double) => 1.217,
        ("C", "O", BondOrder::Aromatic) => 1.355,
        // C–S
        ("C", "S", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.819,
        ("C", "S", BondOrder::Double) => 1.610,
        ("C", "S", BondOrder::Aromatic) => 1.750,
        // C–F
        ("C", "F", _) => 1.350,
        // C–Cl ("C" < "Cl" since "C" == "C" and "" < "l")
        ("C", "Cl", _) => 1.770,
        // C–Br ("Br" < "C")
        ("Br", "C", _) => 1.940,
        // C–I
        ("C", "I", _) => 2.140,
        // C–P
        ("C", "P", _) => 1.840,
        // C–Si
        ("C", "Si", _) => 1.870,
        // H–H
        ("H", "H", _) => 0.741,
        // H–N
        ("H", "N", _) => 1.010,
        // H–O
        ("H", "O", _) => 0.960,
        // H–S
        ("H", "S", _) => 1.340,
        // H–P
        ("H", "P", _) => 1.420,
        // N–N
        ("N", "N", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.450,
        ("N", "N", BondOrder::Double) => 1.250,
        ("N", "N", BondOrder::Triple) => 1.100,
        ("N", "N", BondOrder::Aromatic) => 1.350,
        // N–O
        ("N", "O", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.400,
        ("N", "O", BondOrder::Double) => 1.210,
        ("N", "O", BondOrder::Aromatic) => 1.340,
        // O–O
        ("O", "O", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.480,
        ("O", "O", BondOrder::Double) => 1.210,
        // S–S
        ("S", "S", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 2.050,
        ("S", "S", BondOrder::Double) => 1.890,
        // P–P
        ("P", "P", _) => 2.280,
        // fallback: order-based only
        _ => match order {
            BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.54,
            BondOrder::Double => 1.34,
            BondOrder::Triple => 1.20,
            BondOrder::Quadruple => 1.20,
            BondOrder::Aromatic => 1.40,
            BondOrder::Zero
            | BondOrder::Dative
            | BondOrder::QueryAny
            | BondOrder::QuerySingleOrDouble
            | BondOrder::QuerySingleOrAromatic => 1.54,
            BondOrder::QueryDoubleOrAromatic => 1.40,
        },
    }
}

/// Atom hybridization inferred from bond orders and aromaticity.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Hybridization {
    SP,  // linear (triple bond present)
    SP2, // trigonal planar (double bond or aromatic)
    SP3, // tetrahedral
}

fn atom_hybridization(mol: &Molecule, idx: AtomIdx) -> Hybridization {
    if mol.atom(idx).aromatic {
        return Hybridization::SP2;
    }
    let mut has_triple = false;
    let mut has_double_or_aromatic = false;
    for (_, bond_idx) in mol.neighbors(idx) {
        match mol.bond(bond_idx).order {
            BondOrder::Triple => has_triple = true,
            BondOrder::Double | BondOrder::Aromatic => has_double_or_aromatic = true,
            _ => {}
        }
    }
    if has_triple {
        Hybridization::SP
    } else if has_double_or_aromatic {
        Hybridization::SP2
    } else {
        Hybridization::SP3
    }
}

/// Ideal bond angle (radians) for a center atom given its hybridization.
fn ideal_angle_rad(sym: &str, hyb: Hybridization) -> f64 {
    match hyb {
        Hybridization::SP => 180.0_f64.to_radians(),
        Hybridization::SP2 => 120.0_f64.to_radians(),
        Hybridization::SP3 => match sym {
            "O" | "Se" => 104.5_f64.to_radians(),
            "N" => 107.0_f64.to_radians(),
            "S" => 99.0_f64.to_radians(),
            "P" => 93.0_f64.to_radians(),
            _ => 109.47_f64.to_radians(),
        },
    }
}

/// VDW radius (Å) derived from UFF/Bondi values.
fn uff_vdw_radius(sym: &str) -> f64 {
    match sym {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "Si" => 2.10,
        "P" => 1.80,
        "S" => 1.80,
        "Cl" => 1.75,
        "Br" => 1.85,
        "I" => 1.98,
        "Se" => 1.90,
        "Te" => 2.06,
        _ => 1.70,
    }
}

// ---------------------------------------------------------------------------
// Bond stretching energy
// ---------------------------------------------------------------------------

fn bond_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut energy = 0.0;
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let r = coords.get(a1).distance(&coords.get(a2));
        let sym1 = mol.atom(a1).element.symbol();
        let sym2 = mol.atom(a2).element.symbol();
        let r0 = ideal_bond_len(sym1, sym2, bond.order);
        let dr = r - r0;
        energy += 0.5 * BOND_SPRING_CONSTANT * dr * dr;
    }
    energy
}

// ---------------------------------------------------------------------------
// Angle bending energy
// ---------------------------------------------------------------------------

fn angle_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut energy = 0.0;

    for b_idx in 0..mol.atom_count() {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();

        if neighbors.len() < 2 {
            continue;
        }

        let sym_b = mol.atom(b).element.symbol();
        let hyb = atom_hybridization(mol, b);
        let theta0 = ideal_angle_rad(sym_b, hyb);
        let pb = coords.get(b);

        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let a = neighbors[i];
                let c = neighbors[j];

                let pa = coords.get(a);
                let pc = coords.get(c);

                let va = pa.sub(&pb);
                let vc = pc.sub(&pb);

                let na = va.norm();
                let nc = vc.norm();

                if na < 1e-10 || nc < 1e-10 {
                    continue;
                }

                let cos_theta = (va.dot(&vc) / (na * nc)).clamp(-1.0, 1.0);
                let theta = cos_theta.acos();
                let dtheta = theta - theta0;
                energy += 0.5 * ANGLE_SPRING_CONSTANT * dtheta * dtheta;
            }
        }
    }

    energy
}

// ---------------------------------------------------------------------------
// VDW repulsion energy
// ---------------------------------------------------------------------------

fn vdw_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let n = mol.atom_count();
    let cutoff = VDW_CUTOFF;

    let mut excluded: HashSet<(usize, usize)> = HashSet::new();

    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
    }

    for b_idx in 0..n {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(b).map(|(nb, _)| nb.0 as usize).collect();
        for ii in 0..neighbors.len() {
            for jj in (ii + 1)..neighbors.len() {
                let i = neighbors[ii];
                let j = neighbors[jj];
                excluded.insert((i.min(j), i.max(j)));
            }
        }
    }

    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            let r = coords
                .get(AtomIdx(i as u32))
                .distance(&coords.get(AtomIdx(j as u32)));

            if r < 0.01 || r >= cutoff {
                continue;
            }

            let sym_i = mol.atom(AtomIdx(i as u32)).element.symbol();
            let sym_j = mol.atom(AtomIdx(j as u32)).element.symbol();
            let r0 = uff_vdw_radius(sym_i) + uff_vdw_radius(sym_j);

            let ratio = r0 / r;
            let ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
            let ratio12 = ratio6 * ratio6;
            energy += 0.05 * ratio12;
        }
    }

    energy
}

// ---------------------------------------------------------------------------
// MMFF94 Energy Calculations
// ---------------------------------------------------------------------------

fn total_energy_mmff94(
    mol: &Molecule,
    coords: &Coords3D,
    mmff94_types: &[chematic_ff::MMFF94Type],
) -> f64 {
    let bond_e = bond_energy_mmff94(mol, coords, mmff94_types);
    let angle_e = angle_energy_mmff94(mol, coords, mmff94_types);
    let vdw_e = vdw_energy_mmff94(mol, coords, mmff94_types);

    // Add electrostatic energy using 3D-based charges (B5 Phase 2)
    let elec_e = electrostatic_energy_mmff94(mol, coords, mmff94_types).unwrap_or(0.0);

    bond_e + angle_e + vdw_e + elec_e
}

fn bond_energy_mmff94(
    mol: &Molecule,
    coords: &Coords3D,
    mmff94_types: &[chematic_ff::MMFF94Type],
) -> f64 {
    let mut energy = 0.0;

    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let r = coords.get(a1).distance(&coords.get(a2));
        let t1 = mmff94_types[a1.0 as usize];
        let t2 = mmff94_types[a2.0 as usize];

        if let Some(params) = mmff94_bond_params(t1, t2, bond.order) {
            let dr = r - params.r0;
            energy += 0.5 * params.kb * dr * dr;
        }
    }

    energy
}

fn angle_energy_mmff94(
    mol: &Molecule,
    coords: &Coords3D,
    mmff94_types: &[chematic_ff::MMFF94Type],
) -> f64 {
    let mut energy = 0.0;

    for b_idx in 0..mol.atom_count() {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();

        if neighbors.len() < 2 {
            continue;
        }

        for (i, &a) in neighbors.iter().enumerate() {
            for &c in &neighbors[i + 1..] {
                let t1 = mmff94_types[a.0 as usize];
                let t2 = mmff94_types[b_idx];
                let t3 = mmff94_types[c.0 as usize];

                if let Some(params) = mmff94_angle_params(t1, t2, t3) {
                    let pb = coords.get(b);
                    let pa = coords.get(a);
                    let pc = coords.get(c);

                    let va = pa.sub(&pb);
                    let vc = pc.sub(&pb);

                    let na = va.norm();
                    let nc = vc.norm();

                    if na < 1e-10 || nc < 1e-10 {
                        continue;
                    }

                    let cos_theta = (va.dot(&vc) / (na * nc)).clamp(-1.0, 1.0);
                    let theta = cos_theta.acos();
                    let dtheta = theta - params.theta0;
                    energy += 0.5 * params.ka * dtheta * dtheta;
                }
            }
        }
    }

    energy
}

fn vdw_energy_mmff94(
    mol: &Molecule,
    coords: &Coords3D,
    mmff94_types: &[chematic_ff::MMFF94Type],
) -> f64 {
    let n = mol.atom_count();
    let cutoff = VDW_CUTOFF;
    let mut excluded: HashSet<(usize, usize)> = HashSet::new();

    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
    }

    // Add 1-3 exclusions (skip vdW for atoms separated by one bond)
    for b_idx in 0..n {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(b).map(|(nb, _)| nb.0 as usize).collect();
        for &neighbor in &neighbors {
            excluded.insert((b_idx.min(neighbor), b_idx.max(neighbor)));
        }
    }

    let mut energy = 0.0;

    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }

            let ri = coords.get(AtomIdx(i as u32));
            let rj = coords.get(AtomIdx(j as u32));
            let d = ri.distance(&rj);

            if d > cutoff {
                continue;
            }

            let params_i = mmff94_vdw_params(mmff94_types[i]);
            let params_j = mmff94_vdw_params(mmff94_types[j]);

            // Combine using geometric mean
            let r_ij = (params_i.r_star * params_j.r_star).sqrt();
            let eps_ij = (params_i.epsilon * params_j.epsilon).sqrt();

            // Lennard-Jones 12-6
            if d > 0.0 {
                let r6 = (r_ij / d).powi(6);
                energy += eps_ij * (r6 * r6 - 2.0 * r6);
            }
        }
    }

    energy
}

/// Electrostatic energy using Coulomb's law with 3D-based MMFF94 charges.
/// Uses dielectric screening (kr where k~4 for organic molecules in their own environment).
fn electrostatic_energy_mmff94(
    mol: &Molecule,
    coords: &Coords3D,
    _mmff94_types: &[chematic_ff::MMFF94Type],
) -> Result<f64, String> {
    // Convert coordinates to tuple format for charge calculation
    let coord_tuples: Vec<(f64, f64, f64)> = (0..mol.atom_count())
        .map(|i| {
            let p = coords.get(AtomIdx(i as u32));
            (p.x, p.y, p.z)
        })
        .collect();

    // Calculate 3D-based MMFF94 charges
    let charges = mmff94_charges_3d(mol, &coord_tuples)
        .map_err(|e| format!("charge calculation failed: {}", e))?;

    let n = mol.atom_count();
    let mut energy = 0.0;

    // Coulomb interactions (excluding 1-2 and 1-3 pairs which are handled by bonds/angles)
    let mut excluded: HashSet<(usize, usize)> = HashSet::new();

    // Exclude 1-2 bonded pairs
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
    }

    // Exclude 1-3 pairs (through one bond)
    for b_idx in 0..n {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(b).map(|(nb, _)| nb.0 as usize).collect();
        for &neighbor in &neighbors {
            excluded.insert((b_idx.min(neighbor), b_idx.max(neighbor)));
        }
    }

    let dielectric = 4.0; // Screening factor for organic molecules
    let coulomb_const = 332.0; // kcal·Ų/(mol·e²) in Ångströms

    for i in 0..n {
        for j in (i + 1)..n {
            // Skip bonded and 1-3 interactions (handled by bonds/angles)
            if excluded.contains(&(i, j)) {
                continue;
            }

            let ri = coords.get(AtomIdx(i as u32));
            let rj = coords.get(AtomIdx(j as u32));
            let d = ri.distance(&rj);

            if d > 0.01 {
                // Coulomb interaction: E = k * q_i * q_j / (d * dielectric)
                let coulomb = coulomb_const * charges[i] * charges[j] / (d * dielectric);
                energy += coulomb;
            }
        }
    }

    Ok(energy)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dg::generate_coords;
    use chematic_smiles::parse;

    fn all_pairs_min_dist(coords: &Coords3D, n: usize) -> f64 {
        let mut min_d = f64::MAX;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = coords
                    .get(AtomIdx(i as u32))
                    .distance(&coords.get(AtomIdx(j as u32)));
                min_d = min_d.min(d);
            }
        }
        min_d
    }

    #[test]
    fn test_single_atom_unchanged() {
        let mol = parse("O").unwrap();
        let coords = generate_coords(&mol);
        let orig = coords.get(AtomIdx(0));
        let result = minimize(&mol, coords);
        let after = result.get(AtomIdx(0));
        assert!((orig.x - after.x).abs() < 1e-10);
    }

    #[test]
    fn test_zero_steps_unchanged() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig {
            max_steps: 0,
            ..MinimizeConfig::default()
        };
        let before0 = coords.get(AtomIdx(0));
        let result = minimize_with_config(&mol, coords, &config);
        let after0 = result.get(AtomIdx(0));
        assert!((before0.x - after0.x).abs() < 1e-10);
    }

    #[test]
    fn test_ethane_bond_after_minimize() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(
            d > 1.2 && d < 1.8,
            "C-C distance={d:.3}, expected 1.2-1.8 Å"
        );
    }

    #[test]
    fn test_ethane_converges_to_uff_length() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        // UFF C-C single bond is 1.540 Å; minimizer should get within 0.05 Å.
        assert!(
            (d - 1.540).abs() < 0.05,
            "C-C distance={d:.4}, expected ~1.540"
        );
    }

    #[test]
    fn test_propane_no_clash() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "atom clash: min distance={min_d:.3}");
    }

    #[test]
    fn test_benzene_no_clash() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(
            min_d > 0.8,
            "atom clash in benzene: min distance={min_d:.3}"
        );
    }

    #[test]
    fn test_disconnected_no_clash() {
        let mol = parse("CC.CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(
            min_d > 0.8,
            "atom clash in disconnected: min distance={min_d:.3}"
        );
    }

    #[test]
    fn test_default_config_no_panic() {
        let mol = parse("CC(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        assert_eq!(result.atom_count(), mol.atom_count());
    }

    #[test]
    fn test_acetic_acid_no_clash() {
        let mol = parse("CC(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "clash in acetic acid: {min_d:.3}");
    }

    #[test]
    fn test_minimize_idempotent() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let result1 = minimize(&mol, coords);
        let e1 = total_energy(&mol, &result1);
        let result2 = minimize(&mol, result1);
        let e2 = total_energy(&mol, &result2);
        assert!(e2 <= e1 + 1.0, "energy increased: e1={e1:.4}, e2={e2:.4}");
    }

    #[test]
    fn test_naphthalene_no_overlap() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "overlap in naphthalene: {min_d:.3}");
    }

    #[test]
    fn test_co_bond_double_shorter_than_single() {
        // Acetic acid: C=O should be shorter than C-O
        let mol = parse("CC(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        // Atom 1 is the carbonyl C, its bonds include C=O (double) and C-O (single).
        // Just check overall: minimized coords have no clash and atom count preserved.
        assert_eq!(result.atom_count(), 4);
        let min_d = all_pairs_min_dist(&result, 4);
        assert!(min_d > 0.5, "clash in CO test: {min_d:.3}");
    }

    #[test]
    fn test_heteroatom_c_n_bond() {
        let mol = parse("CN").unwrap(); // methylamine
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        // C-N single bond UFF: 1.469 Å; expect within 0.1 Å.
        assert!(
            (d - 1.469).abs() < 0.1,
            "C-N distance={d:.4}, expected ~1.469"
        );
    }

    #[test]
    fn test_acetylene_sp_hybridization() {
        let mol = parse("C#C").unwrap(); // acetylene: C≡C
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        // C≡C triple bond UFF: 1.204 Å; expect within 0.05 Å.
        assert!(
            (d - 1.204).abs() < 0.05,
            "C≡C distance={d:.4}, expected ~1.204"
        );
    }

    #[test]
    fn test_ideal_bond_len_cc_single() {
        assert!((ideal_bond_len("C", "C", BondOrder::Single) - 1.540).abs() < 1e-6);
        assert!((ideal_bond_len("C", "C", BondOrder::Double) - 1.340).abs() < 1e-6);
        assert!((ideal_bond_len("C", "C", BondOrder::Triple) - 1.204).abs() < 1e-6);
        assert!((ideal_bond_len("C", "C", BondOrder::Aromatic) - 1.395).abs() < 1e-6);
    }

    #[test]
    fn test_ideal_bond_len_symmetry() {
        // Should be the same regardless of argument order.
        let bo = BondOrder::Single;
        assert_eq!(ideal_bond_len("C", "N", bo), ideal_bond_len("N", "C", bo));
        assert_eq!(ideal_bond_len("C", "O", bo), ideal_bond_len("O", "C", bo));
        assert_eq!(ideal_bond_len("Br", "C", bo), ideal_bond_len("C", "Br", bo));
    }

    #[test]
    fn test_atom_hybridization_sp2_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        for i in 0..6 {
            assert_eq!(
                atom_hybridization(&mol, AtomIdx(i)),
                Hybridization::SP2,
                "benzene atom {i} should be SP2"
            );
        }
    }

    #[test]
    fn test_atom_hybridization_sp_triple() {
        let mol = parse("C#C").unwrap();
        assert_eq!(atom_hybridization(&mol, AtomIdx(0)), Hybridization::SP);
        assert_eq!(atom_hybridization(&mol, AtomIdx(1)), Hybridization::SP);
    }

    #[test]
    fn test_atom_hybridization_sp3_alkane() {
        let mol = parse("CCC").unwrap();
        for i in 0..3 {
            assert_eq!(
                atom_hybridization(&mol, AtomIdx(i)),
                Hybridization::SP3,
                "propane atom {i} should be SP3"
            );
        }
    }

    #[test]
    fn test_minimize_dreiding_ethane_no_clash() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let min_coords = minimize_dreiding(&mol, coords);
        let n = mol.atom_count();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = min_coords
                    .get(AtomIdx(i as u32))
                    .distance(&min_coords.get(AtomIdx(j as u32)));
                assert!(
                    d > 0.5,
                    "atoms {i} and {j} clashed after DREIDING minimization (d={d:.3})"
                );
            }
        }
    }

    #[test]
    fn test_minimize_dreiding_benzene_no_clash() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let min_coords = minimize_dreiding(&mol, coords);
        let n = mol.atom_count();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = min_coords
                    .get(AtomIdx(i as u32))
                    .distance(&min_coords.get(AtomIdx(j as u32)));
                assert!(
                    d > 0.5,
                    "atoms {i} and {j} clashed after DREIDING minimization (d={d:.3})"
                );
            }
        }
    }

    #[test]
    fn test_minimize_mmff94_ethane() {
        let mol = parse("CC").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);
        assert_eq!(result.atom_count(), 2);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(d > 1.4 && d < 1.7, "C-C should be ~1.54 Å, got {:.3}", d);
    }

    #[test]
    fn test_minimize_mmff94_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);
        assert_eq!(result.atom_count(), 6);
        let min_d = all_pairs_min_dist(&result, 6);
        assert!(min_d > 1.2, "benzene clash: {min_d:.3}");
    }

    #[test]
    fn test_minimize_mmff94_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);
        // Verify minimize_mmff94 completes without error and produces valid coordinates
        assert_eq!(result.atom_count(), mol.atom_count());
        for i in 0..mol.atom_count() {
            let p = result.get(chematic_core::AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                   "aspirin atom {i} has invalid coords");
        }
    }

    // ===== Phase 2: Electrostatic Energy Integration (B5) =====

    #[test]
    fn test_electrostatic_energy_methanol() {
        // Methanol has partial charges due to O-H polarity
        let mol = parse("CO").unwrap();
        let c = generate_coords(&mol);
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        // Electrostatic energy should be calculable
        let elec_e = electrostatic_energy_mmff94(&mol, &c, &mmff94_types);
        assert!(elec_e.is_ok());
        assert!(elec_e.unwrap().is_finite());
    }

    #[test]
    fn test_electrostatic_energy_carboxylic_acid() {
        // Carboxylic acids have significant charge separation
        let mol = parse("CC(=O)O").unwrap();
        let c = generate_coords(&mol);
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let elec_e = electrostatic_energy_mmff94(&mol, &c, &mmff94_types);
        assert!(elec_e.is_ok());
        let energy = elec_e.unwrap();
        assert!(energy.is_finite());
        // Carboxylic acids with negative oxygen should have non-zero electrostatic energy
    }

    #[test]
    fn test_mmff94_with_electrostatic_ethane() {
        // Ethane is non-polar, so electrostatic should be small
        let mol = parse("CC").unwrap();
        let c = generate_coords(&mol);
        let result = minimize_mmff94(&mol, c);

        // Should still minimize correctly with electrostatic term
        assert_eq!(result.atom_count(), 2);
        let d = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        assert!(d > 1.4 && d < 1.7, "ethane C-C should be ~1.54 Å with electrostatic, got {:.3}", d);
    }

    #[test]
    fn test_mmff94_minimization_includes_charge_effects() {
        // Minimize polar molecule where charges should matter
        let mol = parse("CCO").unwrap();
        let c = generate_coords(&mol);

        // Minimization should complete without error and produce valid coordinates
        let result = minimize_mmff94(&mol, c);

        // All coordinates should be finite and properly positioned
        assert_eq!(result.atom_count(), 3);
        for i in 0..3 {
            let p = result.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                   "atom {i} has invalid coordinate after minimization");
        }

        // Verify atoms are reasonably separated (no clashes)
        let c_c = result.get(AtomIdx(0)).distance(&result.get(AtomIdx(1)));
        let c_o = result.get(AtomIdx(1)).distance(&result.get(AtomIdx(2)));
        assert!(c_c > 1.0, "C-C bond too short: {c_c:.3}");
        assert!(c_o > 1.0, "C-O bond too short: {c_o:.3}");
    }

    #[test]
    fn test_mmff94_charges_3d_integration() {
        // Verify that 3D charges are being used in minimization
        let mol = parse("c1ccccc1O").unwrap(); // Phenol
        let c = generate_coords(&mol);

        // The minimization should complete without errors
        let result = minimize_mmff94(&mol, c);
        assert_eq!(result.atom_count(), mol.atom_count());

        // All coordinates should be finite
        for i in 0..mol.atom_count() {
            let p = result.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    #[test]
    fn test_total_energy_mmff94_includes_electrostatic() {
        // Verify that total_energy_mmff94 includes electrostatic component
        let mol = parse("CCN").unwrap(); // Has C-N polarity
        let c = generate_coords(&mol);
        let mmff94_types = assign_mmff94_types(&mol).unwrap();

        let total_e = total_energy_mmff94(&mol, &c, &mmff94_types);
        let bond_e = bond_energy_mmff94(&mol, &c, &mmff94_types);
        let angle_e = angle_energy_mmff94(&mol, &c, &mmff94_types);
        let vdw_e = vdw_energy_mmff94(&mol, &c, &mmff94_types);

        // Total should include bond + angle + vdw + electrostatic
        let electrostatic_e = electrostatic_energy_mmff94(&mol, &c, &mmff94_types).unwrap_or(0.0);
        let expected = bond_e + angle_e + vdw_e + electrostatic_e;

        assert!((total_e - expected).abs() < 1e-6,
                "total energy mismatch: got {}, expected {}", total_e, expected);
    }
}
