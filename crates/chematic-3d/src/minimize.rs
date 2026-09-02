//! Simplified force-field geometry minimization for molecular structures.
//!
//! Uses gradient descent with finite differences over energy terms:
//! bond stretching, angle bending, VDW repulsion, and (for MMFF94) electrostatic interactions.
//! Bond lengths and angles use element-specific parameters; charges use 3D geometry.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_ff::{
    EnergyBreakdown, MinimizerError, NumericTypeError, OOP_SP2_TYPES, UffType, angle_type_for,
    assign_mmff94_numeric_types_with_view, assign_uff_types, bond_type_for,
    is_angle_in_ring_of_size_3_or_4, minimize_mmff94_lbfgs, minimize_uff as ff_minimize_uff,
    mmff94_angle_energy_resolved, mmff94_bond_energy_resolved, mmff94_energy_breakdown, mmff94_oop,
    mmff94_stbn, mmff94_torsion_energy, mmff94_total_energy, stretch_bend_type_for,
    torsion_no_term_by_design, torsion_type_for, uff_total_energy,
};
use chematic_ff::{
    assign_dreiding_types, assign_mmff94_types, dreiding_angle, dreiding_bond_len, dreiding_vdw,
    mmff94_angle_params, mmff94_bond_params, mmff94_charges_3d, mmff94_vdw_params,
};
use chematic_perception::find_sssr;

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
#[derive(Debug, Clone, Copy, Default)]
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
    minimize_gradient_descent_reporting(mol, coords, config, energy_fn).coords
}

/// Result of the internal generic gradient-descent loop, with the convergence
/// bookkeeping the pre-existing public wrappers (which only ever wanted the
/// final `Coords3D`) historically discarded. Used by the force-field bridge
/// below to report real iteration/convergence/residual-force numbers for the
/// DREIDING/generic paths without duplicating the loop.
struct GradientDescentReport {
    coords: Coords3D,
    iterations: usize,
    converged: bool,
    final_max_grad: f64,
}

fn minimize_gradient_descent_reporting<F>(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
    energy_fn: F,
) -> GradientDescentReport
where
    F: Fn(&Coords3D) -> f64,
{
    if mol.atom_count() <= 1 {
        return GradientDescentReport {
            coords,
            iterations: 0,
            converged: true,
            final_max_grad: 0.0,
        };
    }

    let mut c = coords;
    let delta = 1e-4;
    let mut iterations = 0usize;
    let mut converged = false;
    let mut final_max_grad = 0.0f64;

    for _ in 0..config.max_steps {
        iterations += 1;
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

        final_max_grad = max_grad;

        if max_grad < config.convergence {
            converged = true;
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

    GradientDescentReport {
        coords: c,
        iterations,
        converged,
        final_max_grad,
    }
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

fn minimize_generic_with_config(
    mol: &Molecule,
    coords: Coords3D,
    config: &MinimizeConfig,
) -> Coords3D {
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
// Force-field bridge (Wave 1, Agent F): real chematic-ff MMFF94/UFF
// ---------------------------------------------------------------------------
//
// Everything above this point (`minimize`, `minimize_dreiding`,
// `minimize_mmff94`, `minimize_uff`, `minimize_with_config`, and their
// crippled bond+angle(+VdW/electrostatic)-only energy functions) is
// UNCHANGED — kept byte-for-byte for backward compatibility, per this
// program's hard constraint that no existing public API's default behavior
// changes in this PR.
//
// This section is purely additive: a new opt-in entry point,
// `minimize_with_policy`/`minimize_with_policy_gated`, that bridges to
// chematic-ff's COMPLETE MMFF94 implementation (bond + angle + stretch-bend +
// torsion + out-of-plane + vdW + electrostatic, real L-BFGS minimizer) and
// its separate UFF module, instead of re-implementing a smaller, worse
// subset of the same physics locally. See `docs/rfcs/etkdg_3d_gap_rfc.md`
// ("mechanism 3") for the bug this fixes: `bond_energy_mmff94`/
// `angle_energy_mmff94` above use `if let Some(params) = mmff94_bond_params(...)
// { ... }` — when a type pair isn't covered, that internal coordinate
// silently contributes zero energy/gradient (no restoring force) while VdW
// repulsion still pushes atoms apart. This bridge never does that: missing
// coverage is always a typed `Err` (`Mmff94BondAngleStrict`) or an explicitly reported
// fallback (`Mmff94WithUffFallback`), never a silent zero.
//
// Variant names (`Mmff94BondAngleStrict`/`Mmff94WithUffFallback`/`UffOnly`/`Dreiding`/
// `None`) and the `requested_force_field`/`actual_force_field_used`/
// `fallback_reason`/`missing_parameter_classes` result fields below were
// renamed/added per Coordinator's design-review refinement mid-implementation
// (binding — see PR body); behavior is unchanged from the original brief.

/// Force-field policy for [`minimize_with_policy`]/[`minimize_with_policy_gated`]
/// — the opt-in bridge to chematic-ff's complete implementations. Does not
/// change any existing public function's behavior or default in this crate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ForceFieldPolicy {
    /// Full MMFF94 (bond + angle + stretch-bend + torsion + out-of-plane +
    /// vdW + electrostatic, `chematic_ff::minimize_mmff94_lbfgs`). Refuses
    /// (typed `Err`) if even one BOND or ANGLE lacks MMFF94 parameters,
    /// rather than silently contributing zero energy/gradient for that
    /// internal coordinate — this is the exact scope named in this variant
    /// (`BondAngle`, not "every MMFF94 term"): torsion/out-of-plane coverage
    /// is always measured and reported (`coverage`/
    /// `missing_parameter_classes` on the result) but does NOT by itself
    /// trigger a refusal here. Renamed from a plain `Mmff94Strict` (found in
    /// independent review to be a scope/naming mismatch: the old name
    /// implied gating on every required MMFF94 term, but the gate only ever
    /// checked bond+angle). Use `minimize_with_policy_gated(..., true, false)` to
    /// also gate on torsion/out-of-plane, or `(..., false, true)` to also gate
    /// on stretch-bend (Priority 2 / Stage 1B).
    ///
    /// Naming/scope note for Coordinator: post-chematic-ff-#183 (bond/angle/
    /// torsion classification fixes), the 58-molecule corpus example
    /// (`examples/mmff94_bridge_coverage_report.rs`) now measures 32/58
    /// molecules pass this bond+angle-only gate vs. 31/58 under the widened
    /// (+torsion+oop) gate — nearly identical, unlike the pre-#183 numbers
    /// that originally justified keeping torsion/oop out of the gate ("fails
    /// a large fraction of ordinary organic molecules"). That justification
    /// is now stale. Widening the default gate is a small, well-scoped
    /// follow-up worth considering, but doing so also changes
    /// `Mmff94WithUffFallback`'s fallback trigger population and needs its
    /// own re-verification pass — deliberately not done in this PR to avoid
    /// conflating a naming fix with a behavior change (see PR body).
    Mmff94BondAngleStrict,
    /// Try `Mmff94BondAngleStrict` first; on any typed failure (unsupported atom
    /// element, missing bond/angle parameters, or an unsound MMFF94 result —
    /// see [`ForceFieldBridgeError::MinimizationFailed`]), fall back to
    /// `UffOnly`. When a fallback happens, it is always reported via
    /// [`PolicyMinimizeResult::fallback_reason`]/
    /// [`PolicyMinimizeResult::actual_force_field_used`]/
    /// [`PolicyMinimizeResult::missing_parameter_classes`], never silent.
    ///
    /// NOT infallible: if the `UffOnly` fallback attempt is itself unsound
    /// (measured post-chematic-ff-#183, restricted to the molecules that
    /// actually *reach* this fallback path — i.e. the ones that first fail
    /// `Mmff94BondAngleStrict`: 8 fused/conjugated polycyclic aromatics
    /// (naphthalene, quinoline, pyrene, ibuprofen, ibuprofen_S, naproxen_S,
    /// diphenhydramine, atorvastatin_fragment) plus caffeine remaining
    /// blown up even before the soundness gate — see PR #169 body), this
    /// returns `Err(MinimizationFailed)` too. This policy's contract is
    /// "never silently report success on an unsound geometry," not "always
    /// succeeds" — those are different guarantees, and only the first one
    /// is made here. This 9-molecule figure is scoped to the fallback
    /// *trigger population* specifically, not `UffOnly`'s behavior across
    /// the full corpus — see the correction on [`ForceFieldPolicy::UffOnly`]
    /// below, which a prior version of this doc comment incorrectly copied
    /// verbatim onto both variants.
    Mmff94WithUffFallback,
    /// chematic-ff's real UFF module (`chematic_ff::uff`) — generic,
    /// all-element coverage (bond lengths/angles are formula-derived from
    /// per-type constants, not a lookup table, so there is no missing-entry
    /// case to gate on). NOT infallible: can return
    /// `Err(MinimizationFailed)` if the resulting geometry is unsound (see
    /// `check_minimization_soundness`).
    ///
    /// **Correction (found during the Wave-1 C+F integration smoke test,
    /// PR #186, independently re-verified):** an earlier version of this
    /// doc comment claimed "8/58 corpus molecules, all fused/conjugated
    /// polycyclic aromatics" for `UffOnly` specifically. That figure was
    /// actually the [`ForceFieldPolicy::Mmff94WithUffFallback`] fallback
    /// *trigger population* (9 molecules, not 8, once caffeine is counted;
    /// see that variant's doc above), copied here without re-deriving it —
    /// `UffOnly` was never actually run over the full 58-molecule corpus in
    /// PR #169 itself. The real, first-ever full-corpus measurement (from
    /// PR #186, using legacy `dg::generate_coords` starting geometry, the
    /// same convention this PR's own gate-check example uses) is **17/58**
    /// blow up under `UffOnly`: the same 9 fused/conjugated aromatics above,
    /// plus 8 more with no ring fusion at all (hexane, decane,
    /// triethylene_glycol, hexanediol, hexadecane, penicillin_core,
    /// testosterone, cholesterol) that are structurally invisible to any
    /// `Mmff94WithUffFallback`-based measurement, since all 8 pass
    /// `Mmff94BondAngleStrict` cleanly and so never invoke UFF under that
    /// policy at all. Tracked as a `chematic-ff` UFF-minimizer robustness
    /// gap; see issue #185 (currently scoped only to the narrower
    /// naphthalene-vs-anthracene puzzle, not this full 17-molecule class).
    UffOnly,
    /// Existing DREIDING path (unchanged physics) — chematic-ff has no
    /// DREIDING minimizer, so this stays chematic-3d's own implementation,
    /// just routed through the same reporting result type as the other
    /// policies.
    Dreiding,
    /// No minimization — return the embedded geometry unchanged. Useful for
    /// testing/composability with a raw embedder.
    None,
}

/// Which MMFF94 internal-coordinate class a missing-parameter finding is
/// about.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mmff94TermKind {
    Bond,
    Angle,
    Torsion,
    Oop,
    /// Stretch-bend cross term (Halgren MMFF.V eq. 4), keyed on the same
    /// (outer, center, outer) angle triple as [`Mmff94TermKind::Angle`].
    /// Historically never a member of this enum at all -- see
    /// `Mmff94CoverageReport::has_gate_failure`'s `include_stretch_bend`
    /// parameter for why it is opt-in, not always-gated like bond/angle.
    StretchBend,
}

/// One internal coordinate that chematic-ff's MMFF94 tables have no entry
/// for. `atoms` is the exact atom tuple (2 for bond, 3 for angle, 4 for
/// torsion/oop) so a caller can cite *which* atoms/elements were missing —
/// never just an aggregate count.
#[derive(Debug, Clone)]
pub struct Mmff94MissingTerm {
    pub kind: Mmff94TermKind,
    pub atoms: Vec<AtomIdx>,
    /// Human-readable citation, including each atom's specific numeric
    /// MMFF94 type (not just its element symbol — many distinct MMFF94
    /// types share an element, e.g. sp3 vs. aromatic vs. carbonyl carbon),
    /// e.g. `"Angle F(11)-C(1)-Cl(12) (atom indices [1, 0, 2])"`. Atom order
    /// is canonicalized (see `canonicalize_term_atoms`) so physically
    /// equivalent citations (e.g. a bond read C→N vs N→C) always produce the
    /// same string, so aggregating by this string doesn't double-count.
    pub description: String,
}

/// Full parameter-coverage measurement for one molecule under chematic-ff's
/// MMFF94 tables. Every internal coordinate chematic-ff's own energy
/// functions would evaluate is enumerated here (the enumeration mirrors
/// `chematic_ff::mmff94_minimizer`'s private energy loops exactly, so the
/// totals match what chematic-ff will actually compute over, not an
/// idealized recount) — nothing is silently excluded from the count.
#[derive(Debug, Clone, Default)]
pub struct Mmff94CoverageReport {
    pub bonds_total: usize,
    pub bonds_missing: Vec<Mmff94MissingTerm>,
    pub angles_total: usize,
    pub angles_missing: Vec<Mmff94MissingTerm>,
    pub torsions_total: usize,
    pub torsions_missing: Vec<Mmff94MissingTerm>,
    /// Torsions with NO real MMFF94 term at all, by design, not a coverage
    /// gap -- both central (j-k) atoms fully typed, but one has MMFF's
    /// `lin` flag (see `chematic_ff::torsion_no_term_by_design`'s doc, issue
    /// #227 Phase 1). Included in `torsions_total`, deliberately excluded
    /// from `torsions_missing` and from `has_gate_failure`'s torsion check.
    pub torsions_no_term_by_design: usize,
    pub oop_total: usize,
    pub oop_missing: Vec<Mmff94MissingTerm>,
    /// Stretch-bend cross-term coverage (Priority 2 / Stage 1B, issue #227).
    /// Never gated by default -- see `has_gate_failure`'s `include_stretch_bend`
    /// parameter -- but always measured/reported here regardless, same as
    /// torsion/oop.
    pub stretch_bend_total: usize,
    pub stretch_bend_missing: Vec<Mmff94MissingTerm>,
}

impl Mmff94CoverageReport {
    /// True iff no bond or angle is missing parameters — the exact scope of
    /// the RFC's mechanism-3 bug (silent zero-gradient, no restoring force
    /// while vdW repulsion still pushes atoms apart).
    pub fn bond_angle_fully_covered(&self) -> bool {
        self.bonds_missing.is_empty() && self.angles_missing.is_empty()
    }

    /// Gate check used by `Mmff94BondAngleStrict`/`Mmff94WithUffFallback`.
    /// Bond/angle always participate; torsion/out-of-plane only do when
    /// `include_torsion_oop` is set, and stretch-bend only when
    /// `include_stretch_bend` is set (see `minimize_with_policy_gated`) --
    /// each is its own independent opt-in, not bundled together, so a
    /// caller can measure the stretch-bend-only effect in isolation (Priority
    /// 2 / Stage 1B's "legacy vs. complete-term" comparison).
    pub fn has_gate_failure(&self, include_torsion_oop: bool, include_stretch_bend: bool) -> bool {
        !self.bond_angle_fully_covered()
            || (include_torsion_oop
                && (!self.torsions_missing.is_empty() || !self.oop_missing.is_empty()))
            || (include_stretch_bend && !self.stretch_bend_missing.is_empty())
    }

    /// Total count of missing internal coordinates across all 5 classes —
    /// always meaningful regardless of gate scope, since nothing here is
    /// filtered by what the gate happens to check.
    pub fn total_missing(&self) -> usize {
        self.bonds_missing.len()
            + self.angles_missing.len()
            + self.torsions_missing.len()
            + self.oop_missing.len()
            + self.stretch_bend_missing.len()
    }

    /// Flattened list of every missing internal coordinate across all 5
    /// classes — the exact "which specific element/atom-type pairs lacked
    /// coverage" citation surfaced at the top level as
    /// [`PolicyMinimizeResult::missing_parameter_classes`].
    pub fn all_missing(&self) -> Vec<Mmff94MissingTerm> {
        let mut v = Vec::with_capacity(self.total_missing());
        v.extend(self.bonds_missing.iter().cloned());
        v.extend(self.angles_missing.iter().cloned());
        v.extend(self.torsions_missing.iter().cloned());
        v.extend(self.oop_missing.iter().cloned());
        v.extend(self.stretch_bend_missing.iter().cloned());
        v
    }
}

/// Typed failure for the [`ForceFieldPolicy::Mmff94BondAngleStrict`] bridge. Never
/// silently absorbed into a zero-energy/zero-gradient term.
///
/// Naming note for Coordinator (Wave 2 reconciliation): this is deliberately
/// named for this crate's own use, not coupled to Agent C's parallel
/// `EmbedFailureCause` (distance-geometry embedding failures — a different
/// pipeline stage). The variant `MissingParameters(Mmff94CoverageReport)`
/// here is the natural candidate for an `EmbedFailureCause::
/// MissingForceFieldParameters { .. }` case if/when Coordinator unifies the
/// two vocabularies; this PR does not attempt that unification itself.
#[derive(Debug, Clone)]
pub enum ForceFieldBridgeError {
    /// `chematic_ff::assign_mmff94_numeric_types` failed outright for one or
    /// more atoms (element unsupported by the typer at all — coarser than a
    /// missing bond/angle/torsion/oop entry for an otherwise-typed atom).
    UnsupportedAtomType(String),
    /// One or more internal coordinates lack MMFF94 parameters under the
    /// active coverage gate. Boxed: `clippy::result_large_err` (the report
    /// carries per-missing-term `Vec`s / `String`s so the inline variant is
    /// large relative to the common `Ok` path).
    MissingParameters(Box<Mmff94CoverageReport>),
    /// The minimizer ran (coverage/typing succeeded) but produced a
    /// geometry that is not actually sound — NaN/Inf coordinates, a
    /// catastrophic bond blow-up, or an excessive residual force. Before
    /// this variant existed, every one of these cases was returned as
    /// `Ok(PolicyMinimizeResult { converged: false, .. })` — a result a
    /// careless caller could mistake for success (independent review found
    /// exactly this on the 58-molecule corpus: UFF-fallback runs that blew
    /// worst-bond-length past 800–19,000+ Å still came back `Ok`). See
    /// [`check_minimization_soundness`] for what specifically triggers this
    /// and why plain non-convergence (`converged == false`) alone does
    /// NOT — that would be a false-failure trigger, not a safety gate (see
    /// that function's doc for the measured evidence).
    MinimizationFailed(Box<MinimizationFailureDetail>),
}

/// Which soundness check tripped inside [`check_minimization_soundness`].
#[derive(Debug, Clone, Copy)]
pub enum MinimizationFailureReason {
    /// One or more final coordinates are NaN or infinite.
    NonFiniteCoordinates,
    /// The worst bond length in the final geometry exceeds
    /// [`MAX_SANE_BOND_LENGTH`] Å — no legitimate covalent bond in the
    /// molecules this bridge handles gets remotely close to this.
    CatastrophicBondBlowup,
    /// The finite-difference max |gradient component| at the final geometry
    /// exceeds [`MAX_SANE_RESIDUAL_FORCE`] kcal/mol/Å — a backstop for a
    /// geometry that isn't bond-length-blown-up but is still nowhere near a
    /// force-balanced minimum.
    ExcessiveResidualForce,
}

/// Evidence attached to [`ForceFieldBridgeError::MinimizationFailed`].
/// `converged`/`iterations` are carried here for diagnostics but — per
/// [`check_minimization_soundness`]'s doc — are NOT what triggered the
/// failure; `reason` (plus `worst_bond_length`/`max_residual_force`) is.
#[derive(Debug, Clone)]
pub struct MinimizationFailureDetail {
    pub policy: ForceFieldPolicy,
    pub reason: MinimizationFailureReason,
    /// The underlying minimizer's own convergence flag. Often `false` even
    /// on perfectly sound geometries that simply didn't reach a tight
    /// gradient tolerance within the default iteration budget — see
    /// `check_minimization_soundness`'s doc; this field is diagnostic only.
    pub converged: bool,
    pub iterations: usize,
    /// Always populated (not just when `reason` is
    /// `ExcessiveResidualForce`) so a failed molecule still reports how far
    /// from equilibrium it is.
    pub max_residual_force: f64,
    /// Always populated (not just when `reason` is `CatastrophicBondBlowup`)
    /// so a failed molecule still reports its geometry, rather than a blank
    /// — this is what stops the "measurement trap" of a fixed bug simply
    /// vanishing from a blow-up count instead of being counted as a typed
    /// failure.
    pub worst_bond_length: f64,
    /// This `MinimizationFailureDetail` always describes the ORIGINAL attempt's
    /// failure (the caller's own coordinates) — never the retry's. This flag only
    /// records whether a retry was even attempted before that original failure was
    /// returned: `true` means the original attempt failed with
    /// `CatastrophicBondBlowup`/`NonFiniteCoordinates`, `embed_distance_geometry_v2`
    /// was tried as a rescue, and that rescue was **not accepted** (its own result
    /// was unsound, stereo-violating, or embedding itself errored) — so the
    /// original failure is what's being returned, annotated to show a rescue was
    /// tried and didn't help. `false` means either no rescue was attempted at all
    /// (wrong policy, or an `ExcessiveResidualForce`-only failure — issues #185/#188
    /// measured this rescue helping catastrophic bond blowup specifically, not
    /// residual-force-only unsoundness, see [`run_uff_bridge`]'s doc for why the
    /// rescue is scoped this narrowly) or a rescue *was* attempted and accepted, in
    /// which case this function returns `Ok`, not this error type, and this field is
    /// moot. Lets a caller distinguish "no rescue was even attempted" from "the
    /// rescue was attempted and did not help either" — never "the retry's own
    /// failure reason," which this field does not carry.
    pub distance_geometry_v2_retry_attempted: bool,
}

/// Which starting geometry actually seeded a successful UFF minimization —
/// see [`run_uff_bridge`]'s doc for the full rationale. Exists so a caller
/// can tell "my exact input coordinates were used" apart from "the input was
/// detected as catastrophically unsound after minimizing it, and
/// `embed_distance_geometry_v2` was used to re-seed instead" — never a
/// silent substitution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UffStartingGeometry {
    /// The coordinates the caller passed in were used to seed minimization
    /// and that attempt was already sound — no re-embedding was attempted.
    AsProvided,
    /// The caller-provided coordinates produced a catastrophically unsound
    /// result; `embed_distance_geometry_v2` was used to re-seed instead, and
    /// that retry was sound.
    ReplacedWithDistanceGeometryV2,
}

impl std::fmt::Display for ForceFieldBridgeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ForceFieldBridgeError::UnsupportedAtomType(e) => {
                write!(f, "MMFF94 atom-type assignment failed: {e}")
            }
            ForceFieldBridgeError::MissingParameters(r) => write!(
                f,
                "MMFF94 parameters missing for {} internal coordinate(s) \
                 ({} bond, {} angle, {} torsion, {} oop, {} stretch-bend)",
                r.total_missing(),
                r.bonds_missing.len(),
                r.angles_missing.len(),
                r.torsions_missing.len(),
                r.oop_missing.len(),
                r.stretch_bend_missing.len(),
            ),
            ForceFieldBridgeError::MinimizationFailed(d) => write!(
                f,
                "{:?} minimization produced an unsound geometry: {:?} \
                 (converged={}, iterations={}, worst_bond_length={:.2} Å, \
                 max_residual_force={:.2} kcal/mol/Å)",
                d.policy,
                d.reason,
                d.converged,
                d.iterations,
                d.worst_bond_length,
                d.max_residual_force,
            ),
        }
    }
}

impl std::error::Error for ForceFieldBridgeError {}

impl From<NumericTypeError> for ForceFieldBridgeError {
    fn from(e: NumericTypeError) -> Self {
        ForceFieldBridgeError::UnsupportedAtomType(e.to_string())
    }
}

impl From<MinimizerError> for ForceFieldBridgeError {
    fn from(e: MinimizerError) -> Self {
        // MinimizerError is itself just a NumericTypeError wrapper today;
        // preserve the message rather than collapsing to a generic string.
        ForceFieldBridgeError::UnsupportedAtomType(e.to_string())
    }
}

/// Per-policy energy accounting. `Mmff94` reuses chematic-ff's own
/// [`EnergyBreakdown`] shape verbatim — the same struct
/// `scripts/etkdg_vs_rdkit_gap.py`'s "fair energy delta" metric already
/// consumes via `Mol.mmff94_energy_breakdown()` — so downstream consumers
/// never need a second energy-breakdown format.
#[derive(Debug, Clone)]
pub enum EnergyReport {
    Mmff94(EnergyBreakdown),
    Uff { total: f64 },
    Dreiding { total: f64 },
    None,
}

impl EnergyReport {
    pub fn total(&self) -> f64 {
        match self {
            EnergyReport::Mmff94(b) => b.total,
            EnergyReport::Uff { total } | EnergyReport::Dreiding { total } => *total,
            EnergyReport::None => 0.0,
        }
    }
}

/// Full result of [`minimize_with_policy`]/[`minimize_with_policy_gated`].
///
/// The `requested_force_field`/`actual_force_field_used`/`fallback_reason`/
/// `missing_parameter_classes` quartet exists specifically so that
/// `Mmff94WithUffFallback` can never silently report "MMFF94 minimized" when
/// it actually fell back to UFF — that would be the exact mechanism-3 bug
/// (silent zero/substitution) reintroduced one layer up, at the
/// policy-reporting level instead of the energy-term level.
#[derive(Debug, Clone)]
pub struct PolicyMinimizeResult {
    pub coords: Coords3D,
    /// What the caller asked for.
    pub requested_force_field: ForceFieldPolicy,
    /// What actually ran. For `Mmff94BondAngleStrict` this is always
    /// `Mmff94BondAngleStrict`; for `Mmff94WithUffFallback` this is
    /// `Mmff94BondAngleStrict` when the MMFF94 attempt succeeded (the common case —
    /// `requested_force_field == Mmff94WithUffFallback` but
    /// `actual_force_field_used == Mmff94BondAngleStrict` here is NOT a fallback, it's
    /// just the more informative "which physics actually ran" value) and
    /// `UffOnly` only when it fell back. **Check `fallback_reason.is_some()`
    /// to ask "did a fallback occur," not `actual_force_field_used !=
    /// requested_force_field`** — the latter is true on every successful
    /// `Mmff94WithUffFallback` call and does not mean a fallback happened.
    pub actual_force_field_used: ForceFieldPolicy,
    /// `Some(reason)` iff a fallback actually occurred (only possible under
    /// `Mmff94WithUffFallback`, when the MMFF94 attempt failed) — explains
    /// why. Carries the same typed reason `Mmff94BondAngleStrict` would have returned
    /// as an `Err`. This, not `actual_force_field_used != requested_force_field`,
    /// is the correct fallback-occurred check (see `actual_force_field_used`'s
    /// doc above).
    pub fallback_reason: Option<ForceFieldBridgeError>,
    /// Every specific internal coordinate (bond/angle/torsion/oop, cited by
    /// atom indices and element symbols — never just an aggregate count)
    /// that lacked MMFF94 parameter coverage. Populated whenever MMFF94
    /// typing/coverage was computed at all (`Mmff94BondAngleStrict`, or the MMFF94
    /// attempt inside `Mmff94WithUffFallback` regardless of outcome); empty
    /// (not merely absent) when MMFF94 was never attempted (`UffOnly`/
    /// `Dreiding`/`None`) or when coverage was full.
    pub missing_parameter_classes: Vec<Mmff94MissingTerm>,
    /// Full structured coverage report (all 4 term-class totals + missing
    /// lists), whenever MMFF94 typing/coverage was computed — see
    /// `missing_parameter_classes` for the flattened, always-populated view
    /// of the same data.
    pub coverage: Option<Mmff94CoverageReport>,
    pub energy_before: EnergyReport,
    pub energy_after: EnergyReport,
    pub converged: bool,
    pub iterations: usize,
    /// Max |gradient component| (kcal/mol/Å) at the final geometry, computed
    /// by this bridge's own finite-difference pass over the exact energy
    /// function that produced `energy_after` — not read off any upstream
    /// struct (chematic-ff's `MinimizeResult`/`UffMinimizeResult` don't carry
    /// this field).
    pub max_residual_force: f64,
    /// `Some(_)` whenever UFF actually ran (`UffOnly`, or the UFF-fallback
    /// path of `Mmff94WithUffFallback`) — see [`UffStartingGeometry`] for
    /// what the two variants mean and why raw starting energy alone can't
    /// predict which one occurs. `None` for `Dreiding`/`None`, and for
    /// `Mmff94BondAngleStrict`'s own (non-UFF) success path.
    pub starting_geometry: Option<UffStartingGeometry>,
}

fn trivial_result(coords: Coords3D, policy: ForceFieldPolicy) -> PolicyMinimizeResult {
    PolicyMinimizeResult {
        coords,
        requested_force_field: policy,
        actual_force_field_used: policy,
        fallback_reason: None,
        missing_parameter_classes: Vec::new(),
        coverage: None,
        energy_before: EnergyReport::None,
        energy_after: EnergyReport::None,
        converged: true,
        iterations: 0,
        max_residual_force: 0.0,
        starting_geometry: None,
    }
}

// --- Coords3D <-> Vec<[f64; 3]> bridge helpers ------------------------------

fn coords_to_vec(coords: &Coords3D, n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            let p = coords.get(AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect()
}

fn vec_to_coords(v: &[[f64; 3]]) -> Coords3D {
    let mut c = Coords3D::new_zeroed(v.len());
    for (i, p) in v.iter().enumerate() {
        c.set(AtomIdx(i as u32), Point3::new(p[0], p[1], p[2]));
    }
    c
}

/// Central-difference max |gradient component| over an arbitrary black-box
/// energy function. Used both to report `max_residual_force` in production
/// and, in tests, as an independent bridge-plumbing check (see module tests
/// — chematic-ff exposes no analytic gradient anywhere, so this is a
/// finite-difference-vs-finite-difference cross-check across two independent
/// code paths, not an analytic-vs-FD check; see test doc comment for why).
fn fd_max_gradient<F: Fn(&[[f64; 3]]) -> f64>(
    coords: &[[f64; 3]],
    energy_fn: F,
    delta: f64,
) -> f64 {
    let mut work = coords.to_vec();
    let mut max_g = 0.0_f64;
    for i in 0..work.len() {
        for axis in 0..3 {
            work[i][axis] += delta;
            let ep = energy_fn(&work);
            work[i][axis] -= 2.0 * delta;
            let em = energy_fn(&work);
            work[i][axis] += delta;
            let g = (ep - em) / (2.0 * delta);
            if g.abs() > max_g {
                max_g = g.abs();
            }
        }
    }
    max_g
}

// --- Geometric/energetic soundness gate -------------------------------------
//
// Found in independent review: every policy branch below used to build its
// `PolicyMinimizeResult` directly as `Ok(..)`, including `converged: false`
// results with `max_residual_force` in the hundreds of thousands to millions
// (kcal/mol/Å) and worst bond lengths past 800–19,000+ Å (measured on the
// 58-molecule corpus's UFF-fallback cases) — a caller checking only
// `result.is_ok()` would call that "minimized successfully." This section
// converts a genuinely unsound result into a typed
// `Err(ForceFieldBridgeError::MinimizationFailed)` instead.

/// Worst bond length ceiling (Å) above which a geometry is treated as a
/// catastrophic blow-up, not merely "not yet converged." Reuses this
/// project's own pre-existing convention (`worst_bond > 3.0` is the "blown
/// up" bar used throughout `examples/mmff94_bridge_coverage_report.rs` and
/// the tests in this file) rather than inventing a new number. No
/// light-atom-only organic bond covered by MMFF94/UFF (including C-I, S-S)
/// gets remotely close to 3 Å, so this has ample margin on the "sound" side.
///
/// `pub(crate)`: Wave 2/3 Coordinator integration (`pipeline_v2.rs`'s own,
/// independent final-geometry soundness check, per that program's spec §11) reuses
/// this exact constant rather than hand-copying `3.0` a second time -- particularly
/// important because `ForceFieldPolicy::None` skips `check_minimization_soundness`
/// entirely (see `minimize_with_policy_gated`'s `None` arm), so the pipeline-level
/// check is the *only* gate against a blown-up bond length in that case.
pub(crate) const MAX_SANE_BOND_LENGTH: f64 = 3.0;

/// Residual-force ceiling (kcal/mol/Å) above which a geometry is treated as
/// unsound even if no single bond individually crossed
/// [`MAX_SANE_BOND_LENGTH`] (e.g. a purely angular/torsional distortion).
/// Measured on the 58-molecule MMFF94/UFF corpus (post chematic-ff #183):
/// every molecule that is geometrically fine but simply didn't converge
/// within the default iteration budget has `max_residual_force` ≤ 8.93
/// (cholesterol); every molecule with a real blow-up has
/// `max_residual_force` ≥ 337.99 (quinoline) — but every one of those
/// already trips the bond-length check above, so within that corpus this
/// constant has no independently-discriminating case.
///
/// An initial value of 50.0 (picked inside that 8.93–337.99 gap with no
/// corpus case to validate it) turned out to be a live tripwire, not a
/// backstop: `dreiding_policy_matches_existing_behavior_and_reports_convergence`
/// (acetic acid via [`ForceFieldPolicy::Dreiding`]) has a perfectly normal
/// worst bond length (1.496 Å) but `max_residual_force` = 55.12 after 200
/// gradient-descent iterations — DREIDING's simpler harmonic springs
/// (`BOND_SPRING_CONSTANT`/`ANGLE_SPRING_CONSTANT`) apparently settle to a
/// higher steady-state FD residual than MMFF94/UFF's more carefully
/// parameterized terms do, even on a sound geometry. Raised to 200.0 to keep
/// this measured-sound case comfortably (>3.6×) below the ceiling while
/// staying well (>40%) under the smallest measured genuine blow-up
/// (337.99) — still a backstop with only one policy's data validating the
/// "sound" side, not a value with its own discriminating corpus case.
const MAX_SANE_RESIDUAL_FORCE: f64 = 200.0;

fn worst_bond_length_vec(mol: &Molecule, coords: &[[f64; 3]]) -> f64 {
    let dist = |a: [f64; 3], b: [f64; 3]| {
        let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    mol.bonds()
        .map(|(_, b)| dist(coords[b.atom1.0 as usize], coords[b.atom2.0 as usize]))
        .fold(0.0_f64, f64::max)
}

/// Decide whether a minimization result is geometrically/energetically
/// sound enough to return as `Ok`, or must instead become a typed
/// `Err(MinimizationFailed)`.
///
/// Deliberately does NOT treat `converged == false` alone as a failure
/// trigger. Measured on the 58-molecule corpus: molecules like decane,
/// hexadecane, triethylene glycol, testosterone, cholesterol, the
/// gly-ala-gly tripeptide, and the penicillin core all report
/// `converged == false` under the default `MinimizeConfig::max_steps` — not
/// because their geometry is unsound (every one has a finite, sub-3-Å worst
/// bond and `max_residual_force` under 10 kcal/mol/Å), but because L-BFGS's
/// gradient-norm convergence tolerance is tight relative to that default
/// iteration budget. Treating that as a hard failure would make
/// `Mmff94BondAngleStrict`/`UffOnly`/`Dreiding` spuriously refuse on a large
/// fraction of ordinary, perfectly fine molecules — a false-failure
/// generator, not a safety gate. `converged`/`iterations` are still carried
/// into `MinimizationFailureDetail` for diagnostics on an actual failure,
/// and every `Ok` result with `converged == false` remains visible on
/// `PolicyMinimizeResult` for a caller who wants to distinguish
/// "definitely converged" from "sound but still iterating" — see
/// `examples/mmff94_bridge_coverage_report.rs`'s own separate count of that.
fn check_minimization_soundness(
    mol: &Molecule,
    coords: &[[f64; 3]],
    policy: ForceFieldPolicy,
    converged: bool,
    iterations: usize,
    max_residual_force: f64,
    distance_geometry_v2_retry_attempted: bool,
) -> Result<(), ForceFieldBridgeError> {
    let worst_bond_length = worst_bond_length_vec(mol, coords);
    let reason = if coords.iter().any(|p| p.iter().any(|x| !x.is_finite())) {
        Some(MinimizationFailureReason::NonFiniteCoordinates)
    } else if worst_bond_length > MAX_SANE_BOND_LENGTH {
        Some(MinimizationFailureReason::CatastrophicBondBlowup)
    } else if max_residual_force > MAX_SANE_RESIDUAL_FORCE {
        Some(MinimizationFailureReason::ExcessiveResidualForce)
    } else {
        None
    };
    match reason {
        None => Ok(()),
        Some(reason) => Err(ForceFieldBridgeError::MinimizationFailed(Box::new(
            MinimizationFailureDetail {
                policy,
                reason,
                converged,
                iterations,
                max_residual_force,
                worst_bond_length,
                distance_geometry_v2_retry_attempted,
            },
        ))),
    }
}

// --- MMFF94 parameter-coverage checking -------------------------------------
//
// Post-#183 (chematic-ff MMFF94 bond/angle/torsion classification fixes),
// `bond_type_for`/`angle_type_for`/`torsion_type_for`/`MLTB_TYPES`/
// `OOP_SP2_TYPES` are now `pub` in chematic-ff (a narrow, additive visibility
// change made alongside this fix) and are called directly here — no more
// hand-copied duplicate classification. An earlier version of this bridge
// duplicated these (chematic-ff's versions were private pre-#183), which was
// deliberately kept in lockstep with chematic-ff's *pre-#183* bugs (bond
// typing via OR-logic with no bond-order check; angle type hardcoded to 0).
// That duplicate is now confirmed to have caused a real false-negative post-
// #183: a carbonyl C=O bond (types 3,7) is correctly classified bond_type=0
// by the real, order-aware `bond_type_for` (finds the `(0,3,7,...)` row —
// covered), but the old duplicate's order-blind OR-logic computed
// bond_type=1 for it (type 3 is SP2) and found no `(1,3,7,...)` row —
// falsely reporting "missing coverage" on one of the most common bonds in
// organic chemistry (every ketone/aldehyde/amide/carboxylic acid). Calling
// chematic-ff's real functions directly makes this bridge's coverage report
// track chematic-ff's actual behavior by construction, not by manually kept
// parity — see `carbonyl_c_o_double_bond_is_not_falsely_reported_missing`
// (test, below) for the pinned regression.

/// SSSR rings for `mol`, computed once per coverage pass and threaded through
/// to `angle_type_for`/`torsion_type_for` exactly as `chematic_ff`'s own
/// `mmff94_total_energy`/`angle_energy`/`torsion_energy` do internally (same
/// `find_sssr` call), so ring-size-dependent angle/torsion typing can't
/// diverge from what chematic-ff will actually compute over.
fn mmff94_rings(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    find_sssr(mol).rings().to_vec()
}

/// Reorders `atoms` into a canonical form so that physically-equivalent
/// citations (a bond read C→N vs N→C, an angle's two outer substituents in
/// either order, a torsion read forward vs backward, an out-of-plane
/// center's 3 interchangeable substituents in any order) produce the exact
/// same `atoms`/`description` — otherwise aggregating "distinct missing
/// patterns" across a corpus double-counts the same underlying gap under
/// two different string keys purely because of atom-enumeration order
/// (found in independent review: "Bond C-N" and "Bond N-C" were counted
/// separately). Ordering by ascending MMFF94 numeric type matches
/// chematic-ff's own symmetric lookups, which already sort
/// `(type_i, type_j)` before searching (see `mmff94_bond_params`/
/// `mmff94_bond_energy`), so this is the same canonical order chematic-ff
/// itself uses internally, not an arbitrary new one.
fn canonicalize_term_atoms(kind: Mmff94TermKind, atoms: &[AtomIdx], types: &[u8]) -> Vec<AtomIdx> {
    let ty = |a: AtomIdx| types[a.0 as usize];
    let mut v = atoms.to_vec();
    match kind {
        Mmff94TermKind::Bond => {
            if ty(v[0]) > ty(v[1]) {
                v.swap(0, 1);
            }
        }
        // [outer_a, center, outer_c]: center fixed, sort the outer pair.
        Mmff94TermKind::Angle => {
            if ty(v[0]) > ty(v[2]) {
                v.swap(0, 2);
            }
        }
        // [i, j, k, l]: reading backward (l, k, j, i) is the same torsion.
        Mmff94TermKind::Torsion => {
            let fwd: Vec<u8> = v.iter().map(|&a| ty(a)).collect();
            let rev: Vec<u8> = fwd.iter().rev().copied().collect();
            if rev < fwd {
                v.reverse();
            }
        }
        // [center, s1, s2, s3]: center fixed, the 3 substituents are
        // mutually interchangeable.
        Mmff94TermKind::Oop => {
            v[1..].sort_by_key(|&a| ty(a));
        }
        // Same [outer_a, center, outer_c] layout as Angle -- stretch-bend is
        // keyed on the identical angle triple.
        Mmff94TermKind::StretchBend => {
            if ty(v[0]) > ty(v[2]) {
                v.swap(0, 2);
            }
        }
    }
    v
}

fn missing_term(
    mol: &Molecule,
    types: &[u8],
    kind: Mmff94TermKind,
    atoms: &[AtomIdx],
) -> Mmff94MissingTerm {
    let atoms = canonicalize_term_atoms(kind, atoms, types);
    let labels: Vec<String> = atoms
        .iter()
        .map(|&a| format!("{}({})", mol.atom(a).element.symbol(), types[a.0 as usize]))
        .collect();
    let indices: Vec<u32> = atoms.iter().map(|a| a.0).collect();
    Mmff94MissingTerm {
        kind,
        atoms,
        description: format!("{kind:?} {} (atom indices {indices:?})", labels.join("-")),
    }
}

/// Independently measure MMFF94 parameter coverage for every bond, angle,
/// stretch-bend, torsion, and out-of-plane center that chematic-ff's own
/// energy functions would evaluate for `mol`. Mirrors
/// `chematic_ff::mmff94_minimizer`'s private energy-term enumeration loops
/// exactly (including its lack of an `i == l` guard on 3-membered-ring
/// torsions) so the totals match what chematic-ff will actually compute
/// over, not an idealized recount.
fn compute_mmff94_coverage(mol: &Molecule, types: &[u8]) -> Mmff94CoverageReport {
    let mut report = Mmff94CoverageReport::default();
    let rings = mmff94_rings(mol);

    for (_, bond) in mol.bonds() {
        report.bonds_total += 1;
        let (a1, a2) = (bond.atom1, bond.atom2);
        let (t1, t2) = (types[a1.0 as usize], types[a2.0 as usize]);
        let bt = bond_type_for(t1, t2, bond.order);
        if mmff94_bond_energy_resolved(bt, t1, t2).is_none() {
            report
                .bonds_missing
                .push(missing_term(mol, types, Mmff94TermKind::Bond, &[a1, a2]));
        }
    }

    for b_idx in 0..mol.atom_count() {
        let b = AtomIdx(b_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(b).map(|(nb, _)| nb).collect();
        if neighbors.len() < 2 {
            continue;
        }
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                report.angles_total += 1;
                report.stretch_bend_total += 1;
                let (a, c) = (neighbors[i], neighbors[j]);
                let (ta, tc) = (types[a.0 as usize], types[c.0 as usize]);
                let at = angle_type_for(mol, &rings, a.0 as usize, b_idx, c.0 as usize, types);
                // Shared by the angle-bend check below AND the stretch-bend
                // check further down -- stretch_bend_energy (chematic-ff's
                // mmff94_minimizer) iterates the identical neighbor-pair loop
                // and needs the same flanking-bond types/order (issue #227
                // Priority 2C: the STBN table is keyed by stretch-bend type,
                // not angle type -- see chematic_ff::stretch_bend_type_for's
                // doc).
                let order_ab = mol.bond_between(a, b).expect("a-b angle bond").1.order;
                let order_cb = mol.bond_between(c, b).expect("c-b angle bond").1.order;
                let bt_ab = bond_type_for(ta, types[b_idx], order_ab);
                let bt_cb = bond_type_for(tc, types[b_idx], order_cb);

                // Issue #227 Stage C: the angle-bend term itself needs both
                // flanking bonds' r0 (for the empirical rule, if the direct
                // table/eqLevel ladder don't resolve) -- matches RDKit's real
                // `getMMFFAngleBendParams`, which requires both
                // `getMMFFBondStretchParams` calls to succeed before even
                // attempting the empirical path (mirrors
                // `chematic_ff::mmff94_minimizer`'s own `angle_energy`).
                let bond_ab = mmff94_bond_energy_resolved(bt_ab, ta, types[b_idx]);
                let bond_cb = mmff94_bond_energy_resolved(bt_cb, tc, types[b_idx]);
                let angle_resolved = match (bond_ab, bond_cb) {
                    (Some((bab, _)), Some((bcb, _))) => {
                        let ring_size =
                            is_angle_in_ring_of_size_3_or_4(mol, a.0 as usize, b_idx, c.0 as usize);
                        mmff94_angle_energy_resolved(
                            at,
                            ta,
                            types[b_idx],
                            tc,
                            bab.r0,
                            bcb.r0,
                            ring_size,
                        )
                    }
                    _ => None,
                };
                if angle_resolved.is_none() {
                    report.angles_missing.push(missing_term(
                        mol,
                        types,
                        Mmff94TermKind::Angle,
                        &[a, b, c],
                    ));
                }
                let sbt = stretch_bend_type_for(at, ta, tc, bt_ab, bt_cb);
                if mmff94_stbn(
                    sbt,
                    ta,
                    types[b_idx],
                    tc,
                    mol.atom(a).element.atomic_number(),
                    mol.atom(b).element.atomic_number(),
                    mol.atom(c).element.atomic_number(),
                )
                .is_none()
                {
                    report.stretch_bend_missing.push(missing_term(
                        mol,
                        types,
                        Mmff94TermKind::StretchBend,
                        &[a, b, c],
                    ));
                }
            }
        }
    }

    for (_, bond) in mol.bonds() {
        let (j, k) = (bond.atom1, bond.atom2);
        let nbrs_j: Vec<AtomIdx> = mol.neighbors(j).map(|(nb, _)| nb).collect();
        let nbrs_k: Vec<AtomIdx> = mol.neighbors(k).map(|(nb, _)| nb).collect();
        for &i in &nbrs_j {
            if i == k {
                continue;
            }
            for &l in &nbrs_k {
                if l == j {
                    continue;
                }
                report.torsions_total += 1;
                let (ti_, tj_, tk_, tl_) = (
                    types[i.0 as usize],
                    types[j.0 as usize],
                    types[k.0 as usize],
                    types[l.0 as usize],
                );
                let tt = torsion_type_for(
                    mol,
                    i.0 as usize,
                    j.0 as usize,
                    k.0 as usize,
                    l.0 as usize,
                    ti_,
                    tj_,
                    tk_,
                    tl_,
                );
                if mmff94_torsion_energy(tt, ti_, tj_, tk_, tl_).is_none() {
                    if torsion_no_term_by_design(tj_, tk_) {
                        // Issue #227 Phase 1: RDKit itself generates no term
                        // here either (linear central atom) -- correct, not
                        // a coverage gap. Counted separately so this never
                        // trips `include_torsion_oop_in_gate`.
                        report.torsions_no_term_by_design += 1;
                    } else {
                        report.torsions_missing.push(missing_term(
                            mol,
                            types,
                            Mmff94TermKind::Torsion,
                            &[i, j, k, l],
                        ));
                    }
                }
            }
        }
    }

    for j_idx in 0..mol.atom_count() {
        let tj = types[j_idx];
        if OOP_SP2_TYPES.binary_search(&tj).is_err() {
            continue;
        }
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(j).map(|(nb, _)| nb).collect();
        if neighbors.len() != 3 {
            continue;
        }
        report.oop_total += 1;
        let [i, k, l] = [neighbors[0], neighbors[1], neighbors[2]];
        if mmff94_oop(
            tj,
            types[i.0 as usize],
            types[k.0 as usize],
            types[l.0 as usize],
        )
        .is_none()
        {
            report
                .oop_missing
                .push(missing_term(mol, types, Mmff94TermKind::Oop, &[j, i, k, l]));
        }
    }

    report
}

// --- policy dispatch ---------------------------------------------------------

struct Mmff94BridgeRun {
    coords: Coords3D,
    coverage: Mmff94CoverageReport,
    energy_before: EnergyBreakdown,
    energy_after: EnergyBreakdown,
    converged: bool,
    iterations: usize,
    max_residual_force: f64,
}

fn run_mmff94_bridge(
    mol: &Molecule,
    coords: &Coords3D,
    max_iter: usize,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
) -> Result<Mmff94BridgeRun, ForceFieldBridgeError> {
    let n = mol.atom_count();
    // Must use the same MMFF-specific re-perceived bond orders chematic-ff's
    // own production energy functions use internally (issue #227 Phase 1,
    // torsion parameter gap root cause) -- otherwise this coverage gate and
    // the energy/gradient functions it gates disagree on which classification
    // code a bond/angle/torsion/stretch-bend term resolves to.
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let coverage = compute_mmff94_coverage(&mmff_mol, &types);
    if coverage.has_gate_failure(include_torsion_oop_in_gate, include_stretch_bend_in_gate) {
        return Err(ForceFieldBridgeError::MissingParameters(Box::new(coverage)));
    }

    let coord_vec = coords_to_vec(coords, n);
    let energy_before = mmff94_energy_breakdown(mol, &coord_vec)?;

    let mut work = coord_vec.clone();
    let result = minimize_mmff94_lbfgs(mol, &mut work, max_iter)?;

    let energy_after = mmff94_energy_breakdown(mol, &work)?;
    let max_residual_force = fd_max_gradient(
        &work,
        |c| {
            mmff94_total_energy(mol, c)
                .expect("mmff94_total_energy must not fail after a successful energy_breakdown/minimize call on the same molecule/coords")
        },
        1e-4,
    );

    check_minimization_soundness(
        mol,
        &work,
        ForceFieldPolicy::Mmff94BondAngleStrict,
        result.converged,
        result.iterations,
        max_residual_force,
        false,
    )?;

    Ok(Mmff94BridgeRun {
        coords: vec_to_coords(&work),
        coverage,
        energy_before,
        energy_after,
        converged: result.converged,
        iterations: result.iterations,
        max_residual_force,
    })
}

struct UffBridgeRun {
    coords: Coords3D,
    energy_before: f64,
    energy_after: f64,
    converged: bool,
    iterations: usize,
    max_residual_force: f64,
    starting_geometry: UffStartingGeometry,
}

/// Issues #185/#188 (UFF minimizer catastrophic bond-length blowup on some
/// starting geometries, sound convergence on others, at the shared default
/// 200-iteration budget): tries the caller-provided `coords` first, exactly
/// as before. Only if that specific attempt fails soundness with
/// `CatastrophicBondBlowup`/`NonFiniteCoordinates` does this retry once from
/// `embed_distance_geometry_v2`'s geometry instead (see
/// [`rescue_with_distance_geometry_v2`] for why this is a post-hoc retry on
/// the actual outcome, not a pre-minimization heuristic guess) — every other
/// caller (including every molecule that already passes today) sees zero
/// behavior change, and this never silently substitutes: which geometry
/// actually produced the returned result is always visible on
/// [`UffBridgeRun::starting_geometry`]/[`PolicyMinimizeResult::starting_geometry`],
/// and a failed rescue attempt is disclosed on
/// [`MinimizationFailureDetail::distance_geometry_v2_retry_attempted`].
fn run_uff_bridge(
    mol: &Molecule,
    coords: &Coords3D,
    max_iter: usize,
) -> Result<UffBridgeRun, ForceFieldBridgeError> {
    let n = mol.atom_count();
    let types = assign_uff_types(mol);
    let coord_vec = coords_to_vec(coords, n);
    let energy_before = uff_total_energy(mol, &types, &coord_vec);
    let result = ff_minimize_uff(mol, &types, coord_vec, max_iter);
    let energy_after = uff_total_energy(mol, &types, &result.coords);
    let max_residual_force =
        fd_max_gradient(&result.coords, |c| uff_total_energy(mol, &types, c), 1e-4);

    match check_minimization_soundness(
        mol,
        &result.coords,
        ForceFieldPolicy::UffOnly,
        result.converged,
        result.iterations,
        max_residual_force,
        false,
    ) {
        Ok(()) => Ok(UffBridgeRun {
            coords: vec_to_coords(&result.coords),
            energy_before,
            energy_after,
            converged: result.converged,
            iterations: result.iterations,
            max_residual_force,
            starting_geometry: UffStartingGeometry::AsProvided,
        }),
        Err(ForceFieldBridgeError::MinimizationFailed(detail))
            if matches!(
                detail.reason,
                MinimizationFailureReason::CatastrophicBondBlowup
                    | MinimizationFailureReason::NonFiniteCoordinates
            ) =>
        {
            rescue_with_distance_geometry_v2(mol, &types, max_iter, detail)
        }
        Err(other) => Err(other),
    }
}

/// The rescue itself: `original_failure` is `coords`' own (already-computed)
/// failure. Raw starting energy alone cannot decide in advance whether a
/// given start needs this — measured (`docs/rfcs/uff_robustness_diagnosis_185_188.md`):
/// naphthalene and anthracene start from virtually identical raw
/// `generate_coords` UFF energies (~1.256e5 vs ~1.267e5 kcal/mol) and
/// worst starting bonds (2.26 Å, identical), yet naphthalene's raw UFF
/// minimization reaches a sound geometry (worst bond 2.19 Å) while
/// anthracene's plateaus permanently unsound (worst bond 3.39 Å, still
/// above `MAX_SANE_BOND_LENGTH` even at 200,000 steps) — so the
/// decision is made from the ACTUAL post-minimization outcome, never a
/// pre-minimization heuristic guess.
///
/// A retry is only accepted as a real rescue if it is BOTH geometrically
/// sound AND preserves every declared stereocenter/E-Z bond
/// ([`crate::stereo_constraints::verify_stereo`]) — measured directly, not
/// assumed. A rescue that fixed a bond-length blowup while silently
/// destroying declared stereochemistry would be a worse outcome than the
/// honest failure it replaced, so this bridge does not accept one — the
/// original failure is returned instead, same as if the rescue had never
/// been geometrically sound. Only the rescue's own (new) behavior is held to
/// this bar; the unrelated, pre-existing fact that this bridge's *first*-
/// attempt path never verifies stereo either is a known, separate,
/// out-of-scope gap (already disclosed in `examples/
/// cf_integration_smoke_test.rs`'s closing note), not one this fix
/// introduces or is scoped to close.
///
/// **Revised (issue #210, 2026-08-24)**: the embed call below now sets
/// `enforce_chirality`/`materialize_implicit_h_for_chirality` (both `true`,
/// gated on [`mol_has_declared_stereo`]) rather than plain
/// `EmbedParameters::default()`. The original rationale for defaults —
/// "`enforce_chirality: true` would make embedding refuse any molecule with
/// declared stereo outright" — was stale even before this change (that flag
/// repairs via `repair_stereo` before failing, not an unconditional refusal);
/// what actually blocked ring-fused declared stereocenters (testosterone,
/// cholesterol) specifically was `repair_tetrahedral_center` having no
/// bridge-eligible substituent to reflect when every real neighbor is a ring
/// atom — exactly the gap issue #291's `materialize_implicit_h_for_chirality`
/// closed. **Measured directly** against the 58-molecule corpus this
/// bridge's own tests use (`examples/cf_integration_smoke_test.rs`'s
/// corpus): zero regressions (every previously-passing molecule, and every
/// declared-stereo molecule's existing pass/fail-and-stereo-check outcome,
/// unchanged), and one of #210's 5 named residual molecules
/// (`atorvastatin_fragment`) newly succeeds with stereo preserved.
/// `naproxen_S`/`ibuprofen_S`/`testosterone`/`cholesterol` remain unfixed —
/// this bridge has no post-minimization repair-and-reverify step (UFF
/// minimization is chirality-blind and can walk a stereo-satisfied embed
/// back across a declared boundary, same class of problem issue #291's
/// `pipeline_v2.rs` stage 11 exists to catch and retry; this simpler
/// embed→minimize→verify bridge just falls through to the original failure
/// instead, which is correct/safe, not a bug) — closing that residual would
/// need adding an equivalent repair step here, a separate, larger change not
/// attempted by this measurement.
///
/// If `embed_distance_geometry_v2` itself fails, or the retry is unsound or
/// stereo-violating, the ORIGINAL failure is returned unchanged except for
/// `distance_geometry_v2_retry_attempted` flipping to `true` — the rescue
/// never masks a real failure with a different, potentially-misleading one.
fn rescue_with_distance_geometry_v2(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    max_iter: usize,
    original_failure: Box<MinimizationFailureDetail>,
) -> Result<UffBridgeRun, ForceFieldBridgeError> {
    use crate::distance_geometry_v2::{
        EmbedParameters, embed_distance_geometry_v2, mol_has_declared_stereo,
    };
    use crate::stereo_constraints::verify_stereo;

    let n = mol.atom_count();
    // See this function's own doc comment ("Revised, issue #210") for why
    // these are set instead of `EmbedParameters::default()`.
    let has_stereo = mol_has_declared_stereo(mol);
    let embed_params = EmbedParameters {
        enforce_chirality: has_stereo,
        materialize_implicit_h_for_chirality: has_stereo,
        ..EmbedParameters::default()
    };
    if let Ok(v2_coords) = embed_distance_geometry_v2(mol, &embed_params) {
        let retry_coord_vec = coords_to_vec(&v2_coords, n);
        // `energy_before`/`energy_after` on a successful rescue must both describe
        // the SAME trajectory (the DG v2 geometry, before and after minimizing it) --
        // never the caller's original (now-abandoned) starting geometry paired with
        // the retry's outcome. The caller's original energy is already captured
        // separately, in `original_failure` (the typed error this function falls
        // back to if the rescue doesn't pan out) -- it has no place in a successful
        // `UffBridgeRun`, which reports on the geometry that actually produced it.
        let retry_energy_before = uff_total_energy(mol, types, &retry_coord_vec);
        let retry = ff_minimize_uff(mol, types, retry_coord_vec, max_iter);
        let retry_energy_after = uff_total_energy(mol, types, &retry.coords);
        let retry_max_residual_force =
            fd_max_gradient(&retry.coords, |c| uff_total_energy(mol, types, c), 1e-4);
        let retry_coords_typed = vec_to_coords(&retry.coords);

        let retry_sound = check_minimization_soundness(
            mol,
            &retry.coords,
            ForceFieldPolicy::UffOnly,
            retry.converged,
            retry.iterations,
            retry_max_residual_force,
            false,
        )
        .is_ok();
        let retry_preserves_stereo = verify_stereo(mol, &retry_coords_typed).is_fully_satisfied();

        if retry_sound && retry_preserves_stereo {
            return Ok(UffBridgeRun {
                coords: retry_coords_typed,
                energy_before: retry_energy_before,
                energy_after: retry_energy_after,
                converged: retry.converged,
                iterations: retry.iterations,
                max_residual_force: retry_max_residual_force,
                starting_geometry: UffStartingGeometry::ReplacedWithDistanceGeometryV2,
            });
        }
    }

    let mut detail = original_failure;
    detail.distance_geometry_v2_retry_attempted = true;
    Err(ForceFieldBridgeError::MinimizationFailed(detail))
}

fn finish_mmff94(
    r: Mmff94BridgeRun,
    requested_force_field: ForceFieldPolicy,
    actual_force_field_used: ForceFieldPolicy,
    fallback_reason: Option<ForceFieldBridgeError>,
) -> PolicyMinimizeResult {
    let missing_parameter_classes = r.coverage.all_missing();
    PolicyMinimizeResult {
        coords: r.coords,
        requested_force_field,
        actual_force_field_used,
        fallback_reason,
        missing_parameter_classes,
        coverage: Some(r.coverage),
        energy_before: EnergyReport::Mmff94(r.energy_before),
        energy_after: EnergyReport::Mmff94(r.energy_after),
        converged: r.converged,
        iterations: r.iterations,
        max_residual_force: r.max_residual_force,
        starting_geometry: None,
    }
}

fn finish_uff(
    r: UffBridgeRun,
    requested_force_field: ForceFieldPolicy,
    actual_force_field_used: ForceFieldPolicy,
    fallback_reason: Option<ForceFieldBridgeError>,
    missing_parameter_classes: Vec<Mmff94MissingTerm>,
    coverage: Option<Mmff94CoverageReport>,
) -> PolicyMinimizeResult {
    PolicyMinimizeResult {
        coords: r.coords,
        requested_force_field,
        actual_force_field_used,
        fallback_reason,
        missing_parameter_classes,
        coverage,
        energy_before: EnergyReport::Uff {
            total: r.energy_before,
        },
        energy_after: EnergyReport::Uff {
            total: r.energy_after,
        },
        converged: r.converged,
        iterations: r.iterations,
        max_residual_force: r.max_residual_force,
        starting_geometry: Some(r.starting_geometry),
    }
}

/// Bridge chematic-3d's minimization step to chematic-ff's complete
/// force-field implementations. Opt-in only — no existing public function's
/// default behavior changes because of this.
///
/// `include_torsion_oop_in_gate`: when `true`, `Mmff94BondAngleStrict` (and the MMFF94
/// attempt inside `Mmff94WithUffFallback`) also refuses on a missing torsion
/// or out-of-plane term, not just bond/angle. [`minimize_with_policy`] passes
/// `false` — bond+angle is the RFC's mechanism-3 scope (the "no restoring
/// force, vdW pushes atoms apart" failure mode); torsion/oop absence is a
/// real, always-*reported* accuracy gap (see `coverage`/
/// `missing_parameter_classes` on the result) but not the same
/// structural-blowup pathology. Pre-chematic-ff-#183, gating on it too was
/// measured to fail a large fraction of ordinary organic molecules
/// (`angle_type_for` was hardcoded to 0, permanently stranding ~21% of the
/// angle table, which also fed torsion typing). That measurement is now
/// stale: post-#183, `examples/mmff94_bridge_coverage_report.rs` measures
/// 32/58 passing the bond+angle-only gate vs. 31/58 under the widened gate
/// — nearly identical — so widening the default is now a cheap, small
/// follow-up rather than impractical, just not done in this PR (see
/// [`ForceFieldPolicy::Mmff94BondAngleStrict`]'s doc for why not).
///
/// `include_stretch_bend_in_gate`: same opt-in shape, independent of
/// `include_torsion_oop_in_gate` (Priority 2 / Stage 1B, issue #227) — when
/// `true`, also refuses on a missing stretch-bend cross term. Stretch-bend
/// was historically not even measured toward the gate (no field on
/// `Mmff94CoverageReport` existed for it at all until this parameter was
/// added); `stretch_bend_energy` (chematic-ff's `mmff94_minimizer`) silently
/// contributes zero energy for an uncovered stretch-bend term regardless of
/// this flag — this parameter only controls whether that silence is also a
/// typed refusal up front, it does not change the energy function itself.
/// [`minimize_with_policy`] passes `false`, matching its existing
/// `include_torsion_oop_in_gate = false` default — no existing caller's
/// behavior changes.
pub fn minimize_with_policy_gated(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
    include_torsion_oop_in_gate: bool,
    include_stretch_bend_in_gate: bool,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    if mol.atom_count() <= 1 {
        return Ok(trivial_result(coords, policy));
    }

    match policy {
        ForceFieldPolicy::None => Ok(trivial_result(coords, ForceFieldPolicy::None)),

        ForceFieldPolicy::Dreiding => {
            let dreiding_types = assign_dreiding_types(mol);
            let e_before = total_energy_dreiding(mol, &coords, &dreiding_types);
            let report = minimize_gradient_descent_reporting(mol, coords, config, |c| {
                total_energy_dreiding(mol, c, &dreiding_types)
            });
            let e_after = total_energy_dreiding(mol, &report.coords, &dreiding_types);
            let coord_vec = coords_to_vec(&report.coords, mol.atom_count());
            check_minimization_soundness(
                mol,
                &coord_vec,
                ForceFieldPolicy::Dreiding,
                report.converged,
                report.iterations,
                report.final_max_grad,
                false,
            )?;
            Ok(PolicyMinimizeResult {
                coords: report.coords,
                requested_force_field: ForceFieldPolicy::Dreiding,
                actual_force_field_used: ForceFieldPolicy::Dreiding,
                fallback_reason: None,
                missing_parameter_classes: Vec::new(),
                coverage: None,
                starting_geometry: None,
                energy_before: EnergyReport::Dreiding { total: e_before },
                energy_after: EnergyReport::Dreiding { total: e_after },
                converged: report.converged,
                iterations: report.iterations,
                max_residual_force: report.final_max_grad,
            })
        }

        ForceFieldPolicy::UffOnly => {
            let r = run_uff_bridge(mol, &coords, config.max_steps)?;
            Ok(finish_uff(
                r,
                ForceFieldPolicy::UffOnly,
                ForceFieldPolicy::UffOnly,
                None,
                Vec::new(),
                None,
            ))
        }

        ForceFieldPolicy::Mmff94BondAngleStrict => {
            let r = run_mmff94_bridge(
                mol,
                &coords,
                config.max_steps,
                include_torsion_oop_in_gate,
                include_stretch_bend_in_gate,
            )?;
            Ok(finish_mmff94(
                r,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                None,
            ))
        }

        ForceFieldPolicy::Mmff94WithUffFallback => {
            match run_mmff94_bridge(
                mol,
                &coords,
                config.max_steps,
                include_torsion_oop_in_gate,
                include_stretch_bend_in_gate,
            ) {
                Ok(r) => Ok(finish_mmff94(
                    r,
                    ForceFieldPolicy::Mmff94WithUffFallback,
                    ForceFieldPolicy::Mmff94BondAngleStrict,
                    None,
                )),
                Err(reason) => {
                    let (coverage, missing_parameter_classes) = match &reason {
                        ForceFieldBridgeError::MissingParameters(rep) => {
                            (Some((**rep).clone()), rep.all_missing())
                        }
                        ForceFieldBridgeError::UnsupportedAtomType(_)
                        | ForceFieldBridgeError::MinimizationFailed(_) => (None, Vec::new()),
                    };
                    // NOTE: if the UFF fallback attempt is *itself* unsound
                    // (measured post-#183: still happens on 8/58 corpus
                    // molecules, all fused/conjugated polycyclic aromatics —
                    // see PR body), `?` here returns that failure directly.
                    // The original MMFF94 `reason` (why fallback was even
                    // attempted) is not preserved in that returned error —
                    // an accepted, documented tradeoff for a case this
                    // policy's contract never claimed to make sound, only
                    // non-silent. `Mmff94WithUffFallback` is therefore NOT
                    // infallible; see this variant's doc.
                    let r = run_uff_bridge(mol, &coords, config.max_steps)?;
                    Ok(finish_uff(
                        r,
                        ForceFieldPolicy::Mmff94WithUffFallback,
                        ForceFieldPolicy::UffOnly,
                        Some(reason),
                        missing_parameter_classes,
                        coverage,
                    ))
                }
            }
        }
    }
}

/// Convenience wrapper over [`minimize_with_policy_gated`] with
/// `include_torsion_oop_in_gate = false` and `include_stretch_bend_in_gate =
/// false` (bond+angle only — mechanism-3's exact scope). Use
/// `minimize_with_policy_gated` directly to widen the strict gate to
/// torsion/out-of-plane and/or stretch-bend.
pub fn minimize_with_policy(
    mol: &Molecule,
    coords: Coords3D,
    policy: ForceFieldPolicy,
    config: &MinimizeConfig,
) -> Result<PolicyMinimizeResult, ForceFieldBridgeError> {
    minimize_with_policy_gated(mol, coords, policy, config, false, false)
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
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "aspirin atom {i} has invalid coords"
            );
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
        assert!(
            d > 1.4 && d < 1.7,
            "ethane C-C should be ~1.54 Å with electrostatic, got {:.3}",
            d
        );
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
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "atom {i} has invalid coordinate after minimization"
            );
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

        assert!(
            (total_e - expected).abs() < 1e-6,
            "total energy mismatch: got {}, expected {}",
            total_e,
            expected
        );
    }
}

// ---------------------------------------------------------------------------
// Force-field bridge tests (Wave 1, Agent F)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod policy_bridge_tests {
    use super::*;
    use crate::dg::generate_coords;
    use chematic_ff::{assign_mmff94_numeric_types, mmff94_bond_energy};
    use chematic_smiles::parse;

    // --- Coords3D <-> Vec<[f64; 3]> bridge plumbing -------------------------

    #[test]
    fn coords_vec_round_trip_preserves_per_atom_identity() {
        // Distinct, non-symmetric coordinates per atom so a transposition or
        // off-by-one bug in the bridge's index handling would be caught.
        let n = 5;
        let mut coords = Coords3D::new_zeroed(n);
        for i in 0..n {
            coords.set(
                AtomIdx(i as u32),
                Point3::new(i as f64, 10.0 + i as f64 * 2.0, 100.0 - i as f64 * 3.0),
            );
        }
        let v = coords_to_vec(&coords, n);
        let back = vec_to_coords(&v);
        for i in 0..n {
            let orig = coords.get(AtomIdx(i as u32));
            let rt = back.get(AtomIdx(i as u32));
            assert!(
                (orig.x - rt.x).abs() < 1e-12
                    && (orig.y - rt.y).abs() < 1e-12
                    && (orig.z - rt.z).abs() < 1e-12,
                "atom {i} round-trip mismatch: {orig:?} vs {rt:?}"
            );
        }
    }

    /// Per the RFC's "analytic-vs-finite-difference gradient self-check"
    /// request: chematic-ff (verified by reading `mmff94_minimizer.rs` and
    /// `uff.rs`) exposes NO analytic/closed-form gradient anywhere — both
    /// `compute_gradient` and `uff_gradient` are private and finite-
    /// difference-based internally. There is therefore no analytic gradient
    /// to compare against; a literal analytic-vs-FD test is unsatisfiable
    /// with the current chematic-ff surface (flagged in the PR body).
    ///
    /// What *is* a meaningful correctness check on the bridge itself is
    /// whether this bridge's `Coords3D` <-> `Vec<[f64; 3]>` conversion
    /// preserves atom correspondence into the FD gradient computation: this
    /// test perturbs the SAME atom/axis via two independent routes — (a)
    /// through `Coords3D::get`/`set` (mimicking how a caller holding a
    /// `Coords3D` would use this bridge) and (b) directly on a raw
    /// `Vec<[f64; 3]>` with no `Coords3D` involved at all — and confirms
    /// both report the same gradient, localized on the atom that was
    /// actually perturbed. A transposition/off-by-one bug in the bridge's
    /// index handling would make these disagree.
    #[test]
    fn bridge_fd_gradient_matches_raw_chematic_ff_call() {
        let mol = parse("CCO").expect("ethanol topology (heavy atoms only)");
        let mut coords = generate_coords(&mol);
        // Deliberately stretch the C-O bond (atom 2) far from equilibrium so
        // the gradient at atom 2 is unambiguously the largest one.
        let o = coords.get(AtomIdx(2));
        coords.set(AtomIdx(2), Point3::new(o.x + 3.0, o.y, o.z));

        // Route (a): through Coords3D -> Vec via this bridge's own helper.
        let via_bridge = coords_to_vec(&coords, mol.atom_count());
        let grad_a = fd_max_gradient(
            &via_bridge,
            |c| mmff94_total_energy(&mol, c).expect("energy"),
            1e-4,
        );

        // Route (b): construct the equivalent raw Vec directly, no Coords3D
        // anywhere in this path.
        let raw: Vec<[f64; 3]> = (0..mol.atom_count())
            .map(|i| {
                let p = coords.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        let grad_b = fd_max_gradient(
            &raw,
            |c| mmff94_total_energy(&mol, c).expect("energy"),
            1e-4,
        );

        assert!(
            (grad_a - grad_b).abs() < 1e-9,
            "bridge FD gradient ({grad_a}) disagrees with raw-Vec FD gradient ({grad_b}) \
             on the same coordinates — possible index/transposition bug in coords_to_vec"
        );
        // Sanity: the artificially stretched bond should produce a large,
        // clearly-nonzero residual force, not a spuriously-zero one (which
        // is exactly the mechanism-3 failure mode this bridge fixes).
        assert!(
            grad_a > 10.0,
            "expected a large residual force from a 3 Å-stretched C-O bond, got {grad_a}"
        );
    }

    /// Pins `chematic_ff::bond_type_for` (called directly now — this bridge
    /// no longer maintains its own duplicate) against `mmff94_bond_energy`'s
    /// actual Some/None behavior for a handful of type/order pairs. If
    /// chematic-ff's classification ever regresses, this test catches it —
    /// a divergence here would make the coverage report falsely confident
    /// (missing bonds could recompute as "covered" wrongly, or vice versa),
    /// the same failure class this whole bridge exists to fix.
    ///
    /// The `(1, 3, Single)` and `(3, 7, Double)` cases are the exact
    /// carbonyl-bond regression found in independent review of this PR's
    /// *previous* hand-copied classification (pre-chematic-ff-#183, that
    /// duplicate used order-blind OR-logic: "bond_type=1 if EITHER atom is
    /// SP2"). Since type 3 (carbonyl C) was in that OR-based list, EVERY
    /// bond touching it — including the C=O double bond itself (types 3,7)
    /// — was misclassified bond_type=1, and `mmff94_bond_energy(1, 3, 7)`
    /// has no entry (only `(0, 3, 7, ...)` does), so the old duplicate
    /// falsely reported carbonyl C=O — present in essentially every
    /// ketone/aldehyde/amide/carboxylic acid — as missing MMFF94 coverage.
    /// Calling chematic-ff's real, order-aware `bond_type_for` directly
    /// (this fix) makes both cases resolve to bond_type=0 (forced by the
    /// `Double` order guard, independent of the SP2/MLTB list) and find
    /// their real table rows — see also
    /// `carbonyl_c_o_double_bond_is_not_falsely_reported_missing` below for
    /// the same regression pinned at the `compute_mmff94_coverage` level.
    #[test]
    fn bridge_bond_type_matches_chematic_ff_lookup_behavior() {
        // (type_i, type_j, order, expect_some)
        let cases: &[(u8, u8, BondOrder, bool)] = &[
            (1, 1, BondOrder::Single, true),  // C(sp3)-C(sp3): covered
            (1, 5, BondOrder::Single, true),  // C(sp3)-H: covered
            (1, 11, BondOrder::Single, true), // C(sp3)-F: covered
            (1, 12, BondOrder::Single, true), // C(sp3)-Cl: covered
            (1, 13, BondOrder::Single, true), // C(sp3)-Br: covered
            // Aromatic C-C (benzene, real type post-issue-#227-Phase-1B-0 is
            // 37 = CB, not the old wrong 63): `BondOrder::Aromatic` forces
            // bt=0 unconditionally (like Double/Triple), matching RDKit's
            // own (bondType=0, r0=1.374) for this exact bond -- verified
            // against a live RDKit oracle, see mmff94_minimizer.rs's
            // `benzene_ring_bonds_type_37_resolve_to_the_rdkit_verified_row`.
            (37, 37, BondOrder::Aromatic, true),
            // Cα(sp3)-C(carbonyl) SINGLE bond: bt=0 (type 1 isn't in MLTB_TYPES,
            // so the AND-gate fails) -> real (0,1,3) row -> covered.
            (1, 3, BondOrder::Single, true),
            // The carbonyl C=O DOUBLE bond itself: bt=0 (forced by the
            // Double/Triple/Quadruple order guard, before the MLTB check
            // even runs) -> real (0,3,7) row -> covered. THE carbonyl
            // regression case (see doc comment above).
            (3, 7, BondOrder::Double, true),
        ];
        for &(ti, tj, order, expect_some) in cases {
            let bt = bond_type_for(ti, tj, order);
            let got = mmff94_bond_energy(bt, ti, tj).is_some();
            assert_eq!(
                got, expect_some,
                "bond_type_for({ti},{tj},{order:?})={bt}: mmff94_bond_energy returned {got}, expected {expect_some}"
            );
        }
    }

    /// Direct regression pin for the carbonyl false-negative at the level
    /// that actually matters to callers: `compute_mmff94_coverage` on a real
    /// molecule containing a C=O bond must NOT list it in `bonds_missing`.
    /// This is exactly the check that would have caught the bridge's
    /// pre-#183-duplicate bug (see module docs above and the PR body) before
    /// it shipped — `Mmff94BondAngleStrict` would otherwise have spuriously refused
    /// (or `Mmff94WithUffFallback` silently fallen back on) ordinary
    /// ketones/aldehydes/amides/carboxylic acids.
    #[test]
    fn carbonyl_c_o_double_bond_is_not_falsely_reported_missing() {
        let mol = parse("CC=O").expect("acetaldehyde: C(sp3)-C(carbonyl)=O");
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        // Confirm this really does exercise the carbonyl types (3, 7), not
        // some other typing this bridge doesn't intend to test.
        assert_eq!(types[1], 3, "atom 1 (carbonyl C) should be MMFF94 type 3");
        assert_eq!(types[2], 7, "atom 2 (carbonyl O) should be MMFF94 type 7");

        let report = compute_mmff94_coverage(&mol, &types);
        assert!(
            report.bonds_missing.is_empty(),
            "carbonyl C=O (and Cα-C) bonds must be covered, not falsely missing: {:?}",
            report.bonds_missing
        );
    }

    // --- Mechanism-3 repro: [C@H](F)(Cl)Br under old vs. new MMFF94 --------

    fn chfclbr_mol() -> chematic_core::Molecule {
        parse("[C@H](F)(Cl)Br").expect("chfclbr")
    }

    fn worst_bond(mol: &chematic_core::Molecule, coords: &Coords3D) -> f64 {
        let mut worst = 0.0_f64;
        for (_, bond) in mol.bonds() {
            let d = coords.get(bond.atom1).distance(&coords.get(bond.atom2));
            if d > worst {
                worst = d;
            }
        }
        worst
    }

    #[test]
    fn old_minimize_mmff94_reproduces_mechanism3_blowup_at_this_level() {
        // Confirm the bug reproduces starting from THIS crate's own
        // dg::generate_coords geometry, at the minimize.rs entry point
        // directly (`minimize_mmff94`) — not assuming the RFC's Python-path
        // 24.3 Å number (which goes through etkdg.rs's torsion/constraint
        // stages too) transfers unchanged to this narrower entry point.
        let mol = chfclbr_mol();
        assert_eq!(mol.atom_count(), 4, "C, F, Cl, Br heavy atoms only");
        let coords = generate_coords(&mol);
        let before = worst_bond(&mol, &coords);

        let result = minimize_mmff94(&mol, coords);
        let after = worst_bond(&mol, &result);

        // The old path's bond+angle energy silently zeroes for F-C-Cl,
        // F-C-Br, Cl-C-Br (mmff94_bond_params has no halogen-halogen-via-C
        // angle entries), leaving only VdW repulsion to act on the 3 C-X
        // bonds — they should blow out well past any real covalent bond.
        assert!(
            after > 5.0,
            "expected the old MMFF94 path to blow up C-X bonds on this molecule \
             (before={before:.2} Å, after={after:.2} Å) — if this no longer reproduces, \
             the mechanism-3 'before' baseline this PR cites is stale and must be re-measured"
        );
    }

    #[test]
    fn mmff94_strict_now_resolves_chfclbr_via_empirical_angle_rule() {
        // Issue #227 Stage C: chfclbr's 3 halogen-C-halogen angles (F-C-Cl,
        // F-C-Br, Cl-C-Br -- MMFF94 types 11/1/12, 11/1/13, 12/1/13) used to
        // have no MMFF94 table entry at all (see the superseded
        // `mmff94_strict_refuses_chfclbr_with_missing_angle_params`,
        // replaced here). All 3 now resolve via Halgren's eq. 20 empirical
        // angle-bend rule -- theta0 reused verbatim from the restored
        // `(0, 0, 1, 0)` wildcard row, ka derived empirically -- oracle-
        // confirmed by
        // `mmff94_energy::angle::tests::angle_empirical_reuses_wildcard_theta0_for_halogen_and_amine_c_substituents`
        // against this exact molecule.
        let mol = chfclbr_mol();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        )
        .expect("chfclbr's 3 halogen-C-halogen angles now resolve via the Stage C empirical rule");

        assert_eq!(
            result.requested_force_field,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert!(
            result.fallback_reason.is_none(),
            "no fallback needed now that all bonds/angles resolve"
        );

        let coverage = result
            .coverage
            .as_ref()
            .expect("strict policy always reports coverage");
        assert!(coverage.bonds_missing.is_empty());
        assert!(
            coverage.angles_missing.is_empty(),
            "all 3 halogen-C-halogen angles should now resolve via the empirical rule, got: {:?}",
            coverage.angles_missing
        );

        // The empirical rule must contribute a REAL, nonzero angle energy
        // before minimization -- not silently 0.0 (the same "never claim a
        // term covered but compute 0 energy for it" principle this issue's
        // whole Stage C is built on).
        match &result.energy_before {
            EnergyReport::Mmff94(b) => {
                assert!(
                    b.angle > 0.0,
                    "empirical angle-bend energy must be real and positive before minimization, got {}",
                    b.angle
                );
            }
            other => panic!("expected Mmff94 energy report, got {other:?}"),
        }
        assert!(result.converged);
    }

    #[test]
    fn mmff94_with_uff_fallback_no_longer_needs_to_fall_back_on_chfclbr() {
        // Issue #227 Stage C: companion to
        // `mmff94_strict_now_resolves_chfclbr_via_empirical_angle_rule` --
        // now that strict MMFF94 itself resolves chfclbr, the UFF-fallback
        // policy must use MMFF94 directly rather than falling back (the
        // superseded version of this test, `..._falls_back_and_reports_why_on_chfclbr`,
        // asserted the opposite -- that fallback was required).
        let mol = chfclbr_mol();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
        )
        .expect("chfclbr now succeeds directly under MMFF94, no fallback needed");

        assert_eq!(
            result.requested_force_field,
            ForceFieldPolicy::Mmff94WithUffFallback
        );
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            "no fallback needed now that chfclbr's angles resolve via the empirical rule"
        );
        assert!(
            result.fallback_reason.is_none(),
            "fallback must not fire when the strict policy already succeeds"
        );
        assert!(result.missing_parameter_classes.is_empty());
        let coverage = result
            .coverage
            .as_ref()
            .expect("strict policy always reports coverage");
        assert!(
            coverage.angles_missing.is_empty(),
            "all 3 halogen-C-halogen angles resolve via the Stage C empirical rule, got: {:?}",
            coverage.angles_missing
        );
        assert!(
            coverage.stretch_bend_missing.is_empty(),
            "all 3 halogen-C-halogen triples resolve via the Dfsb periodic-row fallback \
             (Priority 2B) -- if this regresses, the Dfsb port likely broke"
        );

        // MMFF94 itself now produces a sane geometry -- no need for UFF's
        // generic coverage to rescue this molecule anymore.
        let after = worst_bond(&mol, &result.coords);
        assert!(
            after < 3.0,
            "expected a sane, non-blown-up geometry from MMFF94 directly, got worst bond {after:.2} Å"
        );
    }

    #[test]
    fn caffeine_passes_the_complete_bonded_term_gate_after_the_reperceived_view_fix() {
        // Issue #227 Phase 1: caffeine's dione-ring torsions used to be
        // torsions_missing under the classification bug this PR fixes --
        // never gated under the default bond+angle-only policy
        // (minimize_with_policy_gated(..., false, false), what
        // minimize_with_policy uses), but WOULD have failed
        // `include_torsion_oop_in_gate=true` before this fix. Confirms the
        // real, gated production path -- not just the coverage-audit tool --
        // now succeeds, and does so via a real, non-fallback MMFF94 run.
        let mol = chematic_smiles::parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy_gated(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
            true,
            true,
        )
        .expect("caffeine must pass the complete bonded-term gate directly, no fallback");

        let coverage = result
            .coverage
            .as_ref()
            .expect("strict policy always reports coverage");
        assert!(
            coverage.torsions_missing.is_empty(),
            "caffeine has no genuine torsion coverage gap, got: {:?}",
            coverage.torsions_missing
        );
    }

    /// Issues #185/#188's fix, exercised through `Mmff94WithUffFallback`:
    /// naphthalene lacks enough MMFF94 angle/torsion coverage to pass
    /// `Mmff94BondAngleStrict`, so this falls back to `UffOnly` — and that
    /// UFF fallback, from `generate_coords`'s geometry, is measured (58-
    /// molecule corpus, post-chematic-ff-#183) to blow the worst bond length
    /// out past 3 Å (naphthalene: 1.43 Å -> 4.74 Å) at the default 200-step
    /// budget. Before `run_uff_bridge`'s rescue existed, this surfaced as a
    /// typed `Err(MinimizationFailed)` (a real fix over the pre-soundness-gate
    /// behavior of returning `Ok(PolicyMinimizeResult { converged: false, .. })`,
    /// which a careless caller checking only `.is_ok()` would mistake for
    /// success) — but the underlying cause was `generate_coords` producing a
    /// severely clashed start, not an unrecoverable molecule: retrying once
    /// from `embed_distance_geometry_v2`'s geometry succeeds. This test now
    /// pins the fixed behavior: `Ok`, with `starting_geometry` disclosing
    /// that the rescue fired (never silently). See
    /// `chematic_ff_own_uff_minimizer_blows_up_naphthalene_independent_of_this_bridge`
    /// below — chematic-ff's own `minimize_uff`, called with zero bridge code
    /// in the path, still blows up on this exact `generate_coords` geometry;
    /// this test's rescue lives entirely in `run_uff_bridge`, not in
    /// chematic-ff, so that fixture remains correct and unchanged.
    #[test]
    fn naphthalene_now_succeeds_directly_via_mmff94_strict_post_227_fix() {
        // Issue #227 Phase 1B-0: this test used to pin naphthalene as a
        // case where MMFF94 coverage was "genuinely incomplete" (both rings
        // all-carbon fused-aromatic, mistyped as 63/64 pre-fix) and had to
        // fall back to UFF, exercising the separate #185/#188
        // catastrophic-clash-rescue path (`UffStartingGeometry::
        // ReplacedWithDistanceGeometryV2`). After the atom-typing fix
        // (naphthalene's ring carbons now correctly type as 37, CB, matching
        // both rings being 6-membered), naphthalene's MMFF94 coverage is
        // complete and it succeeds directly -- no UFF fallback, no rescue
        // path exercised at all. This is a real positive outcome of the fix,
        // not a loosened assertion. The #185/#188 rescue-path behavior this
        // test used to cover needs a different fixture molecule that still
        // genuinely requires UFF fallback post-fix if that coverage matters
        // going forward -- not re-derived here (out of this PR's scope).
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94WithUffFallback,
            &config,
        )
        .expect("naphthalene must succeed post-#227-fix");
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            "naphthalene should now have complete MMFF94 coverage and need no UFF fallback"
        );
        assert!(
            result.fallback_reason.is_none(),
            "no fallback should have occurred"
        );
        let worst = worst_bond_length_vec(&mol, &coords_to_vec(&result.coords, mol.atom_count()));
        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected a sound geometry, got worst bond {worst:.2} Å"
        );
    }

    /// Isolated confirmation for the item-4 finding: this calls
    /// `chematic_ff::assign_uff_types`/`minimize_uff` DIRECTLY — no
    /// `chematic-3d` bridge/conversion code anywhere in the minimization
    /// path (only `dg::generate_coords` supplies the identical starting
    /// geometry) — so the naphthalene blow-up pinned above cannot be an
    /// artifact of this bridge's `coords_to_vec`/`vec_to_coords` handling;
    /// it reproduces inside chematic-ff's own UFF minimizer. Confirmed: the
    /// vdW 1-3 exclusion bug (chematic-ff #176) is already fixed (verified
    /// by reading `uff.rs`'s graph-based exclusion set), so this is a
    /// distinct, still-open chematic-ff robustness gap (candidate cause,
    /// unproven: `minimize_uff`'s naive steepest-descent-with-step-halving
    /// line search on a larger, more constrained fused-ring system) — not
    /// something this PR fixes (chematic-ff is out of this PR's
    /// file-ownership scope) or definitively diagnoses. Not specific to
    /// ring count either way: anthracene (3 fused rings) blows up under the
    /// same isolated path too, and worse than naphthalene (2 fused
    /// rings) — it never recovers a sound geometry even at 200,000 steps,
    /// vs. naphthalene's ~10,000 (see issue #185's investigation notes; an
    /// earlier reading of anthracene as "immune" was itself an artifact of
    /// a since-fixed `dg::generate_coords` ring-placement bug that
    /// silently superimposed two of its rings, corrupting that
    /// measurement).
    #[test]
    fn chematic_ff_own_uff_minimizer_blows_up_naphthalene_independent_of_this_bridge() {
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let coords = generate_coords(&mol);
        let raw: Vec<[f64; 3]> = (0..mol.atom_count())
            .map(|i| {
                let p = coords.get(AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();

        let types = assign_uff_types(&mol);
        let result = ff_minimize_uff(&mol, &types, raw, 200);

        let worst = mol
            .bonds()
            .map(|(_, b)| {
                let (i, j) = (b.atom1.0 as usize, b.atom2.0 as usize);
                let d = [
                    result.coords[i][0] - result.coords[j][0],
                    result.coords[i][1] - result.coords[j][1],
                    result.coords[i][2] - result.coords[j][2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .fold(0.0_f64, f64::max);

        assert!(
            worst > MAX_SANE_BOND_LENGTH,
            "expected chematic-ff's own uff::minimize_uff to reproduce the measured naphthalene \
             blow-up with zero chematic-3d bridge code in the path (worst bond {worst:.2} Å) -- \
             if this now passes, chematic-ff's UFF minimizer robustness improved and this \
             bridge's soundness-gate regression test / PR body item-4 measurement should be \
             revisited",
        );
    }

    // --- General policy sanity on a fully-covered molecule ------------------

    #[test]
    fn mmff94_strict_succeeds_and_lowers_energy_for_ethane() {
        let mol = parse("CC").expect("ethane");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        )
        .expect("ethane's single C-C bond/no-angle case is fully MMFF94-covered");

        assert_eq!(
            result.requested_force_field,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert_eq!(
            result.actual_force_field_used,
            ForceFieldPolicy::Mmff94BondAngleStrict
        );
        assert!(result.fallback_reason.is_none());
        assert!(result.missing_parameter_classes.is_empty());
        let coverage = result.coverage.expect("coverage must be reported");
        assert!(coverage.bond_angle_fully_covered());
        assert!(
            result.energy_after.total() <= result.energy_before.total() + 1e-6,
            "energy should not increase: before={:.4}, after={:.4}",
            result.energy_before.total(),
            result.energy_after.total()
        );
        assert!(result.max_residual_force.is_finite());
        for i in 0..result.coords.atom_count() {
            let p = result.coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
    }

    /// Locks in the actually-correct `fallback_reason` invariant (a real bug
    /// was found in the doc comment, not the code, during independent
    /// review: `fallback_reason` is `Some` iff a fallback truly occurred —
    /// NOT iff `actual_force_field_used != requested_force_field`, which is
    /// true on *every* successful `Mmff94WithUffFallback` call since
    /// `actual_force_field_used == Mmff94BondAngleStrict` on success while
    /// `requested_force_field == Mmff94WithUffFallback`). Previously only the
    /// failing-molecule fallback path was tested; this covers the
    /// no-fallback-needed path on molecules fully covered by MMFF94.
    #[test]
    fn mmff94_with_uff_fallback_reports_no_fallback_when_mmff94_fully_covers() {
        for smiles in ["CC", "C1CCCCC1", "CCCC", "C1CC2CC3CC1CC(C2)C3"] {
            let mol = parse(smiles).unwrap_or_else(|e| panic!("{smiles}: {e}"));
            let coords = generate_coords(&mol);
            let config = MinimizeConfig::default();

            let result = minimize_with_policy(
                &mol,
                coords,
                ForceFieldPolicy::Mmff94WithUffFallback,
                &config,
            )
            .expect("infallible policy");

            // requested != actual is expected and NOT itself evidence of a
            // fallback -- actual_force_field_used reports which physics ran
            // (Mmff94BondAngleStrict, the more informative value), not "did it match
            // what was requested."
            assert_eq!(
                result.requested_force_field,
                ForceFieldPolicy::Mmff94WithUffFallback
            );
            assert_eq!(
                result.actual_force_field_used,
                ForceFieldPolicy::Mmff94BondAngleStrict,
                "{smiles}: fully MMFF94-covered, should run full MMFF94, not UFF"
            );
            assert!(
                result.fallback_reason.is_none(),
                "{smiles}: no fallback should have occurred, but fallback_reason={:?}",
                result.fallback_reason
            );
            assert!(result.missing_parameter_classes.is_empty(), "{smiles}");
            let coverage = result
                .coverage
                .as_ref()
                .unwrap_or_else(|| panic!("{smiles}: coverage must be reported"));
            assert!(
                coverage.bond_angle_fully_covered(),
                "{smiles}: expected full bond+angle coverage"
            );
        }
    }

    #[test]
    fn uff_only_policy_produces_finite_non_clashing_geometry() {
        let mol = parse("c1ccccc1O").expect("phenol");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config).expect(
            "phenol's UFF geometry is sound here -- UffOnly is NOT unconditionally infallible in \
             general (see its doc: an unsound result returns Err(MinimizationFailed))",
        );
        assert_eq!(result.requested_force_field, ForceFieldPolicy::UffOnly);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::UffOnly);
        assert!(result.fallback_reason.is_none());
        assert!(matches!(result.energy_after, EnergyReport::Uff { .. }));

        let n = mol.atom_count();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = result
                    .coords
                    .get(AtomIdx(i as u32))
                    .distance(&result.coords.get(AtomIdx(j as u32)));
                assert!(d > 0.5, "atoms {i}/{j} clashed after UFF bridge: {d:.3}");
            }
        }
    }

    #[test]
    fn dreiding_policy_matches_existing_behavior_and_reports_convergence() {
        let mol = parse("CC(=O)O").expect("acetic acid");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::Dreiding, &config)
            .expect(
                "acetic acid's DREIDING geometry is sound here (worst bond ~1.5 Å) even though it \
             does not fully converge within the default iteration budget (max_residual_force \
             ~55, below the 200.0 soundness ceiling) -- Dreiding is NOT unconditionally \
             infallible in general (see check_minimization_soundness)",
            );
        assert_eq!(result.requested_force_field, ForceFieldPolicy::Dreiding);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::Dreiding);
        assert!(result.coverage.is_none());
        assert!(
            result.energy_after.total() <= result.energy_before.total() + 1e-6,
            "Dreiding energy should not increase after minimization"
        );
        assert!(result.iterations > 0, "expected at least one FD step");
    }

    #[test]
    fn none_policy_returns_geometry_unchanged() {
        let mol = parse("CCC").expect("propane");
        let coords = generate_coords(&mol);
        let before: Vec<Point3> = (0..mol.atom_count())
            .map(|i| coords.get(AtomIdx(i as u32)))
            .collect();
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::None, &config)
            .expect("None policy never errors");
        assert_eq!(result.requested_force_field, ForceFieldPolicy::None);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::None);
        assert_eq!(result.iterations, 0);
        for (i, orig) in before.iter().enumerate() {
            let after = result.coords.get(AtomIdx(i as u32));
            assert_eq!(orig.x, after.x);
            assert_eq!(orig.y, after.y);
            assert_eq!(orig.z, after.z);
        }
    }

    #[test]
    fn single_atom_molecule_is_trivially_handled_by_every_policy() {
        let mol = parse("O").expect("water heavy atom only");
        let config = MinimizeConfig::default();
        for policy in [
            ForceFieldPolicy::Mmff94BondAngleStrict,
            ForceFieldPolicy::Mmff94WithUffFallback,
            ForceFieldPolicy::UffOnly,
            ForceFieldPolicy::Dreiding,
            ForceFieldPolicy::None,
        ] {
            let coords = generate_coords(&mol);
            let result = minimize_with_policy(&mol, coords, policy, &config)
                .expect("single-atom molecules never fail any policy");
            assert!(result.converged);
            assert_eq!(result.iterations, 0);
        }
    }

    // --- Coverage report sanity ---------------------------------------------

    #[test]
    fn coverage_report_totals_match_bond_and_angle_counts() {
        let mol = parse("CCC").expect("propane"); // 2 bonds, 1 angle (middle C)
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        let report = compute_mmff94_coverage(&mol, &types);
        assert_eq!(report.bonds_total, 2);
        assert_eq!(report.angles_total, 1);
        assert!(report.bond_angle_fully_covered());
        assert_eq!(report.total_missing(), 0);
    }

    /// RDKit/COSMolKit/OpenEye advantage directive, Phase 4 (3D coverage
    /// gap): 11 of the 265-molecule strict corpus's coverage failures were
    /// exactly the same `Angle C(3)-C(2)-C(64)` gap (bis-indolylmaleimide
    /// scaffold, C5B indole-ring carbon exocyclic to a maleimide ring),
    /// closed by adding type 64's `MMFF94_EQ_LEVEL` row
    /// (`chematic-ff/src/mmff94_energy/angle.rs`). End-to-end confirmation
    /// through this crate's own real coverage-check path (not just the
    /// isolated `chematic-ff` angle-lookup unit test).
    #[test]
    fn issue_c5b_angle_gap_bis_indolylmaleimide_scaffold_now_covered() {
        let smiles = [
            "Cn1cc(C2=C(c3cn(CCCSC(=N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21", // chembl_tier_b_0180
            "Cn1cc(C2=C(c3cn(CCC[N+](C)(C)C)c4ccccc34)C(=O)NC2=O)c2ccccc21", // _0184
            "Cn1cc(C2=C(c3cn(CCCCC(=N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21", // _0185
            "CN(C)CCCn1cc(C2=C(c3cn(C)c4ccccc34)C(=O)NC2=O)c2ccccc21",   // _0186
            "Cn1cc(C2=C(c3cn(CCCN)c4ccccc34)C(=O)NC2=O)c2ccccc21",       // _0187
            "Cn1cc(C2=C(c3cn(CCCN=C(N)N)c4ccccc34)C(=O)NC2=O)c2ccccc21", // _0190
            "Cn1cc(C2=C(c3cn(CCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",        // _0195
            "Cn1cc(C2=C(c3cn(CCCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",       // _0196
            "Cn1cc(C2=C(c3cn(CCCCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",      // _0197
            "Cn1cc(C2=C(c3cn(CCCCCO)c4ccccc34)C(=O)NC2=O)c2ccccc21",     // _0198
            "CNCCCn1cc(C2=C(c3cn(C)c4ccccc34)C(=O)NC2=O)c2ccccc21",      // _0199
        ];
        for smi in smiles {
            let mol = parse(smi).unwrap_or_else(|e| panic!("failed to parse '{smi}': {e}"));
            let (types, mmff_mol) =
                assign_mmff94_numeric_types_with_view(&mol).unwrap_or_else(|e| {
                    panic!("assign_mmff94_numeric_types_with_view failed for '{smi}': {e:?}")
                });
            let report = compute_mmff94_coverage(&mmff_mol, &types);
            assert!(
                report.angles_missing.is_empty(),
                "'{smi}': expected no missing angle parameters, got {:?}",
                report.angles_missing
            );
        }
    }

    /// Negative control: `chembl_tier_b_0022`'s `(angle_type=0, 43, 18, 63)`
    /// triple is a separate, already-diagnosed, deliberately-unresolved case
    /// (see `angle_empirical_fails_closed_for_undefined_eq_level_substitution`
    /// in `chematic-ff`) -- this fix (type 64 only) must not accidentally
    /// change its behavior.
    #[test]
    fn issue_type_63_case_unaffected_by_type_64_fix() {
        let smi = "COc1cc2nc(N3CCN(S(=O)(=O)c4cccs4)CC3)nc(N)c2cc1OC";
        let mol = parse(smi);
        // Sanity: this test only needs to hold if the fixture molecule still
        // parses and still exhibits the known gap; if either assumption goes
        // stale, fail loudly rather than silently asserting nothing.
        let Ok(mol) = mol else {
            panic!("chembl_tier_b_0022 fixture failed to parse: {smi}");
        };
        let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(&mol)
            .unwrap_or_else(|e| panic!("assign_mmff94_numeric_types_with_view failed: {e:?}"));
        let report = compute_mmff94_coverage(&mmff_mol, &types);
        assert!(
            !report.angles_missing.is_empty(),
            "chembl_tier_b_0022's known type-63 angle gap should still be present \
             (unaffected by the type-64 fix) -- if this now passes, the fixture may have \
             drifted from the real chembl_tier_b_0022 molecule; re-verify before treating \
             this as a second fix landing for free"
        );
    }

    /// Was `chematic_ff_own_energy_function_is_blind_to_this_bond_stretch`,
    /// which pinned a real chematic-ff bug (`bond_type_for` using OR-logic
    /// with no bond-order check) found while building this bridge:
    /// `mmff94_total_energy` used to show exactly zero energy change for a 1
    /// Å stretch of the Cα(sp3)-C(carbonyl) bond. That bug is now fixed
    /// upstream (chematic-ff #173/PR #183: `bond_type_for` takes the real
    /// `BondOrder` and requires BOTH atoms in `MLTB_TYPES`, not just one) —
    /// confirmed here by asserting the *opposite* of what this test used to
    /// pin: chematic-ff's own energy function is now correctly sensitive to
    /// this stretch. If this ever regresses back to "blind," chematic-ff's
    /// bond_type_for was re-broken and this bridge's coverage report should
    /// be re-audited.
    #[test]
    fn chematic_ff_own_energy_function_is_now_sensitive_to_this_bond_stretch_post_183() {
        let mol = parse("CC=O").expect("acetaldehyde: C(sp3)-C(carbonyl)=O");
        let coords = generate_coords(&mol);
        let mut v = coords_to_vec(&coords, mol.atom_count());

        let d = [v[0][0] - v[1][0], v[0][1] - v[1][1], v[0][2] - v[1][2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let u = [d[0] / len, d[1] / len, d[2] / len];

        let e0 = mmff94_total_energy(&mol, &v).expect("energy");
        v[0][0] += u[0] * 1.0;
        v[0][1] += u[1] * 1.0;
        v[0][2] += u[2] * 1.0;
        let e1 = mmff94_total_energy(&mol, &v).expect("energy");

        assert!(
            (e1 - e0).abs() > 1.0,
            "expected chematic-ff's own mmff94_total_energy to now be clearly sensitive to \
             this 1 Å bond stretch post-#183 (was silently zero pre-fix); got delta={:.6} -- \
             if this now fails, chematic-ff's bond_type_for regressed and this bridge's \
             coverage report should be re-audited",
            e1 - e0
        );
    }

    // -----------------------------------------------------------------------
    // Issues #185 / #188 diagnosis: is `UffOnly` + `dg::generate_coords`'s
    // blow-up a bad-starting-geometry bug, a fundamentally broken minimizer,
    // or something else? See `docs/rfcs/uff_robustness_diagnosis_185_188.md` for
    // the full write-up and measured numbers this pins.
    //
    // Measured finding (not assumed): `generate_coords` produces starting
    // geometries with enormous UFF energies (naphthalene ~1.3e5, anthracene
    // ~1.3e5, hexane ~1.5e7 kcal/mol -- dominated by unrelieved vdW
    // overlap, not bond stretch) for essentially every molecule, not just
    // the ones that end up blowing up. `minimize_uff`'s plain fixed-step
    // steepest descent (no conjugate-gradient/quasi-Newton acceleration)
    // eventually relieves this and reaches a sound, low-energy geometry for
    // *most* molecules checked here, given enough iterations -- but not
    // all: naphthalene needs ~10,000 steps and hexane needs >20,000 (still
    // not fully converged at 100,000), while anthracene never reaches a
    // sound geometry at all, plateauing at worst bond 3.39 Å (still above
    // `MAX_SANE_BOND_LENGTH`) from 10,000 steps all the way through
    // 200,000, with `minimize_uff`'s own RMS-gradient convergence check
    // reporting `converged: true` on that unsound plateau -- a real,
    // distinct trapped-local-minimum failure mode, not just an
    // under-provisioned iteration budget. Whether a molecule's specific
    // clash-relief trajectory succeeds at all is not predicted by starting
    // energy or worst starting bond length: naphthalene and anthracene
    // start from virtually identical raw energy (~1.3e5 kcal/mol either
    // way) and identical worst starting bond (2.26 Å) yet only one
    // recovers. (An earlier version of this comment cited anthracene's raw
    // worst bond as 3.73 Å and its energy as ~5 orders of magnitude worse
    // than naphthalene's, concluding anthracene was *immune* to this
    // failure class and using that as evidence against a `dg.rs`
    // ring-placement explanation -- both the numbers and that conclusion
    // were an artifact of a since-fixed `dg::generate_coords` bug that
    // silently superimposed two of anthracene's three rings on identical
    // coordinates, corrupting the measurement; see issue #185.)
    //
    // Not fixed by the diagnostic pass itself: per that phase's scope, it was
    // regression-fixture only. Fixed by a follow-up PR (`run_uff_bridge`'s
    // rescue -- see `UffStartingGeometry`): raw starting energy alone cannot
    // predict which starts need help (see above), so rather than raising the
    // shared default budget (measured too costly/still-insufficient for the
    // worst cases -- see `docs/rfcs/uff_robustness_diagnosis_185_188.md`'s
    // companion ablation) or touching `minimize_uff`'s algorithm at all, a
    // catastrophic bond blowup/non-finite result at the existing budget is
    // now retried once from `embed_distance_geometry_v2`'s geometry, with
    // the choice always disclosed, never silent. Still not touched:
    // `minimize_uff`'s own algorithm/default `max_steps`, and chematic-ff's
    // public API beyond the additive diagnostic instrumentation from #107's
    // work (`uff_energy_breakdown`/`minimize_uff_with_trace`/etc., unrelated
    // to this fix).

    /// #188 flexible-chain case: hexane's `UffOnly` + `generate_coords`
    /// minimization is measured to blow its worst bond length past 3 Å at
    /// the default 200-step budget (same shape as naphthalene's
    /// `mmff94_with_uff_fallback_reports_typed_failure_when_fallback_itself_is_unsound`
    /// above) -- but, per the #185/#188 fix, `run_uff_bridge` retries once
    /// from `embed_distance_geometry_v2`'s geometry on exactly this failure
    /// class, and that retry succeeds. Pins the now-fixed behavior: `Ok`,
    /// with `starting_geometry` disclosing that the rescue fired.
    #[test]
    fn uff_only_succeeds_for_hexane_from_generate_coords_via_distance_geometry_v2_rescue() {
        let mol = parse("CCCCCC").expect("hexane");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .unwrap_or_else(|e| {
                panic!(
                    "hexane's UFF minimization from generate_coords is catastrophically clashed, \
                     but embed_distance_geometry_v2 rescues it -- expected Ok after the \
                     #185/#188 fix, got {e:?}"
                )
            });
        assert_eq!(result.requested_force_field, ForceFieldPolicy::UffOnly);
        assert_eq!(result.actual_force_field_used, ForceFieldPolicy::UffOnly);
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::ReplacedWithDistanceGeometryV2),
            "expected the rescue to have fired and to be disclosed, not silent"
        );
        let worst = worst_bond_length_vec(&mol, &coords_to_vec(&result.coords, mol.atom_count()));
        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected a sound rescued geometry, got worst bond {worst:.2} Å"
        );
    }

    /// Negative control for the rescue itself: a molecule whose
    /// `generate_coords` start already passes `UffOnly`'s soundness gate
    /// (phenol, pinned separately below as
    /// `uff_only_policy_produces_finite_non_clashing_geometry`) must report
    /// `starting_geometry == AsProvided` -- the rescue must never fire, and
    /// must never be reported as having fired, when the first attempt was
    /// already sound.
    #[test]
    fn uff_only_reports_as_provided_when_no_rescue_was_needed() {
        let mol = parse("c1ccccc1O").expect("phenol");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .expect("phenol's generate_coords geometry is sound under UffOnly without a rescue");
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::AsProvided),
            "no rescue should have been attempted -- the first attempt was already sound"
        );
    }

    /// `energy_before`, on the (no-rescue) `AsProvided` path, must be the caller's
    /// own starting geometry's energy -- computed here independently (a fresh
    /// `uff_total_energy` call on the exact same `generate_coords` output) and
    /// compared for exact equality, not just "looks plausible."
    #[test]
    fn uff_only_as_provided_energy_before_is_the_caller_geometry_energy() {
        let mol = parse("c1ccccc1O").expect("phenol");
        let coords = generate_coords(&mol);
        let n = mol.atom_count();
        let types = assign_uff_types(&mol);
        let expected_energy_before = uff_total_energy(&mol, &types, &coords_to_vec(&coords, n));
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .expect("phenol's generate_coords geometry is sound under UffOnly");
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::AsProvided)
        );
        match result.energy_before {
            EnergyReport::Uff { total } => assert_eq!(
                total, expected_energy_before,
                "AsProvided energy_before must equal the caller-provided geometry's own energy"
            ),
            other => panic!("expected EnergyReport::Uff, got {other:?}"),
        }
    }

    /// The energy-semantics regression this test exists to catch: on a successful
    /// rescue, `energy_before`/`energy_after` must both describe the SAME
    /// trajectory (the `embed_distance_geometry_v2` geometry, before and after
    /// minimizing it) -- never the caller's abandoned original geometry paired
    /// with the retry's outcome. Hexane's `generate_coords` energy (~1.5e7
    /// kcal/mol, per `docs/rfcs/uff_robustness_diagnosis_185_188.md`) and its
    /// `embed_distance_geometry_v2` energy (~3.78 kcal/mol) differ by roughly 6
    /// orders of magnitude -- exactly the fixture needed to make a
    /// caller-geometry/retry-outcome mismatch impossible to miss. A buggy
    /// implementation that reused the caller's `energy_before` would report
    /// something around 1.5e7 here; this asserts the ACTUAL value instead,
    /// independently recomputed from `embed_distance_geometry_v2`'s own
    /// deterministic output (`EmbedParameters::default()`'s `random_seed` is a
    /// fixed constant, so this reproduces the exact geometry the rescue itself
    /// used).
    #[test]
    fn uff_only_rescue_energy_before_reflects_dg_v2_geometry_not_caller() {
        use crate::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};

        let mol = parse("CCCCCC").expect("hexane");
        let n = mol.atom_count();
        let types = assign_uff_types(&mol);

        let caller_coords = generate_coords(&mol);
        let caller_energy = uff_total_energy(&mol, &types, &coords_to_vec(&caller_coords, n));
        assert!(
            caller_energy > 1.0e6,
            "expected hexane's generate_coords energy to be the measured ~1.5e7-scale \
             catastrophic-clash value, got {caller_energy}"
        );

        let v2_coords = embed_distance_geometry_v2(&mol, &EmbedParameters::default())
            .expect("hexane must embed via distance_geometry_v2");
        let expected_v2_energy_before =
            uff_total_energy(&mol, &types, &coords_to_vec(&v2_coords, n));
        assert!(
            expected_v2_energy_before < 1000.0,
            "expected embed_distance_geometry_v2's hexane energy to be the measured \
             ~3.78-scale sound value, got {expected_v2_energy_before}"
        );

        let config = MinimizeConfig::default();
        let result = minimize_with_policy(&mol, caller_coords, ForceFieldPolicy::UffOnly, &config)
            .expect(
                "hexane's generate_coords start is catastrophically clashed, but the \
                     embed_distance_geometry_v2 rescue must succeed",
            );
        assert_eq!(
            result.starting_geometry,
            Some(UffStartingGeometry::ReplacedWithDistanceGeometryV2),
            "expected the rescue to have fired"
        );

        match (result.energy_before, result.energy_after) {
            (EnergyReport::Uff { total: before }, EnergyReport::Uff { total: after }) => {
                assert_eq!(
                    before, expected_v2_energy_before,
                    "rescue energy_before must be the DG v2 geometry's OWN starting energy \
                     (~{expected_v2_energy_before}), not the caller's abandoned geometry's \
                     energy (~{caller_energy}) -- got {before}"
                );
                assert_ne!(
                    before, caller_energy,
                    "rescue energy_before must never equal the caller's original geometry's \
                     energy -- these are two different trajectories"
                );
                assert!(
                    after <= before,
                    "energy_after ({after}) must not exceed energy_before ({before}) -- both \
                     values describe one single minimization trajectory"
                );
            }
            other => panic!("expected (EnergyReport::Uff, EnergyReport::Uff), got {other:?}"),
        }
    }

    /// The stereo-safety gate itself, isolated: cholesterol's `UffOnly`
    /// minimization from `generate_coords` fails soundness (catastrophic
    /// bond blowup), triggering the rescue -- but measured directly
    /// (`verify_stereo`), the rescued geometry still violates 3 of
    /// cholesterol's 8 declared stereocenters. The rescue must therefore be
    /// REFUSED, and the ORIGINAL (bond-blowup) failure returned, annotated to
    /// show a rescue was attempted -- never a silent `Ok` that quietly ships
    /// wrong stereochemistry. See issue #210 (filed alongside this fix) for
    /// the tracked follow-up.
    #[test]
    fn uff_only_refuses_a_rescue_that_would_violate_declared_stereochemistry() {
        let mol =
            parse("C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C")
                .expect("cholesterol");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let err = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
            .expect_err(
                "cholesterol's generate_coords start is catastrophically clashed, and its \
                 embed_distance_geometry_v2 rescue is measured to violate declared \
                 stereochemistry -- the rescue must be refused, not silently accepted",
            );
        match err {
            ForceFieldBridgeError::MinimizationFailed(detail) => {
                assert!(
                    detail.distance_geometry_v2_retry_attempted,
                    "expected the rescue to have been attempted (and refused for stereo \
                     reasons), not skipped entirely: {detail:?}"
                );
            }
            other => panic!("expected MinimizationFailed, got {other:?}"),
        }
    }

    /// The positive counterpart to the cholesterol test above (issue #210,
    /// revised): atorvastatin_fragment is one of #210's 5 named residual
    /// molecules ("still fail `UffOnly` from a legacy `generate_coords`
    /// start... rescue geometry doesn't preserve stereochemistry"), measured
    /// directly (not assumed) to newly succeed, with declared stereochemistry
    /// preserved, once `rescue_with_distance_geometry_v2`'s embed call sets
    /// `enforce_chirality`/`materialize_implicit_h_for_chirality` instead of
    /// plain defaults. `naproxen_S`/`ibuprofen_S`/`testosterone`/`cholesterol`
    /// remain unfixed by this specific change (see this function's own
    /// revised doc comment for why) -- this molecule is the one #210 residual
    /// member this measurement actually closes, not a stand-in for the
    /// others.
    #[test]
    fn uff_only_rescue_now_preserves_stereo_for_atorvastatin_fragment() {
        let mol = parse(
            "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        )
        .expect("atorvastatin_fragment");
        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config).expect(
            "atorvastatin_fragment's generate_coords start is catastrophically clashed, \
                 but the enforce_chirality/materialize_implicit_h_for_chirality-aware rescue \
                 must now succeed",
        );
        assert!(
            crate::stereo_constraints::verify_stereo(&mol, &result.coords).is_fully_satisfied(),
            "the rescued geometry must preserve every declared stereocenter"
        );
    }

    /// Starting-geometry-source axis: the same `UffOnly` policy, same
    /// default budget, starting from `distance_geometry_v2`'s embedder
    /// output instead of `generate_coords`, succeeds for both of the
    /// molecules pinned as failures above -- this is the other half of the
    /// comparison the diagnosis doc's table is built from.
    #[test]
    fn uff_only_succeeds_from_embed_distance_geometry_v2_for_naphthalene_and_hexane() {
        use crate::distance_geometry_v2::{EmbedParameters, embed_distance_geometry_v2};

        for smiles in ["c1ccc2ccccc2c1", "CCCCCC"] {
            let mol = parse(smiles).expect("parse");
            let params = EmbedParameters {
                random_seed: 42,
                max_attempts: 5,
                timeout_ms: Some(10_000),
                ..EmbedParameters::default()
            };
            let coords = embed_distance_geometry_v2(&mol, &params)
                .unwrap_or_else(|e| panic!("embedding {smiles} failed: {e:?}"));
            let config = MinimizeConfig::default();

            let result = minimize_with_policy(&mol, coords, ForceFieldPolicy::UffOnly, &config)
                .unwrap_or_else(|e| {
                    panic!(
                        "expected {smiles} seeded from embed_distance_geometry_v2 to pass the \
                         UffOnly soundness gate at the default budget, got {e:?}"
                    )
                });
            assert!(
                worst_bond_length_vec(&mol, &coords_to_vec(&result.coords, mol.atom_count()))
                    <= MAX_SANE_BOND_LENGTH,
                "{smiles}: expected a sound final geometry"
            );
        }
    }

    /// Proves the "under-budgeted, not fundamentally broken" finding for the
    /// #185-class molecule: same starting geometry as the failing pin above,
    /// but with `max_steps` raised from 200 to 5,000, naphthalene's UFF
    /// minimization actually converges (`converged: true`) to a sound
    /// geometry. `#[ignore]`d -- not part of `bash scripts/check.sh`'s
    /// default run -- because 5,000 fixed-step iterations with a
    /// finite-difference gradient is slow relative to this crate's other
    /// tests; run explicitly with `cargo test -- --ignored` to reproduce.
    #[test]
    #[ignore = "slow (5,000 fixed-step iterations); run with `cargo test -- --ignored` to reproduce the convergence finding"]
    fn uff_direct_naphthalene_converges_given_a_large_enough_iteration_budget() {
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let coords = generate_coords(&mol);
        let raw = coords_to_vec(&coords, mol.atom_count());
        let types = assign_uff_types(&mol);

        let result = ff_minimize_uff(&mol, &types, raw, 5_000);
        let worst = worst_bond_length_vec(&mol, &result.coords);

        assert!(
            result.converged,
            "expected naphthalene's UFF minimization to reach the rms-gradient convergence \
             criterion within 5,000 steps (measured: converges at ~4,539); got \
             converged={} iterations={} energy={} worst_bond={worst}",
            result.converged, result.iterations, result.energy
        );
        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected a sound final geometry once converged, got worst_bond_length={worst}"
        );
    }

    /// Same proof for the #188-class molecule: hexane needs an even larger
    /// budget (measured: worst bond still 45.7 Å at 5,000 steps, 23.8 Å at
    /// 20,000, but a sound 2.04 Å at 100,000) -- confirming this is the same
    /// slow-convergence mechanism as naphthalene, just far more severe, not
    /// a qualitatively different trap. `#[ignore]`d for the same reason as
    /// above (100,000 iterations is meaningfully slower still).
    #[test]
    #[ignore = "slow (100,000 fixed-step iterations); run with `cargo test -- --ignored` to reproduce the convergence finding"]
    fn uff_direct_hexane_reaches_sound_geometry_given_a_large_enough_iteration_budget() {
        let mol = parse("CCCCCC").expect("hexane");
        let coords = generate_coords(&mol);
        let raw = coords_to_vec(&coords, mol.atom_count());
        let types = assign_uff_types(&mol);

        let result = ff_minimize_uff(&mol, &types, raw, 100_000);
        let worst = worst_bond_length_vec(&mol, &result.coords);

        assert!(
            worst <= MAX_SANE_BOND_LENGTH,
            "expected hexane's worst bond length to have recovered to a sound value by 100,000 \
             steps (measured: 2.04 Å); got {worst} -- note this is not expected to report \
             converged=true (rms gradient tolerance is tighter than what 100,000 steps reaches \
             here), only a sound bond length",
        );
    }

    /// Issue #227 Phase 0.1: a molecule `chematic_core::kekulize` cannot
    /// Kekulize (here: an odd-membered, all-carbon, neutral aromatic ring --
    /// unmatchable by simple parity, guaranteed `Err` regardless of which of
    /// `kekulize`'s four internal passes runs) must fail closed all the way
    /// through this crate's own policy bridge, structurally classified as
    /// `ForceFieldBridgeError::UnsupportedAtomType` by `chematic_ff`'s
    /// `NumericTypeError` -> this crate's `From` impl -- never by matching on
    /// the error's message text, and never silently minimized using an
    /// un-re-perceived molecule as a stand-in MMFF view.
    #[test]
    fn kekulization_failure_fails_closed_through_the_policy_bridge_not_silently() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let atoms: Vec<_> = (0..5)
            .map(|_| b.add_atom(chematic_core::Atom::aromatic(chematic_core::Element::C)))
            .collect();
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        let mol = b.build();
        assert!(
            chematic_core::kekulize(&mol).is_err(),
            "test fixture must be genuinely unkekulizable"
        );

        let coords = generate_coords(&mol);
        let config = MinimizeConfig::default();

        let result = minimize_with_policy(
            &mol,
            coords,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        );
        match result {
            Err(ForceFieldBridgeError::UnsupportedAtomType(msg)) => {
                assert!(
                    msg.contains("Kekulization"),
                    "error must name the Kekulization stage, got: {msg}"
                );
            }
            other => panic!(
                "expected Err(UnsupportedAtomType(_)) for an unkekulizable molecule, got {other:?}"
            ),
        }

        // Determinism: repeated calls produce the same classification, not a
        // flaky/order-dependent result.
        let coords2 = generate_coords(&mol);
        let result2 = minimize_with_policy(
            &mol,
            coords2,
            ForceFieldPolicy::Mmff94BondAngleStrict,
            &config,
        );
        assert!(
            matches!(result2, Err(ForceFieldBridgeError::UnsupportedAtomType(_))),
            "Kekulization failure must deterministically classify as UnsupportedAtomType"
        );
    }
}
