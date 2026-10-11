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
                let theta = cos_theta.dacos();
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
                let theta = cos_theta.dacos();
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
                    let theta = cos_theta.dacos();
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

