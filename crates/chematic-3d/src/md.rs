//! Molecular dynamics simulation using Velocity Verlet integration.
//!
//! Provides NVE (microcanonical) and NVT (constant temperature) simulations
//! with optional Berendsen thermostat control. Compatible with WASM.

use std::f64;

use chematic_chem::gasteiger_charges;
use chematic_core::Molecule;
use chematic_ff::{assign_dreiding_types, dreiding_vdw};

use crate::coords::{Coords3D, Point3};

const K_BOLTZMANN: f64 = 0.0019872041; // kcal/(mol·K)
const K_COULOMB: f64 = 332.0637; // kcal·Å/(mol·e²)

/// Thermostat control method.
#[derive(Clone, Debug)]
pub enum Thermostat {
    /// No thermostat (NVE ensemble).
    None,
    /// Berendsen thermostat with coupling constant τ (fs).
    Berendsen { tau_fs: f64 },
}

/// Configuration for molecular dynamics simulation.
#[derive(Clone, Debug)]
pub struct MDConfig {
    /// Timestep in femtoseconds.
    pub timestep_fs: f64,
    /// Number of integration steps.
    pub steps: usize,
    /// Target temperature (K).
    pub temperature_k: f64,
    /// Thermostat control method.
    pub thermostat: Thermostat,
    /// Save trajectory frame every N steps.
    pub save_every: usize,
    /// Include Coulomb interactions.
    pub coulomb: bool,
}

impl Default for MDConfig {
    fn default() -> Self {
        Self {
            timestep_fs: 1.0,
            steps: 100,
            temperature_k: 300.0,
            thermostat: Thermostat::Berendsen { tau_fs: 100.0 },
            save_every: 10,
            coulomb: true,
        }
    }
}

/// A snapshot of the system at a given timestep.
#[derive(Clone, Debug)]
pub struct MDFrame {
    /// Timestep index.
    pub step: usize,
    /// Atomic coordinates (Å).
    pub coords: Coords3D,
    /// Potential energy (kcal/mol).
    pub potential_energy: f64,
    /// Kinetic energy (kcal/mol).
    pub kinetic_energy: f64,
    /// Instantaneous temperature (K).
    pub temperature_k: f64,
}

/// Complete MD trajectory.
#[derive(Clone, Debug)]
pub struct MDTrajectory {
    /// Frames saved at regular intervals.
    pub frames: Vec<MDFrame>,
}

/// Run molecular dynamics simulation.
pub fn run_md(mol: &Molecule, coords: Coords3D, config: &MDConfig) -> MDTrajectory {
    if mol.atom_count() == 0 {
        return MDTrajectory { frames: vec![] };
    }

    let dt = config.timestep_fs / 1000.0;

    // Get atom masses (approximation: use element default)
    let masses = get_atom_masses(mol);

    // Get charges if Coulomb is enabled
    let charges = if config.coulomb {
        gasteiger_charges(mol)
    } else {
        vec![0.0; mol.atom_count()]
    };

    // Initialize velocities from Maxwell-Boltzmann distribution
    let mut prng = crate::prng::Prng::new();
    let mut velocities = initialize_velocities(&masses, config.temperature_k, &mut prng);

    let mut trajectory = MDTrajectory { frames: vec![] };

    let mut current_coords = coords;

    for step in 0..config.steps {
        // Compute forces (finite differences from energy)
        let forces = compute_forces(mol, &current_coords, &charges, config.coulomb);

        // Velocity Verlet: half-step velocity update
        for i in 0..mol.atom_count() {
            let m = masses[i];
            if m > 0.0 {
                velocities[i].x += 0.5 * forces[i].x / m * dt;
                velocities[i].y += 0.5 * forces[i].y / m * dt;
                velocities[i].z += 0.5 * forces[i].z / m * dt;
            }
        }

        // Position update
        for i in 0..mol.atom_count() {
            let m = masses[i];
            let mut p = current_coords.get(chematic_core::AtomIdx(i as u32));
            p.x += velocities[i].x * dt + 0.5 * forces[i].x / m * dt * dt;
            p.y += velocities[i].y * dt + 0.5 * forces[i].y / m * dt * dt;
            p.z += velocities[i].z * dt + 0.5 * forces[i].z / m * dt * dt;
            current_coords.set(chematic_core::AtomIdx(i as u32), p);
        }

        // Recompute forces at new position
        let new_forces = compute_forces(mol, &current_coords, &charges, config.coulomb);

        // Complete velocity update
        for i in 0..mol.atom_count() {
            let m = masses[i];
            if m > 0.0 {
                velocities[i].x += 0.5 * new_forces[i].x / m * dt;
                velocities[i].y += 0.5 * new_forces[i].y / m * dt;
                velocities[i].z += 0.5 * new_forces[i].z / m * dt;
            }
        }

        // Calculate potential energy for frame output
        let current_pot_energy = total_energy(mol, &current_coords, &charges, config.coulomb);

        // Apply thermostat
        let (_kinetic_energy, temperature) =
            compute_kinetic_energy_and_temp(&velocities, &masses, config.temperature_k);

        match config.thermostat {
            Thermostat::None => {}
            Thermostat::Berendsen { tau_fs } => {
                let lambda = if temperature < 1e-6 {
                    1.0
                } else {
                    let tau = tau_fs / 1000.0;
                    (1.0 + (dt / tau) * (config.temperature_k / temperature - 1.0)).sqrt()
                };
                for v in &mut velocities {
                    v.x *= lambda;
                    v.y *= lambda;
                    v.z *= lambda;
                }
            }
        }

        // Save frame
        if (step + 1) % config.save_every == 0 {
            let (ke, temp) =
                compute_kinetic_energy_and_temp(&velocities, &masses, config.temperature_k);
            trajectory.frames.push(MDFrame {
                step: step + 1,
                coords: current_coords.clone(),
                potential_energy: current_pot_energy,
                kinetic_energy: ke,
                temperature_k: temp,
            });
        }
    }

    trajectory
}

fn get_atom_masses(mol: &Molecule) -> Vec<f64> {
    mol.atoms()
        .map(|(_, atom)| {
            let z = atom.element.atomic_number();
            match z {
                1 => 1.008,
                6 => 12.011,
                7 => 14.007,
                8 => 15.999,
                9 => 18.998,
                15 => 30.974,
                16 => 32.06,
                17 => 35.45,
                35 => 79.904,
                53 => 126.90,
                _ => 12.0,
            }
        })
        .collect()
}

fn initialize_velocities(masses: &[f64], temp_k: f64, prng: &mut crate::prng::Prng) -> Vec<Point3> {
    // Maxwell-Boltzmann distribution: v_i ~ sqrt(k_B * T / m)
    // Unit conversion: 1 kcal/mol = 0.01038 amu·Ų/fs²
    // sigma = sqrt(k_B * T * 0.01038 / m) [Å/fs]
    const UNIT_CONVERSION: f64 = 0.01038; // kcal/mol → amu·Ų/fs²

    (0..masses.len())
        .map(|i| {
            let m = masses[i];
            if m > 1e-10 {
                let sigma_sq = K_BOLTZMANN * temp_k * UNIT_CONVERSION / m;
                let sigma = sigma_sq.sqrt();
                Point3::new(
                    gaussian_random(prng) * sigma,
                    gaussian_random(prng) * sigma,
                    gaussian_random(prng) * sigma,
                )
            } else {
                Point3::zero()
            }
        })
        .collect()
}

fn gaussian_random(prng: &mut crate::prng::Prng) -> f64 {
    prng.gaussian_f64()
}

fn compute_forces(
    mol: &Molecule,
    coords: &Coords3D,
    charges: &[f64],
    coulomb: bool,
) -> Vec<Point3> {
    let delta = 1e-4;
    let mut forces = vec![Point3::zero(); mol.atom_count()];

    fn energy_at_delta(
        mol: &Molecule,
        coords: &Coords3D,
        idx: chematic_core::AtomIdx,
        delta: f64,
        axis: impl Fn(&mut Point3, f64),
        charges: &[f64],
        coulomb: bool,
    ) -> f64 {
        let orig = coords.get(idx);
        let mut p = orig;
        axis(&mut p, delta);
        let mut c = coords.clone();
        c.set(idx, p);
        let ep = total_energy(mol, &c, charges, coulomb);

        let mut p = orig;
        axis(&mut p, -delta);
        c.set(idx, p);
        let em = total_energy(mol, &c, charges, coulomb);

        (ep - em) / (2.0 * delta)
    }

    for i in 0..mol.atom_count() {
        let idx = chematic_core::AtomIdx(i as u32);
        forces[i].x = energy_at_delta(mol, coords, idx, delta, |p, d| p.x += d, charges, coulomb);
        forces[i].y = energy_at_delta(mol, coords, idx, delta, |p, d| p.y += d, charges, coulomb);
        forces[i].z = energy_at_delta(mol, coords, idx, delta, |p, d| p.z += d, charges, coulomb);
    }

    forces
}

fn total_energy(mol: &Molecule, coords: &Coords3D, charges: &[f64], coulomb: bool) -> f64 {
    bond_energy(mol, coords)
        + angle_energy(mol, coords)
        + vdw_energy(mol, coords)
        + if coulomb {
            coulomb_energy(coords, charges)
        } else {
            0.0
        }
}

fn bond_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut energy = 0.0;
    let k = 700.0; // Force constant (kcal/mol/Ų)
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let r = coords.get(a1).distance(&coords.get(a2));
        let sym1 = mol.atom(a1).element.symbol();
        let sym2 = mol.atom(a2).element.symbol();
        let ideal = ideal_bond_len(sym1, sym2, bond.order);
        let delta = r - ideal;
        energy += 0.5 * k * delta * delta;
    }
    energy
}

fn angle_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut energy = 0.0;
    let k = 100.0; // Force constant (kcal/mol/rad²)

    for b_idx in 0..mol.atom_count() {
        let b = chematic_core::AtomIdx(b_idx as u32);
        let neighbors: Vec<_> = mol.neighbors(b).map(|(nb, _)| nb).collect();

        if neighbors.len() < 2 {
            continue;
        }

        for (i, &na) in neighbors.iter().enumerate() {
            for &nb in &neighbors[i + 1..] {
                let pa = coords.get(na);
                let pb = coords.get(b);
                let pc = coords.get(nb);

                let v1 = Point3::new(pa.x - pb.x, pa.y - pb.y, pa.z - pb.z);
                let v2 = Point3::new(pc.x - pb.x, pc.y - pb.y, pc.z - pb.z);
                let angle = angle_from_vectors(&v1, &v2);

                let sym = mol.atom(b).element.symbol();
                let hyb = atom_hybridization(mol, b);
                let ideal = ideal_angle_rad(sym, hyb);

                let delta = angle - ideal;
                energy += 0.5 * k * delta * delta;
            }
        }
    }
    energy
}

fn vdw_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut energy = 0.0;
    let n = mol.atom_count();

    // Assign DREIDING types to all atoms
    let types = assign_dreiding_types(mol);

    // Build exclusion list: skip bonded pairs (1-2 and 1-3)
    let mut excluded = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
    }
    // Also exclude 1-3 pairs (atoms separated by one bond)
    let bonds_vec: Vec<_> = mol.bonds().collect();
    for (_, bond1) in &bonds_vec {
        for (_, bond2) in &bonds_vec {
            if bond1.atom2 == bond2.atom1 {
                let i = bond1.atom1.0 as usize;
                let k = bond2.atom2.0 as usize;
                if i != k {
                    excluded.insert((i.min(k), i.max(k)));
                }
            }
        }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            // Skip excluded pairs
            if excluded.contains(&(i, j)) {
                continue;
            }

            let ii = chematic_core::AtomIdx(i as u32);
            let jj = chematic_core::AtomIdx(j as u32);
            let pi = coords.get(ii);
            let pj = coords.get(jj);
            let r = pi.distance(&pj);

            if r < 0.1 {
                continue;
            }

            // Get VDW parameters from DREIDING
            let (r_i, eps_i) = dreiding_vdw(types[i]);
            let (r_j, eps_j) = dreiding_vdw(types[j]);

            // Lorentz-Berthelot combining rules
            let sigma = 0.5 * (r_i + r_j);
            let epsilon = (eps_i * eps_j).sqrt();

            // Lennard-Jones 12-6: E = eps * [(sigma/r)^12 - 2*(sigma/r)^6]
            let sig_r = sigma / r;
            let sig_r6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r;
            let sig_r12 = sig_r6 * sig_r6;

            energy += epsilon * (sig_r12 - 2.0 * sig_r6);
        }
    }
    energy
}

fn ideal_bond_len(sym1: &str, sym2: &str, order: chematic_core::BondOrder) -> f64 {
    use chematic_core::BondOrder;
    let (a, b) = if sym1 <= sym2 {
        (sym1, sym2)
    } else {
        (sym2, sym1)
    };
    match (a, b, order) {
        ("C", "C", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.540,
        ("C", "C", BondOrder::Double) => 1.340,
        ("C", "C", BondOrder::Triple) => 1.204,
        ("C", "C", BondOrder::Aromatic) => 1.395,
        ("C", "H", _) => 1.090,
        ("C", "N", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.469,
        ("C", "N", BondOrder::Double) => 1.279,
        ("C", "N", BondOrder::Triple) => 1.158,
        ("C", "N", BondOrder::Aromatic) => 1.340,
        ("C", "O", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.427,
        ("C", "O", BondOrder::Double) => 1.217,
        ("C", "O", BondOrder::Aromatic) => 1.355,
        ("C", "S", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.819,
        ("H", "N", _) => 1.010,
        ("H", "O", _) => 0.960,
        ("N", "N", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.450,
        ("N", "N", BondOrder::Double) => 1.252,
        ("N", "N", BondOrder::Triple) => 1.098,
        ("N", "O", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.463,
        ("N", "O", BondOrder::Double) => 1.240,
        ("O", "O", BondOrder::Single | BondOrder::Up | BondOrder::Down) => 1.450,
        ("O", "O", BondOrder::Double) => 1.207,
        _ => 1.5,
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Hybridization {
    SP,
    SP2,
    SP3,
}

fn atom_hybridization(mol: &Molecule, idx: chematic_core::AtomIdx) -> Hybridization {
    use chematic_core::BondOrder;
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

fn angle_from_vectors(v1: &Point3, v2: &Point3) -> f64 {
    let dot = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    let mag1 = (v1.x * v1.x + v1.y * v1.y + v1.z * v1.z).sqrt();
    let mag2 = (v2.x * v2.x + v2.y * v2.y + v2.z * v2.z).sqrt();
    if mag1 > 1e-10 && mag2 > 1e-10 {
        let cos_angle = (dot / (mag1 * mag2)).clamp(-1.0, 1.0);
        cos_angle.acos()
    } else {
        0.0
    }
}

fn coulomb_energy(coords: &Coords3D, charges: &[f64]) -> f64 {
    let mut energy = 0.0;
    let n = coords.atom_count();
    for i in 0..n {
        for j in (i + 1)..n {
            let pi = coords.get(chematic_core::AtomIdx(i as u32));
            let pj = coords.get(chematic_core::AtomIdx(j as u32));
            let r = pi.distance(&pj);
            if r > 1e-6 {
                energy += K_COULOMB * charges[i] * charges[j] / r;
            }
        }
    }
    energy
}

fn compute_kinetic_energy_and_temp(
    velocities: &[Point3],
    masses: &[f64],
    _target_temp: f64,
) -> (f64, f64) {
    let mut ke = 0.0;
    for (i, v) in velocities.iter().enumerate() {
        let speed_sq = v.x * v.x + v.y * v.y + v.z * v.z;
        ke += 0.5 * masses[i] * speed_sq;
    }
    let dof = (3 * velocities.len() - 3).max(1) as f64; // 3N - 3 translational
    let temp = 2.0 * ke / (dof * K_BOLTZMANN);
    (ke, temp)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_md_ethane_nve() {
        let mol = parse("CC").expect("ethane");
        let coords = crate::generate_coords(&mol);
        let config = MDConfig {
            timestep_fs: 0.5,
            steps: 50,
            temperature_k: 300.0,
            thermostat: Thermostat::None,
            save_every: 10,
            coulomb: false,
        };
        let traj = run_md(&mol, coords, &config);
        assert!(!traj.frames.is_empty());
        assert!(traj.frames[0].potential_energy > -1e10); // Sanity check
    }

    #[test]
    fn test_md_ethane_nvt() {
        let mol = parse("CC").expect("ethane");
        let coords = crate::generate_coords(&mol);
        let config = MDConfig {
            timestep_fs: 0.5,
            steps: 100,
            temperature_k: 300.0,
            thermostat: Thermostat::Berendsen { tau_fs: 50.0 },
            save_every: 20,
            coulomb: false,
        };
        let traj = run_md(&mol, coords, &config);
        assert!(!traj.frames.is_empty());
        // Check that simulation completes without crashing
        let final_frame = traj.frames.last().unwrap();
        assert!(final_frame.temperature_k > 0.0);
    }

    #[test]
    fn test_berendsen_zero_temperature_no_nan() {
        let mol = parse("C").expect("methane");
        let coords = crate::generate_coords(&mol);
        let config = MDConfig {
            timestep_fs: 1.0,
            steps: 1,
            temperature_k: 300.0,
            thermostat: Thermostat::Berendsen { tau_fs: 100.0 },
            save_every: 1,
            coulomb: false,
        };
        let traj = run_md(&mol, coords, &config);
        assert!(!traj.frames.is_empty());
        let frame = traj.frames.first().unwrap();
        for (i, coord) in frame.coords.points.iter().enumerate() {
            assert!(coord.x.is_finite(), "coord[{}].x is NaN", i);
            assert!(coord.y.is_finite(), "coord[{}].y is NaN", i);
            assert!(coord.z.is_finite(), "coord[{}].z is NaN", i);
        }
    }
}
