//! Simplified force-field geometry minimization for molecular structures.
//!
//! Uses gradient descent with finite differences over three energy terms:
//! bond stretching, angle bending, and VDW repulsion.

use std::collections::HashSet;
use std::f64::consts::PI;

use chematic_core::{AtomIdx, BondOrder, Molecule};

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Configuration for the minimization algorithm.
pub struct MinimizeConfig {
    /// Maximum number of gradient-descent steps.
    pub max_steps: usize,
    /// Base step size for coordinate updates (scaled by max gradient).
    pub step_size: f64,
    /// Convergence threshold: stop when max gradient component < this value.
    pub convergence: f64,
}

impl Default for MinimizeConfig {
    fn default() -> Self {
        Self {
            max_steps: 200,
            step_size: 0.05,
            convergence: 1e-4,
        }
    }
}

/// Minimize molecular geometry using default configuration.
pub fn minimize(mol: &Molecule, coords: Coords3D) -> Coords3D {
    minimize_with_config(mol, coords, &MinimizeConfig::default())
}

/// Minimize molecular geometry using the provided configuration.
pub fn minimize_with_config(mol: &Molecule, coords: Coords3D, config: &MinimizeConfig) -> Coords3D {
    if mol.atom_count() <= 1 {
        return coords;
    }

    let mut c = coords;
    let delta = 1e-4;

    // Central-difference partial derivative of total_energy w.r.t. one
    // component of atom `idx`. `axis` selects x/y/z by mutating the chosen
    // component before each energy evaluation.
    fn partial(
        mol: &Molecule,
        c: &mut Coords3D,
        idx: AtomIdx,
        delta: f64,
        axis: impl Fn(&mut Point3, f64),
    ) -> f64 {
        let orig = c.get(idx);
        let mut p = orig;
        axis(&mut p, delta);
        c.set(idx, p);
        let ep = total_energy(mol, c);
        let mut p = orig;
        axis(&mut p, -delta);
        c.set(idx, p);
        let em = total_energy(mol, c);
        c.set(idx, orig);
        (ep - em) / (2.0 * delta)
    }

    for _ in 0..config.max_steps {
        let mut grad = vec![Point3::zero(); mol.atom_count()];
        let mut max_grad = 0.0f64;

        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            grad[i].x = partial(mol, &mut c, idx, delta, |p, d| p.x += d);
            grad[i].y = partial(mol, &mut c, idx, delta, |p, d| p.y += d);
            grad[i].z = partial(mol, &mut c, idx, delta, |p, d| p.z += d);

            let gmax = grad[i].x.abs().max(grad[i].y.abs()).max(grad[i].z.abs());
            if gmax > max_grad {
                max_grad = gmax;
            }
        }

        if max_grad < config.convergence {
            break;
        }

        // Scale step so the largest gradient component moves `step_size`.
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

// ---------------------------------------------------------------------------
// Total energy
// ---------------------------------------------------------------------------

fn total_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    bond_energy(mol, coords) + angle_energy(mol, coords) + vdw_energy(mol, coords)
}

// ---------------------------------------------------------------------------
// Bond stretching energy
// ---------------------------------------------------------------------------

/// Ideal bond length (Å) based solely on bond order — used for the minimizer.
fn ideal_bond_len_by_order(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1.54,
        BondOrder::Double => 1.34,
        BondOrder::Triple => 1.20,
        BondOrder::Quadruple => 1.20,
        BondOrder::Aromatic => 1.40,
    }
}

fn bond_energy(mol: &Molecule, coords: &Coords3D) -> f64 {
    let mut energy = 0.0;
    for (_, bond) in mol.bonds() {
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let r = coords.get(a1).distance(&coords.get(a2));
        let r0 = ideal_bond_len_by_order(bond.order);
        let dr = r - r0;
        energy += 0.5 * 700.0 * dr * dr;
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
        let deg = neighbors.len();

        if deg < 2 {
            continue;
        }

        // Ideal angle based on degree
        let theta0 = match deg {
            4 => 109.47_f64.to_radians(),
            3 => 120.0_f64.to_radians(),
            2 => PI, // 180°
            _ => 109.47_f64.to_radians(), // > 4 neighbors: use tetrahedral
        };

        let pb = coords.get(b);

        // Iterate unique pairs of neighbors (i < j) to avoid double counting
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
                energy += 0.5 * 100.0 * dtheta * dtheta;
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
    let cutoff = 5.0_f64;

    // Build exclusion set: directly bonded pairs + 1-3 pairs (share a common neighbor)
    let mut excluded: HashSet<(usize, usize)> = HashSet::new();

    // 1-2 exclusions (directly bonded)
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
    }

    // 1-3 exclusions (atoms sharing a common bonded neighbor)
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

            if r < 0.01 {
                continue;
            }
            if r >= cutoff {
                continue;
            }

            let ratio = 2.0 / r;
            let ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
            let ratio12 = ratio6 * ratio6;
            energy += 0.05 * ratio12;
        }
    }

    energy
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
        assert!(d > 1.2 && d < 1.8, "C-C distance={d:.3}, expected 1.2-1.8 Å");
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
        assert!(min_d > 0.8, "atom clash in benzene: min distance={min_d:.3}");
    }

    #[test]
    fn test_disconnected_no_clash() {
        let mol = parse("CC.CC").unwrap();
        let coords = generate_coords(&mol);
        let result = minimize(&mol, coords);
        let min_d = all_pairs_min_dist(&result, mol.atom_count());
        assert!(min_d > 0.8, "atom clash in disconnected: min distance={min_d:.3}");
    }

    #[test]
    fn test_default_config_no_panic() {
        let mol = parse("CC(=O)O").unwrap(); // acetic acid
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
        // Second minimization shouldn't increase energy significantly
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
}
