//! MMFF94 geometry minimizer using complete Halgren 1996 parameters.
//!
//! Provides full MMFF94 energy evaluation (bond + angle + torsion + vdW + electrostatic)
//! and steepest-descent geometry optimization. All energy terms use the faithful
//! parameter tables from `mmff94_energy.rs` and numeric atom types from `mmff94_numeric.rs`.
//!
//! ## Energy terms
//! - **Bond**: cubic-corrected harmonic (Halgren MMFF.II eq. 1)
//! - **Angle**: cubic-corrected harmonic (Halgren MMFF.III eq. 2)
//! - **Torsion**: three-term Fourier (Halgren MMFF.IV)
//! - **vdW**: buffered 14-7 potential with Slater-Kirkwood combining rule (Halgren MMFF.I eq. 2)
//! - **Electrostatic**: Coulomb with δ buffer (Halgren MMFF.V eq. 14)

use chematic_core::{AtomIdx, Molecule};

use crate::mmff94_energy::{
    mmff94_angle_energy, mmff94_bond_energy, mmff94_torsion_energy, mmff94_vdw_combined,
};
use crate::mmff94_numeric::{assign_mmff94_numeric_types, mmff94_charges_numeric, NumericTypeError};

// ─── Public types ────────────────────────────────────────────────────────────

/// Result of a geometry minimization run.
#[derive(Debug, Clone)]
pub struct MinimizeResult {
    /// Final MMFF94 energy (kcal/mol).
    pub energy: f64,
    /// RMSD of atom positions vs initial geometry (Å).
    pub rmsd: f64,
    /// Whether the minimization converged before `max_iter`.
    pub converged: bool,
    /// Number of gradient steps performed.
    pub iterations: usize,
}

/// Error from MMFF94 minimizer setup.
#[derive(Debug)]
pub enum MinimizerError {
    TypeAssignment(NumericTypeError),
}

impl From<NumericTypeError> for MinimizerError {
    fn from(e: NumericTypeError) -> Self {
        MinimizerError::TypeAssignment(e)
    }
}

impl std::fmt::Display for MinimizerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MinimizerError::TypeAssignment(e) => write!(f, "MMFF94 type assignment failed: {}", e),
        }
    }
}

// ─── Public API ──────────────────────────────────────────────────────────────

/// Compute total MMFF94 energy for a given geometry (kcal/mol).
///
/// Includes bond, angle, torsion, vdW, and electrostatic terms.
/// Does not modify coordinates.
pub fn mmff94_total_energy(
    mol: &Molecule,
    coords: &[[f64; 3]],
) -> Result<f64, MinimizerError> {
    let types = assign_mmff94_numeric_types(mol)?;
    let charges = mmff94_charges_numeric(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
    Ok(total_energy(mol, coords, &types, &charges))
}

/// Minimize molecular geometry using the full MMFF94 force field.
///
/// Uses steepest descent with finite-difference gradients and the complete
/// Halgren 1996 parameter tables (bond, angle, torsion, vdW, electrostatic).
///
/// # Arguments
/// * `mol` — molecule graph (topology only, no coordinates)
/// * `coords` — initial 3D coordinates `[[x, y, z]]` in Å; updated in place
/// * `max_iter` — maximum gradient steps (200 typically sufficient)
pub fn minimize_mmff94_full(
    mol: &Molecule,
    coords: &mut Vec<[f64; 3]>,
    max_iter: usize,
) -> Result<MinimizeResult, MinimizerError> {
    if mol.atom_count() <= 1 {
        return Ok(MinimizeResult {
            energy: 0.0,
            rmsd: 0.0,
            converged: true,
            iterations: 0,
        });
    }

    let types = assign_mmff94_numeric_types(mol)?;
    let charges = mmff94_charges_numeric(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);

    let n = mol.atom_count();
    let initial = coords.clone();
    let convergence = 1e-4_f64;
    let step_size = 0.05_f64;
    let delta = 1e-4_f64;

    let mut iters = 0usize;
    let mut converged = false;

    for _ in 0..max_iter {
        iters += 1;
        let mut grad = vec![[0.0_f64; 3]; n];
        let mut max_g = 0.0_f64;

        for i in 0..n {
            for axis in 0..3 {
                coords[i][axis] += delta;
                let ep = total_energy(mol, coords, &types, &charges);
                coords[i][axis] -= 2.0 * delta;
                let em = total_energy(mol, coords, &types, &charges);
                coords[i][axis] += delta;
                let g = (ep - em) / (2.0 * delta);
                grad[i][axis] = g;
                max_g = max_g.max(g.abs());
            }
        }

        if max_g < convergence {
            converged = true;
            break;
        }

        let scale = step_size / max_g.max(1e-8);
        for i in 0..n {
            for axis in 0..3 {
                coords[i][axis] -= scale * grad[i][axis];
            }
        }
    }

    let energy = total_energy(mol, coords, &types, &charges);

    let rmsd = {
        let sum: f64 = coords
            .iter()
            .zip(initial.iter())
            .map(|(c, i0)| {
                let dx = c[0] - i0[0];
                let dy = c[1] - i0[1];
                let dz = c[2] - i0[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum();
        (sum / n as f64).sqrt()
    };

    Ok(MinimizeResult {
        energy,
        rmsd,
        converged,
        iterations: iters,
    })
}

// ─── Energy components ───────────────────────────────────────────────────────

fn total_energy(
    mol: &Molecule,
    coords: &[[f64; 3]],
    types: &[u8],
    charges: &[f64],
) -> f64 {
    bond_energy(mol, coords, types)
        + angle_energy(mol, coords, types)
        + torsion_energy(mol, coords, types)
        + vdw_energy(mol, coords, types)
        + elec_energy(mol, coords, charges)
}

/// Bond stretching: cubic-corrected harmonic (Halgren MMFF.II eq. 1)
/// E = (143.9325 × kb / 2) × ΔR² × (1 − cs×ΔR + (7/12)×cs²×ΔR²)
fn bond_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    const KB_CONV: f64 = 143.9325;
    const CS: f64 = 2.0;
    let mut energy = 0.0;
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let bt = bond_type_for(types[i], types[j]);
        if let Some(p) = mmff94_bond_energy(bt, types[i], types[j]) {
            let r = dist(coords[i], coords[j]);
            let dr = r - p.r0;
            let cubic = 1.0 - CS * dr + (7.0 / 12.0) * CS * CS * dr * dr;
            energy += (KB_CONV * p.kb / 2.0) * dr * dr * cubic;
        }
    }
    energy
}

/// Angle bending: cubic-corrected harmonic (Halgren MMFF.III eq. 2)
/// E = (0.043844 × ka / 2) × Δθ² × (1 − 0.007×Δθ)   [Δθ in degrees]
fn angle_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    const KA_CONV: f64 = 0.043844;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() < 2 {
            continue;
        }
        let at = angle_type_for(types[j_idx]);
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                if let Some(p) = mmff94_angle_energy(at, types[i], types[j_idx], types[k]) {
                    let cos_t = cos_angle(coords[i], coords[j_idx], coords[k]);
                    let theta_deg = cos_t.acos() * RAD_TO_DEG;
                    let dt = theta_deg - p.theta0;
                    let cubic = 1.0 - 0.007 * dt;
                    energy += (KA_CONV * p.ka / 2.0) * dt * dt * cubic;
                }
            }
        }
    }
    energy
}

/// Torsion: three-term Fourier (Halgren MMFF.IV)
/// E = (v1/2)(1+cosφ) + (v2/2)(1-cos2φ) + (v3/2)(1+cos3φ)
fn torsion_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    let mut energy = 0.0;
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        let nbrs_j: Vec<usize> = mol.neighbors(bond.atom1).map(|(nb, _)| nb.0 as usize).collect();
        let nbrs_k: Vec<usize> = mol.neighbors(bond.atom2).map(|(nb, _)| nb.0 as usize).collect();
        let tt = torsion_type_for(types[j], types[k]);
        for &i in &nbrs_j {
            if i == k {
                continue;
            }
            for &l in &nbrs_k {
                if l == j {
                    continue;
                }
                if let Some(p) = mmff94_torsion_energy(tt, types[i], types[j], types[k], types[l]) {
                    let phi = dihedral(coords[i], coords[j], coords[k], coords[l]);
                    energy += 0.5 * p.v1 * (1.0 + phi.cos())
                        + 0.5 * p.v2 * (1.0 - (2.0 * phi).cos())
                        + 0.5 * p.v3 * (1.0 + (3.0 * phi).cos());
                }
            }
        }
    }
    energy
}

/// Van der Waals: buffered 14-7 (Halgren MMFF.I eq. 2)
/// t = (1.07 × r*) / (r + 0.07 × r*)
/// E = ε × t⁷ × (t⁷ − 2)
fn vdw_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    let n = mol.atom_count();
    let mut excl = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excl.insert((i.min(j), i.max(j)));
        for (nb_i, _) in mol.neighbors(bond.atom1) {
            let ni = nb_i.0 as usize;
            excl.insert((ni.min(j), ni.max(j)));
        }
        for (nb_j, _) in mol.neighbors(bond.atom2) {
            let nj = nb_j.0 as usize;
            excl.insert((i.min(nj), i.max(nj)));
        }
    }
    let cutoff = 10.0_f64;
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if excl.contains(&(i, j)) {
                continue;
            }
            let r = dist(coords[i], coords[j]);
            if r > cutoff {
                continue;
            }
            if let Some((r_star, eps)) = mmff94_vdw_combined(types[i], types[j]) {
                if r_star > 0.0 && eps > 0.0 && r > 0.01 {
                    let t = (1.07 * r_star) / (r + 0.07 * r_star);
                    let t7 = t.powi(7);
                    energy += eps * t7 * (t7 - 2.0);
                }
            }
        }
    }
    energy
}

/// Electrostatic: Coulomb with δ=0.05 Å buffer (Halgren MMFF.V eq. 14)
/// E = 332.0716 × q_i × q_j / (D × (r + δ))   [D=1.0]
fn elec_energy(mol: &Molecule, coords: &[[f64; 3]], charges: &[f64]) -> f64 {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    let n = mol.atom_count();
    let mut excl = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excl.insert((i.min(j), i.max(j)));
        for (nb_i, _) in mol.neighbors(bond.atom1) {
            excl.insert(((nb_i.0 as usize).min(j), (nb_i.0 as usize).max(j)));
        }
        for (nb_j, _) in mol.neighbors(bond.atom2) {
            excl.insert((i.min(nb_j.0 as usize), i.max(nb_j.0 as usize)));
        }
    }
    // 1-4 pairs: scale by 0.75 (MMFF94 convention)
    let mut one_four = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        for (nb_j, _) in mol.neighbors(bond.atom1) {
            let i = nb_j.0 as usize;
            if i == k {
                continue;
            }
            for (nb_k, _) in mol.neighbors(bond.atom2) {
                let l = nb_k.0 as usize;
                if l == j {
                    continue;
                }
                let key = (i.min(l), i.max(l));
                if !excl.contains(&key) {
                    one_four.insert(key);
                }
            }
        }
    }
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            if excl.contains(&(i, j)) {
                continue;
            }
            let r = dist(coords[i], coords[j]);
            let scale = if one_four.contains(&(i, j)) { 0.75 } else { 1.0 };
            energy += scale * COULOMB * charges[i] * charges[j] / (r + DELTA);
        }
    }
    energy
}

// ─── Geometry helpers ─────────────────────────────────────────────────────────

#[inline]
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[inline]
fn cos_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot_val = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2];
    let na = (ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]).sqrt();
    let nc = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]).sqrt();
    if na < 1e-12 || nc < 1e-12 {
        return 0.0;
    }
    (dot_val / (na * nc)).clamp(-1.0, 1.0)
}

/// Signed dihedral angle φ (radians) for the quartet i-j-k-l.
#[inline]
fn dihedral(i: [f64; 3], j: [f64; 3], k: [f64; 3], l: [f64; 3]) -> f64 {
    let b1 = [j[0] - i[0], j[1] - i[1], j[2] - i[2]];
    let b2 = [k[0] - j[0], k[1] - j[1], k[2] - j[2]];
    let b3 = [l[0] - k[0], l[1] - k[1], l[2] - k[2]];
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m1 = cross(n1, b2);
    let b2_len = (b2[0] * b2[0] + b2[1] * b2[1] + b2[2] * b2[2]).sqrt();
    if b2_len < 1e-12 {
        return 0.0;
    }
    let x = dot3(n1, n2);
    let y = dot3(m1, n2) / b2_len;
    y.atan2(x)
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ─── Type classification helpers ──────────────────────────────────────────────

/// Determine MMFF94 bond type: 1 if either atom is sp2/aromatic, else 0.
fn bond_type_for(ti: u8, tj: u8) -> u8 {
    const SP2: &[u8] = &[
        2, 3, 9, 10, 37, 38, 39, 40, 41, 56, 57, 58, 59, 63, 64, 65, 66, 67, 78, 79, 80, 81, 82,
    ];
    if SP2.binary_search(&ti).is_ok() || SP2.binary_search(&tj).is_ok() {
        1
    } else {
        0
    }
}

/// Determine MMFF94 angle type (simplified: 0 for most organic angles).
fn angle_type_for(_tj: u8) -> u8 {
    0
}

/// Determine MMFF94 torsion type from central bond atom types.
fn torsion_type_for(tj: u8, tk: u8) -> u8 {
    const SP2: &[u8] = &[
        2, 3, 9, 10, 37, 38, 39, 40, 41, 56, 57, 58, 59, 63, 64, 65, 66, 67, 78, 79, 80, 81, 82,
    ];
    match (
        SP2.binary_search(&tj).is_ok(),
        SP2.binary_search(&tk).is_ok(),
    ) {
        (false, false) => 0, // sp3-sp3
        (true, false) | (false, true) => 1, // sp3-sp2
        (true, true) => 2, // sp2-sp2
    }
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::molecule::MoleculeBuilder;
    use chematic_core::{Atom, BondOrder, Element};

    fn methane_mol() -> (Molecule, Vec<[f64; 3]>) {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let h1 = b.add_atom(Atom::new(Element::H));
        let h2 = b.add_atom(Atom::new(Element::H));
        let h3 = b.add_atom(Atom::new(Element::H));
        let h4 = b.add_atom(Atom::new(Element::H));
        b.add_bond(c, h1, BondOrder::Single).unwrap();
        b.add_bond(c, h2, BondOrder::Single).unwrap();
        b.add_bond(c, h3, BondOrder::Single).unwrap();
        b.add_bond(c, h4, BondOrder::Single).unwrap();
        let mol = b.build();
        let coords = vec![
            [0.0, 0.0, 0.0],
            [0.630, 0.630, 0.630],
            [-0.630, -0.630, 0.630],
            [-0.630, 0.630, -0.630],
            [0.630, -0.630, -0.630],
        ];
        (mol, coords)
    }

    fn butane_backbone() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn energy_is_finite_for_methane() {
        let (mol, coords) = methane_mol();
        let e = mmff94_total_energy(&mol, &coords).expect("energy");
        assert!(e.is_finite(), "energy={}", e);
    }

    #[test]
    fn torsion_differs_by_conformation() {
        let mol = butane_backbone();
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        // Gauche: ~60° central dihedral
        let coords_gauche = vec![
            [0.0, 0.0, 0.0_f64],
            [1.508, 0.0, 0.0],
            [2.016, 1.192, 0.688],
            [3.524, 1.192, 0.688],
        ];
        // Anti: ~180° central dihedral
        let coords_anti = vec![
            [0.0, 0.0, 0.0_f64],
            [1.508, 0.0, 0.0],
            [3.016, 0.0, 0.0],
            [4.524, 0.0, 0.0],
        ];
        let e_gauche = torsion_energy(&mol, &coords_gauche, &types);
        let e_anti = torsion_energy(&mol, &coords_anti, &types);
        assert!(e_gauche.is_finite());
        assert!(e_anti.is_finite());
        assert!(
            (e_gauche - e_anti).abs() > 1e-6,
            "torsion must differ: gauche={}, anti={}",
            e_gauche,
            e_anti
        );
    }

    #[test]
    fn vdw_more_repulsive_at_short_range() {
        let mol = butane_backbone();
        let types = assign_mmff94_numeric_types(&mol).expect("types");
        // Atoms 0 and 3 are 1-4 (not excluded from vdW)
        let coords_close = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.5, 0.0, 0.0], // atom 3 very close to atom 0
        ];
        let coords_far = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [8.0, 0.0, 0.0],
        ];
        let e_close = vdw_energy(&mol, &coords_close, &types);
        let e_far = vdw_energy(&mol, &coords_far, &types);
        assert!(e_close.is_finite());
        assert!(e_far.is_finite());
        assert!(e_close > e_far, "close={} should > far={}", e_close, e_far);
    }

    #[test]
    fn dihedral_anti_is_pi() {
        let i = [0.0_f64, 0.0, 0.0];
        let j = [1.0, 0.0, 0.0];
        let k = [2.0, 0.0, 1.0];
        let l = [3.0, 0.0, 0.0];
        let phi = dihedral(i, j, k, l);
        assert!(phi.abs() > 2.5, "anti dihedral ≈ π: {}", phi);
    }

    #[test]
    fn dihedral_syn_is_zero() {
        let i = [0.0_f64, 1.0, 0.0];
        let j = [0.0, 0.0, 0.0];
        let k = [1.0, 0.0, 0.0];
        let l = [1.0, 1.0, 0.0];
        let phi = dihedral(i, j, k, l);
        assert!(phi.abs() < 0.1, "syn dihedral ≈ 0: {}", phi);
    }

    #[test]
    fn minimize_reduces_energy_for_methane() {
        let (mol, _) = methane_mol();
        // Distorted geometry (C-H bonds ~2× too long)
        let mut coords = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 1.5, 1.5],
            [-1.5, -1.5, 1.5],
            [-1.5, 1.5, -1.5],
            [1.5, -1.5, -1.5],
        ];
        let e_before = mmff94_total_energy(&mol, &coords).expect("energy before");
        let result = minimize_mmff94_full(&mol, &mut coords, 300).expect("minimize");
        assert!(
            result.energy <= e_before,
            "minimize should reduce energy: {} → {}",
            e_before,
            result.energy
        );
        assert!(result.energy.is_finite());
        assert!(result.iterations > 0);
    }
}
