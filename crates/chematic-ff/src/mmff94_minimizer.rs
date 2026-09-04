//! MMFF94 geometry minimizer using complete Halgren 1996 parameters.
//!
//! Provides full MMFF94 energy evaluation (bond + angle + torsion + vdW + electrostatic)
//! and geometry optimization with two algorithms:
//! - **Steepest descent** (`minimize_mmff94_full`) — robust, simple
//! - **L-BFGS** (`minimize_mmff94_lbfgs`) — faster convergence, quasi-Newton
//!
//! ## Energy terms
//! - **Bond**: cubic-corrected harmonic (Halgren MMFF.II eq. 1)
//! - **Angle**: cubic-corrected harmonic (Halgren MMFF.III eq. 2)
//! - **Torsion**: three-term Fourier (Halgren MMFF.IV)
//! - **vdW**: buffered 14-7 potential with Slater-Kirkwood combining rule (Halgren MMFF.I eq. 2)
//! - **Electrostatic**: Coulomb with δ buffer (Halgren MMFF.V eq. 14)

use std::collections::VecDeque;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::find_sssr;

use crate::mmff94_energy::{
    AngleEnergyParams, BondEnergyParams, TorsionEnergyParams, mmff94_angle_energy_resolved,
    mmff94_bond_energy_resolved, mmff94_oop, mmff94_stbn, mmff94_torsion_energy,
    mmff94_vdw_combined,
};
use crate::mmff94_numeric::{
    NumericTypeError, assign_mmff94_numeric_types_with_view, mmff94_charges_numeric,
};

type CoordVec = Vec<[f64; 3]>;
type LbfgsHistory = VecDeque<(CoordVec, CoordVec, f64)>;
type VdwPairs = Vec<(usize, usize)>;
type ElectrostaticPairs = Vec<(usize, usize, f64)>;

#[derive(Clone, Copy)]
struct PreparedBond {
    i: usize,
    j: usize,
    params: BondEnergyParams,
}

#[derive(Clone, Copy)]
struct PreparedAngle {
    i: usize,
    j: usize,
    k: usize,
    params: AngleEnergyParams,
}

#[derive(Clone, Copy)]
struct PreparedTorsion {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    params: TorsionEnergyParams,
}

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

/// Per-term MMFF94 energy breakdown (kcal/mol). Includes all 7 Halgren 1996 energy terms.
#[derive(Debug, Clone, Copy)]
pub struct EnergyBreakdown {
    pub bond: f64,
    pub angle: f64,
    /// Stretch-bend coupling (STRE-BEN, Halgren MMFF.V)
    pub stretch_bend: f64,
    pub torsion: f64,
    /// Out-of-plane bending for sp2 atoms (Halgren MMFF.VI)
    pub oop: f64,
    pub vdw: f64,
    pub electrostatic: f64,
    pub total: f64,
}

/// Reusable MMFF94 energy evaluator for one molecular topology.
///
/// Atom typing, MMFF94's aromatic view, charges, and ring perception are
/// topology-only. Keeping them together avoids repeating that work when a
/// caller evaluates many nearby geometries, as finite-difference gradients do.
#[derive(Clone)]
pub struct Mmff94EnergyModel {
    types: Vec<u8>,
    mmff_mol: Molecule,
    charges: Vec<f64>,
    rings: Vec<Vec<AtomIdx>>,
    vdw_pairs: Vec<(usize, usize)>,
    electrostatic_pairs: Vec<(usize, usize, f64)>,
    bonds: Vec<PreparedBond>,
    angles: Vec<PreparedAngle>,
    torsions: Vec<PreparedTorsion>,
}

impl Mmff94EnergyModel {
    /// Prepare the topology-dependent MMFF94 state once.
    pub fn new(mol: &Molecule) -> Result<Self, MinimizerError> {
        let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
        let charges = mmff94_charges_numeric(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
        let rings = find_sssr(mol).rings().to_vec();
        let (vdw_pairs, electrostatic_pairs) = build_nonbonded_pairs(mol);
        let bonds = build_bond_terms(&mmff_mol, &types);
        let angles = build_angle_terms(&mmff_mol, &types, &rings);
        let torsions = build_torsion_terms(&mmff_mol, &types);
        Ok(Self {
            types,
            mmff_mol,
            charges,
            rings,
            vdw_pairs,
            electrostatic_pairs,
            bonds,
            angles,
            torsions,
        })
    }

    /// Evaluate total MMFF94 energy for coordinates matching the prepared molecule.
    pub fn energy(&self, coords: &[[f64; 3]]) -> f64 {
        self.bond_energy(coords)
            + self.angle_energy(coords)
            + self.stretch_bend_energy(coords)
            + self.torsion_energy(coords)
            + self.oop_energy(coords)
            + self.vdw_energy(coords)
            + self.electrostatic_energy(coords)
    }

    /// Evaluate all MMFF94 energy terms without rebuilding topology state.
    pub fn energy_breakdown(&self, coords: &[[f64; 3]]) -> EnergyBreakdown {
        let b = self.bond_energy(coords);
        let a = self.angle_energy(coords);
        let sb = self.stretch_bend_energy(coords);
        let t = self.torsion_energy(coords);
        let o = self.oop_energy(coords);
        let v = self.vdw_energy(coords);
        let e = self.electrostatic_energy(coords);
        EnergyBreakdown {
            bond: b,
            angle: a,
            stretch_bend: sb,
            torsion: t,
            oop: o,
            vdw: v,
            electrostatic: e,
            total: b + a + sb + t + o + v + e,
        }
    }

    /// Minimize coordinates using the prepared MMFF94 topology state.
    pub fn minimize_lbfgs(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
    ) -> Result<MinimizeResult, MinimizerError> {
        minimize_mmff94_lbfgs_prepared(self, coords, max_iter)
    }

    fn bond_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_bond_energy(coords, &self.bonds)
    }

    fn angle_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_angle_energy(coords, &self.angles)
    }

    fn stretch_bend_energy(&self, coords: &[[f64; 3]]) -> f64 {
        stretch_bend_energy(&self.mmff_mol, coords, &self.types, &self.rings)
    }

    fn torsion_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_torsion_energy(coords, &self.torsions)
    }

    fn oop_energy(&self, coords: &[[f64; 3]]) -> f64 {
        oop_energy(&self.mmff_mol, coords, &self.types)
    }

    fn vdw_energy(&self, coords: &[[f64; 3]]) -> f64 {
        vdw_energy_pairs(coords, &self.types, &self.vdw_pairs)
    }

    fn electrostatic_energy(&self, coords: &[[f64; 3]]) -> f64 {
        electrostatic_energy_pairs(coords, &self.charges, &self.electrostatic_pairs)
    }
}

/// Compute total MMFF94 energy for a given geometry (kcal/mol).
///
/// Includes bond, angle, torsion, vdW, and electrostatic terms.
/// Does not modify coordinates.
pub fn mmff94_total_energy(mol: &Molecule, coords: &[[f64; 3]]) -> Result<f64, MinimizerError> {
    Ok(Mmff94EnergyModel::new(mol)?.energy(coords))
}

/// Scan a torsion dihedral angle i-j-k-l from 0° to 360° in `steps` increments,
/// returning (angle_deg, energy_kcal) pairs. Coordinates are not modified.
///
/// At each step the dihedral is set by rotating atoms past `k` about the j-k bond.
pub fn mmff94_torsion_scan(
    mol: &Molecule,
    coords: &[[f64; 3]],
    atom_i: usize,
    atom_j: usize,
    atom_k: usize,
    atom_l: usize,
    steps: usize,
) -> Result<Vec<(f64, f64)>, MinimizerError> {
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let charges = mmff94_charges_numeric(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
    let ring_set = find_sssr(mol);
    let n = mol.atom_count();
    let steps = steps.max(2);

    let mut results = Vec::with_capacity(steps);

    // Collect atoms on the `l` side of the j-k bond (BFS from k, not crossing j)
    let moving_atoms: Vec<usize> = {
        let mut visited = vec![false; n];
        visited[atom_j] = true;
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(atom_k);
        visited[atom_k] = true;
        let mut group = Vec::new();
        while let Some(cur) = queue.pop_front() {
            group.push(cur);
            for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
                let nbi = nb.0 as usize;
                if !visited[nbi] {
                    visited[nbi] = true;
                    queue.push_back(nbi);
                }
            }
        }
        group
    };

    let mut work = coords.to_vec();

    // Rotate the moving group in `steps` increments of 360°/steps
    let step_rad = 2.0 * std::f64::consts::PI / steps as f64;

    for step in 0..steps {
        let angle_deg = step as f64 * 360.0 / steps as f64;

        if step > 0 {
            // Rotate moving_atoms by step_rad about the j→k axis
            let j = work[atom_j];
            let k = work[atom_k];
            let axis = {
                let d = [k[0] - j[0], k[1] - j[1], k[2] - j[2]];
                let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                if len < 1e-12 {
                    [1.0, 0.0, 0.0]
                } else {
                    [d[0] / len, d[1] / len, d[2] / len]
                }
            };
            let (sin_a, cos_a) = step_rad.sin_cos();
            for &ai in &moving_atoms {
                // Rodrigues' rotation about axis through j
                let p = [work[ai][0] - j[0], work[ai][1] - j[1], work[ai][2] - j[2]];
                let cross = [
                    axis[1] * p[2] - axis[2] * p[1],
                    axis[2] * p[0] - axis[0] * p[2],
                    axis[0] * p[1] - axis[1] * p[0],
                ];
                let dot = axis[0] * p[0] + axis[1] * p[1] + axis[2] * p[2];
                work[ai] = [
                    j[0] + cos_a * p[0] + sin_a * cross[0] + (1.0 - cos_a) * dot * axis[0],
                    j[1] + cos_a * p[1] + sin_a * cross[1] + (1.0 - cos_a) * dot * axis[1],
                    j[2] + cos_a * p[2] + sin_a * cross[2] + (1.0 - cos_a) * dot * axis[2],
                ];
                let _ = (atom_i, atom_l); // suppress unused warnings
            }
        }

        let energy = total_energy(&mmff_mol, &work, &types, &charges, ring_set.rings());
        results.push((angle_deg, energy));
    }

    Ok(results)
}

/// Compute per-term MMFF94 energy breakdown for a given geometry.
pub fn mmff94_energy_breakdown(
    mol: &Molecule,
    coords: &[[f64; 3]],
) -> Result<EnergyBreakdown, MinimizerError> {
    Ok(Mmff94EnergyModel::new(mol)?.energy_breakdown(coords))
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
    coords: &mut [[f64; 3]],
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

    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let charges = mmff94_charges_numeric(mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
    // Ring membership is a topology fact, not a geometry one: compute SSSR
    // once per minimization run rather than once per finite-difference probe
    // (`compute_gradient` alone calls `total_energy` ~6n times per step).
    let ring_set = find_sssr(mol);
    let rings = ring_set.rings();

    let n = mol.atom_count();
    let initial = coords.to_vec();
    let convergence = 1e-4_f64;
    let step_size = 0.05_f64;
    let delta = 1e-4_f64;

    let mut iters = 0usize;
    let mut converged = false;

    for _ in 0..max_iter {
        iters += 1;
        let grad = compute_gradient(&mmff_mol, coords, &types, &charges, rings, delta);
        let max_g = grad
            .iter()
            .flat_map(|v| v.iter())
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max);

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

    let energy = total_energy(&mmff_mol, coords, &types, &charges, rings);

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

/// Minimize molecular geometry using L-BFGS (limited-memory quasi-Newton).
///
/// Typically converges in 2–5× fewer iterations than steepest descent for
/// well-behaved energy surfaces. Falls back to a steepest-descent step when
/// the curvature condition `y·s > 0` is not satisfied.
///
/// Uses finite-difference gradients (δ=1e-4 Å) and backtracking Armijo line search.
pub fn minimize_mmff94_lbfgs(
    mol: &Molecule,
    coords: &mut [[f64; 3]],
    max_iter: usize,
) -> Result<MinimizeResult, MinimizerError> {
    let model = Mmff94EnergyModel::new(mol)?;
    minimize_mmff94_lbfgs_prepared(&model, coords, max_iter)
}

fn minimize_mmff94_lbfgs_prepared(
    model: &Mmff94EnergyModel,
    coords: &mut [[f64; 3]],
    max_iter: usize,
) -> Result<MinimizeResult, MinimizerError> {
    const M: usize = 5; // L-BFGS history size
    const DELTA: f64 = 1e-4; // finite-difference step (Å)
    const CONVERGENCE: f64 = 1e-4; // max |gradient| threshold
    const C_ARMIJO: f64 = 1e-4; // Armijo sufficient-decrease constant
    const TAU: f64 = 0.5; // Armijo backtracking factor

    if model.mmff_mol.atom_count() <= 1 {
        return Ok(MinimizeResult {
            energy: 0.0,
            rmsd: 0.0,
            converged: true,
            iterations: 0,
        });
    }

    let n = model.mmff_mol.atom_count();
    let initial = coords.to_vec();

    // Circular history buffer: (s_k = Δx, y_k = Δg, ρ_k = 1/(y·s))
    let mut history: LbfgsHistory = VecDeque::new();

    let mut g = compute_gradient_prepared(model, coords, DELTA);
    let mut f0 = model.energy(coords);

    let mut iters = 0usize;
    let mut converged = false;

    for _ in 0..max_iter {
        iters += 1;

        // Convergence check
        let max_g = g
            .iter()
            .flat_map(|v| v.iter())
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max);
        if max_g < CONVERGENCE {
            converged = true;
            break;
        }

        // Two-loop L-BFGS recursion → search direction p
        let p = lbfgs_direction(&g, &history);

        // Armijo backtracking line search along p
        let gp: f64 = g.iter().zip(p.iter()).map(|(gi, pi)| dot3(*gi, *pi)).sum();
        let mut alpha = 1.0_f64;
        let (new_coords, f_new) = loop {
            let trial: Vec<[f64; 3]> = coords
                .iter()
                .zip(p.iter())
                .map(|(c, pi)| {
                    [
                        c[0] + alpha * pi[0],
                        c[1] + alpha * pi[1],
                        c[2] + alpha * pi[2],
                    ]
                })
                .collect();
            let f_trial = model.energy(&trial);
            if f_trial <= f0 + C_ARMIJO * alpha * gp {
                break (trial, f_trial);
            }
            alpha *= TAU;
            if alpha < 1e-12 {
                // Line search failed — take a tiny steepest descent step
                let scale = 0.01 / max_g.max(1e-8);
                let trial: Vec<[f64; 3]> = coords
                    .iter()
                    .zip(g.iter())
                    .map(|(c, gi)| {
                        [
                            c[0] - scale * gi[0],
                            c[1] - scale * gi[1],
                            c[2] - scale * gi[2],
                        ]
                    })
                    .collect();
                let f_trial = model.energy(&trial);
                break (trial, f_trial);
            }
        };

        // Compute new gradient
        let g_new = compute_gradient_prepared(model, &new_coords, DELTA);
        // `f_new` is the accepted line-search energy above; do not evaluate
        // the same coordinates a second time after computing the gradient.

        // Compute s = x_new - x, y = g_new - g
        let s: Vec<[f64; 3]> = new_coords
            .iter()
            .zip(coords.iter())
            .map(|(xn, xo)| [xn[0] - xo[0], xn[1] - xo[1], xn[2] - xo[2]])
            .collect();
        let y: Vec<[f64; 3]> = g_new
            .iter()
            .zip(g.iter())
            .map(|(gn, go)| [gn[0] - go[0], gn[1] - go[1], gn[2] - go[2]])
            .collect();
        let ys: f64 = y.iter().zip(s.iter()).map(|(yi, si)| dot3(*yi, *si)).sum();

        // Only store if curvature condition holds
        if ys > 1e-10 {
            if history.len() >= M {
                history.pop_front();
            }
            history.push_back((s, y, 1.0 / ys));
        }

        coords.copy_from_slice(&new_coords);
        g = g_new;
        f0 = f_new;
    }

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
        energy: f0,
        rmsd,
        converged,
        iterations: iters,
    })
}

/// L-BFGS two-loop recursion: compute search direction p = -H_k × g.
fn lbfgs_direction(g: &[[f64; 3]], history: &LbfgsHistory) -> Vec<[f64; 3]> {
    let m = history.len();

    if m == 0 {
        // No history: steepest descent direction
        return g.iter().map(|gi| [-gi[0], -gi[1], -gi[2]]).collect();
    }

    let mut q: Vec<[f64; 3]> = g.to_vec();
    let mut alphas = vec![0.0_f64; m];

    // First loop (backward)
    for i in (0..m).rev() {
        let (s, y, rho) = &history[i];
        let sq: f64 = s.iter().zip(q.iter()).map(|(si, qi)| dot3(*si, *qi)).sum();
        alphas[i] = rho * sq;
        let a = alphas[i];
        for (qi, yi) in q.iter_mut().zip(y.iter()) {
            qi[0] -= a * yi[0];
            qi[1] -= a * yi[1];
            qi[2] -= a * yi[2];
        }
    }

    // Scale by γ = (s_{m-1}·y_{m-1}) / (y_{m-1}·y_{m-1})
    let (s_last, y_last, _) = &history[m - 1];
    let sy: f64 = s_last
        .iter()
        .zip(y_last.iter())
        .map(|(si, yi)| dot3(*si, *yi))
        .sum();
    let yy: f64 = y_last.iter().map(|yi| dot3(*yi, *yi)).sum();
    let gamma = if yy > 1e-20 { sy / yy } else { 1.0 };
    for qi in q.iter_mut() {
        qi[0] *= gamma;
        qi[1] *= gamma;
        qi[2] *= gamma;
    }

    // Second loop (forward)
    for i in 0..m {
        let (s, y, rho) = &history[i];
        let yr: f64 = y.iter().zip(q.iter()).map(|(yi, ri)| dot3(*yi, *ri)).sum();
        let beta = rho * yr;
        let diff = alphas[i] - beta;
        for (qi, si) in q.iter_mut().zip(s.iter()) {
            qi[0] += diff * si[0];
            qi[1] += diff * si[1];
            qi[2] += diff * si[2];
        }
    }

    // p = -H_k g = -q
    q.iter().map(|qi| [-qi[0], -qi[1], -qi[2]]).collect()
}

/// Compute finite-difference gradient: ∂E/∂x_i via central differences.
fn compute_gradient(
    mol: &Molecule,
    coords: &[[f64; 3]],
    types: &[u8],
    charges: &[f64],
    rings: &[Vec<AtomIdx>],
    delta: f64,
) -> Vec<[f64; 3]> {
    let n = coords.len();
    let mut grad = vec![[0.0_f64; 3]; n];
    let mut work = coords.to_vec();
    for i in 0..n {
        for axis in 0..3 {
            work[i][axis] += delta;
            let ep = total_energy(mol, &work, types, charges, rings);
            work[i][axis] -= 2.0 * delta;
            let em = total_energy(mol, &work, types, charges, rings);
            work[i][axis] += delta;
            grad[i][axis] = (ep - em) / (2.0 * delta);
        }
    }
    grad
}

fn compute_gradient_prepared(
    model: &Mmff94EnergyModel,
    coords: &[[f64; 3]],
    delta: f64,
) -> Vec<[f64; 3]> {
    let n = coords.len();
    let mut grad = vec![[0.0_f64; 3]; n];
    let mut work = coords.to_vec();
    for i in 0..n {
        for axis in 0..3 {
            work[i][axis] += delta;
            let ep = model.energy(&work);
            work[i][axis] -= 2.0 * delta;
            let em = model.energy(&work);
            work[i][axis] += delta;
            grad[i][axis] = (ep - em) / (2.0 * delta);
        }
    }
    grad
}

fn build_bond_terms(mol: &Molecule, types: &[u8]) -> Vec<PreparedBond> {
    mol.bonds()
        .filter_map(|(_, bond)| {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            let bt = bond_type_for(types[i], types[j], bond.order);
            mmff94_bond_energy_resolved(bt, types[i], types[j]).map(|(params, _)| PreparedBond {
                i,
                j,
                params,
            })
        })
        .collect()
}

fn build_angle_terms(mol: &Molecule, types: &[u8], rings: &[Vec<AtomIdx>]) -> Vec<PreparedAngle> {
    let mut terms = Vec::new();
    for j_idx in 0..mol.atom_count() {
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let Some((bond_ij, _)) = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx])
                else {
                    continue;
                };
                let Some((bond_kj, _)) = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx])
                else {
                    continue;
                };
                let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                if let Some((params, _)) = mmff94_angle_energy_resolved(
                    at,
                    types[i],
                    types[j_idx],
                    types[k],
                    bond_ij.r0,
                    bond_kj.r0,
                    ring_size,
                ) {
                    terms.push(PreparedAngle {
                        i,
                        j: j_idx,
                        k,
                        params,
                    });
                }
            }
        }
    }
    terms
}

fn build_torsion_terms(mol: &Molecule, types: &[u8]) -> Vec<PreparedTorsion> {
    let mut terms = Vec::new();
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        let nbrs_j: Vec<usize> = mol
            .neighbors(bond.atom1)
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        let nbrs_k: Vec<usize> = mol
            .neighbors(bond.atom2)
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        for &i in &nbrs_j {
            if i == k {
                continue;
            }
            for &l in &nbrs_k {
                if l == j {
                    continue;
                }
                let tt = torsion_type_for(mol, i, j, k, l, types[i], types[j], types[k], types[l]);
                if let Some(params) =
                    mmff94_torsion_energy(tt, types[i], types[j], types[k], types[l])
                {
                    terms.push(PreparedTorsion { i, j, k, l, params });
                }
            }
        }
    }
    terms
}

fn prepared_bond_energy(coords: &[[f64; 3]], terms: &[PreparedBond]) -> f64 {
    const KB_CONV: f64 = 143.9325;
    const CS: f64 = 2.0;
    terms
        .iter()
        .map(|term| {
            let dr = dist(coords[term.i], coords[term.j]) - term.params.r0;
            let cubic = 1.0 - CS * dr + (7.0 / 12.0) * CS * CS * dr * dr;
            (KB_CONV * term.params.kb / 2.0) * dr * dr * cubic
        })
        .sum()
}

fn prepared_angle_energy(coords: &[[f64; 3]], terms: &[PreparedAngle]) -> f64 {
    const KA_CONV: f64 = 0.043844;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    terms
        .iter()
        .map(|term| {
            let dt = cos_angle(coords[term.i], coords[term.j], coords[term.k]).acos() * RAD_TO_DEG
                - term.params.theta0;
            (KA_CONV * term.params.ka / 2.0) * dt * dt * (1.0 - 0.007 * dt)
        })
        .sum()
}

fn prepared_torsion_energy(coords: &[[f64; 3]], terms: &[PreparedTorsion]) -> f64 {
    terms
        .iter()
        .map(|term| {
            let phi = dihedral(
                coords[term.i],
                coords[term.j],
                coords[term.k],
                coords[term.l],
            );
            0.5 * term.params.v1 * (1.0 + phi.cos())
                + 0.5 * term.params.v2 * (1.0 - (2.0 * phi).cos())
                + 0.5 * term.params.v3 * (1.0 + (3.0 * phi).cos())
        })
        .sum()
}

fn build_nonbonded_pairs(mol: &Molecule) -> (VdwPairs, ElectrostaticPairs) {
    let n = mol.atom_count();
    let mut excluded = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        excluded.insert((i.min(j), i.max(j)));
        for (neighbor, _) in mol.neighbors(bond.atom1) {
            let k = neighbor.0 as usize;
            excluded.insert((k.min(j), k.max(j)));
        }
        for (neighbor, _) in mol.neighbors(bond.atom2) {
            let k = neighbor.0 as usize;
            excluded.insert((i.min(k), i.max(k)));
        }
    }

    let mut one_four = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let j = bond.atom1.0 as usize;
        let k = bond.atom2.0 as usize;
        for (neighbor_j, _) in mol.neighbors(bond.atom1) {
            let i = neighbor_j.0 as usize;
            if i == k {
                continue;
            }
            for (neighbor_k, _) in mol.neighbors(bond.atom2) {
                let l = neighbor_k.0 as usize;
                if l == j {
                    continue;
                }
                let pair = (i.min(l), i.max(l));
                if !excluded.contains(&pair) {
                    one_four.insert(pair);
                }
            }
        }
    }

    let mut vdw_pairs = Vec::new();
    let mut electrostatic_pairs = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if excluded.contains(&(i, j)) {
                continue;
            }
            vdw_pairs.push((i, j));
            let scale = if one_four.contains(&(i, j)) {
                0.75
            } else {
                1.0
            };
            electrostatic_pairs.push((i, j, scale));
        }
    }
    (vdw_pairs, electrostatic_pairs)
}

fn vdw_energy_pairs(coords: &[[f64; 3]], types: &[u8], pairs: &[(usize, usize)]) -> f64 {
    let mut energy = 0.0;
    for &(i, j) in pairs {
        let r = dist(coords[i], coords[j]);
        if r > 10.0 {
            continue;
        }
        if let Some((r_star, eps)) = mmff94_vdw_combined(types[i], types[j])
            && r_star > 0.0
            && eps > 0.0
            && r > 0.01
        {
            let t = (1.07 * r_star) / (r + 0.07 * r_star);
            let t7 = t.powi(7);
            energy += eps * t7 * (t7 - 2.0);
        }
    }
    energy
}

fn electrostatic_energy_pairs(
    coords: &[[f64; 3]],
    charges: &[f64],
    pairs: &[(usize, usize, f64)],
) -> f64 {
    const COULOMB: f64 = 332.0716;
    const DELTA: f64 = 0.05;
    pairs
        .iter()
        .map(|&(i, j, scale)| {
            scale * COULOMB * charges[i] * charges[j] / (dist(coords[i], coords[j]) + DELTA)
        })
        .sum()
}

// ─── Energy components ───────────────────────────────────────────────────────

fn total_energy(
    mol: &Molecule,
    coords: &[[f64; 3]],
    types: &[u8],
    charges: &[f64],
    rings: &[Vec<AtomIdx>],
) -> f64 {
    bond_energy(mol, coords, types)
        + angle_energy(mol, coords, types, rings)
        + stretch_bend_energy(mol, coords, types, rings)
        + torsion_energy(mol, coords, types)
        + oop_energy(mol, coords, types)
        + vdw_energy(mol, coords, types)
        + elec_energy(mol, coords, charges)
}

/// Stretch-bend coupling (Halgren MMFF.V eq. 4)
/// E_sb = 2.51210 × (kba_ijk × Δr_ij + kba_kji × Δr_kj) × Δθ   [kcal/mol, Δθ in degrees]
fn stretch_bend_energy(
    mol: &Molecule,
    coords: &[[f64; 3]],
    types: &[u8],
    rings: &[Vec<AtomIdx>],
) -> f64 {
    const CONV: f64 = 2.51210; // md/Å → kcal/(mol·Å·deg)
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    const KB_CONV: f64 = 143.9325;
    const CS: f64 = 2.0;
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() < 2 {
            continue;
        }
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let sbt = stretch_bend_type_for(at, types[i], types[k], bt_ij, bt_kj);
                if let Some((kba_ijk, kba_kji)) = mmff94_stbn(
                    sbt,
                    types[i],
                    types[j_idx],
                    types[k],
                    mol.atom(AtomIdx(i as u32)).element.atomic_number(),
                    mol.atom(AtomIdx(j_idx as u32)).element.atomic_number(),
                    mol.atom(AtomIdx(k as u32)).element.atomic_number(),
                ) {
                    let bond_ij = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx]);
                    let bond_kj = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx]);
                    // Δr_ij
                    let r_ij = dist(coords[i], coords[j_idx]);
                    let dr_ij = bond_ij.map(|(p, _)| r_ij - p.r0).unwrap_or(0.0);
                    // Δr_kj
                    let r_kj = dist(coords[k], coords[j_idx]);
                    let dr_kj = bond_kj.map(|(p, _)| r_kj - p.r0).unwrap_or(0.0);
                    // Δθ in degrees -- both flanking bonds must resolve to feed
                    // r0_ij/r0_jk into the angle empirical rule if it's needed
                    // (matches RDKit's real `getMMFFAngleBendParams`, which
                    // requires both `getMMFFBondStretchParams` calls to succeed
                    // before even attempting the empirical path).
                    let cos_t = cos_angle(coords[i], coords[j_idx], coords[k]);
                    if let (Some((rij, _)), Some((rkj, _))) = (bond_ij, bond_kj) {
                        let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                        if let Some((ap, _)) = mmff94_angle_energy_resolved(
                            at,
                            types[i],
                            types[j_idx],
                            types[k],
                            rij.r0,
                            rkj.r0,
                            ring_size,
                        ) {
                            let dtheta = cos_t.acos() * RAD_TO_DEG - ap.theta0;
                            energy += CONV * (kba_ijk * dr_ij + kba_kji * dr_kj) * dtheta;
                        }
                    }
                    let _ = (KB_CONV, CS); // suppress warnings
                }
            }
        }
    }
    energy
}

/// sp2 atom types eligible for out-of-plane bending (Halgren 1996), used by
/// [`oop_energy`] below. `pub` so downstream coverage-checkers (e.g.
/// `chematic-3d`'s MMFF94 bridge) can mirror exactly which atoms this
/// module's own energy loop would evaluate, instead of hand-copying this
/// list and risking drift.
pub const OOP_SP2_TYPES: &[u8] = &[
    2, 3, 9, 10, 30, 37, 38, 39, 40, 41, 43, 45, 49, 54, 56, 57, 58, 59, 63, 64, 65, 66, 67, 76,
    78, 79, 80, 81, 82,
];

/// Out-of-plane bending for trigonal sp2 centers (Halgren MMFF.VI eq. 6)
/// E_oop = (0.043844 × koop / 2) × χ²  (χ in degrees: Wilson angle of out-of-plane distortion)
fn oop_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8]) -> f64 {
    const CONV: f64 = 0.043844;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        if OOP_SP2_TYPES.binary_search(&types[j_idx]).is_err() {
            continue;
        }
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() != 3 {
            continue; // OOP only for exactly 3 substituents (trigonal)
        }
        let [i, k, l] = [neighbors[0], neighbors[1], neighbors[2]];
        if let Some(koop) = mmff94_oop(types[j_idx], types[i], types[k], types[l]) {
            // Wilson out-of-plane angle: angle between j→l vector and plane (i,j,k)
            let pj = coords[j_idx];
            let pi = coords[i];
            let pk = coords[k];
            let pl = coords[l];
            let rji = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
            let rjk = [pk[0] - pj[0], pk[1] - pj[1], pk[2] - pj[2]];
            let rjl = [pl[0] - pj[0], pl[1] - pj[1], pl[2] - pj[2]];
            let n = cross(rji, rjk); // normal to ijk plane
            let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            let l_len = (rjl[0] * rjl[0] + rjl[1] * rjl[1] + rjl[2] * rjl[2]).sqrt();
            if n_len < 1e-12 || l_len < 1e-12 {
                continue;
            }
            let sin_chi = dot3(n, rjl) / (n_len * l_len);
            let chi_deg = sin_chi.clamp(-1.0, 1.0).asin() * RAD_TO_DEG;
            energy += (CONV * koop / 2.0) * chi_deg * chi_deg;
        }
    }
    energy
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
        let bt = bond_type_for(types[i], types[j], bond.order);
        if let Some((p, _)) = mmff94_bond_energy_resolved(bt, types[i], types[j]) {
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
fn angle_energy(mol: &Molecule, coords: &[[f64; 3]], types: &[u8], rings: &[Vec<AtomIdx>]) -> f64 {
    const KA_CONV: f64 = 0.043844;
    const RAD_TO_DEG: f64 = 180.0 / std::f64::consts::PI;
    let mut energy = 0.0;
    for j_idx in 0..mol.atom_count() {
        let j = AtomIdx(j_idx as u32);
        let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
        if neighbors.len() < 2 {
            continue;
        }
        for (ii, &i) in neighbors.iter().enumerate() {
            for &k in &neighbors[ii + 1..] {
                let at = angle_type_for(mol, rings, i, j_idx, k, types);
                let bt_ij =
                    bond_type_for(types[i], types[j_idx], bond_order_between(mol, i, j_idx));
                let bt_kj =
                    bond_type_for(types[k], types[j_idx], bond_order_between(mol, k, j_idx));
                let r0_ij = mmff94_bond_energy_resolved(bt_ij, types[i], types[j_idx]);
                let r0_kj = mmff94_bond_energy_resolved(bt_kj, types[k], types[j_idx]);
                let (Some((bond_ij, _)), Some((bond_kj, _))) = (r0_ij, r0_kj) else {
                    continue;
                };
                let ring_size = is_angle_in_ring_of_size_3_or_4(mol, i, j_idx, k);
                if let Some((p, _)) = mmff94_angle_energy_resolved(
                    at,
                    types[i],
                    types[j_idx],
                    types[k],
                    bond_ij.r0,
                    bond_kj.r0,
                    ring_size,
                ) {
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
        let nbrs_j: Vec<usize> = mol
            .neighbors(bond.atom1)
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        let nbrs_k: Vec<usize> = mol
            .neighbors(bond.atom2)
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        for &i in &nbrs_j {
            if i == k {
                continue;
            }
            for &l in &nbrs_k {
                if l == j {
                    continue;
                }
                let tt = torsion_type_for(mol, i, j, k, l, types[i], types[j], types[k], types[l]);
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
            if let Some((r_star, eps)) = mmff94_vdw_combined(types[i], types[j])
                && r_star > 0.0
                && eps > 0.0
                && r > 0.01
            {
                let t = (1.07 * r_star) / (r + 0.07 * r_star);
                let t7 = t.powi(7);
                energy += eps * t7 * (t7 - 2.0);
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
            let scale = if one_four.contains(&(i, j)) {
                0.75
            } else {
                1.0
            };
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

/// Atom types with "extended" multiple-bond character (sp2/aromatic/sp
/// centers) — the MMFF94 property table's `mltb ≠ 0` set, i.e. atoms that
/// can be the single-bond end of an MMFF94 "single bond between
/// multiply-bonded atoms" (sbmb) pair (e.g. biphenyl's inter-ring bond, or
/// buta-1,3-diene's central C-C bond).
///
/// Exactly the 16 types that appear in `MMFF94_BOND_ENERGY`'s own
/// `(1, ti, tj, ...)` rows: {2,3,4,9,30,37,39,54,57,58,63,64,67,78,80,81}.
/// Measured, not assumed: a broader candidate (`oop_energy`'s trigonal-only
/// `SP2_TYPES` plus type 4) was tried first and made the 58-molecule corpus
/// harness's bond-miss count *worse* (38 vs. 16) — e.g. amide C(=O)-N
/// (types 3+10) and furan's aromatic C-O (types 37/38+43) both have no `(1,
/// ...)` row for that specific pair, so classifying them bt=1 under the
/// broader set was itself a fresh miss the tighter, evidence-derived set
/// avoids by routing them to their real `(0, ...)` row instead.
pub const MLTB_TYPES: &[u8] = &[2, 3, 4, 9, 30, 37, 39, 54, 57, 58, 63, 64, 67, 78, 80, 81];

/// Real bond order between two bonded atoms, defaulting to `Single` if no
/// bond exists between them (defensive only — every call site here passes
/// atom pairs already known to be bonded via `mol.neighbors`/`mol.bonds`).
fn bond_order_between(mol: &Molecule, a: usize, b: usize) -> BondOrder {
    mol.bond_between(AtomIdx(a as u32), AtomIdx(b as u32))
        .map(|(_, bond)| bond.order)
        .unwrap_or(BondOrder::Single)
}

/// True if ALL of `atoms` belong to a single common SSSR ring of exactly
/// `size` atoms (not just individually each in *some* ring of that size —
/// e.g. an exocyclic substituent atom on a ring must not count).
fn atoms_share_ring_of_size(rings: &[Vec<AtomIdx>], atoms: &[usize], size: usize) -> bool {
    rings
        .iter()
        .any(|ring| ring.len() == size && atoms.iter().all(|&a| ring.contains(&AtomIdx(a as u32))))
}

/// `isAngleInRingOfSize3or4` (`AtomTyper.cpp`, pinned commit): LOCAL bond
/// adjacency, NOT SSSR-based -- deliberately distinct from
/// [`atoms_share_ring_of_size`] above, which IS SSSR-based and feeds angle
/// *type* classification. RDKit's real empirical-rule ring gate and its real
/// angle-type ring gate are two different, independently-defined mechanisms
/// (confirmed by direct source read, issue #227 Stage C). Returns 3 if `i`
/// and `k` are directly bonded (i-j-k-i, a 3-ring), 4 if `i` and `k` share a
/// common neighbor other than `j` (a 4-ring through that neighbor), else 0.
/// `i`/`k` are assumed already bonded to `j` (true for every call site,
/// which only ever iterates real `mol.neighbors(j)` pairs).
///
/// `pub` so downstream coverage-checkers (e.g. `chematic-3d`'s independent
/// MMFF94 coverage measurement) can pass the same `ring_size` into
/// [`crate::mmff94_energy::mmff94_angle_energy_resolved`] this module's own
/// `angle_energy`/`stretch_bend_energy` use, instead of hand-copying this
/// logic and risking drift (same rationale as [`OOP_SP2_TYPES`]).
pub fn is_angle_in_ring_of_size_3_or_4(mol: &Molecule, i: usize, j: usize, k: usize) -> u8 {
    if mol
        .bond_between(AtomIdx(i as u32), AtomIdx(k as u32))
        .is_some()
    {
        return 3;
    }
    let j = j as u32;
    let i_neighbors: std::collections::HashSet<u32> = mol
        .neighbors(AtomIdx(i as u32))
        .map(|(nb, _)| nb.0)
        .filter(|&n| n != j)
        .collect();
    let shares_neighbor = mol
        .neighbors(AtomIdx(k as u32))
        .any(|(nb, _)| nb.0 != j && i_neighbors.contains(&nb.0));
    if shares_neighbor { 4 } else { 0 }
}

/// Determine the MMFF94 bond-type index (Halgren 1996 bond-type index BT).
///
/// BT=1 only for a formally SINGLE, *non-aromatic* bond directly linking two
/// atoms that are BOTH independently "conjugation-capable" ([`MLTB_TYPES`]) —
/// the sbmb special case (e.g. biphenyl's inter-ring bond, r0=1.436 Å).
/// `BondOrder::Aromatic` gets BT=0 unconditionally, same as a real
/// double/triple bond — confirmed against a live RDKit oracle (issue #227
/// Phase 1B-0): benzene's own ring bond (types 37-37, `BondOrder::Aromatic`)
/// resolves to RDKit's `(bondType=0, kb=5.573, r0=1.374)`, not `(1, 5.178,
/// 1.436)` — i.e. the aromatic ring bond and biphenyl's non-aromatic
/// inter-ring bond between the *same two atom types* get *different* BT
/// values, which only the bond order (not the MLTB_TYPES membership alone)
/// can distinguish. A prior version of this function treated `Aromatic`
/// bonds the same as `Single` here and was never checked against a real
/// oracle for that specific claim — the `(1, 63, 63, ...)` row it was
/// "confirmed" against belonged to the *pre-fix*, wrong (63, not 37) atom
/// typing, so the match was coincidental, not a correctness check. A real
/// double/triple bond always gets BT=0 too (its atom types already encode
/// the multiply-bonded hybridization; the `(1, ...)` table rows are the
/// single-bond exception, not an alternate double-bond row), and a bond
/// where either atom is a plain non-conjugating type (e.g. sp3 CR) also
/// always gets BT=0 — confirmed empirically: zero `(1, 1, x, ...)` rows
/// exist in the 493-row bond table for sp3 carbon type 1 with any partner.
pub fn bond_type_for(ti: u8, tj: u8, order: BondOrder) -> u8 {
    if matches!(
        order,
        BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple | BondOrder::Aromatic
    ) {
        return 0;
    }
    if MLTB_TYPES.binary_search(&ti).is_ok() && MLTB_TYPES.binary_search(&tj).is_ok() {
        1
    } else {
        0
    }
}

/// Determine the MMFF94 angle-type index (Halgren 1996, types 0-8).
///
/// `bt_sum` = sum of the two flanking bonds' [`bond_type_for`] indices (0-2).
/// Ring-embedded angles (i,j,k all in one common 3- or 4-membered SSSR ring)
/// get dedicated types; everything else uses `bt_sum` directly:
///
/// | ring   | bt_sum=0 | bt_sum=1 | bt_sum=2 |
/// |--------|----------|----------|----------|
/// | none   | 0        | 1        | 2        |
/// | 3-ring | 3        | 5        | 6        |
/// | 4-ring | 4        | 7        | 8        |
///
/// Matches RDKit's real `getMMFFAngleType` formula (`AtomTyper.cpp:2412-`
/// `2447`, pinned commit — see `scripts/mmff94_provenance/PROVENANCE.md`'s
/// "Stretch-bend" row, Priority 2C addendum): `angleType = ring_size; if
/// bond_type_sum != 0 { angleType += bond_type_sum + ring_size - 2 }`. A
/// prior version of this table gave 3-ring bt_sum=2 -> 8 and 4-ring
/// bt_sum=1 -> 6 / bt_sum=2 -> 7, which disagreed with RDKit's formula —
/// fixed here (issue #227). Measured LATENT on the 265-molecule Wave 1
/// corpus (0/113 reachable ring-embedded angle triples hit these branches),
/// but it is a real, independently-provable formula bug and a needed
/// correct input to [`crate::mmff94_energy::mmff94_stbn`]'s stretch-bend
/// type classification (a wrong angle_type feeds directly into
/// `getMMFFStretchBendType`'s first argument).
///
/// Confirmed empirically against the angle table: `(3, 22, 22, 22)`
/// (all-CR3R, cyclopropane) has θ0≈60°; `(4, 6, 20, 20)` (4-ring) has
/// θ0≈93° — i.e. 3-ring/4-ring are not swapped.
pub fn angle_type_for(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    i: usize,
    j: usize,
    k: usize,
    types: &[u8],
) -> u8 {
    let bt_ij = bond_type_for(types[i], types[j], bond_order_between(mol, i, j));
    let bt_jk = bond_type_for(types[j], types[k], bond_order_between(mol, j, k));
    let bt_sum = bt_ij + bt_jk;

    if atoms_share_ring_of_size(rings, &[i, j, k], 3) {
        return match bt_sum {
            0 => 3,
            1 => 5,
            _ => 6,
        };
    }
    if atoms_share_ring_of_size(rings, &[i, j, k], 4) {
        return match bt_sum {
            0 => 4,
            1 => 7,
            _ => 8,
        };
    }
    bt_sum
}

/// Determine the MMFF94 stretch-bend-type index (Halgren MMFF.V, types
/// 0-11) for angle i-j-k.
///
/// Ported verbatim from RDKit's `getMMFFStretchBendType`
/// (`AtomTyper.cpp:2480-2508`, pinned commit — see
/// `scripts/mmff94_provenance/PROVENANCE.md`'s "Stretch-bend" row,
/// Priority 2C addendum) and the diagnostic's own `resolve_rdkit` (issue
/// #227, `mmff94_stbn_equivalence_diagnostic_227.rs`). Angle types 1, 5, 7
/// each split into two distinct stretch-bend types depending on whether
/// either flanking bond individually has MMFF bond type 1 — [`angle_type_for`]'s
/// `bt_sum` only records the *sum*, discarding exactly the information this
/// split needs, which is why stretch-bend cannot reuse `angle_type` directly
/// as its own [`crate::mmff94_energy::MMFF94_STBN`] table key (issue #227
/// root cause).
///
/// `ta`/`tc` are the outer atom MMFF types (i and k); `bond_type_ij`/
/// `bond_type_jk` are the individual (unswapped) [`bond_type_for`] results
/// for the i-j and j-k bonds. The two bond-type arguments fed into the
/// classification match are canonicalized on `ta <= tc` — this is **not**
/// the same swap rule as [`mmff94_stbn_type_only`](crate::mmff94_energy::mmff94_stbn_type_only)'s
/// own i/k table-lookup canonicalization (that one additionally tie-breaks
/// on bond type when `ta == tc`; see `AtomTyper.cpp:3598-3600`). Do not
/// conflate the two swap rules.
pub fn stretch_bend_type_for(
    angle_type: u8,
    ta: u8,
    tc: u8,
    bond_type_ij: u8,
    bond_type_jk: u8,
) -> u8 {
    let (arg1, arg2) = if ta <= tc {
        (
            bond_type_ij,
            if ta < tc { bond_type_jk } else { bond_type_ij },
        )
    } else {
        (bond_type_jk, bond_type_ij)
    };
    match angle_type {
        1 => {
            if arg1 != 0 || arg1 == arg2 {
                1
            } else {
                2
            }
        }
        2 => 3,
        4 => 4,
        3 => 5,
        5 => {
            if arg1 != 0 || arg1 == arg2 {
                6
            } else {
                7
            }
        }
        6 => 8,
        7 => {
            if arg1 != 0 || arg1 == arg2 {
                9
            } else {
                10
            }
        }
        8 => 11,
        _ => 0,
    }
}

/// Port of RDKit's `isTorsionInRingOfSize4or5` (`AtomTyper.cpp:403-447`,
/// pinned commit — see `scripts/mmff94_provenance/PROVENANCE.md`'s
/// "Torsion" row) and the diagnostic's own `rdkit_ring_size_4_or_5` (issue
/// #227, `mmff94_torsion_equivalence_diagnostic_227.rs`). Purely LOCAL
/// bond-adjacency, NOT SSSR-based: 4-ring iff i-l are directly bonded;
/// 5-ring iff i and l, excluding their ring neighbours j and k
/// respectively, share a common neighbour. Returns 0 if neither.
fn ring_size_4_or_5(mol: &Molecule, i: AtomIdx, j: AtomIdx, k: AtomIdx, l: AtomIdx) -> u8 {
    if mol.bond_between(i, l).is_some() {
        return 4;
    }
    let nbrs_i: Vec<AtomIdx> = mol
        .neighbors(i)
        .map(|(n, _)| n)
        .filter(|&n| n != j)
        .collect();
    let has_common = mol
        .neighbors(l)
        .map(|(n, _)| n)
        .filter(|&n| n != k)
        .any(|n| nbrs_i.contains(&n));
    if has_common { 5 } else { 0 }
}

/// True when RDKit's real torsion resolution generates NO term at all for
/// this i-j-k-l torsion, by design -- not a missing-parameter gap.
///
/// Halgren's empirical-rule cascade (`getMMFFTorsionEmpiricalRuleParams`,
/// `AtomTyper.cpp`, rule (a) per the public transcription cited in
/// `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry) omits the
/// torsion term entirely whenever either central atom (`type_j`/`type_k`,
/// the j-k bond this torsion rotates around) has MMFF's `lin` flag
/// ([`crate::mmff94_numeric_type_registry::Mmff94NumericTypeInfo::linear`],
/// e.g. type 4 CSP/nitrile-or-acetylenic carbon, type 53 `=N=`/cumulated
/// azide nitrogen): rotating around a bond whose other end is a linear
/// (180°) center changes no real geometry, so there is nothing to
/// parameterize. Issue #227 Phase 1: measured as the exact, complete
/// explanation for the only 3 genuine `table_gap` Torsion instances left in
/// the 265-molecule Wave 1 corpus after the classification fix below —
/// oracle-confirmed (`GetMMFFTorsionParams` returns `None` for all 3, and
/// each central atom's registered `linear` flag matches RDKit's own MMFF
/// atom type exactly at that atom). A caller that already omits the term for
/// an unresolved lookup (chematic-ff's own `torsion_energy`, which just adds
/// nothing on `mmff94_torsion_energy(..) == None`) needs no code change for
/// physics; this exists so coverage-reporting/diagnostic callers (the
/// `Mmff94BondAngleStrict` gate under `include_torsion_oop_in_gate`, the
/// Phase 1A audit) can tell "correctly no term" apart from "genuinely
/// missing" instead of counting both as the same failure.
pub fn torsion_no_term_by_design(type_j: u8, type_k: u8) -> bool {
    let linear = |t: u8| {
        crate::mmff94_numeric_type_registry::mmff94_numeric_type_info(t)
            .is_some_and(|info| info.linear)
    };
    linear(type_j) || linear(type_k)
}

/// Determine the MMFF94 torsion-type index (Halgren 1996, types 0-8).
///
/// Ported verbatim from RDKit's `getMMFFTorsionType` (`AtomTyper.cpp:2528-`
/// `2571`, pinned commit — see `scripts/mmff94_provenance/PROVENANCE.md`'s
/// "Torsion" row) and the diagnostic's own `rdkit_torsion_type` (issue #227,
/// `mmff94_torsion_equivalence_diagnostic_227.rs`).
///
/// The base (non-ring) code is the j-k bond's OWN [`bond_type_for`] result —
/// NOT atom-type ([`MLTB_TYPES`]) membership, which is what a prior version
/// of this function used and which disagreed with RDKit's real
/// classification on 100% of the diagnostic's 1,107 routing-bug candidates
/// and 76.3% of ALL 13,530 torsion instances in the 265-molecule corpus (a
/// double/triple/aromatic j-k bond always gets `bond_type_jk=0` from
/// `bond_type_for` regardless of `tj`/`tk`'s own MLTB membership — the prior
/// atom-type-membership formula couldn't see that at all). An
/// empirically-required override bumps the base code to 2 when
/// `bond_type_jk==0 && order_jk==Single && (bond_type_ij==1 ||
/// bond_type_kl==1)` — MMFF.IV page 609's simple condition fails CYGUAN01 in
/// RDKit's own test suite; this corrected condition is RDKit's real code,
/// not derivable from Halgren's paper alone.
///
/// Ring-embedded torsions get dedicated types 4 (4-ring) / 5 (5-ring),
/// determined by [`ring_size_4_or_5`] — a purely LOCAL bond-adjacency check,
/// NOT SSSR-based, replacing this function's prior
/// [`atoms_share_ring_of_size`]-based check entirely (the diagnostic's whole
/// point is that the SSSR-based check was the wrong mechanism here). The
/// 4-ring override additionally requires i-k and j-l to NOT be directly
/// bonded (excludes degenerate/bridged cases); the 5-ring override
/// additionally requires at least one of the four atoms to be MMFF numeric
/// type 1 (`ti==1 || tj==1 || tk==1 || tl==1`) — a condition the prior
/// SSSR-based check had no equivalent of at all.
#[allow(clippy::too_many_arguments)]
pub fn torsion_type_for(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    ti: u8,
    tj: u8,
    tk: u8,
    tl: u8,
) -> u8 {
    let order_ij = bond_order_between(mol, i, j);
    let order_jk = bond_order_between(mol, j, k);
    let order_kl = bond_order_between(mol, k, l);
    let bond_type_ij = bond_type_for(ti, tj, order_ij);
    let bond_type_jk = bond_type_for(tj, tk, order_jk);
    let bond_type_kl = bond_type_for(tk, tl, order_kl);

    let mut torsion_type = bond_type_jk;
    if bond_type_jk == 0
        && order_jk == BondOrder::Single
        && (bond_type_ij == 1 || bond_type_kl == 1)
    {
        torsion_type = 2;
    }

    let (ai, aj, ak, al) = (
        AtomIdx(i as u32),
        AtomIdx(j as u32),
        AtomIdx(k as u32),
        AtomIdx(l as u32),
    );
    let ring_size = ring_size_4_or_5(mol, ai, aj, ak, al);
    if ring_size == 4 && mol.bond_between(ai, ak).is_none() && mol.bond_between(aj, al).is_none() {
        torsion_type = 4;
    } else if ring_size == 5 && (ti == 1 || tj == 1 || tk == 1 || tl == 1) {
        torsion_type = 5;
    }

    torsion_type
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff94_energy::mmff94_stbn_type_only;
    use crate::mmff94_energy::{mmff94_angle_energy, mmff94_bond_energy};
    use crate::mmff94_numeric::{
        assign_mmff94_numeric_types, assign_mmff94_numeric_types_with_view,
    };
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

    #[test]
    fn prepared_energy_model_matches_one_shot_api() {
        let (mol, coords) = methane_mol();
        let model = Mmff94EnergyModel::new(&mol).expect("prepare MMFF94 model");
        let one_shot = mmff94_energy_breakdown(&mol, &coords).expect("one-shot energy");
        let prepared = model.energy_breakdown(&coords);
        assert!((prepared.total - one_shot.total).abs() < 1e-12);
        assert!((model.energy(&coords) - one_shot.total).abs() < 1e-12);
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

    /// First SSSR ring of exactly `size` atoms, in ring (bonded-consecutive)
    /// order, or `None` if the molecule has no ring of that size.
    fn first_ring_of_size(mol: &Molecule, size: usize) -> Option<Vec<usize>> {
        let ring_set = find_sssr(mol);
        ring_set
            .rings()
            .iter()
            .find(|r| r.len() == size)
            .map(|r| r.iter().map(|a| a.0 as usize).collect())
    }

    // ── FF-1: bond_type_for (#173) ─────────────────────────────────────────

    #[test]
    fn bond_type_for_ethene_double_bond_is_type_0_not_1() {
        // Issue #173's primary repro: a real C=C double bond must never be
        // routed to the sbmb bond_type=1 row (r0=1.430, wrong) instead of the
        // correct bond_type=0 row (r0=1.333).
        let mol = chematic_smiles::parse("C=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let (_, bond) = mol.bonds().next().unwrap();
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        assert_eq!(types[i], 2, "vinylic carbon should be type 2");
        assert_eq!(types[j], 2);
        let bt = bond_type_for(types[i], types[j], bond.order);
        assert_eq!(bt, 0, "C=C double bond must classify as bond_type 0");
        let p = mmff94_bond_energy(bt, types[i], types[j]).expect("(0,2,2) row must exist");
        assert!((p.r0 - 1.333).abs() < 1e-6, "r0={} should be 1.333", p.r0);
    }

    #[test]
    fn bond_type_for_acetaldehyde_sp3_carbonyl_single_bond_is_type_0() {
        // Issue #173's secondary repro: CC=O's Cα(sp3)-C(carbonyl) single
        // bond must resolve to the real (0,1,3) row, not silently miss via a
        // wrongly-routed bond_type=1 lookup (which has no (1,1,3) row at all).
        let mol = chematic_smiles::parse("CC=O").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        // atom 0 = CH3 (type 1), atom 1 = C=O carbon (type 3)
        assert_eq!(types[0], 1);
        assert_eq!(types[1], 3);
        let order = bond_order_between(&mol, 0, 1);
        let bt = bond_type_for(types[0], types[1], order);
        assert_eq!(
            bt, 0,
            "sp3-carbonyl single bond must classify as bond_type 0"
        );
        let p = mmff94_bond_energy(bt, types[0], types[1]).expect("(0,1,3) row must exist");
        assert!((p.r0 - 1.492).abs() < 1e-6, "r0={} should be 1.492", p.r0);
    }

    #[test]
    fn ethene_energy_minimum_is_near_1_333_angstrom() {
        // Energy-perturbation reproduction of issue #173's own measurement:
        // mmff94_total_energy("C=C")'s minimum used to sit at r≈1.430 Å; it
        // must now sit near the chemically correct 1.333 Å.
        let mol = chematic_smiles::parse("C=C").unwrap();
        let scan = |r: f64| {
            let coords = vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
            mmff94_total_energy(&mol, &coords).unwrap()
        };
        let e_correct = scan(1.333);
        let e_old_wrong_minimum = scan(1.430);
        let e_far = scan(1.6);
        assert!(
            e_correct < e_old_wrong_minimum,
            "energy at true r0=1.333 ({e_correct}) should be lower than at the old wrong \
             minimum 1.430 ({e_old_wrong_minimum})"
        );
        assert!(
            e_correct < e_far,
            "energy at 1.333 ({e_correct}) should be lower than further out at 1.6 ({e_far})"
        );
    }

    #[test]
    fn acetaldehyde_bond_stretch_now_changes_energy() {
        // Issue #173's exact repro: stretching the Cα-C(carbonyl) bond by a
        // full 1 Å used to produce dE = 0.0 (silent missing parameter).
        let mol = chematic_smiles::parse("CC=O").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        // Build coords: C0 at origin, C1 (carbonyl C) along +x, O at a fixed
        // offset from C1 so the C=O bond itself is untouched by the scan.
        let coords_at = |r: f64| vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0], [r + 1.2, 0.9, 0.0]];
        let e_short = bond_energy(&mol, &coords_at(1.5), &types);
        let e_long = bond_energy(&mol, &coords_at(2.5), &types);
        assert!(
            (e_long - e_short).abs() > 1e-6,
            "bond stretch must now change bond energy: short={e_short} long={e_long} (was \
             identical pre-fix per issue #173)"
        );
    }

    #[test]
    fn bond_gradient_restores_ethene_toward_1_333() {
        // Gradient test using this file's own existing finite-difference
        // pattern (`compute_gradient`, already used by both minimizers).
        let mol = chematic_smiles::parse("C=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap_or_else(|_| vec![0.0; 2]);
        let rings = find_sssr(&mol);
        // Stretched well past the true r0=1.333 minimum.
        let coords = vec![[0.0, 0.0, 0.0_f64], [1.6, 0.0, 0.0]];
        let grad = compute_gradient(&mol, &coords, &types, &charges, rings.rings(), 1e-4);
        // Restoring force on atom 1 should point back toward atom 0 (-x),
        // i.e. dE/dx > 0 at this atom (energy increases as x increases further).
        assert!(
            grad[1][0] > 0.0,
            "gradient on stretched atom should point back toward equilibrium: grad_x={}",
            grad[1][0]
        );
    }

    #[test]
    fn benzene_ring_bonds_type_37_resolve_to_the_rdkit_verified_row() {
        // Issue #227 Phase 1B-0: atom typing fixed to the real Halgren/RDKit
        // numbering (type 37 = CB, not the old wrong 63 = C5A), AND
        // `bond_type_for` fixed to treat `BondOrder::Aromatic` like
        // double/triple (BT=0 unconditionally), not like `Single` run
        // through the MLTB_TYPES-AND check. Both fixes were required
        // together -- fixing only the atom typer would have made every
        // benzene ring bond resolve to the WRONG (1, 37, 37, kb=5.178,
        // r0=1.436) row (biphenyl's *non-aromatic* inter-ring bond value)
        // instead of the correct (0, 37, 37, kb=5.573, r0=1.374) row --
        // still "coverage success," still silently wrong, the same failure
        // shape as the `furan` finding this whole PR exists to close.
        // Cross-checked directly against a live RDKit oracle
        // (`props.GetMMFFBondStretchParams`), not assumed:
        // `AllChem.MMFFGetMoleculeProperties` on benzene returns exactly
        // `(bondType=0, kb=5.573, r0=1.374)` for every ring bond.
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        assert_eq!(mol.bonds().count(), 6);
        for (_, bond) in mol.bonds() {
            let ti = types[bond.atom1.0 as usize];
            let tj = types[bond.atom2.0 as usize];
            assert_eq!(
                ti, 37,
                "benzene carbon must type as 37 (CB), matching RDKit"
            );
            assert_eq!(tj, 37);
            assert_eq!(bond.order, BondOrder::Aromatic);
            let bt = bond_type_for(ti, tj, bond.order);
            assert_eq!(
                bt, 0,
                "aromatic ring bond between two type-37 atoms must resolve to bond_type=0, \
                 matching RDKit's GetMMFFBondStretchParams for benzene"
            );
            let params = mmff94_bond_energy(bt, ti, tj)
                .expect("every benzene ring bond must resolve to a real row, not silently miss");
            assert!(
                (params.r0 - 1.374).abs() < 1e-6,
                "benzene bond r0 must match RDKit's 1.374 A, got {}",
                params.r0
            );
            assert!(
                (params.kb - 5.573).abs() < 1e-6,
                "benzene bond kb must match RDKit's 5.573, got {}",
                params.kb
            );
        }
    }

    #[test]
    fn furan_c_c_bond_no_longer_collides_with_a_nitrogen_row() {
        // The exact regression this whole Phase 1B-0 PR exists to close
        // (issue #227's audit PR #235): before the fix, furan's ring
        // carbons were mistyped as 38/37 (real 38 = NPYD, pyridine-type
        // NITROGEN) and its C-C bond silently resolved to
        // `Some(BondEnergyParams { kb: 5.002, r0: 1.246 })` -- a real table
        // row, but one that belongs to a nitrogen-involved bond, not a
        // furan C-C bond. Pinning the old wrong value here so a future
        // reader can see exactly what this test prevents. After the fix,
        // furan's ring carbons must type as 63/64 (C5A/C5B, both carbon),
        // and their mutual bond must resolve to a genuinely carbon-carbon
        // row (or, honestly, `None` -- anything but a repeat of the old
        // nitrogen-row collision).
        let mol = chematic_smiles::parse("c1ccoc1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        for (_, bond) in mol.bonds() {
            let a1 = mol.atom(bond.atom1);
            let a2 = mol.atom(bond.atom2);
            if a1.element != chematic_core::Element::C || a2.element != chematic_core::Element::C {
                continue;
            }
            let ti = types[bond.atom1.0 as usize];
            let tj = types[bond.atom2.0 as usize];
            assert_ne!(
                (ti, tj),
                (38, 37),
                "must not regress to the old wrong furan C-C typing"
            );
            if let Some(params) = mmff94_bond_energy(bond_type_for(ti, tj, bond.order), ti, tj) {
                assert!(
                    (params.r0 - 1.246).abs() > 1e-6,
                    "furan C-C bond must not resolve to the old nitrogen-row r0=1.246, \
                     types=({ti},{tj})"
                );
            }
        }
    }

    // ── FF-1: torsion_type_for ring types 4/5 (#175) ───────────────────────

    #[test]
    fn torsion_type_for_cyclopentane_ring_is_type_5() {
        let mol = chematic_smiles::parse("C1CCCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 5).expect("cyclopentane must have a 5-ring");
        let (i, j, k, l) = (ring[0], ring[1], ring[2], ring[3]);
        let tt = torsion_type_for(&mol, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_eq!(
            tt, 5,
            "in-ring torsion of a 5-membered ring must classify as type 5"
        );
    }

    #[test]
    fn torsion_type_for_cyclobutane_ring_is_type_4() {
        let mol = chematic_smiles::parse("C1CCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("cyclobutane must have a 4-ring");
        let (i, j, k, l) = (ring[0], ring[1], ring[2], ring[3]);
        let tt = torsion_type_for(&mol, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_eq!(
            tt, 4,
            "in-ring torsion of a 4-membered ring must classify as type 4"
        );
    }

    #[test]
    fn torsion_type_for_ring5_override_requires_a_type_1_atom() {
        // RDKit's real 5-ring override additionally requires at least one of
        // the four torsion atoms to be MMFF numeric type 1 (`ti==1 || tj==1
        // || tk==1 || tl==1`) -- a condition the prior SSSR-based ring check
        // (which this fix replaced entirely) had no equivalent of at all.
        // Same 5-ring topology as `torsion_type_for_cyclopentane_ring_is_type_5`
        // above (where all four real chematic-assigned types happen to be 1,
        // so the condition is trivially satisfied there) but with the four
        // types forced away from 1 here -- purely local ring-adjacency
        // geometry (`ring_size_4_or_5` returning 5) must NOT be sufficient on
        // its own.
        let mol = chematic_smiles::parse("C1CCCC1").unwrap();
        let ring = first_ring_of_size(&mol, 5).expect("cyclopentane must have a 5-ring");
        let (i, j, k, l) = (ring[0], ring[1], ring[2], ring[3]);
        let tt = torsion_type_for(&mol, i, j, k, l, 20, 20, 20, 20);
        assert_ne!(
            tt, 5,
            "5-ring override must not fire when none of the four atoms is MMFF type 1"
        );
    }

    // ── issue #227: torsion_type_for's corrected bond-type-based base
    // classification and empirically-required type-2 override ─────────────

    #[test]
    fn torsion_type_for_aromatic_jk_bond_is_type_0_despite_mltb_membership() {
        // The core of issue #227's torsion bug: a prior version of this
        // function classified the non-ring base case purely from atom-type
        // MLTB_TYPES membership, `(MLTB(tj), MLTB(tk)) -> 0/1/2`, completely
        // ignoring the j-k bond's own real MMFF bond order/type. Benzene's
        // ring carbons are all type 37 (in MLTB_TYPES), so the prior formula
        // would have classified any benzene torsion as type 2 (both
        // "sp2-sp2") -- RDKit's real `getMMFFTorsionType` instead uses
        // `bond_type_for`'s own result as the base code, and an aromatic
        // bond always gets bond_type 0 (same as a real double/triple bond),
        // regardless of the endpoints' own MLTB membership.
        let mol = chematic_smiles::parse("c1ccccc1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        assert_eq!(types[0], 37, "benzene carbon should be type 37");
        assert!(
            MLTB_TYPES.binary_search(&37).is_ok(),
            "type 37 is in MLTB_TYPES -- the prior atom-type-membership formula would wrongly see this as conjugated"
        );
        let tt = torsion_type_for(&mol, 0, 1, 2, 3, types[0], types[1], types[2], types[3]);
        assert_eq!(
            tt, 0,
            "an aromatic j-k bond must classify as torsion type 0 via bond_type_for, not type 2 via MLTB membership"
        );
    }

    #[test]
    fn torsion_type_for_type_2_override_fires_when_flanking_bond_is_sbmb() {
        // Synthetic 4-atom chain i-j-k-l, all single bonds (topology only
        // matters for bond order / ring adjacency here -- the four MMFF
        // types are passed explicitly to isolate the classification
        // formula). ti=tj=2 (both in MLTB_TYPES) makes the i-j bond an sbmb
        // bond (`bond_type_ij` = 1, a single bond directly linking two
        // independently conjugation-capable atoms). tk=tl=1 (not in
        // MLTB_TYPES) makes the j-k bond's own base code 0 (neither a
        // multiple bond nor doubly-conjugated). RDKit's empirically-required
        // override (`bond_type_jk==0 && order_jk==Single && (bond_type_ij==1
        // || bond_type_kl==1)`, needed to pass RDKit's own CYGUAN01
        // regression test, not derivable from Halgren's MMFF.IV page 609
        // formula alone) must bump this to type 2.
        let mut b = MoleculeBuilder::new();
        let i = b.add_atom(Atom::new(Element::C));
        let j = b.add_atom(Atom::new(Element::C));
        let k = b.add_atom(Atom::new(Element::C));
        let l = b.add_atom(Atom::new(Element::C));
        b.add_bond(i, j, BondOrder::Single).unwrap();
        b.add_bond(j, k, BondOrder::Single).unwrap();
        b.add_bond(k, l, BondOrder::Single).unwrap();
        let mol = b.build();
        let tt = torsion_type_for(
            &mol,
            i.0 as usize,
            j.0 as usize,
            k.0 as usize,
            l.0 as usize,
            2,
            2,
            1,
            1,
        );
        assert_eq!(
            tt, 2,
            "base type 0 at j-k with a flanking sbmb (bond_type=1) i-j bond must override to type 2"
        );
    }

    #[test]
    fn torsion_type_for_no_override_when_neither_flank_is_sbmb() {
        // Same topology as the override test above, but with NEITHER flank
        // an sbmb bond (ti=tl=1, not in MLTB_TYPES) -- the override's
        // `bond_type_ij==1 || bond_type_kl==1` condition must not be
        // satisfied trivially just because bond_type_jk==0.
        let mut b = MoleculeBuilder::new();
        let i = b.add_atom(Atom::new(Element::C));
        let j = b.add_atom(Atom::new(Element::C));
        let k = b.add_atom(Atom::new(Element::C));
        let l = b.add_atom(Atom::new(Element::C));
        b.add_bond(i, j, BondOrder::Single).unwrap();
        b.add_bond(j, k, BondOrder::Single).unwrap();
        b.add_bond(k, l, BondOrder::Single).unwrap();
        let mol = b.build();
        let tt = torsion_type_for(
            &mol,
            i.0 as usize,
            j.0 as usize,
            k.0 as usize,
            l.0 as usize,
            1,
            2,
            1,
            1,
        );
        assert_eq!(
            tt, 0,
            "base type 0 must stand when neither flanking bond is an sbmb (bond_type=1) bond"
        );
    }

    #[test]
    fn torsion_type_for_exocyclic_substituent_is_not_ring_typed() {
        // Regression guard for the advisor-flagged failure mode: a torsion
        // whose central j-k bond is a ring bond but whose *terminal* atom is
        // an exocyclic substituent (methylcyclopentane's exocyclic methyl)
        // must NOT be misclassified as a ring torsion.
        let mol = chematic_smiles::parse("CC1CCCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        // atom 0 = exocyclic methyl C, atom 1 = ring C bonded to atom 0.
        let ring = first_ring_of_size(&mol, 5).expect("methylcyclopentane must have a 5-ring");
        assert!(
            ring.contains(&1),
            "atom 1 should be the substituted ring atom"
        );
        let j = 1usize;
        let k = *ring
            .iter()
            .find(|&&a| a != 1 && mol_bonded(&mol, j, a))
            .unwrap();
        let i = 0usize; // exocyclic methyl carbon, NOT a ring member
        let l = *ring
            .iter()
            .find(|&&a| a != j && a != k && mol_bonded(&mol, k, a))
            .unwrap();
        let tt = torsion_type_for(&mol, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_ne!(
            tt, 5,
            "torsion with an exocyclic terminal atom must not be classified as a ring torsion"
        );
    }

    fn mol_bonded(mol: &Molecule, a: usize, b: usize) -> bool {
        mol.bond_between(AtomIdx(a as u32), AtomIdx(b as u32))
            .is_some()
    }

    // ── torsion_no_term_by_design / NoTermByDesign (issue #227 Phase 1) ─────

    #[test]
    fn torsion_no_term_by_design_true_for_linear_central_atom_either_side() {
        // Type 4 (CSP, acetylenic/nitrile C) and 53 (=N=, cumulated
        // azide/diazo N) are the two `linear: true` types this corpus
        // exercises; type 61 is also linear per the registry but unreachable
        // in this corpus.
        assert!(torsion_no_term_by_design(4, 37), "linear on the j side");
        assert!(torsion_no_term_by_design(37, 53), "linear on the k side");
        assert!(torsion_no_term_by_design(4, 53), "linear on both sides");
    }

    #[test]
    fn torsion_no_term_by_design_false_when_neither_central_atom_is_linear() {
        assert!(!torsion_no_term_by_design(1, 1)); // sp3 C - sp3 C
        assert!(!torsion_no_term_by_design(37, 3)); // aromatic C - C=O (caffeine's own case)
    }

    #[test]
    fn torsion_no_term_by_design_unknown_type_fails_closed_to_false() {
        // A type absent from the registry must never be treated as linear
        // (would incorrectly suppress a real coverage-gap report).
        assert!(!torsion_no_term_by_design(250, 1));
    }

    #[test]
    fn nitrile_torsion_has_no_table_row_and_is_flagged_no_term_by_design() {
        // The exact chembl_tier_b_0001 shape (issue #227's 2 real
        // table_gap instances): Ar-Ar-C#N, central k = type 4 (CSP). Oracle-
        // confirmed RDKit's own GetMMFFTorsionParams also returns None here
        // (rdkit==2026.03.4) -- this is genuinely absent from the table on
        // both sides, not a routing bug.
        assert!(mmff94_torsion_energy(1, 37, 37, 4, 42).is_none());
        assert!(torsion_no_term_by_design(37, 4));
    }

    #[test]
    fn cumulated_azide_torsion_has_no_table_row_and_is_flagged_no_term_by_design() {
        // chembl_tier_b_0080's real shape: Ar-N=[N+]=[N-], central k = type
        // 53 (=N=). Oracle-confirmed None as above.
        assert!(mmff94_torsion_energy(0, 37, 9, 53, 47).is_none());
        assert!(torsion_no_term_by_design(9, 53));
    }

    #[test]
    fn caffeine_dione_ring_torsion_resolves_end_to_end_via_the_reperceived_view() {
        // Issue #227 Phase 1's actual root-cause fix, exercised through the
        // real production call shape (assign_mmff94_numeric_types_with_view
        // + torsion_type_for + mmff94_torsion_energy on the SAME reperceived
        // molecule), not just the diagnostic tool. Expected value is the
        // live RDKit oracle's own GetMMFFTorsionParams(mol,4,5,6,7) result
        // for this exact molecule/atom-quadruple (rdkit==2026.03.4):
        // (torsionType=1, 0.0, 2.5, 0.0).
        let m = chematic_smiles::parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        let (types, view) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        let (i, j, k, l) = (4, 5, 6, 7);
        let tt = torsion_type_for(&view, i, j, k, l, types[i], types[j], types[k], types[l]);
        assert_eq!(tt, 1, "must classify as type 1, not type 0");
        let params = mmff94_torsion_energy(tt, types[i], types[j], types[k], types[l])
            .expect("a real table row must resolve once fed the correct bond order");
        assert!((params.v1 - 0.0).abs() < 1e-9);
        assert!((params.v2 - 2.5).abs() < 1e-9, "v2={}", params.v2);
        assert!((params.v3 - 0.0).abs() < 1e-9);
    }

    #[test]
    fn torsion_lookup_is_symmetric_under_reversal_for_the_caffeine_case() {
        let m = chematic_smiles::parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C").unwrap();
        let (types, view) = assign_mmff94_numeric_types_with_view(&m).unwrap();
        let (i, j, k, l) = (4usize, 5usize, 6usize, 7usize);
        let fwd = torsion_type_for(&view, i, j, k, l, types[i], types[j], types[k], types[l]);
        let rev = torsion_type_for(&view, l, k, j, i, types[l], types[k], types[j], types[i]);
        let p_fwd = mmff94_torsion_energy(fwd, types[i], types[j], types[k], types[l]).unwrap();
        let p_rev = mmff94_torsion_energy(rev, types[l], types[k], types[j], types[i]).unwrap();
        assert!((p_fwd.v1 - p_rev.v1).abs() < 1e-12);
        assert!((p_fwd.v2 - p_rev.v2).abs() < 1e-12);
        assert!((p_fwd.v3 - p_rev.v3).abs() < 1e-12);
    }

    #[test]
    fn cyclopentane_ring_torsion_uses_distinct_ring_specific_parameters() {
        // Confirms the classification fix is not inert: the real MMFF94
        // table has a dedicated (5, 1, 1, 1, 1) row distinct from the
        // generic (0, 0, 1, 1, 0) wildcard row old code would have used.
        let generic = mmff94_torsion_energy(0, 1, 1, 1, 1).expect("generic type-0 row");
        let ring5 = mmff94_torsion_energy(5, 1, 1, 1, 1).expect("ring type-5 row");
        assert!(
            (generic.v1, generic.v2, generic.v3) != (ring5.v1, ring5.v2, ring5.v3),
            "ring-specific torsion params must differ from the generic fallback"
        );
    }

    // ── FF-2: angle_type_for (#174) ─────────────────────────────────────────

    #[test]
    fn angle_type_for_cyclopropane_ring_angle_is_type_3() {
        let mol = chematic_smiles::parse("C1CC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 3).expect("cyclopropane must have a 3-ring");
        let (i, j, k) = (ring[0], ring[1], ring[2]);
        let rings = find_sssr(&mol);
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 3,
            "3-ring angle with two sp3 (bond_type 0) flanking bonds must be type 3"
        );
        // The specific (3,1,1,1) row doesn't exist for this crate's current
        // atom typer (which doesn't yet distinguish ring-size-specific sp3
        // carbon types 20/22 from plain type 1) — the type-0 fallback added
        // to `mmff94_angle_energy` must still return a usable value rather
        // than silently dropping the angle term.
        assert!(
            mmff94_angle_energy(at, types[i], types[j], types[k]).is_some(),
            "angle energy lookup must fall back to type 0, not silently miss"
        );
    }

    #[test]
    fn angle_type_for_cyclobutane_ring_angle_is_type_4() {
        let mol = chematic_smiles::parse("C1CCC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("cyclobutane must have a 4-ring");
        let (i, j, k) = (ring[0], ring[1], ring[2]);
        let rings = find_sssr(&mol);
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 4,
            "4-ring angle with two sp3 (bond_type 0) flanking bonds must be type 4"
        );
    }

    // ── FF-2b: angle_type_for's corrected bt_sum=1/2 ring-offset formula
    // (issue #227 Priority 2C) ─────────────────────────────────────────────
    // RDKit's real `getMMFFAngleType`: 3-ring bt_sum=2 -> 6 (was wrongly 8);
    // 4-ring bt_sum=1 -> 7 (was wrongly 6), bt_sum=2 -> 8 (was wrongly 7).
    // Measured LATENT on the 265-molecule corpus (0/113 reachable), so these
    // fixtures are constructed molecules, not corpus-derived.

    #[test]
    fn angle_type_for_radialene_cyclopropane_bt_sum_2_ring_angle_is_type_6_not_8() {
        // [3]radialene (trimethylenecyclopropane): all three cyclopropane
        // ring carbons carry an exocyclic =CH2, so every ring bond is a
        // formal single bond between two MLTB (sp2, type 2) atoms ->
        // bond_type 1 on both flanking bonds of any ring angle -> bt_sum=2.
        let mol = chematic_smiles::parse("C1(=C)C(=C)C1=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 3).expect("radialene core must have a 3-ring");
        let (i, j, k) = (ring[0], ring[1], ring[2]);
        let rings = find_sssr(&mol);
        let bt_ij = bond_type_for(types[i], types[j], bond_order_between(&mol, i, j));
        let bt_jk = bond_type_for(types[j], types[k], bond_order_between(&mol, j, k));
        assert_eq!(
            bt_ij + bt_jk,
            2,
            "both ring-flanking bonds must be bond_type 1"
        );
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 6,
            "3-ring angle with bt_sum=2 must be type 6 (RDKit formula), not the old wrong 8"
        );
    }

    /// Real (molecule-adjacent, not just ring-list-adjacent) neighbors of
    /// `j` that also belong to `ring` -- avoids treating a ring's diagonal
    /// (non-bonded) atom pairs as flanking bonds.
    fn ring_bonded_neighbors(mol: &Molecule, ring: &[usize], j: usize) -> Vec<usize> {
        mol.neighbors(AtomIdx(j as u32))
            .map(|(n, _)| n.0 as usize)
            .filter(|n| ring.contains(n))
            .collect()
    }

    #[test]
    fn angle_type_for_4_ring_bt_sum_1_ring_angle_is_type_7_not_6() {
        // Two adjacent ring carbons carry exocyclic =CH2 (MLTB), the other
        // two are plain sp3 CH2 -- the ring angle centered on one
        // MLTB-substituted carbon has exactly one MLTB-MLTB flanking bond
        // (bond_type 1) and one MLTB-sp3 flanking bond (bond_type 0).
        let mol = chematic_smiles::parse("C1(=C)C(=C)CC1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("must have a 4-ring");
        let rings = find_sssr(&mol);
        let mut found = None;
        for &j in &ring {
            let ring_neighbors = ring_bonded_neighbors(&mol, &ring, j);
            assert_eq!(
                ring_neighbors.len(),
                2,
                "each 4-ring atom has exactly 2 ring neighbors"
            );
            let (i, k) = (ring_neighbors[0], ring_neighbors[1]);
            if !atoms_share_ring_of_size(rings.rings(), &[i, j, k], 4) {
                continue;
            }
            let bt_ij = bond_type_for(types[i], types[j], bond_order_between(&mol, i, j));
            let bt_jk = bond_type_for(types[j], types[k], bond_order_between(&mol, j, k));
            if bt_ij + bt_jk == 1 {
                found = Some((i, j, k));
                break;
            }
        }
        let (i, j, k) = found.expect("must find a bt_sum=1 ring angle in this fixture");
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 7,
            "4-ring angle with bt_sum=1 must be type 7 (RDKit formula), not the old wrong 6"
        );
    }

    #[test]
    fn angle_type_for_4_ring_bt_sum_2_ring_angle_is_type_8_not_7() {
        // Three consecutive ring carbons all carry exocyclic =CH2 (MLTB);
        // the angle centered on the middle one has both flanking ring bonds
        // MLTB-MLTB (bond_type 1 each) -> bt_sum=2.
        let mol = chematic_smiles::parse("C1(=C)C(=C)C(=C)C1").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let ring = first_ring_of_size(&mol, 4).expect("must have a 4-ring");
        let rings = find_sssr(&mol);
        let mut found = None;
        for &j in &ring {
            let ring_neighbors = ring_bonded_neighbors(&mol, &ring, j);
            assert_eq!(
                ring_neighbors.len(),
                2,
                "each 4-ring atom has exactly 2 ring neighbors"
            );
            let (i, k) = (ring_neighbors[0], ring_neighbors[1]);
            if !atoms_share_ring_of_size(rings.rings(), &[i, j, k], 4) {
                continue;
            }
            let bt_ij = bond_type_for(types[i], types[j], bond_order_between(&mol, i, j));
            let bt_jk = bond_type_for(types[j], types[k], bond_order_between(&mol, j, k));
            if bt_ij + bt_jk == 2 {
                found = Some((i, j, k));
                break;
            }
        }
        let (i, j, k) = found.expect("must find a bt_sum=2 ring angle in this fixture");
        let at = angle_type_for(&mol, rings.rings(), i, j, k, &types);
        assert_eq!(
            at, 8,
            "4-ring angle with bt_sum=2 must be type 8 (RDKit formula), not the old wrong 7"
        );
    }

    // ── FF-3: stretch_bend_type_for (issue #227 Priority 2C) ────────────────
    // Reproduces the diagnostic's own self-consistency proof
    // (`mmff94_stbn_equivalence_diagnostic_227.rs`'s top doc comment,
    // ~lines 38-56): `MMFF94_STBN`'s frozen data has exactly one row keyed
    // 5 -- an all-CR3R (type 22, cyclopropane ring carbon) triple -- and its
    // 11 key-4 rows are all CR4R (type 20) triples. A 3-ring/bt_sum=0 angle
    // is angle_type 3, never 5 under `angle_type_for`'s own table, so if the
    // STBN key column really were `angle_type`, the key-5 row would be
    // unreachable garbage. `getMMFFStretchBendType` resolves this: angle
    // type 3 -> stretch-bend type 5 (not 3), and angle type 4 -> stretch-
    // bend type 4 (the one case where the numeral is unchanged).

    #[test]
    fn stretch_bend_type_for_cyclopropane_all_cr3r_is_type_5_not_3() {
        assert_eq!(
            stretch_bend_type_for(3, 22, 22, 0, 0),
            5,
            "angle_type 3 (3-ring, bt_sum=0) must map to stretch-bend type 5"
        );
        // End-to-end: the resulting key 5 must actually resolve the real
        // MMFF94_STBN row for the CR3R triple -- confirming the column
        // really is stretch_bend_type, not angle_type.
        let sbt = stretch_bend_type_for(3, 22, 22, 0, 0);
        assert_eq!(
            mmff94_stbn_type_only(sbt, 22, 22, 22),
            Some((0.0, 0.0)),
            "stretch-bend type 5 must resolve the (5,22,22,22,0.0,0.0) MMFF94_STBN row"
        );
    }

    #[test]
    fn stretch_bend_type_for_cr4r_4_ring_is_type_4() {
        assert_eq!(
            stretch_bend_type_for(4, 20, 20, 0, 0),
            4,
            "angle_type 4 (4-ring, bt_sum=0) must map to stretch-bend type 4 (unchanged numeral)"
        );
    }

    #[test]
    fn angle_type_for_butadiene_sp2_single_bond_is_type_1_and_finds_dedicated_row() {
        // Concrete non-ring demonstration that the fix is not inert: the
        // flanking C1-C2 single bond between two vinylic (type 2) carbons is
        // a real MMFF94 bond_type=1 (sbmb) bond, giving angle_type=1, which
        // has a dedicated (1,2,2,2) row (theta0=121.55) distinct from the
        // old hardcoded angle_type=0's row for the same atom-type triple
        // (theta0=118.043, reachable at angle_type=0 via issue #227 Stage
        // B's eqLevel equivalence ladder -- type 2's own eqLevel table
        // substitutes it to type 1 at Level 4, and a real (0,1,2,1) row
        // exists). Getting `angle_type` wrong here is a silently-DIFFERENT
        // wrong value, not a clean miss -- the same "silent wrong parameter"
        // failure class as the #236 furan collision, which is exactly why
        // `angle_type_for`'s own correctness (asserted below) matters.
        let mol = chematic_smiles::parse("C=CC=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        assert_eq!(types[0], 2);
        assert_eq!(types[1], 2);
        assert_eq!(types[2], 2);
        let rings = find_sssr(&mol);
        let at = angle_type_for(&mol, rings.rings(), 0, 1, 2, &types);
        assert_eq!(at, 1, "C0=C1-C2 angle must classify as angle_type 1");
        let wrong = mmff94_angle_energy(0, 2, 2, 2).expect("angle_type=0 now resolves too");
        let p = mmff94_angle_energy(at, 2, 2, 2).expect("(1,2,2,2) row must exist");
        assert!(
            (p.theta0 - 121.55).abs() < 1e-6,
            "theta0={} should be 121.55",
            p.theta0
        );
        assert!(
            (p.theta0 - wrong.theta0).abs() > 1.0,
            "angle_type=0's wrong theta0={} should differ meaningfully from the correct {}",
            wrong.theta0,
            p.theta0
        );
    }

    #[test]
    fn butadiene_angle_energy_minimum_is_near_correct_theta0() {
        // Energy-perturbation test: with the fix, this angle now has a real
        // restoring force toward its correct 121.55° equilibrium (pre-fix it
        // silently contributed zero energy everywhere — no minimum at all).
        let mol = chematic_smiles::parse("C=CC=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let rings = find_sssr(&mol);
        let scan_energy = |theta_deg: f64| {
            let r = 1.45_f64;
            let half = theta_deg.to_radians() / 2.0;
            // j (atom1) at origin; i (atom0) and k (atom2) placed to form
            // the target angle at j.
            let coords = vec![
                [r * half.cos(), r * half.sin(), 0.0],
                [0.0, 0.0, 0.0],
                [r * half.cos(), -r * half.sin(), 0.0],
                [r * half.cos() + 1.3, -r * half.sin() - 0.9, 0.0],
            ];
            angle_energy(&mol, &coords, &types, rings.rings())
        };
        let e_correct = scan_energy(121.55);
        let e_distorted = scan_energy(100.0);
        assert!(
            e_distorted > e_correct,
            "energy away from the true 121.55° minimum ({e_distorted}) should exceed energy \
             at the minimum ({e_correct})"
        );
    }

    #[test]
    fn angle_gradient_restores_butadiene_toward_theta0() {
        let mol = chematic_smiles::parse("C=CC=C").unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let charges = mmff94_charges_numeric(&mol).unwrap_or_else(|_| vec![0.0; mol.atom_count()]);
        let rings = find_sssr(&mol);
        let r = 1.45_f64;
        let half = 100.0_f64.to_radians() / 2.0; // distorted away from 121.55°
        let coords = vec![
            [r * half.cos(), r * half.sin(), 0.0],
            [0.0, 0.0, 0.0],
            [r * half.cos(), -r * half.sin(), 0.0],
            [r * half.cos() + 1.3, -r * half.sin() - 0.9, 0.0],
        ];
        let grad = compute_gradient(&mol, &coords, &types, &charges, rings.rings(), 1e-4);
        let grad_norm: f64 = grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(
            grad_norm > 1e-6,
            "gradient should be nonzero away from equilibrium angle"
        );
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

    #[test]
    fn lbfgs_reduces_energy_for_methane() {
        let (mol, _) = methane_mol();
        let mut coords = vec![
            [0.0, 0.0, 0.0_f64],
            [1.5, 1.5, 1.5],
            [-1.5, -1.5, 1.5],
            [-1.5, 1.5, -1.5],
            [1.5, -1.5, -1.5],
        ];
        let e_before = mmff94_total_energy(&mol, &coords).expect("energy before");
        let result = minimize_mmff94_lbfgs(&mol, &mut coords, 300).expect("lbfgs");
        assert!(
            result.energy <= e_before,
            "L-BFGS should reduce energy: {} → {}",
            e_before,
            result.energy
        );
        assert!(result.energy.is_finite());
    }

    #[test]
    fn lbfgs_converges_in_fewer_iters_than_sd() {
        let (mol, _) = methane_mol();
        // Moderately distorted — both should converge but L-BFGS faster
        let base_coords = vec![
            [0.0, 0.0, 0.0_f64],
            [1.2, 1.2, 1.2],
            [-1.2, -1.2, 1.2],
            [-1.2, 1.2, -1.2],
            [1.2, -1.2, -1.2],
        ];
        let mut coords_sd = base_coords.clone();
        let mut coords_lbfgs = base_coords;
        let sd = minimize_mmff94_full(&mol, &mut coords_sd, 500).expect("sd");
        let lb = minimize_mmff94_lbfgs(&mol, &mut coords_lbfgs, 500).expect("lbfgs");
        // Both should converge; L-BFGS should need ≤ SD iterations
        assert!(
            lb.iterations <= sd.iterations || lb.converged,
            "L-BFGS iters={} SD iters={}",
            lb.iterations,
            sd.iterations
        );
        assert!(lb.energy.is_finite());
    }

    #[test]
    fn energy_breakdown_sums_to_total() {
        let (mol, coords) = methane_mol();
        let bd = mmff94_energy_breakdown(&mol, &coords).expect("breakdown");
        let sum =
            bd.bond + bd.angle + bd.stretch_bend + bd.torsion + bd.oop + bd.vdw + bd.electrostatic;
        assert!(
            (sum - bd.total).abs() < 1e-10,
            "sum={} total={}",
            sum,
            bd.total
        );
        assert!(bd.total.is_finite());
    }

    /// Frozen 58-molecule corpus, copied from `scripts/etkdg_vs_rdkit_gap.py::CORPUS`
    /// (the same corpus PR #169's `mmff94_bridge_coverage_report` example measured
    /// coverage against). Used here to count MMFF94 parameter-lookup *misses*
    /// per term class (bond/angle/torsion) before vs. after the classification
    /// fixes in this file — the metric that actually matters, since a
    /// "corrected" classification that routes to a still-missing row is a
    /// silent zero-energy regression, not an improvement.
    const CORPUS_58: &[(&str, &str)] = &[
        ("benzene", "c1ccccc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("pyridine", "c1ccncc1"),
        ("furan", "c1ccoc1"),
        ("thiophene", "c1ccsc1"),
        ("adamantane", "C1CC2CC3CC1CC(C2)C3"),
        ("cubane", "C1C2C3C1C4C2C3C4"),
        ("cyclohexane", "C1CCCCC1"),
        ("cyclopentane", "C1CCCC1"),
        ("indole", "c1ccc2[nH]ccc2c1"),
        ("purine", "c1ncc2[nH]cnc2n1"),
        ("quinoline", "c1ccc2ncccc2c1"),
        ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
        ("pyrene", "c1cc2ccc3cccc4ccc(c1)c2c34"),
        ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
        ("butane", "CCCC"),
        ("hexane", "CCCCCC"),
        ("decane", "CCCCCCCCCC"),
        ("triethylene_glycol", "OCCOCCOCCO"),
        ("hexanediol", "OCCCCCCO"),
        ("hexadecane", "CCCCCCCCCCCCCCCC"),
        ("cyclododecane", "C1CCCCCCCCCCC1"),
        ("crown_12_4", "O1CCOCCOCCOCC1"),
        ("cyclooctadecane", "C1CCCCCCCCCCCCCCCCC1"),
        ("l_alanine", "N[C@@H](C)C(=O)O"),
        ("d_alanine", "N[C@H](C)C(=O)O"),
        ("l_serine", "N[C@@H](CO)C(=O)O"),
        ("l_threonine", "C[C@H](O)[C@@H](N)C(=O)O"),
        ("2_butanol_R", "C[C@H](O)CC"),
        ("2_butanol_S", "C[C@@H](O)CC"),
        ("2_chlorobutane_R", "C[C@H](Cl)CC"),
        ("ibuprofen_S", "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"),
        ("naproxen_S", "COc1ccc2cc([C@H](C)C(=O)O)ccc2c1"),
        ("menthol", "C[C@@H]1CC[C@@H](C(C)C)C[C@H]1O"),
        ("chfclbr_R", "[C@H](F)(Cl)Br"),
        ("chfclbr_S", "[C@@H](F)(Cl)Br"),
        ("quaternary_1_R", "[C@](F)(Cl)(Br)I"),
        ("quaternary_1_S", "[C@@](F)(Cl)(Br)I"),
        ("quaternary_2_R", "[C@](C)(N)(O)F"),
        ("quaternary_2_S", "[C@@](C)(N)(O)F"),
        ("but2ene_E", "C/C=C/C"),
        ("but2ene_Z", r"C/C=C\C"),
        ("chloropropene_E", "C(/C=C/C)Cl"),
        ("chloropropene_Z", r"C(/C=C\C)Cl"),
        ("cinnamic_acid_E", "OC(=O)/C=C/c1ccccc1"),
        ("cinnamic_acid_Z", r"OC(=O)/C=C\c1ccccc1"),
        ("pent2ene_E", "CC/C=C/C"),
        ("pent2ene_Z", r"CC/C=C\C"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
        ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
        ("diphenhydramine", "CN(C)CCOC(c1ccccc1)c1ccccc1"),
        (
            "penicillin_core",
            "CC1(C)S[C@@H]2[C@H](NC(=O)C)C(=O)N2[C@H]1C(=O)O",
        ),
        (
            "testosterone",
            "C[C@]12CC[C@H]3[C@@H](CC[C@H]4CCC(=O)C=C34)[C@@H]1CC[C@@H]2O",
        ),
        (
            "cholesterol",
            "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        ),
        (
            "atorvastatin_fragment",
            "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        ),
        ("gly_ala_gly", "NCC(=O)N[C@@H](C)C(=O)NCC(=O)O"),
    ];

    /// Per-term-class MMFF94 parameter-lookup miss counts across [`CORPUS_58`].
    struct CoverageCounts {
        bond_total: usize,
        bond_miss: usize,
        angle_total: usize,
        angle_miss: usize,
        torsion_total: usize,
        torsion_miss: usize,
    }

    fn corpus_coverage() -> CoverageCounts {
        let mut c = CoverageCounts {
            bond_total: 0,
            bond_miss: 0,
            angle_total: 0,
            angle_miss: 0,
            torsion_total: 0,
            torsion_miss: 0,
        };
        for (name, smiles) in CORPUS_58 {
            let mol = chematic_smiles::parse(smiles).unwrap_or_else(|e| {
                panic!("corpus molecule {name} ({smiles}) failed to parse: {e}")
            });
            let types = match assign_mmff94_numeric_types(&mol) {
                Ok(t) => t,
                Err(_) => continue, // unsupported element/type — not this file's concern
            };
            let rings = find_sssr(&mol);

            for (_, bond) in mol.bonds() {
                let i = bond.atom1.0 as usize;
                let j = bond.atom2.0 as usize;
                c.bond_total += 1;
                let bt = bond_type_for(types[i], types[j], bond.order);
                if mmff94_bond_energy(bt, types[i], types[j]).is_none() {
                    c.bond_miss += 1;
                }
            }

            for j_idx in 0..mol.atom_count() {
                let j = AtomIdx(j_idx as u32);
                let neighbors: Vec<usize> = mol.neighbors(j).map(|(nb, _)| nb.0 as usize).collect();
                if neighbors.len() < 2 {
                    continue;
                }
                for (ii, &i) in neighbors.iter().enumerate() {
                    for &k in &neighbors[ii + 1..] {
                        c.angle_total += 1;
                        let at = angle_type_for(&mol, rings.rings(), i, j_idx, k, &types);
                        if mmff94_angle_energy(at, types[i], types[j_idx], types[k]).is_none() {
                            c.angle_miss += 1;
                        }
                    }
                }
            }

            for (_, bond) in mol.bonds() {
                let j = bond.atom1.0 as usize;
                let k = bond.atom2.0 as usize;
                let nbrs_j: Vec<usize> = mol
                    .neighbors(bond.atom1)
                    .map(|(nb, _)| nb.0 as usize)
                    .collect();
                let nbrs_k: Vec<usize> = mol
                    .neighbors(bond.atom2)
                    .map(|(nb, _)| nb.0 as usize)
                    .collect();
                for &i in &nbrs_j {
                    if i == k {
                        continue;
                    }
                    for &l in &nbrs_k {
                        if l == j {
                            continue;
                        }
                        c.torsion_total += 1;
                        let tt = torsion_type_for(
                            &mol, i, j, k, l, types[i], types[j], types[k], types[l],
                        );
                        if mmff94_torsion_energy(tt, types[i], types[j], types[k], types[l])
                            .is_none()
                        {
                            c.torsion_miss += 1;
                        }
                    }
                }
            }
        }
        c
    }

    #[test]
    fn corpus_58_parameter_coverage_does_not_regress() {
        // Numbers on THIS 58-molecule corpus (CORPUS_58, distinct from the
        // 265-molecule Wave 1 corpus used elsewhere in issue #227) have
        // drifted across several since-merged fixes; re-measured fresh
        // rather than trusted from stale history:
        //   original bond_type_for/angle_type_for/torsion_type_for:
        //     bond 129/585 miss, angle 292/737 miss, torsion 305/841 miss
        //   immediately pre-this-fix (post atom-type-parity + stretch-bend
        //   fixes, i.e. this branch's parent commit):
        //     bond   1/585 miss, angle  24/737 miss, torsion  22/841 miss
        //   post this fix (issue #227 torsion classification: bond-type-jk-
        //   based base case + local ring-4/5 override, replacing the old
        //   atom-type-membership + SSSR-based formula):
        //     bond   1/585 miss, angle  24/737 miss, torsion   4/841 miss
        // Torsion coverage DOES improve now (unlike the original bond/angle-
        // only fix this comment used to describe): the corrected
        // classification routes several of this corpus's torsions to a
        // table row chematic's own existing fallback chain can reach that
        // the old, wrong classification code could not.
        let c = corpus_coverage();
        eprintln!(
            "bond: {}/{} miss, angle: {}/{} miss, torsion: {}/{} miss",
            c.bond_miss, c.bond_total, c.angle_miss, c.angle_total, c.torsion_miss, c.torsion_total
        );
        // Loose tripwires (see MEMORY.md's note on avoiding over-tuned
        // pseudo-gates), not tight equality: a future change may legitimately
        // shift these by a few rows (e.g. an atom-typer improvement), but
        // must never regress back toward the pre-fix counts above.
        assert!(
            c.bond_miss <= 45,
            "bond misses regressed toward pre-fix (129): {}",
            c.bond_miss
        );
        assert!(
            c.angle_miss <= 290,
            "angle misses regressed toward pre-fix (292): {}",
            c.angle_miss
        );
        assert!(
            c.torsion_miss <= 20,
            "torsion misses regressed toward pre-this-fix (22): {}",
            c.torsion_miss
        );
    }

    #[test]
    fn energy_breakdown_bond_term_positive_for_distorted() {
        let (mol, _) = methane_mol();
        // Very distorted C-H bonds → high bond energy
        let stretched = vec![
            [0.0, 0.0, 0.0_f64],
            [2.0, 2.0, 2.0],
            [-2.0, -2.0, 2.0],
            [-2.0, 2.0, -2.0],
            [2.0, -2.0, -2.0],
        ];
        let bd = mmff94_energy_breakdown(&mol, &stretched).expect("breakdown");
        assert!(
            bd.bond > 0.0,
            "stretched bond energy should be positive: {}",
            bd.bond
        );
    }
}
