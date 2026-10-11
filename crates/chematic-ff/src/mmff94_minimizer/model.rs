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
    /// Central atom of a linear MMFF type: `143.9325 ka (1 + cos θ)`.
    linear: bool,
}

#[derive(Clone, Copy)]
struct PreparedStretchBend {
    i: usize,
    j: usize,
    k: usize,
    kba_ij: f64,
    kba_kj: f64,
    r0_ij: f64,
    r0_kj: f64,
    theta0: f64,
}

#[derive(Clone, Copy)]
struct PreparedTorsion {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    params: TorsionEnergyParams,
}

#[derive(Clone, Copy)]
struct PreparedOop {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    koop: f64,
}

#[derive(Clone, Copy)]
struct PreparedVdwPair {
    i: usize,
    j: usize,
    r_star: f64,
    epsilon: f64,
}

#[derive(Clone, Copy)]
struct PreparedElectrostaticPair {
    i: usize,
    j: usize,
    scale_charge_product: f64,
}

// ─── Public types ────────────────────────────────────────────────────────────

/// Result of a geometry minimization run.
#[derive(Debug, Clone)]
pub struct MinimizeResult {
    /// Final MMFF94 energy (kcal/mol).
    pub energy: f64,
    /// RMSD of atom positions vs initial geometry (Å).
    pub rmsd: f64,
    /// Whether the optimizer's gradient convergence check passed.
    pub converged: bool,
    /// Number of optimizer loop iterations performed.
    pub iterations: usize,
    /// Why this run stopped; non-convergence is not always an iteration cap.
    pub termination: Mmff94TerminationReason,
}

/// Terminal condition of an MMFF94 minimization run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mmff94TerminationReason {
    /// The gradient fell below the configured threshold.
    GradientConverged,
    /// The configured iteration count was exhausted, including a zero budget.
    IterationLimit,
    /// Backtracking and steepest-descent recovery could not produce a finite,
    /// descending proposal accepted by the caller's geometry predicate.
    /// Includes non-finite input states. The last accepted coordinates remain.
    ConstraintRejectedFallback,
}

/// Error from MMFF94 minimizer setup.
#[derive(Debug)]
pub enum MinimizerError {
    TypeAssignment(NumericTypeError),
    ChargeCalculation(NumericTypeError),
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
            MinimizerError::ChargeCalculation(e) => {
                write!(f, "MMFF94 charge calculation failed: {}", e)
            }
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

/// Options for building an MMFF94 energy model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Mmff94Options {
    /// Leave out van der Waals and electrostatic pairs between disconnected
    /// fragments (a salt's ions, a solvent molecule), as RDKit's
    /// `MMFFGetMoleculeForceField` does by default
    /// (`ignoreInterfragInteractions=True`). Off by default: chematic's
    /// minimizer keeps the fragments of a salt from overlapping through them.
    pub ignore_interfragment_interactions: bool,
}

/// Connected-component id of every atom.
fn fragment_ids(mol: &Molecule) -> Vec<usize> {
    let n = mol.atom_count();
    let mut id = vec![usize::MAX; n];
    let mut next = 0;
    for start in 0..n {
        if id[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        id[start] = next;
        while let Some(a) = stack.pop() {
            for (nb, _) in mol.neighbors(AtomIdx(a as u32)) {
                let b = nb.0 as usize;
                if id[b] == usize::MAX {
                    id[b] = next;
                    stack.push(b);
                }
            }
        }
        next += 1;
    }
    id
}

/// Reusable MMFF94 energy evaluator for one molecular topology.
///
/// Atom typing, MMFF94's aromatic view, charges, and ring perception are
/// topology-only. Keeping them together avoids repeating that work when a
/// caller evaluates many nearby geometries, as finite-difference gradients do.
#[derive(Clone)]
pub struct Mmff94EnergyModel {
    mmff_mol: Molecule,
    vdw_pairs: VdwPairs,
    electrostatic_pairs: ElectrostaticPairs,
    bonds: Vec<PreparedBond>,
    angles: Vec<PreparedAngle>,
    stretch_bends: Vec<PreparedStretchBend>,
    torsions: Vec<PreparedTorsion>,
    oops: Vec<PreparedOop>,
}

impl Mmff94EnergyModel {
    /// Prepare the topology-dependent MMFF94 state once.
    pub fn new(mol: &Molecule) -> Result<Self, MinimizerError> {
        Self::new_with_options(mol, Mmff94Options::default())
    }

    /// [`Self::new`] with [`Mmff94Options`].
    pub fn new_with_options(
        mol: &Molecule,
        options: Mmff94Options,
    ) -> Result<Self, MinimizerError> {
        let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
        let charges = mmff94_charges_numeric(mol).map_err(MinimizerError::ChargeCalculation)?;
        let rings = find_sssr(mol).rings().to_vec();
        let (mut vdw_pairs, mut electrostatic_pairs) = build_nonbonded_pairs(mol, &types, &charges);
        if options.ignore_interfragment_interactions {
            let frag = fragment_ids(mol);
            vdw_pairs.retain(|p| frag[p.i] == frag[p.j]);
            electrostatic_pairs.retain(|p| frag[p.i] == frag[p.j]);
        }
        let bonds = build_bond_terms(&mmff_mol, &types);
        let angles = build_angle_terms(&mmff_mol, &types, &rings);
        let stretch_bends = build_stretch_bend_terms(&mmff_mol, &types, &rings);
        let torsions = build_torsion_terms(&mmff_mol, &types);
        let oops = build_oop_terms(&mmff_mol, &types);
        Ok(Self {
            mmff_mol,
            vdw_pairs,
            electrostatic_pairs,
            bonds,
            angles,
            stretch_bends,
            torsions,
            oops,
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

    fn energy_cutoff_nonbonded(&self, coords: &[[f64; 3]]) -> f64 {
        let (vdw, electrostatic) = self.cutoff_nonbonded_energy(coords);
        self.bond_energy(coords)
            + self.angle_energy(coords)
            + self.stretch_bend_energy(coords)
            + self.torsion_energy(coords)
            + self.oop_energy(coords)
            + vdw
            + electrostatic
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

    /// Analytic gradient for the prepared non-bonded terms.
    ///
    /// This is intentionally exposed as a bounded building block: bonded,
    /// torsional, out-of-plane, and stretch-bend terms still have separate
    /// gradient gates in the minimizer.
    pub fn nonbonded_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        prepared_nonbonded_gradient(coords, &self.vdw_pairs, &self.electrostatic_pairs)
    }

    /// Evaluate the prepared non-bonded terms with the bounded 10 Å cutoff.
    ///
    /// This is opt-in because the compatibility-default MMFF94 evaluator
    /// retains its historical all-pair electrostatic path.  The returned
    /// tuple is `(vdw, electrostatic)` in kcal/mol; both components use the
    /// same coordinate-dependent cutoff boundary.
    pub fn cutoff_nonbonded_energy(&self, coords: &[[f64; 3]]) -> (f64, f64) {
        (
            vdw_energy_cutoff_pairs(coords, &self.vdw_pairs),
            electrostatic_energy_cutoff_pairs(coords, &self.electrostatic_pairs),
        )
    }

    /// Analytic gradient for the bounded-cutoff non-bonded terms.
    ///
    /// This pairs [`Self::cutoff_nonbonded_energy`] with the same strict
    /// `distance <= 10 Å` inclusion rule.  It is experimental and does not
    /// alter [`Self::nonbonded_gradient`].
    pub fn cutoff_nonbonded_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        prepared_cutoff_nonbonded_gradient(coords, &self.vdw_pairs, &self.electrostatic_pairs)
    }

    /// Analytic gradient for the prepared bond and angle terms.
    ///
    /// Torsion, out-of-plane, and stretch-bend terms are deliberately kept
    /// outside this bounded building block until their singular geometries
    /// have dedicated soundness gates.
    pub fn bond_angle_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        prepared_bond_angle_gradient(coords, &self.bonds, &self.angles)
    }

    /// Analytic gradient for the prepared stretch-bend coupling terms.
    pub fn stretch_bend_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        prepared_stretch_bend_gradient(coords, &self.stretch_bends)
    }

    /// Analytic gradient for prepared torsion terms.
    ///
    /// Collinear or zero-area dihedrals are skipped because their signed
    /// angle has no stable Cartesian derivative.  The existing minimizer
    /// continues to use its finite-difference gradient until the remaining
    /// out-of-plane and singular-geometry gates are complete.
    pub fn torsion_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        prepared_torsion_gradient(coords, &self.torsions)
    }

    /// Analytic gradient for prepared, non-singular Wilson out-of-plane terms.
    pub fn oop_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        prepared_oop_gradient(coords, &self.oops)
    }

    /// Combine the currently soundness-gated analytic term gradients.
    ///
    /// The prepared non-bonded pair set is intentionally fixed for this
    /// no-cutoff evaluator, matching [`Self::energy`]. The 3D pipeline uses
    /// this exact energy/gradient pair; the cutoff variant remains opt-in.
    pub fn bounded_analytic_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.bounded_analytic_gradient_with_nonbonded(false, coords)
    }

    /// Analytic gradient with the bounded 10 Å coordinate-dependent
    /// non-bonded neighbor lists enabled.
    ///
    /// This is an experimental opt-in path. It keeps the compatibility
    /// default in [`Self::bounded_analytic_gradient`] unchanged while making
    /// the cutoff list usable by the analytic minimizer as one consistent
    /// energy/gradient pair.
    pub fn bounded_analytic_gradient_cutoff(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.bounded_analytic_gradient_with_nonbonded(true, coords)
    }

    fn bounded_analytic_gradient_with_nonbonded(
        &self,
        use_cutoff_nonbonded: bool,
        coords: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let nonbonded = if use_cutoff_nonbonded {
            self.cutoff_nonbonded_gradient(coords)
        } else {
            self.nonbonded_gradient(coords)
        };
        let components = [
            self.bond_angle_gradient(coords),
            self.stretch_bend_gradient(coords),
            self.torsion_gradient(coords),
            self.oop_gradient(coords),
            nonbonded,
        ];
        let mut gradient = vec![[0.0; 3]; coords.len()];
        for component in components {
            for (total, part) in gradient.iter_mut().zip(component) {
                for axis in 0..3 {
                    total[axis] += part[axis];
                }
            }
        }
        gradient
    }

    /// Minimize coordinates using the prepared MMFF94 topology state.
    pub fn minimize_lbfgs(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
    ) -> Result<MinimizeResult, MinimizerError> {
        minimize_mmff94_lbfgs_prepared(self, coords, max_iter)
    }

    /// L-BFGS path using the bounded prepared analytic gradient.
    ///
    /// Singular derivative geometries fail closed. The finite-difference
    /// [`Self::minimize_lbfgs`] method remains available as a reference path;
    /// the 3D MMFF94 bridge uses this analytic path after its corpus-level
    /// gradient, geometry, convergence, and timeout gates.
    pub fn minimize_lbfgs_bounded_analytic(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
    ) -> Result<MinimizeResult, MinimizerError> {
        minimize_mmff94_lbfgs_prepared_with_mode(self, coords, max_iter, true, false, |_| true)
    }

    /// Analytic L-BFGS with a caller-supplied acceptance constraint.
    ///
    /// The predicate is evaluated for each energy-decreasing line-search
    /// proposal. Rejected proposals are backtracked; if no bounded step can
    /// satisfy the predicate, minimization stops at the last accepted
    /// geometry with `converged = false`. The initial coordinates remain the
    /// caller's responsibility and are never replaced by a rejected proposal.
    pub fn minimize_lbfgs_bounded_analytic_with_constraint<F>(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
        accept: F,
    ) -> Result<MinimizeResult, MinimizerError>
    where
        F: Fn(&[[f64; 3]]) -> bool,
    {
        minimize_mmff94_lbfgs_prepared_with_mode(self, coords, max_iter, true, false, accept)
    }

    /// Dense BFGS with RDKit's line search and gradient scaling
    /// (`BFGSOpt::minimize` as called by `ForceField::minimize`).
    ///
    /// With [`Mmff94Convergence::Rdkit`] this is RDKit's minimizer: from the
    /// same coordinates it follows RDKit's trajectory, so `MMFFOptimizeMolecule
    /// (maxIters=n)` and this method with `max_iter = n` reach the same
    /// geometry (to rounding). With [`Mmff94Convergence::MaxGradient`] it
    /// runs the same algorithm to chematic's absolute residual-force test,
    /// never accepts an energy increase and bounds each line search's first
    /// trial by [`BFGS_MAX_COORDINATE_STEP`]. The inverse Hessian is a dense
    /// `3N x 3N` matrix; above [`BFGS_DENSE_MAX_ATOMS`] atoms the call
    /// delegates to [`Self::minimize_lbfgs_bounded_analytic`] (MaxGradient)
    /// so memory stays bounded.
    pub fn minimize_bfgs(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
        convergence: Mmff94Convergence,
    ) -> Result<MinimizeResult, MinimizerError> {
        minimize_mmff94_bfgs_prepared(self, coords, max_iter, convergence, |_| true)
    }

    /// [`Self::minimize_bfgs`] with a caller-supplied acceptance constraint.
    ///
    /// A line-search proposal the predicate rejects is backtracked like one
    /// that fails the sufficient-decrease test. If no admissible step is
    /// left, the last accepted coordinates are kept and the run reports
    /// [`Mmff94TerminationReason::ConstraintRejectedFallback`].
    pub fn minimize_bfgs_with_constraint<F>(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
        convergence: Mmff94Convergence,
        accept: F,
    ) -> Result<MinimizeResult, MinimizerError>
    where
        F: Fn(&[[f64; 3]]) -> bool,
    {
        minimize_mmff94_bfgs_prepared(self, coords, max_iter, convergence, accept)
    }

    /// Experimental analytic L-BFGS using the bounded 10 Å coordinate-
    /// dependent non-bonded neighbor lists.
    ///
    /// The cutoff energy and gradient are switched together, so line search
    /// and curvature updates see one consistent objective. The established
    /// finite-difference and non-cutoff analytic methods remain unchanged.
    pub fn minimize_lbfgs_bounded_analytic_cutoff(
        &self,
        coords: &mut [[f64; 3]],
        max_iter: usize,
    ) -> Result<MinimizeResult, MinimizerError> {
        minimize_mmff94_lbfgs_prepared_with_mode(self, coords, max_iter, true, true, |_| true)
    }

    fn bond_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_bond_energy(coords, &self.bonds)
    }

    fn angle_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_angle_energy(coords, &self.angles)
    }

    fn stretch_bend_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_stretch_bend_energy(coords, &self.stretch_bends)
    }

    fn torsion_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_torsion_energy(coords, &self.torsions)
    }

    fn oop_energy(&self, coords: &[[f64; 3]]) -> f64 {
        prepared_oop_energy(coords, &self.oops)
    }

    fn vdw_energy(&self, coords: &[[f64; 3]]) -> f64 {
        vdw_energy_pairs(coords, &self.vdw_pairs)
    }

    fn electrostatic_energy(&self, coords: &[[f64; 3]]) -> f64 {
        electrostatic_energy_pairs(coords, &self.electrostatic_pairs)
    }
}

/// Prepare only the topology state consumed by the RDKit-order force-field
/// builder. The general [`Mmff94EnergyModel`] also constructs its own angle,
/// stretch-bend and non-bonded tables; rebuilding those here was redundant
/// because [`crate::rdkit_mmff`] must create them in RDKit contribution order.
pub(crate) fn rdkit_prepared_parts(
    mol: &Molecule,
) -> Result<(crate::rdkit_mmff::PreparedParts, Vec<u8>, Molecule), MinimizerError> {
    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let bonds = build_bond_terms(&mmff_mol, &types);
    let torsions = build_torsion_terms(&mmff_mol, &types);
    let oops = build_oop_terms(&mmff_mol, &types);
    let parts = crate::rdkit_mmff::PreparedParts {
        bonds: bonds
            .iter()
            .map(|b| (b.i, b.j, b.params.r0, b.params.kb))
            .collect(),
        oops: oops.iter().map(|o| (o.j, o.koop)).collect(),
        torsions: torsions
            .iter()
            .map(|t| (t.i, t.j, t.k, t.l, t.params.v1, t.params.v2, t.params.v3))
            .collect(),
    };
    Ok((parts, types, mmff_mol))
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
    let charges = mmff94_charges_numeric(mol).map_err(MinimizerError::ChargeCalculation)?;
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
            let (sin_a, cos_a) = step_rad.dsin_cos();
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

/// [`mmff94_energy_breakdown`] with [`Mmff94Options`] (e.g. RDKit's default
/// of no interactions between fragments).
pub fn mmff94_energy_breakdown_with_options(
    mol: &Molecule,
    coords: &[[f64; 3]],
    options: Mmff94Options,
) -> Result<EnergyBreakdown, MinimizerError> {
    Ok(Mmff94EnergyModel::new_with_options(mol, options)?.energy_breakdown(coords))
}
