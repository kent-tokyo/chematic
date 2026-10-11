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
            termination: Mmff94TerminationReason::GradientConverged,
        });
    }

    let (types, mmff_mol) = assign_mmff94_numeric_types_with_view(mol)?;
    let charges = mmff94_charges_numeric(mol).map_err(MinimizerError::ChargeCalculation)?;
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
        termination: if converged {
            Mmff94TerminationReason::GradientConverged
        } else {
            Mmff94TerminationReason::IterationLimit
        },
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
    minimize_mmff94_lbfgs_prepared_with_mode(model, coords, max_iter, false, false, |_| true)
}

fn minimize_mmff94_lbfgs_prepared_with_mode<F>(
    model: &Mmff94EnergyModel,
    coords: &mut [[f64; 3]],
    max_iter: usize,
    use_analytic_gradient: bool,
    use_cutoff_nonbonded: bool,
    accept: F,
) -> Result<MinimizeResult, MinimizerError>
where
    F: Fn(&[[f64; 3]]) -> bool,
{
    // The 265-molecule A6 gate shows that retaining 80 correction pairs reaches
    // substantially more stationary points than the historical five-pair
    // window without a measurable throughput penalty. Keep this bounded: each
    // entry stores two coordinate-sized vectors and is discarded at run end.
    const CONVERGENCE: f64 = 1e-4; // max |gradient| threshold

    if model.mmff_mol.atom_count() <= 1 {
        return Ok(MinimizeResult {
            energy: 0.0,
            rmsd: 0.0,
            converged: true,
            iterations: 0,
            termination: Mmff94TerminationReason::GradientConverged,
        });
    }

    let initial = coords.to_vec();

    // Circular history buffer: (s_k = Δx, y_k = Δg, ρ_k = 1/(y·s))
    let mut history: LbfgsHistory = VecDeque::new();

    let mut g = lbfgs_gradient(model, coords, use_analytic_gradient, use_cutoff_nonbonded);
    let mut f0 = lbfgs_energy(model, coords, use_cutoff_nonbonded);

    let mut iters = 0usize;
    let mut termination = Mmff94TerminationReason::IterationLimit;

    'iterations: for _ in 0..max_iter {
        iters += 1;

        // Convergence check
        let max_g = finite_max_gradient(&g);
        if !f0.is_finite() || !max_g.is_finite() {
            termination = Mmff94TerminationReason::ConstraintRejectedFallback;
            break;
        }
        if max_g < CONVERGENCE {
            break;
        }

        // Two-loop L-BFGS recursion → search direction p. Numerical noise in
        // a nearly singular curvature update can occasionally produce a
        // non-descent or non-finite direction. Feeding that direction to the
        // Armijo loop causes repeated rejected objective evaluations (and can
        // consume a pipeline timeout without making progress). Fall back to
        // steepest descent, which is always a valid local descent direction
        // for a finite, non-zero gradient.
        let (p, gp) = descent_direction(&g, &mut history);
        let step = lbfgs_line_search(
            model,
            coords,
            &g,
            &p,
            f0,
            gp,
            max_g,
            use_cutoff_nonbonded,
            &mut history,
            &accept,
        );
        let Some((new_coords, f_new)) = step else {
            // No admissible bounded step exists along either the L-BFGS or
            // tiny steepest-descent direction. Preserve the last accepted
            // coordinates and report ordinary non-convergence; never cross a
            // caller-defined geometry boundary just to keep iterating.
            termination = Mmff94TerminationReason::ConstraintRejectedFallback;
            break 'iterations;
        };

        // Compute new gradient
        let g_new = lbfgs_gradient(
            model,
            &new_coords,
            use_analytic_gradient,
            use_cutoff_nonbonded,
        );
        // `f_new` is the accepted line-search energy above; do not evaluate
        // the same coordinates a second time after computing the gradient.

        // Compute s = x_new - x, y = g_new - g
        update_lbfgs_history(&mut history, coords, &new_coords, &g, &g_new);

        coords.copy_from_slice(&new_coords);
        g = g_new;
        f0 = f_new;
    }

    // The final permitted step can reach the tolerance. Report the accepted
    // state instead of requiring another loop iteration to observe it.
    let converged = f0.is_finite() && finite_max_gradient(&g) < CONVERGENCE;
    if converged {
        termination = Mmff94TerminationReason::GradientConverged;
    }

    Ok(MinimizeResult {
        energy: f0,
        rmsd: coordinate_rmsd(coords, &initial),
        converged,
        iterations: iters,
        termination,
    })
}

fn lbfgs_gradient(
    model: &Mmff94EnergyModel,
    coords: &[[f64; 3]],
    analytic: bool,
    cutoff: bool,
) -> CoordVec {
    if !analytic {
        return compute_gradient_prepared(model, coords, 1e-4);
    }
    if cutoff {
        model.bounded_analytic_gradient_cutoff(coords)
    } else {
        model.bounded_analytic_gradient(coords)
    }
}

fn lbfgs_energy(model: &Mmff94EnergyModel, coords: &[[f64; 3]], cutoff: bool) -> f64 {
    if cutoff {
        model.energy_cutoff_nonbonded(coords)
    } else {
        model.energy(coords)
    }
}

fn finite_max_gradient(gradient: &[[f64; 3]]) -> f64 {
    gradient.iter().flatten().fold(0.0, |maximum, &value| {
        if value.is_finite() {
            maximum.max(value.abs())
        } else {
            f64::INFINITY
        }
    })
}

fn descent_direction(g: &[[f64; 3]], history: &mut LbfgsHistory) -> (CoordVec, f64) {
    let mut direction = lbfgs_direction(g, history);
    let mut slope: f64 = g
        .iter()
        .zip(&direction)
        .map(|(gi, pi)| dot3(*gi, *pi))
        .sum();
    if !slope.is_finite() || slope >= 0.0 {
        history.clear();
        direction = g.iter().map(|gi| [-gi[0], -gi[1], -gi[2]]).collect();
        slope = -g.iter().map(|gi| dot3(*gi, *gi)).sum::<f64>();
    }
    (direction, slope)
}

#[allow(clippy::too_many_arguments)]
fn lbfgs_line_search<F>(
    model: &Mmff94EnergyModel,
    coords: &[[f64; 3]],
    gradient: &[[f64; 3]],
    direction: &[[f64; 3]],
    energy: f64,
    slope: f64,
    max_gradient: f64,
    cutoff: bool,
    history: &mut LbfgsHistory,
    accept: &F,
) -> Option<(CoordVec, f64)>
where
    F: Fn(&[[f64; 3]]) -> bool,
{
    let objective = |candidate: &[[f64; 3]]| lbfgs_energy(model, candidate, cutoff);
    armijo_backtracking(
        coords,
        direction,
        (energy, slope),
        1.0,
        &objective,
        accept,
    )
    .or_else(|| {
        // A stale curvature history or a geometry constraint can exhaust the
        // quasi-Newton search. Restart with a bounded steepest-descent step
        // and backtrack it too: a fixed displacement can overshoot near a
        // stationary point.
        history.clear();
        let fallback_direction = gradient
            .iter()
            .map(|g| [-g[0], -g[1], -g[2]])
            .collect::<CoordVec>();
        let fallback_slope = -gradient.iter().map(|g| dot3(*g, *g)).sum::<f64>();
        armijo_backtracking(
            coords,
            &fallback_direction,
            (energy, fallback_slope),
            0.01 / max_gradient,
            &objective,
            accept,
        )
    })
}

/// Backtrack until the candidate decreases the objective and satisfies the
/// caller's geometry constraint. Exhaustion preserves the last accepted
/// state; no untested fixed-size step is committed.
fn armijo_backtracking<E, A>(
    coords: &[[f64; 3]],
    direction: &[[f64; 3]],
    objective_and_slope: (f64, f64),
    initial_step: f64,
    energy: &E,
    accept: &A,
) -> Option<(CoordVec, f64)>
where
    E: Fn(&[[f64; 3]]) -> f64,
    A: Fn(&[[f64; 3]]) -> bool,
{
    let (objective, slope) = objective_and_slope;
    if !objective.is_finite()
        || !slope.is_finite()
        || slope >= 0.0
        || !initial_step.is_finite()
        || initial_step <= 0.0
    {
        return None;
    }
    let mut alpha = initial_step;
    for _ in 0..40 {
        let trial = coords
            .iter()
            .zip(direction)
            .map(|(c, p)| [c[0] + alpha * p[0], c[1] + alpha * p[1], c[2] + alpha * p[2]])
            .collect::<CoordVec>();
        if trial.iter().flatten().all(|value| value.is_finite()) {
            let value = energy(&trial);
            if value.is_finite() && value <= objective + 1e-4 * alpha * slope && accept(&trial) {
                return Some((trial, value));
            }
        }
        alpha *= 0.5;
    }
    None
}

fn update_lbfgs_history(
    history: &mut LbfgsHistory,
    old_coords: &[[f64; 3]],
    new_coords: &[[f64; 3]],
    old_gradient: &[[f64; 3]],
    new_gradient: &[[f64; 3]],
) {
    let s = vector_difference(new_coords, old_coords);
    let y = vector_difference(new_gradient, old_gradient);
    let ys = y.iter().zip(&s).map(|(yi, si)| dot3(*yi, *si)).sum::<f64>();
    if ys <= 1e-10 {
        return;
    }
    if history.len() >= LBFGS_HISTORY_SIZE {
        history.pop_front();
    }
    history.push_back((s, y, 1.0 / ys));
}

fn vector_difference(left: &[[f64; 3]], right: &[[f64; 3]]) -> CoordVec {
    left.iter()
        .zip(right)
        .map(|(l, r)| [l[0] - r[0], l[1] - r[1], l[2] - r[2]])
        .collect()
}

fn coordinate_rmsd(coords: &[[f64; 3]], initial: &[[f64; 3]]) -> f64 {
    let sum = coords
        .iter()
        .zip(initial)
        .map(|(c, start)| {
            let delta = [c[0] - start[0], c[1] - start[1], c[2] - start[2]];
            dot3(delta, delta)
        })
        .sum::<f64>();
    (sum / coords.len() as f64).sqrt()
}

/// Largest molecule (atoms, hydrogens included) for which
/// [`Mmff94EnergyModel::minimize_bfgs`] keeps a dense inverse Hessian
/// (1,500 x 1,500 doubles, 18 MB); larger molecules use L-BFGS.
pub const BFGS_DENSE_MAX_ATOMS: usize = 500;

/// Longest move (Å) along any coordinate that one line search of
/// [`Mmff94EnergyModel::minimize_bfgs`] starts from under
/// [`Mmff94Convergence::MaxGradient`]. RDKit's first steps can move an atom
/// several Å, which on the A6 corpus carried strained stereocentres across
/// their plane; 0.2-0.5 Å give the same successes and nearly the same
/// converged rows, 0.3 is used.
pub const BFGS_MAX_COORDINATE_STEP: f64 = 0.3;

/// Convergence test for [`Mmff94EnergyModel::minimize_bfgs`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Mmff94Convergence {
    /// Every gradient component below this value (kcal/mol/Å). chematic's
    /// L-BFGS and the 3D pipeline use 1e-4.
    MaxGradient(f64),
    /// RDKit's `BFGSOpt` tests: the step test (`TOLX` = 1.2e-7, relative
    /// to coordinates) or `max |g_i| max(|x_i|, 1) / max(|E| s, 1) <
    /// grad_tol`, with `g` RDKit's scaled gradient and `s` its scale factor.
    /// `MMFFOptimizeMolecule` passes `grad_tol = 1e-4`. A converged run can
    /// keep residual forces well above 1e-4 kcal/mol/Å.
    Rdkit {
        /// RDKit's `forceTol` (1e-4 by default).
        grad_tol: f64,
    },
}

impl Default for Mmff94Convergence {
    fn default() -> Self {
        Mmff94Convergence::MaxGradient(1e-4)
    }
}

/// RDKit `ForceFieldsHelper::calcGradient`: scales the gradient by 0.1 and,
/// when a component still exceeds 10, halves the factor until none does.
/// Returns (RDKit's returned scale, largest unscaled component).
fn rdkit_scaled_gradient(model: &Mmff94EnergyModel, x: &[[f64; 3]], out: &mut [f64]) -> (f64, f64) {
    let g = model.bounded_analytic_gradient(x);
    let mut raw_max = 0.0_f64;
    let mut scaled_max = 0.0_f64;
    for (o, v) in out.iter_mut().zip(g.iter().flat_map(|v| v.iter())) {
        raw_max = raw_max.max(v.abs());
        *o = v * 0.1;
        scaled_max = scaled_max.max(o.abs());
    }
    let mut scale = 0.1;
    if scaled_max > 10.0 {
        while scaled_max * scale > 10.0 {
            scale *= 0.5;
        }
        for o in out.iter_mut() {
            *o *= scale;
        }
    }
    (scale, raw_max)
}

enum BfgsLineSearch {
    /// Sufficient decrease met by an admissible point.
    Accepted(f64),
    /// RDKit's `resCode = 1`: the step fell below its position-scaled
    /// minimum. Carries the last trial and whether it was admissible.
    StepTooSmall(f64, bool),
    /// The direction is not a descent direction.
    BadDirection,
}

/// RDKit `BFGSOpt::linearSearch` (Numerical Recipes `lnsrch`): Armijo
/// backtracking with quadratic then cubic interpolation. A point the caller's
/// predicate rejects is treated as failing the decrease test and halves the
/// step without entering the interpolation.
#[allow(clippy::too_many_arguments)]
fn bfgs_line_search<F>(
    model: &Mmff94EnergyModel,
    x3: &[[f64; 3]],
    f0: f64,
    g: &[f64],
    dir: &mut [f64],
    max_step: f64,
    max_component: Option<f64>,
    out3: &mut [[f64; 3]],
    accept: &F,
) -> BfgsLineSearch
where
    F: Fn(&[[f64; 3]]) -> bool,
{
    const FUNCTOL: f64 = 1e-4;
    const MOVETOL: f64 = 1e-7;
    const MAX_ITER: usize = 1000;
    let x = x3.as_flattened();
    let len = dir.iter().map(|d| d * d).sum::<f64>().sqrt();
    if len > max_step {
        for d in dir.iter_mut() {
            *d *= max_step / len;
        }
    }
    if let Some(cap) = max_component {
        let big = dir.iter().fold(0.0_f64, |m, d| m.max(d.abs()));
        if big > cap {
            for d in dir.iter_mut() {
                *d *= cap / big;
            }
        }
    }
    let slope: f64 = dir.iter().zip(g).map(|(d, gi)| d * gi).sum();
    // Also rejects a NaN slope.
    if slope.partial_cmp(&0.0) != Some(std::cmp::Ordering::Less) {
        return BfgsLineSearch::BadDirection;
    }
    let test = dir
        .iter()
        .zip(x)
        .map(|(d, xi)| d.abs() / xi.abs().max(1.0))
        .fold(0.0_f64, f64::max);
    let lambda_min = MOVETOL / test;
    let mut lambda = 1.0_f64;
    let mut lambda2 = 0.0_f64;
    let mut f2 = 0.0_f64;
    let mut have_previous = false;
    let mut last = (f0, false);
    for _ in 0..MAX_ITER {
        if lambda < lambda_min {
            return BfgsLineSearch::StepTooSmall(last.0, last.1);
        }
        for ((o, xi), d) in out3.as_flattened_mut().iter_mut().zip(x).zip(dir.iter()) {
            *o = xi + lambda * d;
        }
        let f = model.energy(out3);
        let admissible = f.is_finite() && accept(out3);
        last = (f, admissible);
        if admissible && f - f0 <= FUNCTOL * lambda * slope {
            return BfgsLineSearch::Accepted(f);
        }
        let next = if !admissible {
            0.5 * lambda
        } else if !have_previous {
            -slope / (2.0 * (f - f0 - slope))
        } else {
            let rhs1 = f - f0 - lambda * slope;
            let rhs2 = f2 - f0 - lambda2 * slope;
            let a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
            let b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2))
                / (lambda - lambda2);
            let t = if a == 0.0 {
                -slope / (2.0 * b)
            } else {
                let disc = b * b - 3.0 * a * slope;
                if disc < 0.0 {
                    0.5 * lambda
                } else if b <= 0.0 {
                    (-b + disc.sqrt()) / (3.0 * a)
                } else {
                    -slope / (b + disc.sqrt())
                }
            };
            t.min(0.5 * lambda)
        };
        if admissible {
            lambda2 = lambda;
            f2 = f;
            have_previous = true;
        }
        // NaN-safe: a non-finite interpolant falls back to the 0.1 floor.
        lambda = if next.is_nan() {
            0.1 * lambda
        } else {
            next.max(0.1 * lambda)
        };
    }
    out3.copy_from_slice(x3);
    BfgsLineSearch::StepTooSmall(f0, false)
}

fn minimize_mmff94_bfgs_prepared<F>(
    model: &Mmff94EnergyModel,
    coords: &mut [[f64; 3]],
    max_iter: usize,
    convergence: Mmff94Convergence,
    accept: F,
) -> Result<MinimizeResult, MinimizerError>
where
    F: Fn(&[[f64; 3]]) -> bool,
{
    let n = model.mmff_mol.atom_count();
    if n <= 1 {
        return Ok(MinimizeResult {
            energy: 0.0,
            rmsd: 0.0,
            converged: true,
            iterations: 0,
            termination: Mmff94TerminationReason::GradientConverged,
        });
    }
    if n > BFGS_DENSE_MAX_ATOMS {
        return minimize_mmff94_lbfgs_prepared_with_mode(
            model, coords, max_iter, true, false, accept,
        );
    }
    let initial = coords.to_vec();
    let mut state = BfgsState::new(model, coords, convergence);
    let mut iters = 0usize;
    let mut outcome = if state.initial_gradient_converged(convergence) {
        BfgsIteration::Converged
    } else {
        BfgsIteration::Continue
    };
    while outcome == BfgsIteration::Continue && iters < max_iter {
        iters += 1;
        outcome = state.iterate(model, convergence, &accept);
    }
    coords.copy_from_slice(&state.coords);
    let converged = outcome == BfgsIteration::Converged;
    let termination = match outcome {
        BfgsIteration::Converged => Mmff94TerminationReason::GradientConverged,
        BfgsIteration::Rejected => Mmff94TerminationReason::ConstraintRejectedFallback,
        BfgsIteration::Continue => Mmff94TerminationReason::IterationLimit,
    };
    Ok(MinimizeResult {
        energy: state.energy,
        rmsd: coordinate_rmsd(coords, &initial),
        converged,
        iterations: iters,
        termination,
    })
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BfgsIteration {
    Continue,
    Converged,
    Rejected,
}

struct BfgsState {
    dim: usize,
    coords: CoordVec,
    gradient: Vec<f64>,
    gradient_delta: Vec<f64>,
    hessian_gradient: Vec<f64>,
    direction: Vec<f64>,
    trial: CoordVec,
    inverse_hessian: Vec<f64>,
    hessian_is_identity: bool,
    energy: f64,
    raw_gradient_max: f64,
    max_step: f64,
    max_component: Option<f64>,
}

impl BfgsState {
    fn new(model: &Mmff94EnergyModel, coords: &[[f64; 3]], convergence: Mmff94Convergence) -> Self {
        let dim = 3 * coords.len();
        let mut state = Self {
            dim,
            coords: coords.to_vec(),
            gradient: vec![0.0; dim],
            gradient_delta: vec![0.0; dim],
            hessian_gradient: vec![0.0; dim],
            direction: vec![0.0; dim],
            trial: vec![[0.0; 3]; coords.len()],
            inverse_hessian: vec![0.0; dim * dim],
            hessian_is_identity: true,
            energy: model.energy(coords),
            raw_gradient_max: 0.0,
            max_step: 0.0,
            max_component: matches!(convergence, Mmff94Convergence::MaxGradient(_))
                .then_some(BFGS_MAX_COORDINATE_STEP),
        };
        state.reset_inverse_hessian();
        (_, state.raw_gradient_max) =
            rdkit_scaled_gradient(model, &state.coords, &mut state.gradient);
        state.set_steepest_descent_direction();
        let norm = state
            .coords
            .as_flattened()
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        state.max_step = 100.0 * norm.max(dim as f64);
        state
    }

    fn iterate<F>(
        &mut self,
        model: &Mmff94EnergyModel,
        convergence: Mmff94Convergence,
        accept: &F,
    ) -> BfgsIteration
    where
        F: Fn(&[[f64; 3]]) -> bool,
    {
        const TOLX: f64 = 1.2e-7;
        let search = bfgs_line_search(
            model,
            &self.coords,
            self.energy,
            &self.gradient,
            &mut self.direction,
            self.max_step,
            self.max_component,
            &mut self.trial,
            accept,
        );
        let new_energy = match search {
            BfgsLineSearch::Accepted(energy) => energy,
            BfgsLineSearch::StepTooSmall(energy, true)
                if matches!(convergence, Mmff94Convergence::Rdkit { .. }) =>
            {
                energy
            }
            BfgsLineSearch::BadDirection | BfgsLineSearch::StepTooSmall(..)
                if !self.hessian_is_identity =>
            {
                self.restart();
                return BfgsIteration::Continue;
            }
            BfgsLineSearch::StepTooSmall(..)
                if matches!(convergence, Mmff94Convergence::Rdkit { .. }) =>
            {
                return BfgsIteration::Converged;
            }
            BfgsLineSearch::BadDirection | BfgsLineSearch::StepTooSmall(..) => {
                return BfgsIteration::Rejected;
            }
        };
        self.energy = new_energy;
        let step_test = self.accept_trial();
        if matches!(convergence, Mmff94Convergence::Rdkit { .. }) && step_test < TOLX {
            return BfgsIteration::Converged;
        }
        let (scale, raw) = rdkit_scaled_gradient(model, &self.coords, &mut self.gradient);
        self.raw_gradient_max = raw;
        if self.gradient_converged(convergence, scale) {
            return BfgsIteration::Converged;
        }
        self.update_inverse_hessian();
        self.update_direction();
        BfgsIteration::Continue
    }

    fn accept_trial(&mut self) -> f64 {
        let mut step_test = 0.0_f64;
        let coords = self.coords.as_flattened_mut();
        let trial = self.trial.as_flattened();
        for i in 0..self.dim {
            self.direction[i] = trial[i] - coords[i];
            coords[i] = trial[i];
            step_test = step_test.max(self.direction[i].abs() / coords[i].abs().max(1.0));
            self.gradient_delta[i] = self.gradient[i];
        }
        step_test
    }

    fn gradient_converged(&self, convergence: Mmff94Convergence, scale: f64) -> bool {
        match convergence {
            Mmff94Convergence::MaxGradient(tolerance) => self.raw_gradient_max < tolerance,
            Mmff94Convergence::Rdkit { grad_tol } => {
                let denominator = (self.energy.abs() * scale).max(1.0);
                self.gradient
                    .iter()
                    .zip(self.coords.as_flattened())
                    .map(|(g, x)| g.abs() * x.abs().max(1.0))
                    .fold(0.0_f64, f64::max)
                    / denominator
                    < grad_tol
            }
        }
    }

    fn initial_gradient_converged(&self, convergence: Mmff94Convergence) -> bool {
        matches!(convergence, Mmff94Convergence::MaxGradient(t) if self.raw_gradient_max < t)
    }

    fn update_inverse_hessian(&mut self) {
        const EPS: f64 = 3e-8;
        for i in 0..self.dim {
            self.gradient_delta[i] = self.gradient[i] - self.gradient_delta[i];
        }
        let (mut fac, mut fae, mut sum_dg, mut sum_step) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..self.dim {
            let row = &self.inverse_hessian[i * self.dim..(i + 1) * self.dim];
            self.hessian_gradient[i] = row
                .iter()
                .zip(&self.gradient_delta)
                .map(|(h, d)| h * d)
                .sum();
            fac += self.gradient_delta[i] * self.direction[i];
            fae += self.gradient_delta[i] * self.hessian_gradient[i];
            sum_dg += self.gradient_delta[i] * self.gradient_delta[i];
            sum_step += self.direction[i] * self.direction[i];
        }
        if fac <= (EPS * sum_dg * sum_step).sqrt() {
            return;
        }
        let inv_fac = 1.0 / fac;
        let inv_fae = 1.0 / fae;
        for i in 0..self.dim {
            self.gradient_delta[i] =
                inv_fac * self.direction[i] - inv_fae * self.hessian_gradient[i];
        }
        for i in 0..self.dim {
            let (step_i, hdg_i, delta_i) = (
                inv_fac * self.direction[i],
                inv_fae * self.hessian_gradient[i],
                fae * self.gradient_delta[i],
            );
            for j in i..self.dim {
                let value = self.inverse_hessian[i * self.dim + j]
                    + step_i * self.direction[j]
                    - hdg_i * self.hessian_gradient[j]
                    + delta_i * self.gradient_delta[j];
                self.inverse_hessian[i * self.dim + j] = value;
                self.inverse_hessian[j * self.dim + i] = value;
            }
        }
        self.hessian_is_identity = false;
    }

    fn update_direction(&mut self) {
        for i in 0..self.dim {
            let row = &self.inverse_hessian[i * self.dim..(i + 1) * self.dim];
            self.direction[i] = -row
                .iter()
                .zip(&self.gradient)
                .map(|(h, g)| h * g)
                .sum::<f64>();
        }
    }

    fn restart(&mut self) {
        self.reset_inverse_hessian();
        self.hessian_is_identity = true;
        self.set_steepest_descent_direction();
    }

    fn reset_inverse_hessian(&mut self) {
        self.inverse_hessian.fill(0.0);
        for i in 0..self.dim {
            self.inverse_hessian[i * self.dim + i] = 1.0;
        }
    }

    fn set_steepest_descent_direction(&mut self) {
        for (direction, gradient) in self.direction.iter_mut().zip(&self.gradient) {
            *direction = -*gradient;
        }
    }
}

/// L-BFGS two-loop recursion: compute search direction p = -H_k × g.
fn lbfgs_direction(g: &[[f64; 3]], history: &LbfgsHistory) -> Vec<[f64; 3]> {
    let m = history.len();

    if m == 0 {
        // No history: steepest descent direction
        return g.iter().map(|gi| [-gi[0], -gi[1], -gi[2]]).collect();
    }

    let mut q: Vec<[f64; 3]> = g.to_vec();
    let mut alphas = [0.0_f64; LBFGS_HISTORY_SIZE];

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

    // p = -H_k g = -q. Reuse q's allocation; this runs once per optimizer
    // iteration and avoids allocating/copying a second coordinate-sized
    // direction buffer without changing operation order.
    for qi in &mut q {
        qi[0] = -qi[0];
        qi[1] = -qi[1];
        qi[2] = -qi[2];
    }
    q
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
    if n < 16 {
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
        return grad;
    }

    grad.par_iter_mut().enumerate().for_each(|(i, atom_grad)| {
        let mut work = coords.to_vec();
        for (axis, component) in atom_grad.iter_mut().enumerate() {
            work[i][axis] += delta;
            let ep = total_energy(mol, &work, types, charges, rings);
            work[i][axis] -= 2.0 * delta;
            let em = total_energy(mol, &work, types, charges, rings);
            work[i][axis] += delta;
            *component = (ep - em) / (2.0 * delta);
        }
    });
    grad
}

fn compute_gradient_prepared(
    model: &Mmff94EnergyModel,
    coords: &[[f64; 3]],
    delta: f64,
) -> Vec<[f64; 3]> {
    let n = coords.len();
    let mut grad = vec![[0.0_f64; 3]; n];
    // Large molecules benefit from independent atom probes. Keep small
    // molecules on the allocation-light sequential path to avoid rayon
    // scheduling and one coordinate-buffer clone per worker.
    if n < 16 {
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
        return grad;
    }

    grad.par_iter_mut().enumerate().for_each(|(i, atom_grad)| {
        let mut work = coords.to_vec();
        for (axis, component) in atom_grad.iter_mut().enumerate() {
            work[i][axis] += delta;
            let ep = model.energy(&work);
            work[i][axis] -= 2.0 * delta;
            let em = model.energy(&work);
            work[i][axis] += delta;
            *component = (ep - em) / (2.0 * delta);
        }
    });
    grad
}
