//! RDKit's BFGS minimizer (`Code/Numerics/Optimizer/BFGSOpt.h`) and the
//! gradient scaling of `ForceFieldsHelper::calcGradient`, shared by the
//! RDKit-exact force fields ([`crate::rdkit_mmff`], [`crate::rdkit_uff`]) and
//! usable by any other port of an RDKit minimization (e.g. ETKDG's
//! distance-geometry and 3D force fields).

/// `std::max(a, b)`: `b` only when `a < b`.
fn smax(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

/// `ForceFieldsHelper::calcGradient`'s scaling of a freshly computed raw
/// gradient: multiply by 0.1, then, if any component still exceeds 10 in
/// magnitude, keep halving the scale until it does not and apply it again.
/// Returns the final scale factor (the value `calcGradient` returns).
pub fn scale_gradient(grad: &mut [f64]) -> f64 {
    let mut max_grad = -1e8_f64;
    let mut grad_scale = 0.1;
    for g in grad.iter_mut() {
        *g *= grad_scale;
        if g.abs() > max_grad {
            max_grad = g.abs();
        }
    }
    if max_grad > 10.0 {
        while max_grad * grad_scale > 10.0 {
            grad_scale *= 0.5;
        }
        for g in grad.iter_mut() {
            *g *= grad_scale;
        }
    }
    grad_scale
}

/// RDKit's `BFGSOpt::minimize` (RDKit 2026.03, Numerical Recipes
/// `dfpmin`) with `BFGSOpt::linearSearch`, constant for constant and
/// operation for operation, so a port of an RDKit force field minimized
/// with it follows RDKit's trajectory bit for bit.
///
/// * `pos`: flat coordinates (`3 * n_atoms`), updated in place.
/// * `grad_tol`: RDKit's `forceTol` (`ForceField::minimize` default 1e-4).
/// * `max_its`: maximum iterations (`maxIts`).
/// * `func`: energy at a position.
/// * `grad_fn`: fills the gradient at a position and returns the scale
///   applied to it. `ForceField::minimize` passes RDKit's
///   `ForceFieldsHelper::calcGradient`, i.e. the raw gradient followed by
///   [`scale_gradient`]; the returned scale feeds the convergence test.
///
/// Returns RDKit's status: 0 converged, 1 `max_its` reached, -1 when the
/// line search fails (RDKit throws "bad direction in linearSearch" there).
///
/// `ForceField::minimize(maxIts, forceTol, energyTol)` is this call with
/// `grad_fn = |p, g| { g.fill(0.0); raw_grad(p, g); scale_gradient(g) }`,
/// after RDKit's early return of 0 when the force field has no
/// contributions. (`energyTol` is unused by RDKit's BFGS.)
pub fn bfgs_minimize(
    pos: &mut [f64],
    grad_tol: f64,
    max_its: u32,
    func: impl Fn(&[f64]) -> f64,
    grad_fn: impl Fn(&[f64], &mut [f64]) -> f64,
) -> i32 {
    bfgs_minimize_detailed(pos, grad_tol, max_its, func, grad_fn).status
}

/// Outcome of [`bfgs_minimize_detailed`].
#[derive(Clone, Debug, PartialEq)]
pub struct BfgsOutcome {
    /// RDKit's status: 0 converged, 1 `max_its` reached, -1 line-search
    /// failure.
    pub status: i32,
    /// The point of the last energy evaluation when it is not the returned
    /// position. RDKit's `ForceField` caches interatomic distances per
    /// evaluation and its `calcEnergy()` (no positions argument) reuses that
    /// cache, so after such a minimization RDKit's reported energy takes
    /// its distances from this point and everything else from the final
    /// position. This happens only when the last line search gives up
    /// (step below `MOVETOL`): RDKit then restores the previous point and
    /// reports convergence.
    pub stale_distance_point: Option<Vec<f64>>,
}

/// [`bfgs_minimize`], also reporting whether RDKit's distance cache is
/// stale afterwards (see [`BfgsOutcome::stale_distance_point`]).
pub fn bfgs_minimize_detailed(
    pos: &mut [f64],
    grad_tol: f64,
    max_its: u32,
    func: impl Fn(&[f64]) -> f64,
    grad_fn: impl Fn(&[f64], &mut [f64]) -> f64,
) -> BfgsOutcome {
    let done = |status: i32, stale: Option<Vec<f64>>| BfgsOutcome {
        status,
        stale_distance_point: stale,
    };
    const EPS: f64 = 3e-8;
    const TOLX: f64 = 4.0 * EPS;
    const MAXSTEP: f64 = 100.0;
    let dim = pos.len();
    let mut grad = vec![0.0; dim];
    let mut dgrad = vec![0.0; dim];
    let mut hess_dgrad = vec![0.0; dim];
    let mut xi = vec![0.0; dim];
    let mut inv_hessian = vec![0.0; dim * dim];
    let mut new_pos = vec![0.0; dim];
    let mut fp = func(pos);
    grad_fn(pos, &mut grad);
    let mut sum = 0.0;
    for i in 0..dim {
        inv_hessian[i * dim + i] = 1.0;
        xi[i] = -grad[i];
        sum += pos[i] * pos[i];
    }
    let max_step = MAXSTEP * smax(sum.sqrt(), dim as f64);
    for _iter in 1..=max_its {
        let (status, func_val, stale) =
            linear_search(pos, fp, &grad, &mut xi, &mut new_pos, &func, max_step);
        if status < 0 {
            // RDKit's CHECK_INVARIANT("bad direction in linearSearch").
            return done(-1, None);
        }
        fp = func_val;
        let mut test = 0.0_f64;
        for i in 0..dim {
            xi[i] = new_pos[i] - pos[i];
            pos[i] = new_pos[i];
            let temp = xi[i].abs() / smax(pos[i].abs(), 1.0);
            if temp > test {
                test = temp;
            }
            dgrad[i] = grad[i];
        }
        if test < TOLX {
            return done(0, stale);
        }
        let grad_scale = grad_fn(pos, &mut grad);
        test = 0.0;
        let term = smax(func_val * grad_scale, 1.0);
        for i in 0..dim {
            let temp = grad[i].abs() * smax(pos[i].abs(), 1.0);
            test = smax(test, temp);
            dgrad[i] = grad[i] - dgrad[i];
        }
        test /= term;
        if test < grad_tol {
            return done(0, None);
        }
        let (mut fac, mut fae, mut sum_dgrad, mut sum_xi) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..dim {
            let row = &inv_hessian[i * dim..(i + 1) * dim];
            let mut h = 0.0;
            for (iv, dg) in row.iter().zip(&dgrad) {
                h += iv * dg;
            }
            hess_dgrad[i] = h;
            fac += dgrad[i] * xi[i];
            fae += dgrad[i] * hess_dgrad[i];
            sum_dgrad += dgrad[i] * dgrad[i];
            sum_xi += xi[i] * xi[i];
        }
        if fac > (EPS * sum_dgrad * sum_xi).sqrt() {
            fac = 1.0 / fac;
            let fad = 1.0 / fae;
            for i in 0..dim {
                dgrad[i] = fac * xi[i] - fad * hess_dgrad[i];
            }
            for i in 0..dim {
                let pxi = fac * xi[i];
                let hdgi = fad * hess_dgrad[i];
                let dgi = fae * dgrad[i];
                for j in i..dim {
                    inv_hessian[i * dim + j] += pxi * xi[j] - hdgi * hess_dgrad[j] + dgi * dgrad[j];
                    inv_hessian[j * dim + i] = inv_hessian[i * dim + j];
                }
            }
        }
        for i in 0..dim {
            let row = &inv_hessian[i * dim..(i + 1) * dim];
            let mut p = 0.0;
            for (iv, g) in row.iter().zip(&grad) {
                p -= iv * g;
            }
            xi[i] = p;
        }
    }
    done(1, None)
}

/// `BFGSOpt::linearSearch`: (resCode, newVal, last trial point when the
/// search gave up and restored `old_pt`).
fn linear_search(
    old_pt: &[f64],
    old_val: f64,
    grad: &[f64],
    dir: &mut [f64],
    new_pt: &mut [f64],
    func: &impl Fn(&[f64]) -> f64,
    max_step: f64,
) -> (i32, f64, Option<Vec<f64>>) {
    const FUNCTOL: f64 = 1e-4;
    const MOVETOL: f64 = 1e-7;
    const MAX_ITER_LINEAR_SEARCH: u32 = 1000;
    let dim = old_pt.len();
    let mut new_val = 0.0;
    let mut sum = 0.0;
    for d in dir.iter() {
        sum += d * d;
    }
    sum = sum.sqrt();
    if sum > max_step {
        for d in dir.iter_mut() {
            *d *= max_step / sum;
        }
    }
    let mut slope = 0.0;
    for i in 0..dim {
        slope += dir[i] * grad[i];
    }
    if slope >= 0.0 {
        return (-1, new_val, None);
    }
    let mut test = 0.0;
    for i in 0..dim {
        let temp = dir[i].abs() / smax(old_pt[i].abs(), 1.0);
        if temp > test {
            test = temp;
        }
    }
    let lambda_min = MOVETOL / test;
    let mut lambda = 1.0_f64;
    let mut lambda2 = 0.0_f64;
    let mut val2 = 0.0_f64;
    let mut tmp_lambda;
    let mut it = 0;
    while it < MAX_ITER_LINEAR_SEARCH {
        if lambda < lambda_min {
            // RDKit breaks out and restores the old point (keeping the last
            // trial's value as newVal).
            let trial = (it > 0).then(|| new_pt.to_vec());
            new_pt.copy_from_slice(old_pt);
            return (1, new_val, trial);
        }
        for i in 0..dim {
            new_pt[i] = old_pt[i] + lambda * dir[i];
        }
        new_val = func(new_pt);
        if new_val - old_val <= FUNCTOL * lambda * slope {
            return (0, new_val, None);
        }
        if it == 0 {
            tmp_lambda = -slope / (2.0 * (new_val - old_val - slope));
        } else {
            let rhs1 = new_val - old_val - lambda * slope;
            let rhs2 = val2 - old_val - lambda2 * slope;
            let a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / (lambda - lambda2);
            let b = (-lambda2 * rhs1 / (lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2))
                / (lambda - lambda2);
            if a == 0.0 {
                tmp_lambda = -slope / (2.0 * b);
            } else {
                let disc = b * b - 3.0 * a * slope;
                if disc < 0.0 {
                    tmp_lambda = 0.5 * lambda;
                } else if b <= 0.0 {
                    tmp_lambda = (-b + disc.sqrt()) / (3.0 * a);
                } else {
                    tmp_lambda = -slope / (b + disc.sqrt());
                }
            }
            if tmp_lambda > 0.5 * lambda {
                tmp_lambda = 0.5 * lambda;
            }
        }
        lambda2 = lambda;
        val2 = new_val;
        lambda = smax(tmp_lambda, 0.1 * lambda);
        it += 1;
    }
    new_pt.copy_from_slice(old_pt);
    (-1, new_val, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_converges() {
        let mut pos = vec![3.0, -2.0];
        let status = bfgs_minimize(
            &mut pos,
            1e-4,
            200,
            |p| (p[0] - 1.0).powi(2) + 2.0 * (p[1] + 0.5).powi(2),
            |p, g| {
                g[0] = 2.0 * (p[0] - 1.0);
                g[1] = 4.0 * (p[1] + 0.5);
                1.0
            },
        );
        assert_eq!(status, 0);
        assert!((pos[0] - 1.0).abs() < 1e-3 && (pos[1] + 0.5).abs() < 1e-3);
    }

    /// RDKit's line search gives up below `MOVETOL`, restores the start
    /// point and reports convergence; the last trial point is reported as
    /// the stale distance-cache point.
    #[test]
    fn line_search_give_up_restores_point() {
        let start = vec![0.5];
        let mut pos = start.clone();
        let out = bfgs_minimize_detailed(
            &mut pos,
            1e-4,
            200,
            |p| if p[0] == 0.5 { 1.0 } else { 2.0 },
            |_, g| {
                g[0] = 1.0;
                1.0
            },
        );
        assert_eq!(out.status, 0);
        assert_eq!(pos, start);
        let stale = out.stale_distance_point.expect("trial point");
        assert!(stale[0] < 0.5 && stale[0] > 0.5 - 1e-6);
    }
}
