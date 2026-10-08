//! `BFGSOpt::minimize` / `BFGSOpt::linearSearch` (RDKit 2026.03,
//! Numerical Recipes `dfpmin`/`lnsrch`) operation for operation.

const FUNCTOL: f64 = 1e-4;
const MOVETOL: f64 = 1e-7;
const EPS: f64 = 3e-8;
const TOLX: f64 = 4.0 * EPS;
const MAXSTEP: f64 = 100.0;

/// RDKit's `CHECK_INVARIANT(status >= 0, "bad direction in linearSearch")`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BadDirection;

#[allow(clippy::too_many_arguments)]
fn linear_search(
    old_pt: &[f64],
    old_val: f64,
    grad: &[f64],
    dir: &mut [f64],
    new_pt: &mut [f64],
    func: &mut impl FnMut(&[f64]) -> f64,
    max_step: f64,
) -> (i32, f64) {
    const MAX_ITER_LINEAR_SEARCH: u32 = 1000;
    let dim = old_pt.len();
    let mut new_val = 0.0;
    let mut res_code = -1;
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
        return (res_code, new_val);
    }
    let mut test = 0.0;
    for i in 0..dim {
        let temp = dir[i].abs() / old_pt[i].abs().max(1.0);
        if temp > test {
            test = temp;
        }
    }
    let lambda_min = MOVETOL / test;
    let mut lambda: f64 = 1.0;
    let mut lambda2 = 0.0;
    let mut val2 = 0.0;
    let mut tmp_lambda;
    let mut it = 0;
    while it < MAX_ITER_LINEAR_SEARCH {
        if lambda < lambda_min {
            res_code = 1;
            break;
        }
        for i in 0..dim {
            new_pt[i] = old_pt[i] + lambda * dir[i];
        }
        new_val = func(new_pt);
        if new_val - old_val <= FUNCTOL * lambda * slope {
            return (0, new_val);
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
        lambda = tmp_lambda.max(0.1 * lambda);
        it += 1;
    }
    new_pt.copy_from_slice(old_pt);
    (res_code, new_val)
}

/// `BFGSOpt::minimize(dim, pos, gradTol, ..., maxIts)`; `grad_fn` fills the
/// gradient and returns its scale factor (RDKit's `calcGradient` functor).
/// Returns 0 (converged) or 1 (more iterations needed).
pub(crate) fn minimize(
    pos: &mut [f64],
    grad_tol: f64,
    max_its: u32,
    mut func: impl FnMut(&[f64]) -> f64,
    mut grad_fn: impl FnMut(&[f64], &mut [f64]) -> f64,
) -> Result<i32, BadDirection> {
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
    let max_step = MAXSTEP * sum.sqrt().max(dim as f64);
    for _iter in 1..=max_its {
        let (status, func_val) =
            linear_search(pos, fp, &grad, &mut xi, &mut new_pos, &mut func, max_step);
        if status < 0 {
            return Err(BadDirection);
        }
        fp = func_val;
        let mut test = 0.0;
        for i in 0..dim {
            xi[i] = new_pos[i] - pos[i];
            pos[i] = new_pos[i];
            let temp = xi[i].abs() / pos[i].abs().max(1.0);
            if temp > test {
                test = temp;
            }
            dgrad[i] = grad[i];
        }
        if test < TOLX {
            return Ok(0);
        }
        let grad_scale = grad_fn(pos, &mut grad);
        test = 0.0;
        let term = (func_val * grad_scale).max(1.0);
        for i in 0..dim {
            let temp = grad[i].abs() * pos[i].abs().max(1.0);
            test = f64::max(test, temp);
            dgrad[i] = grad[i] - dgrad[i];
        }
        test /= term;
        if test < grad_tol {
            return Ok(0);
        }
        let (mut fac, mut fae, mut sum_dgrad, mut sum_xi) = (0.0, 0.0, 0.0, 0.0);
        for i in 0..dim {
            let row = &inv_hessian[i * dim..(i + 1) * dim];
            let mut h = 0.0;
            for j in 0..dim {
                h += row[j] * dgrad[j];
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
            for j in 0..dim {
                p -= row[j] * grad[j];
            }
            xi[i] = p;
        }
    }
    Ok(1)
}

/// RDKit's `ForceFieldsHelper::calcGradient` scaling (in place); returns the
/// scale factor.
pub(crate) fn scale_gradient(grad: &mut [f64]) -> f64 {
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
