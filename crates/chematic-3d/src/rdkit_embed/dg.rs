//! RDKit `DistGeom`: bounds matrix, triangle smoothing, random distance
//! matrices, metric-matrix initial coordinates and the distance-geometry
//! error function (`DistViolation`, `ChiralViolation`, `FourthDim`).

use super::rng::MinstdRand;

/// `DistGeom::BoundsMatrix`: upper bounds in the upper triangle, lower
/// bounds in the lower triangle of a row-major `n x n` matrix.
#[derive(Clone, Debug)]
pub struct BoundsMatrix {
    n: usize,
    data: Vec<f64>,
}

impl BoundsMatrix {
    /// All zeros (`BoundsMatrix(N)`).
    pub fn new(n: usize) -> Self {
        BoundsMatrix {
            n,
            data: vec![0.0; n * n],
        }
    }

    pub fn num_rows(&self) -> usize {
        self.n
    }

    #[inline]
    pub(crate) fn val(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    #[inline]
    pub(crate) fn set_val(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.n + j] = v;
    }

    pub fn upper(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.val(i, j)
        } else {
            self.val(j, i)
        }
    }

    pub fn lower(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.val(j, i)
        } else {
            self.val(i, j)
        }
    }

    pub fn set_upper(&mut self, i: usize, j: usize, v: f64) {
        if i < j {
            self.set_val(i, j, v)
        } else {
            self.set_val(j, i, v)
        }
    }

    pub fn set_lower(&mut self, i: usize, j: usize, v: f64) {
        if i < j {
            self.set_val(j, i, v)
        } else {
            self.set_val(i, j, v)
        }
    }

    #[allow(dead_code)] // used by the bounds-matrix builder (next step)
    pub(crate) fn set_upper_if_better(&mut self, i: usize, j: usize, v: f64) {
        if v < self.upper(i, j) && v > self.lower(i, j) {
            self.set_upper(i, j, v);
        }
    }

    #[allow(dead_code)]
    pub(crate) fn set_lower_if_better(&mut self, i: usize, j: usize, v: f64) {
        if v > self.lower(i, j) && v < self.upper(i, j) {
            self.set_lower(i, j, v);
        }
    }

    /// Build from a full `n x n` row-major matrix in RDKit's layout.
    pub fn from_raw(n: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), n * n);
        BoundsMatrix { n, data }
    }

    pub fn raw(&self) -> &[f64] {
        &self.data
    }
}

/// `DistGeom::triangleSmoothBounds(boundsMat, tol)`.
pub fn triangle_smooth_bounds(m: &mut BoundsMatrix, tol: f64) -> bool {
    let npt = m.n;
    for k in 0..npt {
        for i in 0..npt.saturating_sub(1) {
            if i == k {
                continue;
            }
            let (ii, ik) = if i > k { (k, i) } else { (i, k) };
            let uik = m.val(ii, ik);
            let lik = m.val(ik, ii);
            for j in i + 1..npt {
                if j == k {
                    continue;
                }
                let (jj, jk) = if j > k { (k, j) } else { (j, k) };
                let ukj = m.val(jj, jk);
                let sum_uik_ukj = uik + ukj;
                if m.val(i, j) > sum_uik_ukj {
                    m.set_val(i, j, sum_uik_ukj);
                }
                let diff_lik_ujk = lik - ukj;
                let diff_ljk_uik = m.val(jk, jj) - uik;
                if m.val(j, i) < diff_lik_ujk {
                    m.set_val(j, i, diff_lik_ujk);
                } else if m.val(j, i) < diff_ljk_uik {
                    m.set_val(j, i, diff_ljk_uik);
                }
                let l_bound = m.val(j, i);
                let u_bound = m.val(i, j);
                if tol > 0.0
                    && (l_bound - u_bound) / l_bound > 0.0
                    && (l_bound - u_bound) / l_bound < tol
                {
                    m.set_val(i, j, l_bound);
                } else if l_bound - u_bound > 0.0 {
                    return false;
                }
            }
        }
    }
    true
}

/// Packed lower-triangle symmetric matrix (`RDNumeric::SymmMatrix`).
#[derive(Clone, Debug)]
pub(crate) struct SymmMatrix {
    pub n: usize,
    pub data: Vec<f64>,
}

impl SymmMatrix {
    pub(crate) fn new(n: usize) -> Self {
        SymmMatrix {
            n,
            data: vec![0.0; n * (n + 1) / 2],
        }
    }

    #[inline]
    pub(crate) fn get(&self, i: usize, j: usize) -> f64 {
        if i >= j {
            self.data[i * (i + 1) / 2 + j]
        } else {
            self.data[j * (j + 1) / 2 + i]
        }
    }

    #[inline]
    pub(crate) fn set(&mut self, i: usize, j: usize, v: f64) {
        if i >= j {
            self.data[i * (i + 1) / 2 + j] = v;
        } else {
            self.data[j * (j + 1) / 2 + i] = v;
        }
    }
}

/// `DistGeom::pickRandomDistMat`.
pub(crate) fn pick_random_dist_mat(
    m: &BoundsMatrix,
    dist: &mut SymmMatrix,
    rng: &mut MinstdRand,
) -> f64 {
    let npt = m.n;
    let mut largest = -1.0;
    for i in 1..npt {
        let id = i * (i + 1) / 2;
        for j in 0..i {
            let ub = m.upper(i, j);
            let lb = m.lower(i, j);
            let rval = rng.next_f64();
            let d = lb + rval * (ub - lb);
            dist.data[id + j] = d;
            if d > largest {
                largest = d;
            }
        }
    }
    largest
}

/// `RDNumeric::multiply(SymmMatrix, Vector, Vector)`.
fn symm_mul(a: &SymmMatrix, x: &[f64], y: &mut [f64]) {
    let n = a.n;
    for i in 0..n {
        y[i] = 0.0;
        let mut id = i * (i + 1) / 2;
        for j in 0..=i {
            y[i] += a.data[id] * x[j];
            id += 1;
        }
        id -= 1;
        for j in i + 1..n {
            id += j;
            y[i] += a.data[id] * x[j];
        }
    }
}

fn norm_l2(v: &[f64]) -> f64 {
    let mut r = 0.0;
    for x in v {
        r += x * x;
    }
    r.sqrt()
}

/// `Vector::normalize`; `false` for a (near) zero vector (RDKit throws).
fn normalize(v: &mut [f64]) -> bool {
    let l = norm_l2(v);
    if l < 1e-16 {
        return false;
    }
    for x in v.iter_mut() {
        *x /= l;
    }
    true
}

/// `RDNumeric::EigenSolvers::powerEigenSolver`. Writes eigenvalues and
/// eigenvectors (rows) for the components it converged on.
pub(crate) fn power_eigen_solver(
    num_eig: usize,
    mat: &mut SymmMatrix,
    eig_vals: &mut [f64],
    eig_vecs: &mut [Vec<f64>],
    seed: i32,
) -> Result<bool, ()> {
    const MAX_ITERATIONS: u32 = 1000;
    const TOLERANCE: f64 = 0.001;
    const HUGE_EIGVAL: f64 = 1.0e10;
    const TINY_EIGVAL: f64 = 1.0e-10;
    let n = mat.n;
    let mut v = vec![0.0; n];
    let mut z = vec![0.0; n];
    // seed <= 0 would use clock(): never the case for a seeded embedding of
    // a real molecule; treat as seed 1 deterministically.
    let mut seed: i32 = if seed <= 0 { 1 } else { seed };
    let mut converged = false;
    for ei in 0..num_eig {
        let mut eig_val = -HUGE_EIGVAL;
        seed = seed.wrapping_add(ei as i32);
        // Vector::setToRandom(seed): its own generator, seeded with
        // (unsigned) seed, or clock()+1 for seed 0.
        let mut vrng = MinstdRand::new(if seed as u32 > 0 { seed as u32 } else { 1 });
        for x in v.iter_mut() {
            *x = vrng.next_f64();
        }
        if !normalize(&mut v) {
            return Err(());
        }
        converged = false;
        for _iter in 0..MAX_ITERATIONS {
            symm_mul(mat, &v, &mut z);
            let prev_val = eig_val;
            // largestAbsValId
            let mut res = -1.0;
            let mut eval_id = n;
            for (i, x) in z.iter().enumerate() {
                if x.abs() > res {
                    res = x.abs();
                    eval_id = i;
                }
            }
            eig_val = z[eval_id];
            if eig_val.abs() < TINY_EIGVAL {
                break;
            }
            v.copy_from_slice(&z);
            for x in v.iter_mut() {
                *x /= eig_val;
            }
            if (eig_val - prev_val).abs() < TOLERANCE {
                converged = true;
                break;
            }
        }
        if !converged {
            break;
        }
        if !normalize(&mut v) {
            return Err(());
        }
        eig_vecs[ei].copy_from_slice(&v);
        eig_vals[ei] = eig_val;
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                mat.data[id + j] -= eig_val * v[i] * v[j];
            }
        }
    }
    Ok(converged)
}

const EIGVAL_TOL: f64 = 0.001;

/// `DistGeom::computeInitialCoords` into `pos` (`n * dim`, row-major).
pub(crate) fn compute_initial_coords(
    dist: &SymmMatrix,
    pos: &mut [f64],
    dim: usize,
    rng: &mut MinstdRand,
    rand_neg_eig: bool,
    num_zero_fail: u32,
) -> Result<bool, ()> {
    let n = dist.n;
    let mut sq = SymmMatrix::new(n);
    let mut t = SymmMatrix::new(n);
    let mut eig_vecs = vec![vec![0.0; n]; dim];
    let mut eig_vals = vec![0.0; dim];
    let mut sum_sq_d2 = 0.0;
    for i in 0..dist.data.len() {
        sq.data[i] = dist.data[i] * dist.data[i];
        sum_sq_d2 += sq.data[i];
    }
    sum_sq_d2 /= (n as u32).wrapping_mul(n as u32) as f64;
    let mut sq_d0i = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            sq_d0i[i] += sq.get(i, j);
        }
        sq_d0i[i] /= n as f64;
        sq_d0i[i] -= sum_sq_d2;
        if sq_d0i[i] < EIGVAL_TOL && n > 3 {
            return Ok(false);
        }
    }
    for i in 0..n {
        for j in 0..=i {
            let val = 0.5 * (sq_d0i[i] + sq_d0i[j] - sq.get(i, j));
            t.set(i, j, val);
        }
    }
    let n_eigs = dim.min(n);
    power_eigen_solver(
        n_eigs,
        &mut t,
        &mut eig_vals,
        &mut eig_vecs,
        (sum_sq_d2 * n as f64) as i32,
    )?;
    let mut found_neg = false;
    let mut zero_eigs = 0;
    for e in eig_vals.iter_mut() {
        if *e > EIGVAL_TOL {
            *e = e.sqrt();
        } else if e.abs() < EIGVAL_TOL {
            *e = 0.0;
            zero_eigs += 1;
        } else {
            found_neg = true;
        }
    }
    if found_neg && !rand_neg_eig {
        return Ok(false);
    }
    if zero_eigs >= num_zero_fail && n > 3 {
        return Ok(false);
    }
    for i in 0..n {
        for j in 0..dim {
            pos[i * dim + j] = if eig_vals[j] >= 0.0 {
                eig_vals[j] * eig_vecs[j][i]
            } else {
                1.0 - 2.0 * rng.next_f64()
            };
        }
    }
    Ok(true)
}

/// `DistGeom::ChiralSet`.
#[derive(Clone, Debug)]
pub struct ChiralSet {
    pub idx: [usize; 5],
    pub vol_lower: f64,
    pub vol_upper: f64,
    pub in_fused_small_rings: bool,
}

/// `DistGeom::calcChiralVolume` on flat positions of dimension `dim`.
pub(crate) fn chiral_volume(
    i1: usize,
    i2: usize,
    i3: usize,
    i4: usize,
    pos: &[f64],
    dim: usize,
) -> f64 {
    let p = |a: usize, k: usize| pos[a * dim + k];
    let v1 = [
        p(i1, 0) - p(i4, 0),
        p(i1, 1) - p(i4, 1),
        p(i1, 2) - p(i4, 2),
    ];
    let v2 = [
        p(i2, 0) - p(i4, 0),
        p(i2, 1) - p(i4, 1),
        p(i2, 2) - p(i4, 2),
    ];
    let v3 = [
        p(i3, 0) - p(i4, 0),
        p(i3, 1) - p(i4, 1),
        p(i3, 2) - p(i4, 2),
    ];
    let c = cross(v2, v3);
    dot(v1, c)
}

#[inline]
pub(crate) fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        -a[0] * b[2] + a[2] * b[0],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
pub(crate) fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

struct DistContrib {
    i: usize,
    j: usize,
    ub2: f64,
    lb2: f64,
    weight: f64,
}

struct ChiralContrib {
    idx: [usize; 4],
    vol_upper: f64,
    vol_lower: f64,
    weight: f64,
}

/// The force field of `DistGeom::constructForceField`.
pub(crate) struct DgField {
    dim: usize,
    dist: Vec<DistContrib>,
    chiral: Vec<ChiralContrib>,
    fourth: Vec<(usize, f64)>,
    /// Contrib order: dist (if any), chiral (if any), fourth (if any).
    has_dist: bool,
    has_chiral: bool,
    has_fourth: bool,
}

impl DgField {
    pub(crate) fn new(
        m: &BoundsMatrix,
        dim: usize,
        csets: &[ChiralSet],
        weight_chiral: f64,
        weight_fourth: f64,
        basin_size_tol: f64,
    ) -> Self {
        let n = m.n;
        let mut dist = Vec::new();
        for i in 1..n {
            for j in 0..i {
                let l = m.lower(i, j);
                let u = m.upper(i, j);
                if u - l <= basin_size_tol {
                    dist.push(DistContrib {
                        i,
                        j,
                        ub2: u * u,
                        lb2: l * l,
                        weight: 1.0,
                    });
                }
            }
        }
        let mut chiral = Vec::new();
        if weight_chiral > 1.0e-8 {
            for c in csets {
                chiral.push(ChiralContrib {
                    idx: [c.idx[1], c.idx[2], c.idx[3], c.idx[4]],
                    vol_upper: c.vol_upper,
                    vol_lower: c.vol_lower,
                    weight: weight_chiral,
                });
            }
        }
        let mut fourth = Vec::new();
        if dim == 4 && weight_fourth > 1.0e-8 {
            for i in 0..n {
                fourth.push((i, weight_fourth));
            }
        }
        DgField {
            dim,
            has_dist: !dist.is_empty(),
            has_chiral: !chiral.is_empty(),
            has_fourth: !fourth.is_empty(),
            dist,
            chiral,
            fourth,
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        !(self.has_dist || self.has_chiral || self.has_fourth)
    }

    fn dist_energy(&self, pos: &[f64]) -> f64 {
        let dim = self.dim;
        let mut accum = 0.0;
        for c in &self.dist {
            let d2 = dist2(c.i, c.j, pos, dim);
            let mut val = 0.0;
            if d2 > c.ub2 {
                val = (d2 / c.ub2) - 1.0;
            } else if d2 < c.lb2 {
                val = ((2.0 * c.lb2) / (c.lb2 + d2)) - 1.0;
            }
            if val > 0.0 {
                accum += c.weight * val * val;
            }
        }
        accum
    }

    fn chiral_energy(&self, pos: &[f64]) -> f64 {
        let mut res = 0.0;
        for c in &self.chiral {
            let vol = chiral_volume(c.idx[0], c.idx[1], c.idx[2], c.idx[3], pos, self.dim);
            if vol < c.vol_lower {
                res += c.weight * (vol - c.vol_lower) * (vol - c.vol_lower);
            } else if vol > c.vol_upper {
                res += c.weight * (vol - c.vol_upper) * (vol - c.vol_upper);
            }
        }
        res
    }

    fn fourth_energy(&self, pos: &[f64]) -> f64 {
        let mut res = 0.0;
        for &(i, w) in &self.fourth {
            let pid = i * 4 + 3;
            res += w * pos[pid] * pos[pid];
        }
        res
    }

    /// `ForceField::calcEnergy(pos)`.
    pub(crate) fn energy(&self, pos: &[f64]) -> f64 {
        let mut res = 0.0;
        if self.has_dist {
            res += self.dist_energy(pos);
        }
        if self.has_chiral {
            res += self.chiral_energy(pos);
        }
        if self.has_fourth {
            res += self.fourth_energy(pos);
        }
        res
    }

    /// `ForceField::calcGrad(pos, grad)` (accumulating; caller zeroes).
    pub(crate) fn grad(&self, pos: &[f64], grad: &mut [f64]) {
        let dim = self.dim;
        if self.has_dist {
            for c in &self.dist {
                let d2 = dist2(c.i, c.j, pos, dim);
                let d;
                let pre;
                if d2 > c.ub2 {
                    d = d2.sqrt();
                    pre = 4.0 * (((d * d) / c.ub2) - 1.0) * (d / c.ub2);
                } else if d2 < c.lb2 {
                    d = d2.sqrt();
                    let l2d2 = d2 + c.lb2;
                    pre = 8.0 * c.lb2 * d * (1.0 - 2.0 * c.lb2 / l2d2) / (l2d2 * l2d2);
                } else {
                    continue;
                }
                for k in 0..dim {
                    let p1 = dim * c.i + k;
                    let p2 = dim * c.j + k;
                    let g = if d > 0.0 {
                        c.weight * pre * (pos[p1] - pos[p2]) / d
                    } else {
                        c.weight * pre * (pos[p1] - pos[p2])
                    };
                    grad[p1] += g;
                    grad[p2] -= g;
                }
            }
        }
        if self.has_chiral {
            for c in &self.chiral {
                let [i1, i2, i3, i4] = c.idx;
                let p = |a: usize, k: usize| pos[a * dim + k];
                let v1 = [
                    p(i1, 0) - p(i4, 0),
                    p(i1, 1) - p(i4, 1),
                    p(i1, 2) - p(i4, 2),
                ];
                let v2 = [
                    p(i2, 0) - p(i4, 0),
                    p(i2, 1) - p(i4, 1),
                    p(i2, 2) - p(i4, 2),
                ];
                let v3 = [
                    p(i3, 0) - p(i4, 0),
                    p(i3, 1) - p(i4, 1),
                    p(i3, 2) - p(i4, 2),
                ];
                let vol = dot(v1, cross(v2, v3));
                let pre = if vol < c.vol_lower {
                    c.weight * (vol - c.vol_lower)
                } else if vol > c.vol_upper {
                    c.weight * (vol - c.vol_upper)
                } else {
                    continue;
                };
                let (x, y, zz) = (0, 1, 2);
                grad[dim * i1] += pre * (v2[y] * v3[zz] - v3[y] * v2[zz]);
                grad[dim * i1 + 1] += pre * (v3[x] * v2[zz] - v2[x] * v3[zz]);
                grad[dim * i1 + 2] += pre * (v2[x] * v3[y] - v3[x] * v2[y]);
                grad[dim * i2] += pre * (v3[y] * v1[zz] - v3[zz] * v1[y]);
                grad[dim * i2 + 1] += pre * (v3[zz] * v1[x] - v3[x] * v1[zz]);
                grad[dim * i2 + 2] += pre * (v3[x] * v1[y] - v3[y] * v1[x]);
                grad[dim * i3] += pre * (v2[zz] * v1[y] - v2[y] * v1[zz]);
                grad[dim * i3 + 1] += pre * (v2[x] * v1[zz] - v2[zz] * v1[x]);
                grad[dim * i3 + 2] += pre * (v2[y] * v1[x] - v2[x] * v1[y]);
                grad[dim * i4] += pre
                    * (p(i1, 2) * (p(i2, 1) - p(i3, 1))
                        + p(i2, 2) * (p(i3, 1) - p(i1, 1))
                        + p(i3, 2) * (p(i1, 1) - p(i2, 1)));
                grad[dim * i4 + 1] += pre
                    * (p(i1, 0) * (p(i2, 2) - p(i3, 2))
                        + p(i2, 0) * (p(i3, 2) - p(i1, 2))
                        + p(i3, 0) * (p(i1, 2) - p(i2, 2)));
                grad[dim * i4 + 2] += pre
                    * (p(i1, 1) * (p(i2, 0) - p(i3, 0))
                        + p(i2, 1) * (p(i3, 0) - p(i1, 0))
                        + p(i3, 1) * (p(i1, 0) - p(i2, 0)));
            }
        }
        if self.has_fourth {
            for &(i, w) in &self.fourth {
                let pid = i * 4 + 3;
                grad[pid] += w * pos[pid];
            }
        }
    }

    /// `ForceField::minimize(maxIts, forceTol)`.
    pub(crate) fn minimize(
        &self,
        pos: &mut [f64],
        max_its: u32,
        force_tol: f64,
    ) -> Result<i32, super::bfgs::BadDirection> {
        if self.is_empty() {
            return Ok(0);
        }
        super::bfgs::minimize(
            pos,
            force_tol,
            max_its,
            |p| self.energy(p),
            |p, g| {
                g.iter_mut().for_each(|x| *x = 0.0);
                self.grad(p, g);
                super::bfgs::scale_gradient(g)
            },
        )
    }
}

#[inline]
fn dist2(i: usize, j: usize, pos: &[f64], dim: usize) -> f64 {
    let mut d2 = 0.0;
    for k in 0..dim {
        let d = pos[dim * i + k] - pos[dim * j + k];
        d2 += d * d;
    }
    d2
}
