//! `DGeomHelpers::EmbeddingOps::embedPoints` and its checks.

use super::dg::{
    BoundsMatrix, ChiralSet, DgField, SymmMatrix, chiral_volume, compute_initial_coords, cross,
    dot, pick_random_dist_mat,
};
use super::rng::MinstdRand;

const ERROR_TOL: f64 = 0.00001;
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;
const MIN_TETRAHEDRAL_CHIRAL_VOL: f64 = 0.50;
const TETRAHEDRAL_CENTERINVOLUME_TOL: f64 = 0.30;

/// The subset of `EmbedParameters` the numeric pipeline reads.
#[derive(Clone, Debug)]
pub struct EmbedParams {
    pub max_iterations: u32,
    pub random_seed: i32,
    pub rand_neg_eig: bool,
    pub num_zero_fail: u32,
    pub optimizer_force_tol: f64,
    pub enforce_chirality: bool,
    pub basin_thresh: f64,
}

impl Default for EmbedParams {
    fn default() -> Self {
        EmbedParams {
            max_iterations: 0,
            random_seed: 42,
            rand_neg_eig: true,
            num_zero_fail: 1,
            optimizer_force_tol: 1e-3,
            enforce_chirality: true,
            basin_thresh: 5.0,
        }
    }
}

/// Everything `embedPoints` needs about one fragment.
pub struct EmbedArgs<'a> {
    pub mmat: &'a BoundsMatrix,
    pub chiral_centers: &'a [ChiralSet],
    pub tetrahedral_centers: &'a [ChiralSet],
    /// (neighbor, atom, double-bond partner).
    pub double_bond_ends: &'a [(usize, usize, usize)],
    /// ([controlling atoms], +1 trans / -1 cis).
    pub stereo_double_bonds: &'a [([usize; 4], i32)],
    /// The (ET)(K)DG stage: minimizes the 3D positions in place and
    /// returns the planarity verdict; `None` for plain DG.
    #[allow(clippy::type_complexity)]
    pub exp_torsions: Option<&'a dyn Fn(&mut [[f64; 3]]) -> Result<bool, EmbedError>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EmbedError {
    /// RDKit throws an invariant violation ("bad direction in linearSearch").
    BadDirection,
    /// RDKit throws ("Cannot normalize a zero length vector").
    ZeroLengthVector,
}

/// Map RDKit's minimizer status, turning its thrown invariant into an error.
pub(crate) fn check_status(status: i32) -> Result<i32, EmbedError> {
    if status < 0 {
        Err(EmbedError::BadDirection)
    } else {
        Ok(status)
    }
}

fn p3(pos: &[f64], dim: usize, i: usize) -> [f64; 3] {
    [pos[i * dim], pos[i * dim + 1], pos[i * dim + 2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn length(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn normalized(a: [f64; 3]) -> Result<[f64; 3], EmbedError> {
    let l = length(a);
    if l < 1e-16 {
        return Err(EmbedError::ZeroLengthVector);
    }
    Ok([a[0] / l, a[1] / l, a[2] / l])
}

fn volume_test(c: &ChiralSet, pos: &[f64], dim: usize) -> Result<bool, EmbedError> {
    let p0 = p3(pos, dim, c.idx[0]);
    let v1 = normalized(sub(p0, p3(pos, dim, c.idx[1])))?;
    let v2 = normalized(sub(p0, p3(pos, dim, c.idx[2])))?;
    let v3 = normalized(sub(p0, p3(pos, dim, c.idx[3])))?;
    let v4 = normalized(sub(p0, p3(pos, dim, c.idx[4])))?;
    let vol_scale = if c.in_fused_small_rings { 0.25 } else { 1.0 };
    let lim = vol_scale * MIN_TETRAHEDRAL_CHIRAL_VOL;
    if dot(cross(v1, v2), v3).abs() < lim {
        return Ok(false);
    }
    if dot(cross(v1, v2), v4).abs() < lim {
        return Ok(false);
    }
    if dot(cross(v1, v3), v4).abs() < lim {
        return Ok(false);
    }
    Ok(dot(cross(v2, v3), v4).abs() >= lim)
}

fn same_side(
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    p0: [f64; 3],
    tol: f64,
) -> bool {
    let normal = cross(sub(v2, v1), sub(v3, v1));
    let d1 = dot(normal, sub(v4, v1));
    let d2 = dot(normal, sub(p0, v1));
    if d1.abs() < tol || d2.abs() < tol {
        return false;
    }
    !((d1 < 0.0) ^ (d2 < 0.0))
}

fn center_in_volume(c: &ChiralSet, pos: &[f64], dim: usize, tol: f64) -> bool {
    if c.idx[0] == c.idx[4] {
        return true;
    }
    let [p0, p1, p2, p3_, p4] = c.idx.map(|i| p3(pos, dim, i));
    same_side(p1, p2, p3_, p4, p0, tol)
        && same_side(p2, p3_, p4, p1, p0, tol)
        && same_side(p3_, p4, p1, p2, p0, tol)
        && same_side(p4, p1, p2, p3_, p0, tol)
}

fn bounds_fulfilled(atoms: &[usize], mmat: &BoundsMatrix, pos: &[f64], dim: usize) -> bool {
    for i in 0..atoms.len() - 1 {
        for j in i + 1..atoms.len() {
            let (a1, a2) = (atoms[i], atoms[j]);
            let d2 = length(sub(p3(pos, dim, a1), p3(pos, dim, a2)));
            let lb = mmat.lower(a1, a2);
            let ub = mmat.upper(a1, a2);
            if (d2 < lb && (d2 - lb).abs() > 0.1 * ub) || (d2 > ub && (d2 - ub).abs() > 0.1 * ub) {
                return false;
            }
        }
    }
    true
}

fn have_opposite_sign(a: f64, b: f64) -> bool {
    a.is_sign_negative() ^ b.is_sign_negative()
}

fn check_chiral_centers(pos: &[f64], dim: usize, args: &EmbedArgs) -> bool {
    for c in args.chiral_centers {
        let vol = chiral_volume(c.idx[1], c.idx[2], c.idx[3], c.idx[4], pos, dim);
        let lb = c.vol_lower;
        let ub = c.vol_upper;
        if (lb > 0.0 && vol < lb && (vol / lb < 0.8 || have_opposite_sign(vol, lb)))
            || (ub < 0.0 && vol > ub && (vol / ub < 0.8 || have_opposite_sign(vol, ub)))
        {
            return false;
        }
    }
    true
}

fn first_minimization(
    pos: &mut [f64],
    dim: usize,
    args: &EmbedArgs,
    p: &EmbedParams,
) -> Result<bool, EmbedError> {
    let field = DgField::new(
        args.mmat,
        dim,
        args.chiral_centers,
        1.0,
        0.1,
        p.basin_thresh,
    );
    if field.energy(pos) > ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            need_more = check_status(field.minimize(pos, 400, p.optimizer_force_tol))?;
        }
    }
    let local_e = field.energy(pos);
    let n = pos.len() / dim;
    Ok(local_e / (n as f64) < MAX_MINIMIZED_E_PER_ATOM)
}

fn minimize_fourth_dimension(
    pos: &mut [f64],
    dim: usize,
    args: &EmbedArgs,
    p: &EmbedParams,
) -> Result<bool, EmbedError> {
    let field = DgField::new(
        args.mmat,
        dim,
        args.chiral_centers,
        0.2,
        1.0,
        p.basin_thresh,
    );
    if field.energy(pos) > ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            need_more = check_status(field.minimize(pos, 200, p.optimizer_force_tol))?;
        }
    }
    Ok(true)
}

fn double_bond_geometry_checks(
    pos: &[f64],
    dim: usize,
    args: &EmbedArgs,
) -> Result<bool, EmbedError> {
    const LINEAR_TOL: f64 = 1e-3;
    for &(a0, a1, a2) in args.double_bond_ends {
        let p0 = p3(pos, dim, a0);
        let p1 = p3(pos, dim, a1);
        let p2 = p3(pos, dim, a2);
        let v1 = normalized(sub(p1, p0))?;
        let v2 = normalized(sub(p1, p2))?;
        if dot(v1, v2) + 1.0 < LINEAR_TOL {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `RDGeom::computeDihedralAngle` (`Point3D::angleTo` of the normals).
fn dihedral_angle(p1: [f64; 3], p2: [f64; 3], p3_: [f64; 3], p4: [f64; 3]) -> f64 {
    let beg_end = sub(p3_, p2);
    let beg_nbr = sub(p1, p2);
    let crs1 = cross(beg_nbr, beg_end);
    let end_nbr = sub(p4, p3_);
    let crs2 = cross(end_nbr, beg_end);
    let lsq = dot(crs1, crs1) * dot(crs2, crs2);
    let mut d = dot(crs1, crs2);
    d /= lsq.sqrt();
    if d <= -1.0 {
        return std::f64::consts::PI;
    }
    if d >= 1.0 {
        return 0.0;
    }
    d.acos()
}

fn double_bond_stereo_checks(pos: &[f64], dim: usize, args: &EmbedArgs) -> bool {
    for (atoms, sign) in args.stereo_double_bonds {
        let [a0, a1, a2, a3] = atoms.map(|i| p3(pos, dim, i));
        let d = dihedral_angle(a0, a1, a2, a3);
        if (d - std::f64::consts::FRAC_PI_2) * (*sign as f64) < 0.0 {
            return false;
        }
    }
    true
}

fn final_chiral_checks(pos: &[f64], dim: usize, args: &EmbedArgs) -> bool {
    if !check_chiral_centers(pos, dim, args) {
        return false;
    }
    let mut atoms = std::collections::BTreeSet::new();
    for c in args.chiral_centers {
        if c.idx[0] != c.idx[4] {
            atoms.extend(c.idx);
        }
    }
    let atoms: Vec<usize> = atoms.into_iter().collect();
    if !atoms.is_empty() && !bounds_fulfilled(&atoms, args.mmat, pos, dim) {
        return false;
    }
    args.chiral_centers
        .iter()
        .all(|c| center_in_volume(c, pos, dim, 0.1))
}

/// `EmbeddingOps::embedPoints`: positions (`n * dim`) on success.
pub fn embed_points(
    n: usize,
    args: &EmbedArgs,
    params: &EmbedParams,
    seed: i32,
) -> Result<Option<Vec<f64>>, EmbedError> {
    let dim = if args.chiral_centers.is_empty() { 3 } else { 4 };
    let max_iterations = if params.max_iterations == 0 {
        10 * n as u32
    } else {
        params.max_iterations
    };
    let mut dist = SymmMatrix::new(n);
    let mut rng = MinstdRand::new(seed as u32);
    let mut pos = vec![0.0; n * dim];
    let mut got = false;
    let mut iter = 0;
    while !got && iter < max_iterations {
        iter += 1;
        pick_random_dist_mat(args.mmat, &mut dist, &mut rng);
        got = compute_initial_coords(
            &dist,
            &mut pos,
            dim,
            &mut rng,
            params.rand_neg_eig,
            params.num_zero_fail,
        )
        .map_err(|_| EmbedError::ZeroLengthVector)?;
        if !got {
            continue;
        }
        got = first_minimization(&mut pos, dim, args, params)?;
        if got {
            for t in args.tetrahedral_centers {
                if !volume_test(t, &pos, dim)?
                    || !center_in_volume(t, &pos, dim, TETRAHEDRAL_CENTERINVOLUME_TOL)
                {
                    got = false;
                    break;
                }
            }
        }
        if got && params.enforce_chirality && !args.chiral_centers.is_empty() {
            got = check_chiral_centers(&pos, dim, args);
        }
        if got && !args.chiral_centers.is_empty() {
            got = minimize_fourth_dimension(&mut pos, dim, args, params)?;
        }
        if got && let Some(etk) = args.exp_torsions {
            let mut p3s: Vec<[f64; 3]> = (0..n).map(|i| p3(&pos, dim, i)).collect();
            got = etk(&mut p3s)?;
            for (i, p) in p3s.iter().enumerate() {
                pos[i * dim..i * dim + 3].copy_from_slice(p);
            }
        }
        if got {
            got = double_bond_geometry_checks(&pos, dim, args)?;
        }
        if params.enforce_chirality && got {
            if !args.chiral_centers.is_empty() {
                got = final_chiral_checks(&pos, dim, args);
            }
            if got && !args.stereo_double_bonds.is_empty() {
                got = double_bond_stereo_checks(&pos, dim, args);
            }
        }
    }
    if !got {
        return Ok(None);
    }
    Ok(Some((0..n).flat_map(|i| p3(&pos, dim, i)).collect()))
}
