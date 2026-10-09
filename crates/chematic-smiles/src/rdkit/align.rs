//! RDKit's point and molecule alignment (RDKit 2026.03.1):
//! `Numerics/Alignment/AlignPoints.cpp` (quaternion fit with Jacobi
//! diagonalisation), `Geometry/Transform3D` and `GraphMol/MolAlign`
//! (`AlignMol`, `GetBestRMS`, `CalcRMS`), with RDKit's operation order so
//! results agree bit for bit.

use super::RdkitSmilesError;
use super::mol::{BondType, Mol};
use super::substruct::substruct_matches;

const TOLERANCE: f64 = 1.0e-6;

type Point = [f64; 3];

/// `RDGeom::Transform3D`: a row-major 4x4 matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Transform3D(pub [f64; 16]);

impl Transform3D {
    fn identity() -> Self {
        let mut d = [0.0; 16];
        for i in 0..4 {
            d[i * 5] = 1.0;
        }
        Self(d)
    }

    /// `TransformPoint`.
    pub fn transform_point(&self, p: Point) -> Point {
        let d = &self.0;
        [
            d[0] * p[0] + d[1] * p[1] + d[2] * p[2] + d[3],
            d[4] * p[0] + d[5] * p[1] + d[6] * p[2] + d[7],
            d[8] * p[0] + d[9] * p[1] + d[10] * p[2] + d[11],
        ]
    }

    fn set_translation(&mut self, m: Point) {
        self.0[3] = m[0];
        self.0[7] = m[1];
        self.0[11] = m[2];
        self.0[15] = 1.0;
    }

    fn set_rotation_from_quaternion(&mut self, q: [f64; 4]) {
        let q00 = q[0] * q[0];
        let q11 = q[1] * q[1];
        let q22 = q[2] * q[2];
        let q33 = q[3] * q[3];
        let sum_sq = q00 + q11 + q22 + q33;
        let q01 = 2.0 * q[0] * q[1];
        let q02 = 2.0 * q[0] * q[2];
        let q03 = 2.0 * q[0] * q[3];
        let q12 = 2.0 * q[1] * q[2];
        let q13 = 2.0 * q[1] * q[3];
        let q23 = 2.0 * q[2] * q[3];
        let d = &mut self.0;
        d[0] = (q00 + q11 - q22 - q33) / sum_sq;
        d[1] = (q12 + q03) / sum_sq;
        d[2] = (q13 - q02) / sum_sq;
        d[4] = (q12 - q03) / sum_sq;
        d[5] = (q00 - q11 + q22 - q33) / sum_sq;
        d[6] = (q23 + q01) / sum_sq;
        d[8] = (q13 + q02) / sum_sq;
        d[9] = (q23 - q01) / sum_sq;
        d[10] = (q00 - q11 - q22 + q33) / sum_sq;
    }

    fn reflect(&mut self) {
        for i in 0..3 {
            for j in 0..3 {
                self.0[i * 4 + j] *= -1.0;
            }
        }
    }

    /// The matrix as rows.
    pub fn rows(&self) -> [[f64; 4]; 4] {
        let d = &self.0;
        [
            [d[0], d[1], d[2], d[3]],
            [d[4], d[5], d[6], d[7]],
            [d[8], d[9], d[10], d[11]],
            [d[12], d[13], d[14], d[15]],
        ]
    }
}

fn length_sq(p: Point) -> f64 {
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
}

/// `_weightedSumOfPoints`.
fn weighted_sum_of_points(points: &[Point], weights: Option<&[f64]>) -> Point {
    let mut res = [0.0; 3];
    for (i, p) in points.iter().enumerate() {
        let mut t = *p;
        if let Some(w) = weights {
            for v in &mut t {
                *v *= w[i];
            }
        }
        for k in 0..3 {
            res[k] += t[k];
        }
    }
    res
}

/// `_weightedSumOfLenSq`.
fn weighted_sum_of_len_sq(points: &[Point], weights: Option<&[f64]>) -> f64 {
    let mut res = 0.0;
    for (i, p) in points.iter().enumerate() {
        let mut l = length_sq(*p);
        if let Some(w) = weights {
            l *= w[i];
        }
        res += l;
    }
    res
}

/// `jacobi`: eigenvalues (ascending) and eigenvectors (columns) of `quad`.
fn jacobi(quad: &mut [[f64; 4]; 4], max_iter: u32) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut vals = [0.0; 4];
    let mut vecs = [[0.0; 4]; 4];
    for j in 0..4 {
        for row in &mut vecs {
            row[j] = 0.0;
        }
        vecs[j][j] = 1.0;
        vals[j] = quad[j][j];
    }
    'iterate: for _ in 0..max_iter {
        let mut diag_norm = 0.0;
        let mut off_diag_norm = 0.0;
        for j in 0..4 {
            diag_norm += f64::abs(vals[j]);
            for i in 0..j {
                off_diag_norm += f64::abs(quad[i][j]);
            }
        }
        if diag_norm.abs() > 1.0e-16 && (off_diag_norm / diag_norm) <= TOLERANCE {
            break 'iterate;
        }
        for j in 1..4 {
            for i in 0..j {
                let b = quad[i][j];
                if b.abs() > 0.0 {
                    let dma = vals[j] - vals[i];
                    let t = if (dma.abs() + b.abs()) <= dma.abs() {
                        b / dma
                    } else {
                        let q = 0.5 * dma / b;
                        let t = 1.0 / (q.abs() + (1.0 + q * q).sqrt());
                        if q < 0.0 { -t } else { t }
                    };
                    let c = 1.0 / (t * t + 1.0).sqrt();
                    let s = t * c;
                    quad[i][j] = 0.0;
                    for k in 0..i {
                        let atemp = c * quad[k][i] - s * quad[k][j];
                        quad[k][j] = s * quad[k][i] + c * quad[k][j];
                        quad[k][i] = atemp;
                    }
                    for k in i + 1..j {
                        let atemp = c * quad[i][k] - s * quad[k][j];
                        quad[k][j] = s * quad[i][k] + c * quad[k][j];
                        quad[i][k] = atemp;
                    }
                    for k in j + 1..4 {
                        let atemp = c * quad[i][k] - s * quad[j][k];
                        quad[j][k] = s * quad[i][k] + c * quad[j][k];
                        quad[i][k] = atemp;
                    }
                    for row in vecs.iter_mut() {
                        let vtemp = c * row[i] - s * row[j];
                        row[j] = s * row[i] + c * row[j];
                        row[i] = vtemp;
                    }
                    let dtemp = c * c * vals[i] + s * s * vals[j] - 2.0 * c * s * b;
                    vals[j] = s * s * vals[i] + c * c * vals[j] + 2.0 * c * s * b;
                    vals[i] = dtemp;
                }
            }
        }
    }
    for j in 0..3 {
        let mut k = j;
        let mut dtemp = vals[k];
        for (i, &v) in vals.iter().enumerate().skip(j + 1) {
            if v < dtemp {
                k = i;
                dtemp = v;
            }
        }
        if k > j {
            vals[k] = vals[j];
            vals[j] = dtemp;
            for row in vecs.iter_mut() {
                row.swap(j, k);
            }
        }
    }
    (vals, vecs)
}

/// `RDNumeric::Alignments::AlignPoints`: the sum of squared residuals of
/// the best fit of `probe` onto `reference`, and the transform achieving it.
pub fn align_points(
    reference: &[Point],
    probe: &[Point],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<(f64, Transform3D), RdkitSmilesError> {
    let npt = reference.len();
    if npt != probe.len() || weights.is_some_and(|w| w.len() != npt) {
        return Err(RdkitSmilesError::Unsupported(
            "Mismatch in number of points".into(),
        ));
    }
    let mut trans = Transform3D::identity();
    let wts_sum = match weights {
        Some(w) => {
            let mut res = 0.0;
            for &x in w {
                if x <= 0.0 {
                    return Err(RdkitSmilesError::Unsupported(
                        "Negative weight specified for a point".into(),
                    ));
                }
                res += x;
            }
            res
        }
        None => npt as f64,
    };
    let mut rpt_sum = weighted_sum_of_points(reference, weights);
    let ppt_sum = weighted_sum_of_points(probe, weights);
    let rpt_sum_len_sq = weighted_sum_of_len_sq(reference, weights);
    let ppt_sum_len_sq = weighted_sum_of_len_sq(probe, weights);

    // _computeCovarianceMat
    let mut cov = [[0.0f64; 3]; 3];
    for i in 0..npt {
        let r = reference[i];
        let p = probe[i];
        let w = weights.map_or(1.0, |w| w[i]);
        for a in 0..3 {
            for b in 0..3 {
                cov[a][b] += w * p[a] * r[b];
            }
        }
    }
    if reflect {
        for v in &mut rpt_sum {
            *v *= -1.0;
        }
        for row in &mut cov {
            for v in row {
                *v = -*v;
            }
        }
    }

    // _covertCovMatToQuad
    let temp = ppt_sum[0] / wts_sum;
    let pxrx = cov[0][0] - temp * rpt_sum[0];
    let pxry = cov[0][1] - temp * rpt_sum[1];
    let pxrz = cov[0][2] - temp * rpt_sum[2];
    let temp = ppt_sum[1] / wts_sum;
    let pyrx = cov[1][0] - temp * rpt_sum[0];
    let pyry = cov[1][1] - temp * rpt_sum[1];
    let pyrz = cov[1][2] - temp * rpt_sum[2];
    let temp = ppt_sum[2] / wts_sum;
    let pzrx = cov[2][0] - temp * rpt_sum[0];
    let pzry = cov[2][1] - temp * rpt_sum[1];
    let pzrz = cov[2][2] - temp * rpt_sum[2];
    let mut quad = [[0.0f64; 4]; 4];
    quad[0][0] = -2.0 * (pxrx + pyry + pzrz);
    quad[1][1] = -2.0 * (pxrx - pyry - pzrz);
    quad[2][2] = -2.0 * (pyry - pzrz - pxrx);
    quad[3][3] = -2.0 * (pzrz - pxrx - pyry);
    quad[0][1] = 2.0 * (pyrz - pzry);
    quad[1][0] = quad[0][1];
    quad[0][2] = 2.0 * (pzrx - pxrz);
    quad[2][0] = quad[0][2];
    quad[0][3] = 2.0 * (pxry - pyrx);
    quad[3][0] = quad[0][3];
    quad[1][2] = -2.0 * (pxry + pyrx);
    quad[2][1] = quad[1][2];
    quad[1][3] = -2.0 * (pzrx + pxrz);
    quad[3][1] = quad[1][3];
    quad[2][3] = -2.0 * (pyrz + pzry);
    quad[3][2] = quad[2][3];

    let (vals, vecs) = jacobi(&mut quad, max_iterations);
    let quater = [vecs[0][0], vecs[1][0], vecs[2][0], vecs[3][0]];
    trans.set_rotation_from_quaternion(quater);
    if reflect {
        trans.reflect();
    }
    let mut ssr = vals[0] - (length_sq(ppt_sum) + length_sq(rpt_sum)) / wts_sum
        + rpt_sum_len_sq
        + ppt_sum_len_sq;
    if ssr < 0.0 && ssr.abs() < TOLERANCE {
        ssr = 0.0;
    }
    if reflect {
        for v in &mut rpt_sum {
            *v *= -1.0;
        }
    }
    let moved = trans.transform_point(ppt_sum);
    let mut mv = rpt_sum;
    for k in 0..3 {
        mv[k] -= moved[k];
    }
    for v in &mut mv {
        *v /= wts_sum;
    }
    trans.set_translation(mv);
    Ok((ssr, trans))
}

/// An atom map: `(probe atom, reference atom)` pairs.
pub type AtomMap = Vec<(usize, usize)>;

fn check_coords(mol: &Mol, coords: &[Point], what: &str) -> Result<(), RdkitSmilesError> {
    if coords.len() != mol.atoms.len() {
        return Err(RdkitSmilesError::Unsupported(format!(
            "{what} has {} positions for {} atoms",
            coords.len(),
            mol.atoms.len()
        )));
    }
    Ok(())
}

/// `alignConfsOnAtomMap`: mean squared deviation after the best fit.
fn align_confs_on_atom_map(
    probe: &[Point],
    reference: &[Point],
    map: &[(usize, usize)],
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<(f64, Transform3D), RdkitSmilesError> {
    let mut prb = Vec::with_capacity(map.len());
    let mut refp = Vec::with_capacity(map.len());
    for &(p, r) in map {
        let (Some(&pp), Some(&rr)) = (probe.get(p), reference.get(r)) else {
            return Err(RdkitSmilesError::Unsupported(
                "atom map index out of range".into(),
            ));
        };
        prb.push(pp);
        refp.push(rr);
    }
    let (ssr, trans) = align_points(&refp, &prb, weights, reflect, max_iterations)?;
    Ok((ssr / prb.len() as f64, trans))
}

/// `calcMSDInternal`: mean squared deviation without moving the probe.
fn calc_msd(
    probe: &[Point],
    reference: &[Point],
    map: &[(usize, usize)],
    weights: Option<&[f64]>,
) -> f64 {
    let mut ssr = 0.0;
    for (i, &(p, r)) in map.iter().enumerate() {
        let (pp, rr) = (probe[p], reference[r]);
        let d = [pp[0] - rr[0], pp[1] - rr[1], pp[2] - rr[2]];
        ssr += weights.map_or(1.0, |w| w[i]) * length_sq(d);
    }
    ssr / map.len() as f64
}

/// The probe molecule as `getAllMatchesPrbRef` matches it: per atom
/// whether its formal charge is ignored, per bond whether it matches any
/// single or double bond (`details::symmetrizeTerminalAtoms`).
fn symmetrized_terminal_groups(mol: &Mol) -> (Vec<bool>, Vec<bool>) {
    // [O,N;D1;$([O,N;D1]-[*]=[O,N;D1]),$([O,N;D1]=[*]-[O,N;D1])]~[*]
    let terminal = |a: usize| matches!(mol.atoms[a].anum, 7 | 8) && mol.degree(a) == 1;
    let mut neutral = vec![false; mol.atoms.len()];
    let mut either = vec![false; mol.bonds.len()];
    for x in 0..mol.atoms.len() {
        if !terminal(x) {
            continue;
        }
        let b = mol.atom_bonds[x][0];
        let y = mol.bonds[b].other(x);
        let other_end = |want_first: BondType, want_second: BondType| {
            mol.bonds[b].bt == want_first
                && mol.atom_bonds[y].iter().any(|&b2| {
                    let z = mol.bonds[b2].other(y);
                    b2 != b && z != x && mol.bonds[b2].bt == want_second && terminal(z)
                })
        };
        if other_end(BondType::Single, BondType::Double)
            || other_end(BondType::Double, BondType::Single)
        {
            neutral[x] = true;
            either[b] = true;
        }
    }
    (neutral, either)
}

/// `Atom::Match(what)` of a plain query atom (`ignore_charge`: its formal
/// charge was set to 0).
fn atom_match(q: &Mol, qa: usize, m: &Mol, ma: usize, ignore_charge: bool) -> bool {
    let (a, b) = (&q.atoms[qa], &m.atoms[ma]);
    if a.anum != b.anum {
        return false;
    }
    if a.anum == 0 {
        return !(a.isotope != 0 && b.isotope != 0 && a.isotope != b.isotope);
    }
    let charge = if ignore_charge { 0 } else { a.charge };
    !((charge != 0 && charge != b.charge)
        || (a.isotope != 0 && a.isotope != b.isotope)
        || (a.radicals != 0 && a.radicals != b.radicals))
}

/// All matches of the probe (as the query) in the reference, as
/// `getAllMatchesPrbRef` (`uniquify=false`) enumerates them; here both are
/// the same molecule.
fn all_matches(
    mol: &Mol,
    max_matches: usize,
    symmetrize: bool,
) -> Result<Vec<AtomMap>, RdkitSmilesError> {
    let (neutral, either) = if symmetrize {
        symmetrized_terminal_groups(mol)
    } else {
        (vec![false; mol.atoms.len()], vec![false; mol.bonds.len()])
    };
    let vc = |q: usize, t: usize| atom_match(mol, q, mol, t, neutral[q]);
    let ec = |qb: usize, tb: usize| {
        let (q, t) = (&mol.bonds[qb], &mol.bonds[tb]);
        if either[qb] {
            return matches!(t.bt, BondType::Single | BondType::Double);
        }
        if q.bt != t.bt {
            return false;
        }
        if q.bt == BondType::Dative {
            // The direction must match too.
            return atom_match(mol, q.begin, mol, t.begin, neutral[q.begin])
                && atom_match(mol, q.end, mol, t.end, neutral[q.end]);
        }
        true
    };
    let matches = substruct_matches(mol, mol, &vc, &ec, max_matches);
    if matches.is_empty() {
        return Err(RdkitSmilesError::Unsupported(
            "No sub-structure match found between the reference and probe mol".into(),
        ));
    }
    Ok(matches)
}

/// Result of an RDKit-style alignment.
#[derive(Clone, Debug, PartialEq)]
pub struct RdkitAlignment {
    /// The RMSD.
    pub rmsd: f64,
    /// The transform taking the probe onto the reference.
    pub transform: Transform3D,
    /// The atom map used, `(probe atom, reference atom)`.
    pub atom_map: AtomMap,
}

/// `getAlignmentTransform` / `AlignMol`.
pub(crate) fn align_mol(
    mol: &Mol,
    probe: &[Point],
    reference: &[Point],
    atom_map: Option<&[(usize, usize)]>,
    weights: Option<&[f64]>,
    reflect: bool,
    max_iterations: u32,
) -> Result<RdkitAlignment, RdkitSmilesError> {
    check_coords(mol, probe, "probe")?;
    check_coords(mol, reference, "reference")?;
    let map: AtomMap = match atom_map {
        Some(m) => m.to_vec(),
        None => {
            // `SubstructMatch(refMol, prbMol, match)`: the first match.
            let mut first = all_matches_plain(mol, 1)?;
            first.swap_remove(0)
        }
    };
    let (msd, transform) =
        align_confs_on_atom_map(probe, reference, &map, weights, reflect, max_iterations)?;
    Ok(RdkitAlignment {
        rmsd: msd.sqrt(),
        transform,
        atom_map: map,
    })
}

fn all_matches_plain(mol: &Mol, max_matches: usize) -> Result<Vec<AtomMap>, RdkitSmilesError> {
    all_matches(mol, max_matches, false).map_err(|_| {
        RdkitSmilesError::Unsupported(
            "No sub-structure match found between the probe and query mol".into(),
        )
    })
}

/// `getBestRMSInternal` over `matches`: with `align`, the best fit per
/// match (`GetBestRMS`), otherwise the deviation in place (`CalcRMS`).
fn best_rms_internal(
    probe: &[Point],
    reference: &[Point],
    matches: &[AtomMap],
    weights: Option<&[f64]>,
    align: bool,
) -> Result<RdkitAlignment, RdkitSmilesError> {
    let mut msd_best = f64::MAX;
    let mut best = 0;
    let mut best_trans = Transform3D::identity();
    for (k, m) in matches.iter().enumerate() {
        if weights.is_some_and(|w| w.len() != m.len()) {
            return Err(RdkitSmilesError::Unsupported(
                "Mismatch in number of weights".into(),
            ));
        }
        let (msd, trans) = if align {
            align_confs_on_atom_map(probe, reference, m, weights, false, 50)?
        } else {
            (
                calc_msd(probe, reference, m, weights),
                Transform3D::identity(),
            )
        };
        if msd < msd_best {
            msd_best = msd;
            best = k;
            best_trans = trans;
        }
    }
    Ok(RdkitAlignment {
        rmsd: msd_best.sqrt(),
        transform: best_trans,
        atom_map: matches[best].clone(),
    })
}

/// `GetBestRMS(prbMol, refMol, maxMatches, symmetrizeConjugatedTerminalGroups,
/// weights)` (one thread) for two conformers of `mol`; also returns the
/// best transform and match (`GetBestAlignmentTransform`).
pub(crate) fn best_rms(
    mol: &Mol,
    probe: &[Point],
    reference: &[Point],
    max_matches: usize,
    symmetrize: bool,
    weights: Option<&[f64]>,
) -> Result<RdkitAlignment, RdkitSmilesError> {
    check_coords(mol, probe, "probe")?;
    check_coords(mol, reference, "reference")?;
    let matches = all_matches(mol, max_matches, symmetrize)?;
    best_rms_internal(probe, reference, &matches, weights, true)
}

/// `CalcRMS(prbMol, refMol, maxMatches, symmetrizeConjugatedTerminalGroups,
/// weights)`: the smallest in-place RMSD over all matches.
pub(crate) fn calc_rms(
    mol: &Mol,
    probe: &[Point],
    reference: &[Point],
    max_matches: usize,
    symmetrize: bool,
    weights: Option<&[f64]>,
) -> Result<f64, RdkitSmilesError> {
    check_coords(mol, probe, "probe")?;
    check_coords(mol, reference, "reference")?;
    let matches = all_matches(mol, max_matches, symmetrize)?;
    Ok(best_rms_internal(probe, reference, &matches, weights, false)?.rmsd)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aligning_a_rotated_copy_recovers_it() {
        let reference = [
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [1.5, 1.2, 0.3],
            [-0.4, 0.9, -1.1],
        ];
        // Rotate 90 degrees about z and translate.
        let probe: Vec<Point> = reference
            .iter()
            .map(|p| [-p[1] + 2.0, p[0] - 1.0, p[2] + 0.5])
            .collect();
        let (ssr, t) = align_points(&reference, &probe, None, false, 50).unwrap();
        assert!(ssr.abs() < 1e-12);
        for (p, r) in probe.iter().zip(&reference) {
            let q = t.transform_point(*p);
            for k in 0..3 {
                assert!((q[k] - r[k]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn best_rms_tries_symmetric_matches() {
        // Propane: swapping the terminal carbons is a second match.
        let mol = crate::parse("CCC").unwrap();
        let reference = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 1.4, 0.0]];
        let probe = [[2.0, 1.4, 0.0], [1.5, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let calc = crate::rdkit_calc_rms(&mol, &probe, &reference, 1_000_000, true, None).unwrap();
        assert_eq!(calc, 0.0);
        let best = crate::rdkit_best_rms(&mol, &probe, &reference, 1_000_000, true, None).unwrap();
        assert!(best.rmsd < 1e-6);
        assert_eq!(best.atom_map.len(), 3);
    }
}
