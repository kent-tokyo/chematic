//! ETKDG's experimental-torsion stage: `CrystalFF::getExperimentalTorsions`,
//! `DistGeom::construct3DForceField` (torsion M6, UFF inversion, distance
//! and angle constraints), `EmbeddingOps::minimizeWithExpTorsions` and the
//! planarity check (`construct3DImproperForceField`).

use std::sync::OnceLock;

use chematic_smarts::QueryMolecule;
use chematic_smiles::RdkitMolView;

use super::bounds::BondsAndAngles;
use super::dg::BoundsMatrix;
use super::substruct::{Target, substruct_match};
use super::torsion_prefs::{
    TORSION_PREFERENCES_MACROCYCLES, TORSION_PREFERENCES_SMALL_RINGS, TORSION_PREFERENCES_V1,
    TORSION_PREFERENCES_V2,
};

const MIN_MACROCYCLE_SIZE: usize = 9;
const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_FORCE_CONSTANT: f64 = 100.0;
const ERROR_TOL: f64 = 0.00001;

/// One parsed line of a torsion table (`ExpTorsionAngle`).
pub(crate) struct ExpTorsionAngle {
    pub smarts: String,
    pub pattern: QueryMolecule,
    /// query atom indices of map numbers 1..4
    pub idx: [usize; 4],
    pub signs: [i32; 6],
    pub v: [f64; 6],
}

fn parse_table(data: &str) -> Vec<ExpTorsionAngle> {
    let mut out = Vec::new();
    for line in data.split('\n') {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut tok = line.split(' ').filter(|t| !t.is_empty());
        let smarts = tok.next().expect("smarts").to_string();
        let mut signs = [0i32; 6];
        let mut v = [0f64; 6];
        for i in 0..6 {
            signs[i] = tok.next().expect("sign").parse().expect("int");
            v[i] = tok.next().expect("V").parse().expect("double");
        }
        let pattern = chematic_smarts::parse_smarts(&smarts)
            .unwrap_or_else(|e| panic!("torsion SMARTS {smarts}: {e:?}"));
        let mut idx = [0usize; 4];
        for (i, a) in pattern.atoms.iter().enumerate() {
            if let Some(n) = a.atom_map
                && (1..5).contains(&n)
            {
                idx[usize::from(n) - 1] = i;
            }
        }
        out.push(ExpTorsionAngle {
            smarts,
            pattern,
            idx,
            signs,
            v,
        });
    }
    out
}

/// `ExpTorsionAngleCollection::getParams(version, smallRings, macrocycles)`.
fn params(version: u32, small_rings: bool, macrocycles: bool) -> &'static [ExpTorsionAngle] {
    static CACHE: OnceLock<[OnceLock<Vec<ExpTorsionAngle>>; 8]> = OnceLock::new();
    let cache = CACHE.get_or_init(Default::default);
    let key = (if version == 1 { 0 } else { 4 })
        + usize::from(small_rings) * 2
        + usize::from(macrocycles);
    cache[key].get_or_init(|| {
        let mut s = String::from(if version == 1 {
            TORSION_PREFERENCES_V1
        } else {
            TORSION_PREFERENCES_V2
        });
        if small_rings {
            s.push_str(TORSION_PREFERENCES_SMALL_RINGS);
        }
        if macrocycles {
            s.push_str(TORSION_PREFERENCES_MACROCYCLES);
        }
        parse_table(&s)
    })
}

/// `CrystalFFDetails` (the parts used without coordMap / CPCI).
#[derive(Clone, Debug, Default)]
pub struct CrystalFfDetails {
    pub exp_torsion_atoms: Vec<[usize; 4]>,
    /// (signs, force constants)
    pub exp_torsion_angles: Vec<([i32; 6], [f64; 6])>,
    /// [nbr0, centre, nbr1, nbr2, atomic number of centre, bound to sp2 O]
    pub improper_atoms: Vec<[usize; 6]>,
    pub bonds_angles: BondsAndAngles,
    pub bounds_mat_force_scaling: f64,
}

/// `CrystalFF::getExperimentalTorsions` (no constrained atoms).
pub fn get_experimental_torsions(
    v: &RdkitMolView,
    use_exp_torsions: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_basic_knowledge: bool,
    version: u32,
) -> CrystalFfDetails {
    let nb = v.num_bonds();
    let na = v.num_atoms();
    let mut d = CrystalFfDetails {
        bounds_mat_force_scaling: 1.0,
        ..Default::default()
    };
    let mut excluded = vec![false; nb];
    let br = &v.bond_rings;
    for i in 0..br.len() {
        let mut rs1 = vec![false; nb];
        for &b in &br[i] {
            rs1[b] = true;
        }
        for j in i + 1..br.len() {
            if br[i].len() >= MIN_MACROCYCLE_SIZE && br[j].len() >= MIN_MACROCYCLE_SIZE {
                continue;
            }
            let mut n_common = 0;
            for &b in &br[j] {
                if rs1[b] {
                    n_common += 1;
                    if n_common > 1 {
                        break;
                    }
                }
            }
            if n_common > 1 {
                if br[i].len() < MIN_MACROCYCLE_SIZE {
                    for &b in &br[i] {
                        excluded[b] = true;
                    }
                }
                if br[j].len() < MIN_MACROCYCLE_SIZE {
                    for &b in &br[j] {
                        excluded[b] = true;
                    }
                }
            }
        }
    }
    let mut num_bond_rings = vec![0usize; nb];
    for r in br {
        for &b in r {
            num_bond_rings[b] += 1;
        }
    }
    let mut done_bonds = vec![false; nb];
    if use_exp_torsions {
        let target = Target::new(v);
        for p in params(version, use_small_ring_torsions, use_macrocycle_torsions) {
            for m in substruct_match(&target, &p.pattern) {
                let (a1, a2, a3, a4) = (m[p.idx[0]], m[p.idx[1]], m[p.idx[2]], m[p.idx[3]]);
                let bid2 = v.bond_between(a2, a3).expect("bond between central atoms");
                if excluded[bid2] || num_bond_rings[bid2] > 3 {
                    done_bonds[bid2] = true;
                }
                if !done_bonds[bid2] {
                    done_bonds[bid2] = true;
                    d.exp_torsion_atoms.push([a1, a2, a3, a4]);
                    d.exp_torsion_angles.push((p.signs, p.v));
                }
            }
            let _ = &p.smarts;
        }
    }
    if use_basic_knowledge {
        for aid2 in 0..na {
            let at = &v.atoms[aid2];
            let z = at.atomic_num;
            if (z == 6 || z == 7 || z == 8) && at.hybridization == 3 && v.degree(aid2) == 3 {
                let mut atoms = [usize::MAX; 4];
                atoms[1] = aid2;
                let mut i = 0;
                let mut bound_sp2_o = false;
                for x in v.neighbors(aid2) {
                    atoms[i] = x;
                    if !bound_sp2_o {
                        bound_sp2_o =
                            z == 6 && v.atoms[x].atomic_num == 8 && v.atoms[x].hybridization == 3;
                    }
                    if i == 0 {
                        i += 1;
                    }
                    i += 1;
                }
                d.improper_atoms.push([
                    atoms[0],
                    atoms[1],
                    atoms[2],
                    atoms[3],
                    z as usize,
                    usize::from(bound_sp2_o),
                ]);
            }
        }
        for ring in &v.atom_rings {
            let r = ring.len();
            if !(4..=6).contains(&r) {
                continue;
            }
            for i in 0..r {
                let (a1, a2, a3, a4) = (
                    ring[i],
                    ring[(i + 1) % r],
                    ring[(i + 2) % r],
                    ring[(i + 3) % r],
                );
                let bid2 = v.bond_between(a2, a3).expect("ring bond");
                let sp2 = |a: usize| v.atoms[a].hybridization == 3;
                if !done_bonds[bid2] && sp2(a1) && sp2(a2) && sp2(a3) && sp2(a4) {
                    done_bonds[bid2] = true;
                    d.exp_torsion_atoms.push([a1, a2, a3, a4]);
                    let mut signs = [1i32; 6];
                    signs[1] = -1;
                    let mut fc = [0f64; 6];
                    fc[1] = 100.0;
                    d.exp_torsion_angles.push((signs, fc));
                }
            }
        }
    }
    d
}

// ---- force-field terms ------------------------------------------------------

type P3 = [f64; 3];

#[inline]
fn pt(pos: &[f64], i: usize) -> P3 {
    [pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]]
}
#[inline]
fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        -a[0] * b[2] + a[2] * b[0],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn len_sq(a: P3) -> f64 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}
#[inline]
fn length(a: P3) -> f64 {
    len_sq(a).sqrt()
}
#[inline]
fn div(a: P3, s: f64) -> P3 {
    [a[0] / s, a[1] / s, a[2] / s]
}
#[inline]
fn is_double_zero(x: f64) -> bool {
    x < 1.0e-10 && x > -1.0e-10
}

struct Torsion {
    idx: [usize; 4],
    fc: [f64; 6],
    signs: [i32; 6],
}

struct Inversion {
    idx: [usize; 4],
    c0: f64,
    c1: f64,
    c2: f64,
    force_constant: f64,
}

struct DistC {
    i: usize,
    j: usize,
    min_len: f64,
    max_len: f64,
    fc: f64,
}

struct AngleC {
    idx: [usize; 3],
    min_angle: f64,
    max_angle: f64,
    fc: f64,
}

enum Contrib {
    Torsions(Vec<Torsion>),
    Inversions(Vec<Inversion>),
    Distances(Vec<DistC>),
    Angles(Vec<AngleC>),
}

/// A 3D `ForceFields::ForceField` with the ETKDG contribution kinds.
pub(crate) struct Etk3dField {
    contribs: Vec<Contrib>,
}

fn torsion_cos_phi(i: P3, j: P3, k: P3, l: P3) -> f64 {
    let r1 = sub(i, j);
    let r2 = sub(k, j);
    let r3 = sub(j, k);
    let r4 = sub(l, k);
    let t1 = cross(r1, r2);
    let t2 = cross(r3, r4);
    let t1_len = length(t1);
    let t2_len = length(t2);
    if is_double_zero(t1_len) || is_double_zero(t2_len) {
        return 0.0;
    }
    (dot(t1, t2) / (t1_len * t2_len)).clamp(-1.0, 1.0)
}

fn torsion_energy_m6(fc: &[f64; 6], s: &[i32; 6], cos_phi: f64) -> f64 {
    let c2 = cos_phi * cos_phi;
    let c3 = cos_phi * c2;
    let c4 = cos_phi * c3;
    let c5 = cos_phi * c4;
    let c6 = cos_phi * c5;
    let cos2 = 2.0 * c2 - 1.0;
    let cos3 = 4.0 * c3 - 3.0 * cos_phi;
    let cos4 = 8.0 * c4 - 8.0 * c2 + 1.0;
    let cos5 = 16.0 * c5 - 20.0 * c3 + 5.0 * cos_phi;
    let cos6 = 32.0 * c6 - 48.0 * c4 + 18.0 * c2 - 1.0;
    fc[0] * (1.0 + f64::from(s[0]) * cos_phi)
        + fc[1] * (1.0 + f64::from(s[1]) * cos2)
        + fc[2] * (1.0 + f64::from(s[2]) * cos3)
        + fc[3] * (1.0 + f64::from(s[3]) * cos4)
        + fc[4] * (1.0 + f64::from(s[4]) * cos5)
        + fc[5] * (1.0 + f64::from(s[5]) * cos6)
}

/// `UFF::Utils::calculateCosY`.
fn calculate_cos_y(i: P3, j: P3, k: P3, l: P3) -> f64 {
    const ZERO_TOL: f64 = 1.0e-16;
    let rji = sub(i, j);
    let rjk = sub(k, j);
    let rjl = sub(l, j);
    let l2ji = len_sq(rji);
    let l2jk = len_sq(rjk);
    let l2jl = len_sq(rjl);
    if l2ji < ZERO_TOL || l2jk < ZERO_TOL || l2jl < ZERO_TOL {
        return 0.0;
    }
    let n = div(cross(rji, rjk), l2ji.sqrt() * l2jk.sqrt());
    let l2n = len_sq(n);
    if l2n < ZERO_TOL {
        return 0.0;
    }
    dot(n, rjl) / (l2jl.sqrt() * l2n.sqrt())
}

fn normalize(a: P3) -> P3 {
    let l = length(a);
    // RDKit throws below 1e-16; callers guard against zero lengths.
    div(a, l)
}

const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

impl Etk3dField {
    pub(crate) fn energy(&self, pos: &[f64]) -> f64 {
        let mut res = 0.0;
        for c in &self.contribs {
            res += match c {
                Contrib::Torsions(v) => {
                    let mut acc = 0.0;
                    for t in v {
                        let [a, b, cc, dd] = t.idx.map(|i| pt(pos, i));
                        acc += torsion_energy_m6(&t.fc, &t.signs, torsion_cos_phi(a, b, cc, dd));
                    }
                    acc
                }
                Contrib::Inversions(v) => {
                    let mut acc = 0.0;
                    for t in v {
                        let [p1, p2, p3, p4] = t.idx.map(|i| pt(pos, i));
                        let cos_y = calculate_cos_y(p1, p2, p3, p4);
                        let sin_y_sq = 1.0 - cos_y * cos_y;
                        let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
                        let cos2w = 2.0 * sin_y * sin_y - 1.0;
                        acc += t.force_constant * (t.c0 + t.c1 * sin_y + t.c2 * cos2w);
                    }
                    acc
                }
                Contrib::Distances(v) => {
                    let mut acc = 0.0;
                    for c in v {
                        let d2 = len_sq(sub(pt(pos, c.i.min(c.j)), pt(pos, c.i.max(c.j))));
                        let diff = if d2 < c.min_len * c.min_len {
                            c.min_len - d2.sqrt()
                        } else if d2 > c.max_len * c.max_len {
                            d2.sqrt() - c.max_len
                        } else {
                            continue;
                        };
                        acc += 0.5 * c.fc * diff * diff;
                    }
                    acc
                }
                Contrib::Angles(v) => {
                    let mut acc = 0.0;
                    for c in v {
                        let angle = angle_deg(pos, c.idx);
                        let term = angle_term(angle, c);
                        acc += c.fc * term * term;
                    }
                    acc
                }
            };
        }
        res
    }

    pub(crate) fn grad(&self, pos: &[f64], grad: &mut [f64]) {
        for c in &self.contribs {
            match c {
                Contrib::Torsions(v) => {
                    for t in v {
                        let [ip, jp, kp, lp] = t.idx.map(|i| pt(pos, i));
                        let r = [sub(ip, jp), sub(kp, jp), sub(jp, kp), sub(lp, kp)];
                        let mut tt = [cross(r[0], r[1]), cross(r[2], r[3])];
                        let d = [length(tt[0]), length(tt[1])];
                        if is_double_zero(d[0]) || is_double_zero(d[1]) {
                            // RDKit returns from the whole contribution here.
                            break;
                        }
                        tt[0] = div(tt[0], d[0]);
                        tt[1] = div(tt[1], d[1]);
                        let cos_phi = dot(tt[0], tt[1]).clamp(-1.0, 1.0);
                        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
                        let sin_phi = if sin_phi_sq > 0.0 {
                            sin_phi_sq.sqrt()
                        } else {
                            0.0
                        };
                        let c2 = cos_phi * cos_phi;
                        let c3 = cos_phi * c2;
                        let c4 = cos_phi * c3;
                        let c5 = cos_phi * c4;
                        let fc = &t.fc;
                        let s = t.signs.map(f64::from);
                        let de_dphi = -fc[0] * s[0] * sin_phi
                            - 2.0 * fc[1] * s[1] * (2.0 * cos_phi * sin_phi)
                            - 3.0 * fc[2] * s[2] * (4.0 * c2 * sin_phi - sin_phi)
                            - 4.0 * fc[3] * s[3] * (8.0 * c3 * sin_phi - 4.0 * cos_phi * sin_phi)
                            - 5.0
                                * fc[4]
                                * s[4]
                                * (16.0 * c4 * sin_phi - 12.0 * c2 * sin_phi + sin_phi)
                            - 6.0
                                * fc[4]
                                * s[4]
                                * (32.0 * c5 * sin_phi - 32.0 * c3 * sin_phi + 6.0 * sin_phi);
                        let sin_term = -de_dphi
                            * (if is_double_zero(sin_phi) {
                                1.0 / cos_phi
                            } else {
                                1.0 / sin_phi
                            });
                        torsion_grad(&r, &tt, &d, t.idx, grad, sin_term, cos_phi);
                    }
                }
                Contrib::Inversions(v) => {
                    for t in v {
                        let [p1, p2, p3, p4] = t.idx.map(|i| pt(pos, i));
                        let mut rji = sub(p1, p2);
                        let mut rjk = sub(p3, p2);
                        let mut rjl = sub(p4, p2);
                        let dji = length(rji);
                        let djk = length(rjk);
                        let djl = length(rjl);
                        if is_double_zero(dji) || is_double_zero(djk) || is_double_zero(djl) {
                            break;
                        }
                        rji = normalize(rji);
                        rjk = normalize(rjk);
                        rjl = normalize(rjl);
                        let n = normalize(cross([-rji[0], -rji[1], -rji[2]], rjk));
                        let cos_y = dot(n, rjl).clamp(-1.0, 1.0);
                        let sin_y_sq = 1.0 - cos_y * cos_y;
                        let sin_y = f64::max(sin_y_sq.sqrt(), 1.0e-8);
                        let cos_theta = dot(rji, rjk).clamp(-1.0, 1.0);
                        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
                        let sin_theta = f64::max(sin_theta_sq.sqrt(), 1.0e-8);
                        let de_dw = -t.force_constant * (t.c1 * cos_y - 4.0 * t.c2 * cos_y * sin_y);
                        let t1 = cross(rjl, rjk);
                        let t2 = cross(rji, rjl);
                        let t3 = cross(rjk, rji);
                        let term1 = sin_y * sin_theta;
                        let term2 = cos_y / (sin_y * sin_theta_sq);
                        let mut tg1 = [0.0; 3];
                        let mut tg3 = [0.0; 3];
                        let mut tg4 = [0.0; 3];
                        for k in 0..3 {
                            tg1[k] = (t1[k] / term1 - (rji[k] - rjk[k] * cos_theta) * term2) / dji;
                            tg3[k] = (t2[k] / term1 - (rjk[k] - rji[k] * cos_theta) * term2) / djk;
                            tg4[k] = (t3[k] / term1 - rjl[k] * cos_y / sin_y) / djl;
                        }
                        for k in 0..3 {
                            grad[3 * t.idx[0] + k] += de_dw * tg1[k];
                            grad[3 * t.idx[1] + k] += -de_dw * (tg1[k] + tg3[k] + tg4[k]);
                            grad[3 * t.idx[2] + k] += de_dw * tg3[k];
                            grad[3 * t.idx[3] + k] += de_dw * tg4[k];
                        }
                    }
                }
                Contrib::Distances(v) => {
                    for c in v {
                        let d2 = len_sq(sub(pt(pos, c.i.min(c.j)), pt(pos, c.i.max(c.j))));
                        let (dist, mut pre);
                        if d2 < c.min_len * c.min_len {
                            dist = d2.sqrt();
                            pre = dist - c.min_len;
                        } else if d2 > c.max_len * c.max_len {
                            dist = d2.sqrt();
                            pre = dist - c.max_len;
                        } else {
                            continue;
                        }
                        pre *= c.fc;
                        pre /= f64::max(1.0e-8, dist);
                        for k in 0..3 {
                            let g = pre * (pos[3 * c.i + k] - pos[3 * c.j + k]);
                            grad[3 * c.i + k] += g;
                            grad[3 * c.j + k] -= g;
                        }
                    }
                }
                Contrib::Angles(v) => {
                    for c in v {
                        let [p1, p2, p3] = c.idx.map(|i| pt(pos, i));
                        let r = [sub(p1, p2), sub(p3, p2)];
                        let rl = [
                            f64::max(1.0e-5, len_sq(r[0])),
                            f64::max(1.0e-5, len_sq(r[1])),
                        ];
                        let cos_t = (dot(r[0], r[1]) / (rl[0] * rl[1]).sqrt()).clamp(-1.0, 1.0);
                        let angle = RAD2DEG * cos_t.acos();
                        let term = angle_term(angle, c);
                        let de_dtheta = 2.0 * RAD2DEG * c.fc * term;
                        let rp = cross(r[1], r[0]);
                        let prefactor = de_dtheta / f64::max(1.0e-5, length(rp));
                        let t = [-prefactor / rl[0], prefactor / rl[1]];
                        let d0 = cross(r[0], rp).map(|x| x * t[0]);
                        let d2 = cross(r[1], rp).map(|x| x * t[1]);
                        let d1 = [-d0[0] - d2[0], -d0[1] - d2[1], -d0[2] - d2[2]];
                        for (g, dd) in c.idx.iter().zip([d0, d1, d2]) {
                            grad[3 * g] += dd[0];
                            grad[3 * g + 1] += dd[1];
                            grad[3 * g + 2] += dd[2];
                        }
                    }
                }
            }
        }
    }

    pub(crate) fn minimize(&self, pos: &mut [f64], max_its: u32, force_tol: f64) -> i32 {
        if self.contribs.is_empty() {
            return 0;
        }
        chematic_ff::rdkit_bfgs::bfgs_minimize(
            pos,
            force_tol,
            max_its,
            |p| self.energy(p),
            |p, g| {
                g.fill(0.0);
                self.grad(p, g);
                chematic_ff::rdkit_bfgs::scale_gradient(g)
            },
        )
    }
}

fn angle_deg(pos: &[f64], idx: [usize; 3]) -> f64 {
    let [p1, p2, p3] = idx.map(|i| pt(pos, i));
    let r = [sub(p1, p2), sub(p3, p2)];
    let rl = [
        f64::max(1.0e-5, len_sq(r[0])),
        f64::max(1.0e-5, len_sq(r[1])),
    ];
    let cos_t = (dot(r[0], r[1]) / (rl[0] * rl[1]).sqrt()).clamp(-1.0, 1.0);
    RAD2DEG * cos_t.acos()
}

fn angle_term(angle: f64, c: &AngleC) -> f64 {
    if angle < c.min_angle {
        angle - c.min_angle
    } else if angle > c.max_angle {
        angle - c.max_angle
    } else {
        0.0
    }
}

/// `MMFF::Utils::calcTorsionGrad`.
fn torsion_grad(
    r: &[P3; 4],
    t: &[P3; 2],
    d: &[f64; 2],
    idx: [usize; 4],
    grad: &mut [f64],
    sin_term: f64,
    cos_phi: f64,
) {
    let (x, y, z) = (0, 1, 2);
    let dc = [
        1.0 / d[0] * (t[1][x] - cos_phi * t[0][x]),
        1.0 / d[0] * (t[1][y] - cos_phi * t[0][y]),
        1.0 / d[0] * (t[1][z] - cos_phi * t[0][z]),
        1.0 / d[1] * (t[0][x] - cos_phi * t[1][x]),
        1.0 / d[1] * (t[0][y] - cos_phi * t[1][y]),
        1.0 / d[1] * (t[0][z] - cos_phi * t[1][z]),
    ];
    let g = |k: usize, c: usize| 3 * idx[k] + c;
    grad[g(0, 0)] += sin_term * (dc[2] * r[1][y] - dc[1] * r[1][z]);
    grad[g(0, 1)] += sin_term * (dc[0] * r[1][z] - dc[2] * r[1][x]);
    grad[g(0, 2)] += sin_term * (dc[1] * r[1][x] - dc[0] * r[1][y]);

    grad[g(1, 0)] += sin_term
        * (dc[1] * (r[1][z] - r[0][z])
            + dc[2] * (r[0][y] - r[1][y])
            + dc[4] * (-r[3][z])
            + dc[5] * (r[3][y]));
    grad[g(1, 1)] += sin_term
        * (dc[0] * (r[0][z] - r[1][z])
            + dc[2] * (r[1][x] - r[0][x])
            + dc[3] * (r[3][z])
            + dc[5] * (-r[3][x]));
    grad[g(1, 2)] += sin_term
        * (dc[0] * (r[1][y] - r[0][y])
            + dc[1] * (r[0][x] - r[1][x])
            + dc[3] * (-r[3][y])
            + dc[4] * (r[3][x]));

    grad[g(2, 0)] += sin_term
        * (dc[1] * (r[0][z])
            + dc[2] * (-r[0][y])
            + dc[4] * (r[3][z] - r[2][z])
            + dc[5] * (r[2][y] - r[3][y]));
    grad[g(2, 1)] += sin_term
        * (dc[0] * (-r[0][z])
            + dc[2] * (r[0][x])
            + dc[3] * (r[2][z] - r[3][z])
            + dc[5] * (r[3][x] - r[2][x]));
    grad[g(2, 2)] += sin_term
        * (dc[0] * (r[0][y])
            + dc[1] * (-r[0][x])
            + dc[3] * (r[3][y] - r[2][y])
            + dc[4] * (r[2][x] - r[3][x]));

    grad[g(3, 0)] += sin_term * (dc[4] * r[2][z] - dc[5] * r[2][y]);
    grad[g(3, 1)] += sin_term * (dc[5] * r[2][x] - dc[3] * r[2][z]);
    grad[g(3, 2)] += sin_term * (dc[3] * r[2][y] - dc[4] * r[2][x]);
}

/// UFF `calcInversionCoefficientsAndForceConstant` for C/N/O centres.
fn inversion(idx: [usize; 4], z: usize, c_bound_to_o: bool, scale: f64) -> Inversion {
    debug_assert!(matches!(z, 6..=8));
    let res = (if c_bound_to_o { 50.0 } else { 6.0 }) / 3.0;
    Inversion {
        idx,
        c0: 1.0,
        c1: -1.0,
        c2: 0.0,
        force_constant: res * scale,
    }
}

/// `addImproperTorsionTerms`.
fn improper_terms(improper: &[[usize; 6]], scale: f64, constrained: &mut [bool]) -> Vec<Inversion> {
    let mut out = Vec::new();
    for ia in improper {
        for i in 0..3 {
            let n = match i {
                0 => [0, 1, 2, 3],
                1 => [0, 1, 3, 2],
                _ => [2, 1, 3, 0],
            };
            out.push(inversion(
                [ia[n[0]], ia[n[1]], ia[n[2]], ia[n[3]]],
                ia[4],
                ia[5] != 0,
                scale,
            ));
            constrained[ia[n[1]]] = true;
        }
    }
    out
}

/// `DistGeom::construct3DForceField` (no CPCI).
pub(crate) fn construct_3d_force_field(
    m: &BoundsMatrix,
    pos: &[f64],
    d: &CrystalFfDetails,
) -> Etk3dField {
    let n = m.num_rows();
    let mut atom_pairs = vec![false; n * n];
    let mut improper_constrained = vec![false; n];
    let mut contribs = Vec::new();

    let mut tors = Vec::new();
    for (t, a) in d.exp_torsion_atoms.iter().enumerate() {
        let (i, l) = (a[0], a[3]);
        atom_pairs[i.min(l) * n + i.max(l)] = true;
        let (signs, fc) = d.exp_torsion_angles[t];
        tors.push(Torsion { idx: *a, fc, signs });
    }
    if !tors.is_empty() {
        contribs.push(Contrib::Torsions(tors));
    }
    let inv = improper_terms(&d.improper_atoms, 10.0, &mut improper_constrained);
    if !inv.is_empty() {
        contribs.push(Contrib::Inversions(inv));
    }
    let mut d12 = Vec::new();
    for &(i, j) in &d.bonds_angles.bonds {
        atom_pairs[i.min(j) * n + i.max(j)] = true;
        let dd = length(sub(pt(pos, i), pt(pos, j)));
        d12.push(DistC {
            i,
            j,
            min_len: dd - KNOWN_DIST_TOL,
            max_len: dd + KNOWN_DIST_TOL,
            fc: KNOWN_DIST_FORCE_CONSTANT,
        });
    }
    if !d12.is_empty() {
        contribs.push(Contrib::Distances(d12));
    }
    let mut d13 = Vec::new();
    let mut a13 = Vec::new();
    for a in &d.bonds_angles.angles {
        let (i, j, k) = (a[0], a[1], a[2]);
        atom_pairs[i.min(k) * n + i.max(k)] = true;
        if a[3] != 0 {
            a13.push(AngleC {
                idx: [i, j, k],
                min_angle: 179.0,
                max_angle: 180.0,
                fc: 1.0,
            });
        } else if improper_constrained[j] {
            d13.push(DistC {
                i,
                j: k,
                min_len: m.lower(i, k),
                max_len: m.upper(i, k),
                fc: KNOWN_DIST_FORCE_CONSTANT,
            });
        } else {
            let dd = length(sub(pt(pos, i), pt(pos, k)));
            d13.push(DistC {
                i,
                j: k,
                min_len: dd - KNOWN_DIST_TOL,
                max_len: dd + KNOWN_DIST_TOL,
                fc: KNOWN_DIST_FORCE_CONSTANT,
            });
        }
    }
    if !a13.is_empty() {
        contribs.push(Contrib::Angles(a13));
    }
    if !d13.is_empty() {
        contribs.push(Contrib::Distances(d13));
    }
    let mut lr = Vec::new();
    for i in 1..n {
        for j in 0..i {
            if !atom_pairs[j * n + i] {
                let fdist = d.bounds_mat_force_scaling * 10.0;
                lr.push(DistC {
                    i,
                    j,
                    min_len: m.lower(i, j),
                    max_len: m.upper(i, j),
                    fc: fdist,
                });
            }
        }
    }
    if !lr.is_empty() {
        contribs.push(Contrib::Distances(lr));
    }
    Etk3dField { contribs }
}

/// `DistGeom::construct3DImproperForceField`.
pub(crate) fn construct_3d_improper_force_field(n: usize, d: &CrystalFfDetails) -> Etk3dField {
    let mut constrained = vec![false; n];
    let mut contribs = Vec::new();
    let inv = improper_terms(&d.improper_atoms, 10.0, &mut constrained);
    if !inv.is_empty() {
        contribs.push(Contrib::Inversions(inv));
    }
    let mut ang = Vec::new();
    for a in &d.bonds_angles.angles {
        if a[3] != 0 {
            ang.push(AngleC {
                idx: [a[0], a[1], a[2]],
                min_angle: 179.0,
                max_angle: 180.0,
                fc: 10.0,
            });
        }
    }
    if !ang.is_empty() {
        contribs.push(Contrib::Angles(ang));
    }
    Etk3dField { contribs }
}

/// `EmbeddingOps::minimizeWithExpTorsions` (ETKDG: useBasicKnowledge).
/// Returns RDKit's planarity verdict, or `None` when RDKit's minimizer throws.
pub(crate) fn minimize_with_exp_torsions(
    pos3: &mut [f64],
    m: &BoundsMatrix,
    d: &CrystalFfDetails,
    force_tol: f64,
    use_basic_knowledge: bool,
) -> Option<bool> {
    if !use_basic_knowledge {
        unimplemented!("plain ETDG (constructPlain3DForceField) is not ported");
    }
    let field = construct_3d_force_field(m, pos3, d);
    if field.energy(pos3) > ERROR_TOL && field.minimize(pos3, 300, force_tol) < 0 {
        return None;
    }
    let field2 = construct_3d_improper_force_field(m.num_rows(), d);
    let planarity_tolerance = 0.7;
    Some(field2.energy(pos3) <= d.improper_atoms.len() as f64 * planarity_tolerance)
}

#[cfg(test)]
mod numerical_contract_tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
            "{actual} != {expected}"
        );
    }

    fn check_gradient(field: &Etk3dField, pos: &[f64]) {
        let mut grad = vec![0.0; pos.len()];
        field.grad(pos, &mut grad);
        for i in 0..pos.len() {
            let mut plus = pos.to_vec();
            let mut minus = pos.to_vec();
            plus[i] += 1e-6;
            minus[i] -= 1e-6;
            let numerical = (field.energy(&plus) - field.energy(&minus)) / 2e-6;
            close(grad[i], numerical, 2e-5);
        }
        for axis in 0..3 {
            close(grad.iter().skip(axis).step_by(3).sum(), 0.0, 1e-8);
        }
    }

    #[test]
    fn interval_angle_restraints_have_expected_energy_and_gradient() {
        let field = Etk3dField {
            contribs: vec![Contrib::Angles(vec![AngleC {
                idx: [0, 1, 2],
                min_angle: 60.0,
                max_angle: 120.0,
                fc: 2.0,
            }])],
        };
        for (degrees, expected) in [
            (30.0_f64, 1800.0),
            (60.0, 0.0),
            (90.0, 0.0),
            (120.0, 0.0),
            (150.0, 1800.0),
        ] {
            let theta = degrees.to_radians();
            let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, theta.cos(), theta.sin(), 0.0];
            close(field.energy(&pos), expected, 1e-10);
            if degrees != 60.0 && degrees != 120.0 {
                check_gradient(&field, &pos);
            } else {
                // Central differences cross the piecewise restraint boundary.
                let mut grad = [0.0; 9];
                field.grad(&pos, &mut grad);
                assert!(grad.iter().all(|v| v.abs() < 1e-8));
            }
            let translated: Vec<_> = pos
                .chunks_exact(3)
                .flat_map(|p| [p[0] + 3.0, p[1] - 2.0, p[2] + 4.0])
                .collect();
            close(field.energy(&translated), expected, 1e-10);
            let rotated: Vec<_> = pos
                .chunks_exact(3)
                .flat_map(|p| [-p[1], p[0], p[2]])
                .collect();
            close(field.energy(&rotated), expected, 1e-10);
        }
        // Coincident/collinear coordinates must remain finite at the guarded denominator.
        for pos in [[0.0; 9], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0]] {
            assert!(field.energy(&pos).is_finite());
            let mut grad = [0.0; 9];
            field.grad(&pos, &mut grad);
            assert!(grad.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn distance_interval_is_symmetric_and_its_gradient_matches_energy() {
        for (i, j) in [(0, 1), (1, 0)] {
            let field = Etk3dField {
                contribs: vec![Contrib::Distances(vec![DistC {
                    i,
                    j,
                    min_len: 1.0,
                    max_len: 2.0,
                    fc: 4.0,
                }])],
            };
            for (distance, expected) in [(0.5, 0.5), (1.5, 0.0), (3.0, 2.0)] {
                let pos = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
                close(field.energy(&pos), expected, 1e-12);
                check_gradient(&field, &pos);
            }
            let mut coincident = [0.0; 6];
            field.grad(&[0.0; 6], &mut coincident);
            assert_eq!(coincident, [0.0; 6]);
        }
    }

    #[test]
    fn empty_force_field_is_already_minimized() {
        let field = Etk3dField { contribs: vec![] };
        let mut pos = [1.0, 2.0, 3.0];
        assert_eq!(field.energy(&pos), 0.0);
        assert_eq!(field.minimize(&mut pos, 100, 1e-4), 0);
        assert_eq!(pos, [1.0, 2.0, 3.0]);
        let mut grad = [0.0; 3];
        field.grad(&pos, &mut grad);
        assert_eq!(grad, [0.0; 3]);
    }

    #[test]
    fn inversion_normal_and_degenerate_geometry_contracts() {
        let field = Etk3dField {
            contribs: vec![Contrib::Inversions(vec![Inversion {
                idx: [0, 1, 2, 3],
                c0: 1.0,
                c1: -1.0,
                c2: 0.0,
                force_constant: 6.0,
            }])],
        };
        let planar = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, -1.0, 0.0];
        close(field.energy(&planar), 0.0, 1e-12);
        check_gradient(&field, &planar);
        let nonplanar = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, -1.0, 0.5];
        check_gradient(&field, &nonplanar);
        let orthogonal = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        close(field.energy(&orthogonal), 6.0, 1e-12);
        for pos in [
            [0.0; 12],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        ] {
            assert!(field.energy(&pos).is_finite());
            let mut grad = [0.0; 12];
            field.grad(&pos, &mut grad);
            assert!(grad.iter().all(|v| v.is_finite()));
        }
    }
}
