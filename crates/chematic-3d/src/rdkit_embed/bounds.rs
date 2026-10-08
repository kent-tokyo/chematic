//! `DGeomHelpers::setTopolBounds` (RDKit 2026.03 `BoundsMatrixBuilder.cpp`)
//! on an [`RdkitMolView`], plus `collectBondsAndAngles`.

use std::collections::HashSet;
use std::f64::consts::PI as M_PI;

use chematic_smiles::RdkitMolView;

use super::dg::BoundsMatrix;
use super::stdsort::std_sort;

const DIST12_DELTA: f64 = 0.01;
const DIST13_TOL: f64 = 0.04;
const GEN_DIST_TOL: f64 = 0.06;
const DIST15_TOL: f64 = 0.08;
const VDW_SCALE_15: f64 = 0.7;
const H_BOND_LENGTH: f64 = 1.8;
const MAX_UPPER: f64 = 1000.0;
const MIN_MACROCYCLE_RING_SIZE: usize = 9;

// RDKit enum values used by the view.
const SINGLE: u8 = 1;
const DOUBLE: u8 = 2;
const TRIPLE: u8 = 3;
const HYB_SP: u8 = 2;
const HYB_SP2: u8 = 3;
const HYB_SP3: u8 = 4;
const HYB_SP3D: u8 = 5;
const HYB_SP3D2: u8 = 6;
const STEREO_ANY: u8 = 1;
const STEREO_Z: u8 = 2;
const STEREO_E: u8 = 3;

/// RDKit throws an invariant violation while building the bounds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundsInvariant(pub &'static str);

#[derive(Clone, Copy, PartialEq, Eq)]
enum Path14Type {
    Cis,
    Trans,
    Other,
}

struct Path14 {
    bid1: usize,
    bid2: usize,
    bid3: usize,
    ty: Path14Type,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DistType {
    D12,
    D13,
    D14,
}

/// Packed symmetric matrix (`RDNumeric::SymmMatrix`) with an initial value.
struct Symm<T: Copy> {
    data: Vec<T>,
}

impl<T: Copy> Symm<T> {
    fn new(n: usize, v: T) -> Self {
        Symm {
            data: vec![v; n * (n + 1) / 2],
        }
    }
    fn idx(i: usize, j: usize) -> usize {
        if i >= j {
            i * (i + 1) / 2 + j
        } else {
            j * (j + 1) / 2 + i
        }
    }
    fn get(&self, i: usize, j: usize) -> T {
        self.data[Self::idx(i, j)]
    }
    fn set(&mut self, i: usize, j: usize, v: T) {
        self.data[Self::idx(i, j)] = v;
    }
}

struct ComputedData {
    bond_lengths: Vec<f64>,
    bond_adj: Symm<i64>,
    bond_angles: Symm<f64>,
    paths14: Vec<Path14>,
    cis_paths: HashSet<u64>,
    trans_paths: HashSet<u64>,
    set15_atoms: Vec<bool>,
    visited12: Vec<bool>,
    visited13: Vec<bool>,
    visited14: Vec<bool>,
}

impl ComputedData {
    fn new(na: usize, nb: usize) -> Self {
        ComputedData {
            bond_lengths: vec![0.0; nb],
            bond_adj: Symm::new(nb, -1),
            bond_angles: Symm::new(nb, -1.0),
            paths14: Vec::new(),
            cis_paths: HashSet::new(),
            trans_paths: HashSet::new(),
            set15_atoms: vec![false; na * na],
            visited12: vec![false; na * na],
            visited13: vec![false; na * na],
            visited14: vec![false; na * na],
        }
    }

    fn visited_bound(&self, pid: usize, max: DistType) -> bool {
        (max >= DistType::D12 && self.visited12[pid])
            || (max >= DistType::D13 && self.visited13[pid])
            || (max >= DistType::D14 && self.visited14[pid])
    }

    fn add_path(set: &mut HashSet<u64>, nb: u64, b1: usize, b2: usize, b3: usize) {
        set.insert(b1 as u64 * nb * nb + b2 as u64 * nb + b3 as u64);
        set.insert(b3 as u64 * nb * nb + b2 as u64 * nb + b1 as u64);
    }
}

/// `initBoundsMat(mmat, 0.0, 1000.0)`.
pub fn init_bounds_mat(m: &mut BoundsMatrix) {
    let n = m.num_rows();
    for i in 1..n {
        for j in 0..i {
            m.set_upper(i, j, 1000.0);
            m.set_lower(i, j, 0.0);
        }
    }
}

fn check_and_set_bounds(
    i: usize,
    j: usize,
    lb: f64,
    ub: f64,
    m: &mut BoundsMatrix,
) -> Result<(), BoundsInvariant> {
    let clb = m.lower(i, j);
    let cub = m.upper(i, j);
    if !(ub > lb) {
        return Err(BoundsInvariant("upper bound not greater than lower bound"));
    }
    if !(lb > DIST12_DELTA || clb > DIST12_DELTA) {
        return Err(BoundsInvariant("bad lower bound"));
    }
    if clb <= DIST12_DELTA {
        m.set_lower(i, j, lb);
    } else if lb < clb && lb > DIST12_DELTA {
        m.set_lower(i, j, lb);
    }
    if cub >= MAX_UPPER {
        m.set_upper(i, j, ub);
    } else if ub > cub && ub < MAX_UPPER {
        m.set_upper(i, j, ub);
    }
    Ok(())
}

fn compute13_dist(d1: f64, d2: f64, angle: f64) -> f64 {
    let res = d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * angle.cos();
    res.sqrt()
}

fn compute14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

fn compute14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// `RDGeom::compute14Dist3D` (rotation about the x axis by `tor`).
fn compute14_dist_3d(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64, tor: f64) -> f64 {
    let p1 = [d1 * ang12.cos(), d1 * ang12.sin(), 0.0];
    let p4 = [d2 - d3 * ang23.cos(), d3 * ang23.sin(), 0.0];
    let (c, s) = (tor.cos(), tor.sin());
    // Transform3D::TransformPoint with the X-axis rotation matrix.
    let x = 1.0 * p4[0] + 0.0 * p4[1] + 0.0 * p4[2] + 0.0;
    let y = 0.0 * p4[0] + c * p4[1] + (-s) * p4[2] + 0.0;
    let z = 0.0 * p4[0] + s * p4[1] + c * p4[2] + 0.0;
    let d = [x - p1[0], y - p1[1], z - p1[2]];
    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
}

fn pid(a: usize, b: usize, na: usize) -> usize {
    a.min(b) * na + a.max(b)
}

fn set12_bounds(
    mol: &RdkitMolView,
    labels: &[String],
    m: &mut BoundsMatrix,
    acc: &mut ComputedData,
) {
    let na = mol.num_atoms();
    let mut squish = vec![false; na];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.is_conjugated
            && (mol.atoms[b.begin].atomic_num > 10 || mol.atoms[b.end].atomic_num > 10)
            && mol.is_bond_in_ring_of_size(bi, 5)
        {
            squish[b.begin] = true;
            squish[b.end] = true;
        }
    }
    for (bi, b) in mol.bonds.iter().enumerate() {
        let (beg, end) = (b.begin, b.end);
        let order = b.bond_type_as_double;
        let rest = if order > 0.0 {
            chematic_ff::rdkit_uff::rdkit_uff_bond_rest_length(order, &labels[beg], &labels[end])
        } else {
            None
        };
        if let Some(bl) = rest {
            let extra = if squish[beg] || squish[end] { 0.2 } else { 0.0 };
            acc.bond_lengths[bi] = bl;
            m.set_upper(beg, end, bl + extra + DIST12_DELTA);
            m.set_lower(beg, end, bl - extra - DIST12_DELTA);
        } else {
            let vw1 = mol.atoms[beg].rvdw;
            let vw2 = mol.atoms[end].rvdw;
            let bl = (vw1 + vw2) / 2.0;
            acc.bond_lengths[bi] = bl;
            m.set_upper(beg, end, 1.5 * bl);
            m.set_lower(beg, end, 0.5 * bl);
        }
        acc.visited12[pid(beg, end, na)] = true;
    }
}

fn is_larger_sp2_atom(mol: &RdkitMolView, a: usize) -> bool {
    mol.atoms[a].atomic_num > 13
        && mol.atoms[a].hybridization == HYB_SP2
        && mol.num_atom_rings(a) > 0
}

fn set13_helper(
    aid1: usize,
    aid: usize,
    aid3: usize,
    angle: f64,
    acc: &ComputedData,
    m: &mut BoundsMatrix,
    mol: &RdkitMolView,
) -> Result<(), BoundsInvariant> {
    let bid1 = mol.bond_between(aid1, aid).expect("bond");
    let bid2 = mol.bond_between(aid, aid3).expect("bond");
    let mut dl = compute13_dist(acc.bond_lengths[bid1], acc.bond_lengths[bid2], angle);
    let mut dist_tol = DIST13_TOL;
    if is_larger_sp2_atom(mol, aid1) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom(mol, aid) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom(mol, aid3) {
        dist_tol *= 2.0;
    }
    let du = dl + dist_tol;
    dl -= dist_tol;
    check_and_set_bounds(aid1, aid3, dl, du, m)
}

fn ring_angle(hyb: u8, ring_size: usize) -> f64 {
    if (hyb == HYB_SP2 && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        M_PI * (1.0 - 2.0 / ring_size as f64)
    } else if hyb == HYB_SP3 {
        if ring_size == 5 {
            104.0 * M_PI / 180.0
        } else {
            109.5 * M_PI / 180.0
        }
    } else if hyb == HYB_SP3D {
        105.0 * M_PI / 180.0
    } else if hyb == HYB_SP3D2 {
        90.0 * M_PI / 180.0
    } else {
        120.0 * M_PI / 180.0
    }
}

fn set13_bounds(
    mol: &RdkitMolView,
    m: &mut BoundsMatrix,
    acc: &mut ComputedData,
) -> Result<(), BoundsInvariant> {
    let npt = mol.num_atoms();
    let mut atom_rings = mol.atom_rings.clone();
    std_sort(&mut atom_rings, &|a: &Vec<usize>, b: &Vec<usize>| {
        a.len() < b.len()
    });
    let mut visited = vec![0i64; npt];
    let mut angle_taken = vec![0.0f64; npt];
    let nb = mol.num_bonds();
    let mut done_paths = vec![false; nb * nb];
    for ring in &atom_rings {
        let r = ring.len();
        let mut aid1 = ring[r - 1];
        for i in 0..r {
            let aid2 = ring[i];
            let aid3 = if i == r - 1 { ring[0] } else { ring[i + 1] };
            let bid1 = mol
                .bond_between(aid1, aid2)
                .ok_or(BoundsInvariant("no bond found"))?;
            let bid2 = mol
                .bond_between(aid2, aid3)
                .ok_or(BoundsInvariant("no bond found"))?;
            let id1 = nb * bid1 + bid2;
            let id2 = nb * bid2 + bid1;
            let p = pid(aid1, aid3, npt);
            if !done_paths[id1] && !done_paths[id2] {
                let angle = ring_angle(mol.atoms[aid2].hybridization, r);
                if !acc.visited_bound(p, DistType::D12) {
                    set13_helper(aid1, aid2, aid3, angle, acc, m, mol)?;
                    acc.visited13[p] = true;
                }
                acc.bond_angles.set(bid1, bid2, angle);
                acc.bond_adj.set(bid1, bid2, aid2 as i64);
                visited[aid2] += 1;
                angle_taken[aid2] += angle;
                done_paths[id1] = true;
                done_paths[id2] = true;
            }
            aid1 = aid2;
        }
    }
    for aid2 in 0..npt {
        let deg = mol.degree(aid2);
        let n13 = (deg * deg.saturating_sub(1) / 2) as i64;
        if n13 == visited[aid2] {
            continue;
        }
        let ahyb = mol.atoms[aid2].hybridization;
        let bonds = &mol.atom_bonds[aid2];
        if visited[aid2] >= 1 {
            for (k1, &bid1) in bonds.iter().enumerate() {
                let aid1 = mol.other_atom(bid1, aid2);
                for &bid2 in &bonds[..k1] {
                    let aid3 = mol.other_atom(bid2, aid2);
                    if acc.bond_angles.get(bid1, bid2) < 0.0 {
                        let angle = if ahyb == HYB_SP2 {
                            (2.0 * M_PI - angle_taken[aid2]) / ((n13 - visited[aid2]) as u32 as f64)
                        } else if ahyb == HYB_SP3 {
                            if mol.is_atom_in_ring_of_size(aid2, 3) {
                                116.0 * M_PI / 180.0
                            } else if mol.is_atom_in_ring_of_size(aid2, 4) {
                                112.0 * M_PI / 180.0
                            } else {
                                109.5 * M_PI / 180.0
                            }
                        } else if deg == 5 {
                            105.0 * M_PI / 180.0
                        } else if deg == 6 {
                            135.0 * M_PI / 180.0
                        } else {
                            120.0 * M_PI / 180.0
                        };
                        let p = pid(aid1, aid3, npt);
                        if !acc.visited_bound(p, DistType::D12) {
                            set13_helper(aid1, aid2, aid3, angle, acc, m, mol)?;
                            acc.visited13[p] = true;
                        }
                        acc.bond_angles.set(bid1, bid2, angle);
                        acc.bond_adj.set(bid1, bid2, aid2 as i64);
                        angle_taken[aid2] += angle;
                        visited[aid2] += 1;
                    }
                }
            }
        } else {
            for (k1, &bid1) in bonds.iter().enumerate() {
                let aid1 = mol.other_atom(bid1, aid2);
                for &bid2 in &bonds[..k1] {
                    let aid3 = mol.other_atom(bid2, aid2);
                    let angle = if ahyb == HYB_SP {
                        M_PI
                    } else if ahyb == HYB_SP2 {
                        2.0 * M_PI / 3.0
                    } else if ahyb == HYB_SP3 {
                        109.5 * M_PI / 180.0
                    } else if ahyb == HYB_SP3D {
                        105.0 * M_PI / 180.0
                    } else if ahyb == HYB_SP3D2 {
                        135.0 * M_PI / 180.0
                    } else {
                        120.0 * M_PI / 180.0
                    };
                    let p = pid(aid1, aid3, npt);
                    if !acc.visited_bound(p, DistType::D12) {
                        if deg <= 4 {
                            set13_helper(aid1, aid2, aid3, angle, acc, m, mol)?;
                        } else {
                            let dmax = acc.bond_lengths[bid1] + acc.bond_lengths[bid2];
                            check_and_set_bounds(aid1, aid3, 1.0, dmax * 1.2, m)?;
                        }
                        acc.visited13[p] = true;
                    }
                    acc.bond_angles.set(bid1, bid2, angle);
                    acc.bond_adj.set(bid1, bid2, aid2 as i64);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                }
            }
        }
    }
    Ok(())
}

/// `_getAtomStereo`.
fn atom_stereo(mol: &RdkitMolView, b: usize, aid1: usize, aid4: usize) -> u8 {
    let bond = &mol.bonds[b];
    let mut st = bond.stereo;
    if st > STEREO_ANY
        && bond.stereo_atoms.len() >= 2
        && ((bond.stereo_atoms[0] != aid1) ^ (bond.stereo_atoms[1] != aid4))
    {
        // Z <-> E (CIS/TRANS do not occur in the view).
        st = match st {
            STEREO_Z => STEREO_E,
            STEREO_E => STEREO_Z,
            s => s,
        };
    }
    st
}

struct Ctx<'a> {
    mol: &'a RdkitMolView,
    dmat: &'a [f64],
    na: usize,
    nb: u64,
}

impl Ctx<'_> {
    fn d(&self, a: usize, b: usize) -> f64 {
        self.dmat[a.max(b) * self.na + a.min(b)]
    }
}

/// Atoms (aid1, atm2, atm3, aid4) of the path bid1-bid2-bid3.
fn path_atoms(
    cx: &Ctx,
    acc: &ComputedData,
    b1: usize,
    b2: usize,
    b3: usize,
) -> (usize, usize, usize, usize) {
    let atm2 = acc.bond_adj.get(b1, b2) as usize;
    let atm3 = acc.bond_adj.get(b2, b3) as usize;
    (
        cx.mol.other_atom(b1, atm2),
        atm2,
        atm3,
        cx.mol.other_atom(b3, atm3),
    )
}

fn path_geom(
    acc: &ComputedData,
    b1: usize,
    b2: usize,
    b3: usize,
) -> Result<(f64, f64, f64, f64, f64), BoundsInvariant> {
    let ba12 = acc.bond_angles.get(b1, b2);
    let ba23 = acc.bond_angles.get(b2, b3);
    if !(ba12 > 0.0) || !(ba23 > 0.0) {
        return Err(BoundsInvariant("bad bond angle"));
    }
    Ok((
        acc.bond_lengths[b1],
        acc.bond_lengths[b2],
        acc.bond_lengths[b3],
        ba12,
        ba23,
    ))
}

fn set_in_ring14(
    cx: &Ctx,
    b1: usize,
    b2: usize,
    b3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    ring_size: usize,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let (aid1, atm2, atm3, aid4) = path_atoms(cx, acc, b1, b2, b3);
    let ahyb2 = mol.atoms[atm2].hybridization;
    let ahyb3 = mol.atoms[atm3].hybridization;
    let p = pid(aid1, aid4, cx.na);
    if acc.visited_bound(p, DistType::D13) {
        return Ok(());
    }
    if cx.d(aid1, aid4) < 2.9 {
        return Ok(());
    }
    let (bl1, bl2, bl3, ba12, ba23) = path_geom(acc, b1, b2, b3)?;
    let stype = atom_stereo(mol, b2, aid1, aid4);
    let mut prefer_cis = false;
    let mut prefer_trans = false;
    if ring_size <= 8 && ahyb2 == HYB_SP2 && ahyb3 == HYB_SP2 && stype != STEREO_E {
        if mol.num_bond_rings(b2) > 1 {
            if mol.num_bond_rings(b1) == 1 && mol.num_bond_rings(b3) == 1 {
                for br in &mol.bond_rings {
                    if br.contains(&b1) {
                        if br.contains(&b3) {
                            prefer_cis = true;
                        }
                        break;
                    }
                }
            }
        } else {
            prefer_cis = true;
        }
    } else if stype == STEREO_Z {
        prefer_cis = true;
    } else if stype == STEREO_E {
        prefer_trans = true;
    }
    let ty = if prefer_cis {
        ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
        Path14Type::Cis
    } else if prefer_trans {
        ComputedData::add_path(&mut acc.trans_paths, cx.nb, b1, b2, b3);
        Path14Type::Trans
    } else {
        Path14Type::Other
    };
    acc.paths14.push(Path14 {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        ty,
    });
    let (dl, du);
    if prefer_cis {
        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        du = dl + 2.0 * GEN_DIST_TOL;
    } else if prefer_trans {
        dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        du = dl + 2.0 * GEN_DIST_TOL;
    } else {
        let mut l = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut u = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if u < l {
            std::mem::swap(&mut u, &mut l);
        }
        if (u - l).abs() < DIST12_DELTA {
            l -= GEN_DIST_TOL;
            u += GEN_DIST_TOL;
        }
        dl = l;
        du = u;
    }
    acc.visited14[p] = true;
    check_and_set_bounds(aid1, aid4, dl, du, m)
}

fn set_two_in_same_ring14(
    cx: &Ctx,
    b1: usize,
    b2: usize,
    b3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    macrocycle: bool,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let (aid1, atm2, atm3, aid4) = path_atoms(cx, acc, b1, b2, b3);
    let p = pid(aid1, aid4, cx.na);
    if acc.visited_bound(p, DistType::D13) {
        return Ok(());
    }
    if cx.d(aid1, aid4) < 2.9 {
        return Ok(());
    }
    if mol.bond_between(aid1, atm3).is_some() || mol.bond_between(aid4, atm2).is_some() {
        return Ok(());
    }
    let (bl1, bl2, bl3, ba12, ba23) = path_geom(acc, b1, b2, b3)?;
    let (mut dl, mut du);
    let ty;
    let amide = |bb1: usize, bb3: usize, a1: usize, a2: usize, a3: usize, a4: usize| {
        let n = |a: usize| mol.atoms[a].atomic_num;
        n(a1) != 1
            && n(a3) == 6
            && mol.bonds[bb3].bond_type == DOUBLE
            && (n(a4) == 8 || n(a4) == 7)
            && mol.bonds[bb1].bond_type == SINGLE
            && (n(a2) == 8 || n(a2) == 7)
    };
    if macrocycle {
        if amide(b1, b3, aid1, atm2, atm3, aid4) || amide(b3, b1, aid4, atm3, atm2, aid1) {
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
            ty = Path14Type::Cis;
            du = dl;
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        } else {
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
            if du < dl {
                std::mem::swap(&mut du, &mut dl);
            }
            if (du - dl).abs() < DIST12_DELTA {
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            }
            ty = Path14Type::Other;
        }
    } else if mol.atoms[atm2].hybridization == HYB_SP2 && mol.atoms[atm3].hybridization == HYB_SP2 {
        dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        du = dl;
        dl -= GEN_DIST_TOL;
        du += GEN_DIST_TOL;
        ty = Path14Type::Trans;
        ComputedData::add_path(&mut acc.trans_paths, cx.nb, b1, b2, b3);
    } else {
        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut du, &mut dl);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        ty = Path14Type::Other;
    }
    check_and_set_bounds(aid1, aid4, dl, du, m)?;
    acc.paths14.push(Path14 {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        ty,
    });
    acc.visited14[p] = true;
    Ok(())
}

fn check_h2_nx3h1_ox2(mol: &RdkitMolView, a: usize) -> bool {
    let at = &mol.atoms[a];
    (at.atomic_num == 6 && at.total_num_hs_with_neighbors == 2)
        || (at.atomic_num == 8 && at.total_num_hs_with_neighbors == 0)
        || (at.atomic_num == 7 && mol.degree(a) == 3 && at.total_num_hs_with_neighbors == 1)
}

#[allow(dead_code)]
fn check_nh_ch_ch_nh(mol: &RdkitMolView, a1: usize, a2: usize, a3: usize, a4: usize) -> bool {
    mol.atoms[a1].atomic_num != 1
        && mol.atoms[a4].atomic_num != 1
        && check_h2_nx3h1_ox2(mol, a2)
        && check_h2_nx3h1_ox2(mol, a3)
}

fn check_amide_ester14(
    mol: &RdkitMolView,
    bnd1: usize,
    bnd3: usize,
    a2: usize,
    a3: usize,
    a4: usize,
) -> bool {
    let n = |a: usize| mol.atoms[a].atomic_num;
    n(a3) == 6
        && mol.bonds[bnd3].bond_type == DOUBLE
        && (n(a4) == 8 || n(a4) == 7)
        && mol.bonds[bnd1].bond_type == SINGLE
        && (n(a2) == 8 || (n(a2) == 7 && mol.atoms[a2].total_num_hs_with_neighbors == 1))
}

fn check_macrocycle_all_in_same_ring_amide_ester14(
    mol: &RdkitMolView,
    a1: usize,
    a2: usize,
    a3: usize,
    a4: usize,
) -> bool {
    let n = |a: usize| mol.atoms[a].atomic_num;
    if n(a3) != 6 {
        return false;
    }
    if (n(a2) == 7 || n(a2) == 8) && mol.degree(a2) == 3 && mol.degree(a3) == 3 {
        for nbr in mol.neighbors(a2) {
            if nbr != a1 && nbr != a3 {
                let b = mol.bond_between(a2, nbr).expect("bond");
                if (n(nbr) != 6 && n(nbr) != 1) || mol.bonds[b].bond_type != SINGLE {
                    return false;
                }
                break;
            }
        }
        for nbr in mol.neighbors(a3) {
            if nbr != a2 && nbr != a4 {
                let b = mol.bond_between(a3, nbr).expect("bond");
                if n(nbr) != 8 || mol.bonds[b].bond_type != DOUBLE {
                    return false;
                }
                break;
            }
        }
        return true;
    }
    false
}

fn is_carbonyl(mol: &RdkitMolView, a: usize) -> bool {
    if mol.atoms[a].atomic_num == 6 && mol.degree(a) > 2 {
        for nbr in mol.neighbors(a) {
            let n = mol.atoms[nbr].atomic_num;
            if (n == 8 || n == 7)
                && mol.bonds[mol.bond_between(a, nbr).expect("bond")].bond_type == DOUBLE
            {
                return true;
            }
        }
    }
    false
}

fn check_amide_ester15(mol: &RdkitMolView, bnd1: usize, bnd3: usize, a2: usize, a3: usize) -> bool {
    let n2 = mol.atoms[a2].atomic_num;
    (n2 == 8 || (n2 == 7 && mol.atoms[a2].total_num_hs_with_neighbors == 1))
        && mol.bonds[bnd1].bond_type == SINGLE
        && mol.atoms[a3].atomic_num == 6
        && mol.bonds[bnd3].bond_type == SINGLE
        && is_carbonyl(mol, a3)
}

#[allow(unused_assignments)]
/// `_setChain14Bounds` (`macrocycle == false`) and
/// `_setMacrocycleAllInSameRing14Bounds` (`macrocycle == true`).
fn set_chain14(
    cx: &Ctx,
    b1: usize,
    b2: usize,
    b3: usize,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    force_trans_amides: bool,
    macrocycle: bool,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let (aid1, atm2, atm3, aid4) = path_atoms(cx, acc, b1, b2, b3);
    let p = pid(aid1, aid4, cx.na);
    if acc.visited_bound(p, DistType::D13) {
        return Ok(());
    }
    let (bl1, bl2, bl3, ba12, ba23) = path_geom(acc, b1, b2, b3)?;
    let mut set_the_bound = true;
    let (mut dl, mut du) = (0.0, 0.0);
    let ty;
    let nb = cx.nb;
    let bt = |b: usize| mol.bonds[b].bond_type;
    let an = |a: usize| mol.atoms[a].atomic_num;
    let sec_amide_h = |h: usize, n: usize| {
        an(h) == 1
            && an(n) == 7
            && mol.degree(n) == 3
            && mol.atoms[n].total_num_hs_with_neighbors == 1
    };
    match bt(b2) {
        DOUBLE => {
            if bt(b1) == DOUBLE || bt(b3) == DOUBLE {
                dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                du = dl + 2.0 * GEN_DIST_TOL;
                ty = Path14Type::Cis;
                ComputedData::add_path(&mut acc.cis_paths, nb, b1, b2, b3);
            } else if mol.bonds[b2].stereo > STEREO_ANY {
                let st = atom_stereo(mol, b2, aid1, aid4);
                if st == STEREO_Z {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    du = dl + 2.0 * GEN_DIST_TOL;
                    ty = Path14Type::Cis;
                    ComputedData::add_path(&mut acc.cis_paths, nb, b1, b2, b3);
                } else {
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    dl = du;
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                    ty = Path14Type::Trans;
                    ComputedData::add_path(&mut acc.trans_paths, nb, b1, b2, b3);
                }
            } else {
                dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                if (du - dl).abs() < DIST12_DELTA {
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                }
                ty = Path14Type::Other;
            }
        }
        SINGLE => {
            if an(atm2) == 16 && an(atm3) == 16 && mol.degree(atm2) == 2 && mol.degree(atm3) == 2 {
                dl = compute14_dist_3d(bl1, bl2, bl3, ba12, ba23, M_PI / 2.0) - GEN_DIST_TOL;
                du = dl + 2.0 * GEN_DIST_TOL;
                ty = Path14Type::Other;
            } else if macrocycle
                && (check_macrocycle_all_in_same_ring_amide_ester14(mol, aid1, atm2, atm3, aid4)
                    || check_macrocycle_all_in_same_ring_amide_ester14(mol, aid4, atm3, atm2, aid1))
            {
                dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23) + 0.1;
                ty = Path14Type::Trans;
                ComputedData::add_path(&mut acc.trans_paths, nb, b1, b2, b3);
                du = dl;
                dl -= GEN_DIST_TOL;
                du += GEN_DIST_TOL;
            } else if !macrocycle
                && (check_amide_ester14(mol, b1, b3, atm2, atm3, aid4)
                    || check_amide_ester14(mol, b3, b1, atm3, atm2, aid1))
            {
                if force_trans_amides {
                    if sec_amide_h(aid1, atm2) || sec_amide_h(aid4, atm3) {
                        dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        ty = Path14Type::Trans;
                        ComputedData::add_path(&mut acc.trans_paths, nb, b1, b2, b3);
                    } else {
                        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        ty = Path14Type::Cis;
                        ComputedData::add_path(&mut acc.cis_paths, nb, b1, b2, b3);
                    }
                    du = dl;
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                } else {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    ty = Path14Type::Other;
                }
            } else if check_amide_ester15(mol, b1, b3, atm2, atm3)
                || check_amide_ester15(mol, b3, b1, atm3, atm2)
            {
                if macrocycle {
                    if an(atm2) == 7
                        && mol.degree(atm2) == 3
                        && an(aid1) == 1
                        && mol.atoms[atm2].total_num_hs_with_neighbors == 1
                    {
                        set_the_bound = false;
                        ty = Path14Type::Other;
                    } else {
                        dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        ty = Path14Type::Trans;
                        ComputedData::add_path(&mut acc.trans_paths, nb, b1, b2, b3);
                    }
                    du = dl;
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                } else if force_trans_amides {
                    if sec_amide_h(aid1, atm2) || sec_amide_h(aid4, atm3) {
                        dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        ty = Path14Type::Cis;
                        ComputedData::add_path(&mut acc.cis_paths, nb, b1, b2, b3);
                    } else {
                        dl = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        ty = Path14Type::Trans;
                        ComputedData::add_path(&mut acc.trans_paths, nb, b1, b2, b3);
                    }
                    du = dl;
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                } else {
                    dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                    du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    ty = Path14Type::Other;
                }
            } else {
                dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                ty = Path14Type::Other;
            }
        }
        _ => {
            dl = compute14_dist_cis(bl1, bl2, bl3, ba12, ba23);
            du = compute14_dist_trans(bl1, bl2, bl3, ba12, ba23);
            ty = Path14Type::Other;
        }
    }
    if set_the_bound {
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        check_and_set_bounds(aid1, aid4, dl, du, m)?;
        acc.paths14.push(Path14 {
            bid1: b1,
            bid2: b2,
            bid3: b3,
            ty,
        });
        acc.visited14[p] = true;
    }
    Ok(())
}

/// `_record14Path`.
fn record14_path(cx: &Ctx, b1: usize, b2: usize, b3: usize, acc: &mut ComputedData) {
    let atm2 = acc.bond_adj.get(b1, b2) as usize;
    let atm3 = acc.bond_adj.get(b2, b3) as usize;
    let ty = if cx.mol.atoms[atm2].hybridization == HYB_SP2
        && cx.mol.atoms[atm3].hybridization == HYB_SP2
    {
        ComputedData::add_path(&mut acc.cis_paths, cx.nb, b1, b2, b3);
        Path14Type::Cis
    } else {
        Path14Type::Other
    };
    acc.paths14.push(Path14 {
        bid1: b1,
        bid2: b2,
        bid3: b3,
        ty,
    });
}

fn set14_bounds(
    cx: &Ctx,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
    macro14: bool,
    force_trans_amides: bool,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let nb = cx.nb;
    let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
    let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut done_paths: HashSet<u64> = HashSet::new();
    for bring in &mol.bond_rings {
        let r = bring.len();
        if r < 3 {
            continue;
        }
        let mut bid1 = bring[r - 1];
        for i in 0..r {
            let bid2 = bring[i];
            let bid3 = bring[(i + 1) % r];
            let (b1, b2, b3) = (bid1 as u64, bid2 as u64, bid3 as u64);
            ring_bond_pairs.insert(b1 * nb + b2);
            ring_bond_pairs.insert(b2 * nb + b1);
            done_paths.insert(b1 * nb * nb + b2 * nb + b3);
            done_paths.insert(b3 * nb * nb + b2 * nb + b1);
            if r > 5 {
                if macro14 && r >= MIN_MACROCYCLE_RING_SIZE {
                    set_chain14(cx, bid1, bid2, bid3, acc, m, force_trans_amides, true)?;
                    bid_is_macrocycle.insert(bid2);
                } else {
                    set_in_ring14(cx, bid1, bid2, bid3, acc, m, r)?;
                }
            } else {
                record14_path(cx, bid1, bid2, bid3, acc);
            }
            bid1 = bid2;
        }
    }
    for bid2 in 0..mol.num_bonds() {
        let aid2 = mol.bonds[bid2].begin;
        let aid3 = mol.bonds[bid2].end;
        for &bid1 in &mol.atom_bonds[aid2] {
            if bid1 == bid2 {
                continue;
            }
            for &bid3 in &mol.atom_bonds[aid3] {
                if bid3 == bid2 {
                    continue;
                }
                let (b1, b2, b3) = (bid1 as u64, bid2 as u64, bid3 as u64);
                let id1 = b1 * nb * nb + b2 * nb + b3;
                let id2 = b3 * nb * nb + b2 * nb + b1;
                if done_paths.contains(&id1) || done_paths.contains(&id2) {
                    continue;
                }
                let pairs = [b1 * nb + b2, b2 * nb + b1, b2 * nb + b3, b3 * nb + b2];
                if pairs.iter().any(|p| ring_bond_pairs.contains(p)) {
                    if macro14 && bid_is_macrocycle.contains(&bid2) {
                        set_two_in_same_ring14(cx, bid1, bid2, bid3, acc, m, true)?;
                    } else {
                        set_two_in_same_ring14(cx, bid1, bid2, bid3, acc, m, false)?;
                    }
                } else if (mol.num_bond_rings(bid1) > 0 && mol.num_bond_rings(bid2) > 0)
                    || (mol.num_bond_rings(bid2) > 0 && mol.num_bond_rings(bid3) > 0)
                {
                    set_in_ring14(cx, bid1, bid2, bid3, acc, m, 0)?;
                } else if mol.num_bond_rings(bid2) > 0 {
                    set_in_ring14(cx, bid1, bid2, bid3, acc, m, 0)?;
                } else {
                    set_chain14(cx, bid1, bid2, bid3, acc, m, force_trans_amides, false)?;
                }
            }
        }
    }
    Ok(())
}

fn c15_cis_cis(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() - d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 + a23).cos()) / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 - cval.acos())
}

fn c15_cis_trans(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() - d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 + a23).cos()) / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 + cval.acos())
}

fn c15_trans_trans(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() + d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 - a23).cos()) / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 + cval.acos())
}

fn c15_trans_cis(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = d2 - d3 * a23.cos() - d1 * a12.cos();
    let dy14 = d3 * a23.sin() + d1 * a12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * a23.cos() + d1 * (a12 - a23).cos()) / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 - cval.acos())
}

#[allow(unused_assignments)]
fn set15_helper(
    cx: &Ctx,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    ty: Path14Type,
    acc: &mut ComputedData,
    m: &mut BoundsMatrix,
) -> Result<(), BoundsInvariant> {
    let mol = cx.mol;
    let na = cx.na;
    let nb = mol.num_bonds();
    let aid2 = acc.bond_adj.get(bid1, bid2) as usize;
    let aid1 = mol.other_atom(bid1, aid2);
    let aid3 = acc.bond_adj.get(bid2, bid3) as usize;
    let aid4 = mol.other_atom(bid3, aid3);
    let d1 = acc.bond_lengths[bid1];
    let d2 = acc.bond_lengths[bid2];
    let d3 = acc.bond_lengths[bid3];
    let ang12 = acc.bond_angles.get(bid1, bid2);
    let ang23 = acc.bond_angles.get(bid2, bid3);
    for i in 0..nb {
        let mut du = -1.0;
        let mut dl = 0.0;
        if acc.bond_adj.get(bid3, i) == aid4 as i64 {
            let aid5 = mol.other_atom(i, aid4);
            let p = pid(aid1, aid5, na);
            if acc.visited_bound(p, DistType::D14) {
                return Ok(());
            }
            if cx.d(aid1, aid5) < 3.9 {
                continue;
            }
            if aid1 != aid5
                && (m.lower(aid1, aid5) < DIST12_DELTA
                    || acc.set15_atoms[aid1 * na + aid5]
                    || acc.set15_atoms[aid5 * na + aid1])
            {
                let d4 = acc.bond_lengths[i];
                let ang34 = acc.bond_angles.get(bid3, i);
                let path_id =
                    bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + i as u64;
                let cis = acc.cis_paths.contains(&path_id);
                let trans = acc.trans_paths.contains(&path_id);
                match ty {
                    Path14Type::Cis => {
                        if cis {
                            dl = c15_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                            du = dl + DIST15_TOL;
                            dl -= DIST15_TOL;
                        } else if trans {
                            dl = c15_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                            du = dl + DIST15_TOL;
                            dl -= DIST15_TOL;
                        } else {
                            dl = c15_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                            du = c15_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                        }
                    }
                    Path14Type::Trans => {
                        if cis {
                            dl = c15_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                            du = dl + DIST15_TOL;
                            dl -= DIST15_TOL;
                        } else if trans {
                            dl = c15_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                            du = dl + DIST15_TOL;
                            dl -= DIST15_TOL;
                        } else {
                            dl = c15_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL;
                            du = c15_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL;
                        }
                    }
                    Path14Type::Other => {
                        if cis {
                            dl = c15_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                            du = c15_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                        } else if trans {
                            dl = c15_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
                            du = c15_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;
                        } else {
                            let vw1 = mol.atoms[aid1].rvdw;
                            let vw5 = mol.atoms[aid5].rvdw;
                            dl = VDW_SCALE_15 * (vw1 + vw5);
                        }
                    }
                }
                if du < 0.0 {
                    du = MAX_UPPER;
                }
                check_and_set_bounds(aid1, aid5, dl, du, m)?;
                acc.set15_atoms[aid1 * na + aid5] = true;
                acc.set15_atoms[aid5 * na + aid1] = true;
            }
        }
    }
    Ok(())
}

fn set_lower_bound_vdw(mol: &RdkitMolView, m: &mut BoundsMatrix, dmat: &[f64]) {
    let npt = mol.num_atoms();
    let mut h_donor = vec![false; npt];
    let mut acceptor = vec![false; npt];
    let is_acc = |a: usize| matches!(mol.atoms[a].atomic_num, 7 | 8);
    let is_h_donor = |a: usize| {
        mol.atoms[a].atomic_num == 1
            && mol
                .neighbors(a)
                .any(|n| matches!(mol.atoms[n].atomic_num, 7 | 8))
    };
    for i in 1..npt {
        let vw1 = mol.atoms[i].rvdw;
        if is_h_donor(i) {
            h_donor[i] = true;
        }
        if is_acc(i) {
            acceptor[i] = true;
        }
        for j in 0..i {
            let vw2 = mol.atoms[j].rvdw;
            if m.lower(i, j) < DIST12_DELTA {
                if (h_donor[i] && acceptor[j]) || (acceptor[i] && h_donor[j]) {
                    m.set_lower(i, j, H_BOND_LENGTH);
                } else if dmat[i * npt + j] == 4.0 {
                    m.set_lower(i, j, VDW_SCALE_15 * (vw1 + vw2));
                } else if dmat[i * npt + j] == 5.0 {
                    m.set_lower(
                        i,
                        j,
                        (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2),
                    );
                } else {
                    m.set_lower(i, j, vw1 + vw2);
                }
            }
        }
    }
}

/// `MolOps::getDistanceMat` (topological distances; 1e8 between fragments).
fn distance_mat(mol: &RdkitMolView) -> Vec<f64> {
    let n = mol.num_atoms();
    let mut d = vec![1e8; n * n];
    for s in 0..n {
        d[s * n + s] = 0.0;
        let mut q = std::collections::VecDeque::from([s]);
        while let Some(a) = q.pop_front() {
            for nb in mol.neighbors(a) {
                if d[s * n + nb] == 1e8 {
                    d[s * n + nb] = d[s * n + a] + 1.0;
                    q.push_back(nb);
                }
            }
        }
    }
    d
}

/// Output of [`set_topol_bounds`] besides the matrix: `collectBondsAndAngles`.
#[derive(Clone, Debug, Default)]
pub struct BondsAndAngles {
    pub bonds: Vec<(usize, usize)>,
    /// (aid1, aid2 (centre), aid3, flag for linear triple/cumulene)
    pub angles: Vec<[usize; 4]>,
}

/// `collectBondsAndAngles`.
pub fn collect_bonds_and_angles(mol: &RdkitMolView) -> BondsAndAngles {
    let mut out = BondsAndAngles::default();
    let nb = mol.num_bonds();
    for i in 0..nb {
        let bi = &mol.bonds[i];
        out.bonds.push((bi.begin, bi.end));
        for j in i + 1..nb {
            let bj = &mol.bonds[j];
            let (a11, a12, a21, a22) = (bi.begin, bi.end, bj.begin, bj.end);
            if a11 != a21 && a11 != a22 && a12 != a21 && a12 != a22 {
                continue;
            }
            let mut t = [0usize; 4];
            if a12 == a21 {
                t = [a11, a12, a22, 0];
            } else if a12 == a22 {
                t = [a11, a12, a21, 0];
            } else if a11 == a21 {
                t = [a12, a11, a22, 0];
            } else if a11 == a22 {
                t = [a12, a11, a21, 0];
            }
            if bi.bond_type == TRIPLE || bj.bond_type == TRIPLE {
                t[3] = 1;
            } else if bi.bond_type == DOUBLE && bj.bond_type == DOUBLE && mol.degree(t[1]) == 2 {
                t[3] = 1;
            }
            out.angles.push(t);
        }
    }
    out
}

/// `DGeomHelpers::setTopolBounds(mol, mmat, set15bounds, scaleVDW=false,
/// useMacrocycle14config, forceTransAmides)` on an initialized matrix.
pub fn set_topol_bounds(
    mol: &RdkitMolView,
    uff_labels: &[String],
    m: &mut BoundsMatrix,
    set15: bool,
    use_macrocycle14config: bool,
    force_trans_amides: bool,
) -> Result<(), BoundsInvariant> {
    let na = mol.num_atoms();
    let nb = mol.num_bonds();
    let mut acc = ComputedData::new(na, nb);
    let dmat = distance_mat(mol);
    set12_bounds(mol, uff_labels, m, &mut acc);
    set13_bounds(mol, m, &mut acc)?;
    let cx = Ctx {
        mol,
        dmat: &dmat,
        na,
        nb: nb as u64,
    };
    set14_bounds(&cx, &mut acc, m, use_macrocycle14config, force_trans_amides)?;
    if set15 {
        let paths: Vec<(usize, usize, usize, Path14Type)> = acc
            .paths14
            .iter()
            .map(|p| (p.bid1, p.bid2, p.bid3, p.ty))
            .collect();
        for (b1, b2, b3, ty) in paths {
            set15_helper(&cx, b1, b2, b3, ty, &mut acc, m)?;
            set15_helper(&cx, b3, b2, b1, ty, &mut acc, m)?;
        }
    }
    set_lower_bound_vdw(mol, m, &dmat);
    Ok(())
}
