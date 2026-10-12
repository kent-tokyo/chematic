//! `DGeomHelpers::setTopolBounds` (RDKit 2026.09.1 `BoundsMatrixBuilder.cpp`)
//! on an [`RdkitMolView`], plus `collectBondsAndAngles`.

use std::collections::{HashMap, HashSet};
use std::f64::consts::PI as M_PI;

use chematic_smiles::RdkitMolView;

use super::dg::BoundsMatrix;
use super::stdsort::std_sort;

const DIST12_DELTA: f64 = 0.01;
// RDKit 2026.09.1 derives 1-3 distance tolerances from an angular range
// instead of applying a fixed distance-space delta. This keeps the bounds
// physically consistent when either adjacent bond is itself flexible.
const ANGLE_DELTA: f64 = 0.035;
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
    forced: bool,
}

/// A 1-4 amide/ester configuration that RDKit 2026.09.1 also feeds into
/// ETKDG's basic-knowledge torsion stage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ForcedTorsion {
    pub(crate) atoms: [usize; 4],
    pub(crate) cis: bool,
}

#[derive(Clone, Copy)]
struct Collected14Bounds {
    lower: f64,
    upper: f64,
    aid1: usize,
    aid4: usize,
}

fn merge_14_bounds(mut bounds: Vec<Collected14Bounds>) -> Collected14Bounds {
    debug_assert!(!bounds.is_empty());
    bounds.sort_by(|left, right| left.lower.total_cmp(&right.lower));
    let mut current = bounds[0];
    let mut component_upper = current.upper;
    let mut result_lower = None;
    for bound in bounds.iter().copied().skip(1) {
        if bound.lower <= current.upper {
            current.lower = bound.lower;
            current.upper = current.upper.min(bound.upper);
        } else {
            result_lower.get_or_insert(current.lower);
            current.lower = bound.lower;
            current.upper = if bound.lower <= component_upper {
                component_upper.min(bound.upper)
            } else {
                bound.upper
            };
        }
        component_upper = component_upper.max(bound.upper);
    }
    Collected14Bounds {
        lower: result_lower.unwrap_or(current.lower),
        upper: current.upper,
        aid1: bounds[0].aid1,
        aid4: bounds[0].aid4,
    }
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
    // Match the contraction produced for RDGeom::compute13Dist by the
    // compilers used for the pinned RDKit wheels. Keeping the subtraction as
    // one fused operation minimizes last-bit drift between the two kernels.
    let res = (-2.0 * d1 * d2).mul_add(angle.cos(), d1 * d1 + d2 * d2);
    res.sqrt()
}

fn compute14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = (-d1).mul_add(ang12.cos(), (-d3).mul_add(ang23.cos(), d2));
    let dy = (-d1).mul_add(ang12.sin(), d3 * ang23.sin());
    dx.mul_add(dx, dy * dy).sqrt()
}

fn compute14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = (-d1).mul_add(ang12.cos(), (-d3).mul_add(ang23.cos(), d2));
    let dy = d1.mul_add(ang12.sin(), d3 * ang23.sin());
    dx.mul_add(dx, dy * dy).sqrt()
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
    d[0].mul_add(d[0], d[1].mul_add(d[1], d[2] * d[2])).sqrt()
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
            // RDKit first marks every endpoint of a conjugated bond to a
            // larger heteroatom in a five-membered ring, then widens only a
            // bond whose *two* endpoints are marked. Using `||` here widens
            // the next carbon-carbon bond as well and lets that error spread
            // through triangle smoothing (notably for thiophene).
            let extra = if squish[beg] && squish[end] { 0.2 } else { 0.0 };
            acc.bond_lengths[bi] = bl;
            m.set_upper(beg, end, bl + extra + DIST12_DELTA);
            m.set_lower(beg, end, bl - extra - DIST12_DELTA);
        } else {
            // RDKit 2026.09 changed its UFF fallback for ordinary covalent
            // bonds from a broad van-der-Waals interval to covalent radii
            // with Pauling's bond-order correction. Keep the broad interval
            // only for bond types outside RDKit's ordinary single..aromatic
            // range.
            let unusual_bond_type = b.bond_type < 1 || b.bond_type > 12;
            let (bl, upper_scale, lower_scale) = if unusual_bond_type {
                let vw1 = mol.atoms[beg].rvdw;
                let vw2 = mol.atoms[end].rvdw;
                ((vw1 + vw2) / 2.0, 1.5, 0.75)
            } else {
                let rcov1 =
                    chematic_smiles::rdkit::rdkit_covalent_radius(mol.atoms[beg].atomic_num);
                let rcov2 =
                    chematic_smiles::rdkit::rdkit_covalent_radius(mol.atoms[end].atomic_num);
                let radii_sum = rcov1 + rcov2;
                let bl = (-0.1332 * radii_sum).mul_add(order.ln(), radii_sum);
                (bl, 1.1, 0.9)
            };
            acc.bond_lengths[bi] = bl;
            m.set_upper(beg, end, upper_scale * bl);
            m.set_lower(beg, end, lower_scale * bl);
        }
        acc.visited12[pid(beg, end, na)] = true;
    }
}

fn is_larger_sp2_atom(mol: &RdkitMolView, a: usize) -> bool {
    mol.atoms[a].atomic_num > 13
        && mol.atoms[a].hybridization == HYB_SP2
        && mol.num_atom_rings(a) > 0
}

fn scaled_angle_tolerance(mol: &RdkitMolView, atoms: [usize; 3]) -> f64 {
    atoms.into_iter().fold(ANGLE_DELTA, |tolerance, atom| {
        if is_larger_sp2_atom(mol, atom) {
            tolerance * 2.0
        } else {
            tolerance
        }
    })
}

fn angle_range(angle: f64, tolerance: f64) -> (f64, f64) {
    if angle + tolerance >= M_PI {
        (M_PI - tolerance, M_PI)
    } else if angle - tolerance <= 0.0 {
        (0.0, tolerance)
    } else {
        (angle - tolerance, angle + tolerance)
    }
}

fn set13_helper(
    aid1: usize,
    aid: usize,
    aid3: usize,
    angle: f64,
    m: &mut BoundsMatrix,
    mol: &RdkitMolView,
) -> Result<(), BoundsInvariant> {
    let tolerance = scaled_angle_tolerance(mol, [aid1, aid, aid3]);
    let (angle_lower, angle_upper) = angle_range(angle, tolerance);
    let du = compute13_dist(m.upper(aid1, aid), m.upper(aid, aid3), angle_upper);
    let dl = compute13_dist(m.lower(aid1, aid), m.lower(aid, aid3), angle_lower);
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

fn set_ring_13_bounds(
    mol: &RdkitMolView,
    ring: &[usize],
    m: &mut BoundsMatrix,
    acc: &mut ComputedData,
    visited: &mut [i64],
    angle_taken: &mut [f64],
    done_paths: &mut [bool],
) -> Result<(), BoundsInvariant> {
    let bond_count = mol.num_bonds();
    let mut aid1 = ring[ring.len() - 1];
    for position in 0..ring.len() {
        let aid2 = ring[position];
        let aid3 = ring[(position + 1) % ring.len()];
        let bid1 = mol
            .bond_between(aid1, aid2)
            .ok_or(BoundsInvariant("no bond found"))?;
        let bid2 = mol
            .bond_between(aid2, aid3)
            .ok_or(BoundsInvariant("no bond found"))?;
        let paths = [bond_count * bid1 + bid2, bond_count * bid2 + bid1];
        if !done_paths[paths[0]] && !done_paths[paths[1]] {
            let angle = ring_angle(mol.atoms[aid2].hybridization, ring.len());
            let pair = pid(aid1, aid3, mol.num_atoms());
            if !acc.visited_bound(pair, DistType::D12) {
                set13_helper(aid1, aid2, aid3, angle, m, mol)?;
                acc.visited13[pair] = true;
            }
            acc.bond_angles.set(bid1, bid2, angle);
            acc.bond_adj.set(bid1, bid2, aid2 as i64);
            visited[aid2] += 1;
            angle_taken[aid2] += angle;
            done_paths[paths[0]] = true;
            done_paths[paths[1]] = true;
        }
        aid1 = aid2;
    }
    Ok(())
}

fn remaining_partial_angle(
    mol: &RdkitMolView,
    center: usize,
    n13: i64,
    visited: i64,
    angle_taken: f64,
) -> f64 {
    let hybridization = mol.atoms[center].hybridization;
    if hybridization == HYB_SP2 {
        (2.0 * M_PI - angle_taken) / ((n13 - visited) as u32 as f64)
    } else if hybridization == HYB_SP3 {
        if mol.is_atom_in_ring_of_size(center, 3) {
            116.0 * M_PI / 180.0
        } else if mol.is_atom_in_ring_of_size(center, 4) {
            112.0 * M_PI / 180.0
        } else {
            109.5 * M_PI / 180.0
        }
    } else if mol.degree(center) == 5 {
        105.0 * M_PI / 180.0
    } else if mol.degree(center) == 6 {
        135.0 * M_PI / 180.0
    } else {
        120.0 * M_PI / 180.0
    }
}

fn unvisited_angle(hybridization: u8) -> f64 {
    match hybridization {
        HYB_SP => M_PI,
        HYB_SP2 => 2.0 * M_PI / 3.0,
        HYB_SP3 => 109.5 * M_PI / 180.0,
        HYB_SP3D => 105.0 * M_PI / 180.0,
        HYB_SP3D2 => 135.0 * M_PI / 180.0,
        _ => 120.0 * M_PI / 180.0,
    }
}

#[allow(clippy::too_many_arguments)]
fn set_atom_13_bounds(
    mol: &RdkitMolView,
    center: usize,
    m: &mut BoundsMatrix,
    acc: &mut ComputedData,
    visited: &mut [i64],
    angle_taken: &mut [f64],
) -> Result<(), BoundsInvariant> {
    let degree = mol.degree(center);
    let n13 = (degree * degree.saturating_sub(1) / 2) as i64;
    if n13 == visited[center] {
        return Ok(());
    }
    let had_ring_angle = visited[center] > 0;
    let bonds = &mol.atom_bonds[center];
    for (k1, &bid1) in bonds.iter().enumerate() {
        let aid1 = mol.other_atom(bid1, center);
        for &bid2 in &bonds[..k1] {
            if had_ring_angle && acc.bond_angles.get(bid1, bid2) >= 0.0 {
                continue;
            }
            let aid3 = mol.other_atom(bid2, center);
            let angle = if had_ring_angle {
                remaining_partial_angle(mol, center, n13, visited[center], angle_taken[center])
            } else {
                unvisited_angle(mol.atoms[center].hybridization)
            };
            let pair = pid(aid1, aid3, mol.num_atoms());
            if !acc.visited_bound(pair, DistType::D12) {
                if had_ring_angle || degree <= 4 {
                    set13_helper(aid1, center, aid3, angle, m, mol)?;
                } else {
                    let dmax = acc.bond_lengths[bid1] + acc.bond_lengths[bid2];
                    check_and_set_bounds(aid1, aid3, 1.0, dmax * 1.2, m)?;
                }
                acc.visited13[pair] = true;
            }
            acc.bond_angles.set(bid1, bid2, angle);
            acc.bond_adj.set(bid1, bid2, center as i64);
            angle_taken[center] += angle;
            visited[center] += 1;
        }
    }
    Ok(())
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
    let mut done_paths = vec![false; mol.num_bonds() * mol.num_bonds()];
    for ring in &atom_rings {
        set_ring_13_bounds(
            mol,
            ring,
            m,
            acc,
            &mut visited,
            &mut angle_taken,
            &mut done_paths,
        )?;
    }
    for center in 0..npt {
        set_atom_13_bounds(mol, center, m, acc, &mut visited, &mut angle_taken)?;
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

fn path_angles(
    acc: &ComputedData,
    b1: usize,
    b2: usize,
    b3: usize,
) -> Result<(f64, f64), BoundsInvariant> {
    let ba12 = acc.bond_angles.get(b1, b2);
    let ba23 = acc.bond_angles.get(b2, b3);
    if !(ba12 > 0.0) || !(ba23 > 0.0) {
        return Err(BoundsInvariant("bad bond angle"));
    }
    Ok((ba12, ba23))
}

mod bounds_14;
use bounds_14::set14_bounds;
fn c15_cis_cis(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = (-d1).mul_add(a12.cos(), (-d3).mul_add(a23.cos(), d2));
    let dy14 = (-d1).mul_add(a12.sin(), d3 * a23.sin());
    let d14 = dx14.mul_add(dx14, dy14 * dy14).sqrt();
    let numerator = d1.mul_add((a12 + a23).cos(), (-d2).mul_add(a23.cos(), d3));
    let cval = (numerator / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 - cval.acos())
}

fn c15_cis_trans(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = (-d1).mul_add(a12.cos(), (-d3).mul_add(a23.cos(), d2));
    let dy14 = (-d1).mul_add(a12.sin(), d3 * a23.sin());
    let d14 = dx14.mul_add(dx14, dy14 * dy14).sqrt();
    let numerator = d1.mul_add((a12 + a23).cos(), (-d2).mul_add(a23.cos(), d3));
    let cval = (numerator / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 + cval.acos())
}

fn c15_trans_trans(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = (-d1).mul_add(a12.cos(), (-d3).mul_add(a23.cos(), d2));
    let dy14 = d1.mul_add(a12.sin(), d3 * a23.sin());
    let d14 = dx14.mul_add(dx14, dy14 * dy14).sqrt();
    let numerator = d1.mul_add((a12 - a23).cos(), (-d2).mul_add(a23.cos(), d3));
    let cval = (numerator / d14).clamp(-1.0, 1.0);
    compute13_dist(d14, d4, a34 + cval.acos())
}

fn c15_trans_cis(d1: f64, d2: f64, d3: f64, d4: f64, a12: f64, a23: f64, a34: f64) -> f64 {
    let dx14 = (-d1).mul_add(a12.cos(), (-d3).mul_add(a23.cos(), d2));
    let dy14 = d1.mul_add(a12.sin(), d3 * a23.sin());
    let d14 = dx14.mul_add(dx14, dy14 * dy14).sqrt();
    let numerator = d1.mul_add((a12 - a23).cos(), (-d2).mul_add(a23.cos(), d3));
    let cval = (numerator / d14).clamp(-1.0, 1.0);
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
pub(crate) fn set_topol_bounds_with_forced_torsions(
    mol: &RdkitMolView,
    uff_labels: &[String],
    m: &mut BoundsMatrix,
    set15: bool,
    use_macrocycle14config: bool,
    force_trans_amides: bool,
) -> Result<Vec<ForcedTorsion>, BoundsInvariant> {
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
    Ok(acc
        .paths14
        .iter()
        .filter_map(|path| {
            if !path.forced {
                return None;
            }
            let (aid1, aid2, aid3, aid4) = path_atoms(&cx, &acc, path.bid1, path.bid2, path.bid3);
            Some(ForcedTorsion {
                atoms: [aid1, aid2, aid3, aid4],
                cis: path.ty == Path14Type::Cis,
            })
        })
        .collect())
}

/// Build RDKit-compatible topology bounds while preserving the historical
/// public return type. The ETKDG embedding pipeline uses the internal variant
/// above to retain RDKit 2026.09.1's forced amide/ester torsion metadata.
pub fn set_topol_bounds(
    mol: &RdkitMolView,
    uff_labels: &[String],
    m: &mut BoundsMatrix,
    set15: bool,
    use_macrocycle14config: bool,
    force_trans_amides: bool,
) -> Result<(), BoundsInvariant> {
    set_topol_bounds_with_forced_torsions(
        mol,
        uff_labels,
        m,
        set15,
        use_macrocycle14config,
        force_trans_amides,
    )
    .map(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn propane_one_three_bounds_match_rdkit_2026_09_1() {
        let molecule =
            chematic_chem::add_hydrogens(&chematic_smiles::parse("CCC").expect("valid propane"));
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit view");
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&molecule);
        let mut bounds = BoundsMatrix::new(view.num_atoms());
        init_bounds_mat(&mut bounds);
        set_topol_bounds(&view, &labels, &mut bounds, true, true, true).expect("bounds build");

        assert!((bounds.lower(0, 2) - 2.425_702_290_291).abs() < 1.0e-12);
        assert!((bounds.upper(0, 2) - 2.519_525_675_977).abs() < 1.0e-12);
    }

    #[test]
    fn butane_one_four_bounds_match_rdkit_2026_09_1() {
        let molecule =
            chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCC").expect("valid butane"));
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit view");
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&molecule);
        let mut bounds = BoundsMatrix::new(view.num_atoms());
        init_bounds_mat(&mut bounds);
        set_topol_bounds(&view, &labels, &mut bounds, true, true, true).expect("bounds build");

        assert!((bounds.lower(0, 3) - 2.408_255_053_938_461).abs() < 1.0e-12);
        assert!((bounds.upper(0, 3) - 3.875_383_240_895_435_4).abs() < 1.0e-12);
    }

    #[test]
    fn thiophene_squish_does_not_widen_adjacent_carbon_bonds() {
        let molecule = chematic_chem::add_hydrogens(
            &chematic_smiles::parse("c1ccsc1").expect("valid thiophene"),
        );
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit view");
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&molecule);
        let mut bounds = BoundsMatrix::new(view.num_atoms());
        init_bounds_mat(&mut bounds);
        set_topol_bounds(&view, &labels, &mut bounds, true, true, true).expect("bounds build");

        // RDKit 2026.09.1 GetMoleculeBoundsMatrix(ETKDGv3 parameters).
        assert!((bounds.lower(0, 4) - 1.369_256_405_400_789).abs() < 1.0e-12);
        assert!((bounds.upper(0, 4) - 1.389_256_405_400_789).abs() < 1.0e-12);
        assert!((bounds.lower(0, 2) - 2.186_996_595_462_591_4).abs() < 1.0e-12);
        assert!((bounds.upper(0, 2) - 2.276_098_883_787_998).abs() < 1.0e-12);
    }

    #[test]
    fn cubane_uses_the_rdkit_2026_09_relevant_ring_family() {
        let molecule = chematic_chem::add_hydrogens(
            &chematic_smiles::parse("C1C2C3C1C4C2C3C4").expect("valid cubane"),
        );
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit view");

        assert_eq!(view.atom_rings.len(), 7);
        assert_eq!(view.atom_rings[0], vec![0, 1, 2, 3]);
        assert_eq!(view.atom_rings[6], vec![2, 3, 4, 7, 6]);
    }

    #[test]
    fn macrolactam_exocyclic_oxygen_bounds_match_rdkit_2026_09_1() {
        let molecule = chematic_chem::add_hydrogens(
            &chematic_smiles::parse("O=C1CCCCCCCCCCN1").expect("valid macrolactam"),
        );
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit view");
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&molecule);
        let mut bounds = BoundsMatrix::new(view.num_atoms());
        init_bounds_mat(&mut bounds);
        set_topol_bounds(&view, &labels, &mut bounds, true, true, true).expect("bounds build");

        assert!((bounds.upper(0, 11) - 3.692_529_513_413).abs() < 1.0e-12);
        assert!((bounds.lower(0, 10) - 2.275).abs() < 1.0e-12);
        assert_eq!(bounds.upper(0, 10), MAX_UPPER);
    }

    #[test]
    fn acetamide_exports_forced_cn_torsions_for_etkdg_basic_knowledge() {
        let molecule = chematic_chem::add_hydrogens(
            &chematic_smiles::parse("CC(=O)NC").expect("valid acetamide"),
        );
        let view = chematic_smiles::rdkit_mol_view(&molecule).expect("RDKit view");
        let labels = chematic_ff::rdkit_uff::rdkit_uff_atom_labels(&molecule);
        let mut bounds = BoundsMatrix::new(view.num_atoms());
        init_bounds_mat(&mut bounds);
        let forced =
            set_topol_bounds_with_forced_torsions(&view, &labels, &mut bounds, true, true, true)
                .expect("bounds build");

        assert!(forced.iter().any(|torsion| {
            let [_, center1, center2, _] = torsion.atoms;
            !torsion.cis
                && [
                    view.atoms[center1].atomic_num,
                    view.atoms[center2].atomic_num,
                ]
                .into_iter()
                .collect::<HashSet<_>>()
                    == HashSet::from([6, 7])
        }));
    }

    #[test]
    fn angle_range_stays_inside_zero_to_pi() {
        assert_eq!(angle_range(0.01, 0.035), (0.0, 0.035));
        assert_eq!(angle_range(M_PI - 0.01, 0.035), (M_PI - 0.035, M_PI));
    }

    fn collected(lower: f64, upper: f64) -> Collected14Bounds {
        Collected14Bounds {
            lower,
            upper,
            aid1: 2,
            aid4: 7,
        }
    }

    #[test]
    fn overlapping_one_four_bounds_are_intersected() {
        let merged = merge_14_bounds(vec![collected(1.0, 4.0), collected(2.0, 3.0)]);
        assert_eq!((merged.lower, merged.upper), (2.0, 3.0));
        assert_eq!((merged.aid1, merged.aid4), (2, 7));
    }

    #[test]
    fn disjoint_one_four_bounds_are_spanned() {
        let merged = merge_14_bounds(vec![collected(1.0, 2.0), collected(3.0, 4.0)]);
        assert_eq!((merged.lower, merged.upper), (1.0, 4.0));
    }

    #[test]
    fn chained_one_four_bounds_keep_maximal_overlap_components() {
        let merged = merge_14_bounds(vec![
            collected(1.0, 5.0),
            collected(2.0, 3.0),
            collected(4.0, 6.0),
        ]);
        assert_eq!((merged.lower, merged.upper), (2.0, 5.0));
    }
}
