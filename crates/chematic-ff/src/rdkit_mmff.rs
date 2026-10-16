//! RDKit's MMFF94 force field and minimizer, operation for operation.
//!
//! [`RdkitMmffField`] holds the terms RDKit 2026.03's
//! `MMFF::constructForceField` builds (bond stretch, angle bend,
//! stretch-bend, out-of-plane, torsion and a combined non-bonded
//! contribution), in RDKit's order and with RDKit's energy and gradient
//! expressions (`Code/ForceField/MMFF/*.cpp`), so energies and gradients
//! reproduce RDKit within platform floating-point tolerance (and are
//! bit-identical on the recorded Linux comparison lane).
//! [`RdkitMmffField::minimize`] is
//! `ForceField::minimize` (`BFGSOpt::minimize` with its line search and
//! RDKit's gradient scaling), so `MMFFOptimizeMolecule` lands on the same
//! coordinates. The MMFF parameters come from chematic's MMFF94 typing.

// RDKit's spelling (`x * -1.0`, `-1. * k`) is kept: those are exact
// negations, so rewriting them would not change any result.
#![allow(clippy::neg_multiply, clippy::type_complexity)]

use chematic_core::{AtomIdx, BondOrder, Molecule};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cell::RefCell;

mod dielectric;
pub use dielectric::MmffDielectricError;

use crate::mmff94_minimizer::{MinimizerError, rdkit_prepared_parts};
use crate::mmff94_numeric::mmff94_charges_rdkit_order;

const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;
const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

pub(crate) fn clip_to_one(x: f64) -> f64 {
    x.clamp(-1.0, 1.0)
}

/// `std::max(a, b)`: `b` only when `a < b` (keeps a NaN `a`).
pub(crate) fn smax(a: f64, b: f64) -> f64 {
    if a < b { b } else { a }
}

pub(crate) fn is_double_zero(x: f64) -> bool {
    x < 1.0e-10 && x > -1.0e-10
}

/// Bonded parameters of chematic's prepared MMFF94 model.
pub(crate) struct PreparedParts {
    /// (i, j, r0, kb)
    pub bonds: Vec<(usize, usize, f64, f64)>,
    /// (central atom, koop), one entry per out-of-plane term
    pub oops: Vec<(usize, f64)>,
    /// (i, j, k, l, v1, v2, v3)
    pub torsions: Vec<(usize, usize, usize, usize, f64, f64, f64)>,
}

// ─── Point3D arithmetic as RDKit's `RDGeom::Point3D` ─────────────────────────

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct P3 {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

impl P3 {
    pub(crate) fn at(pos: &[f64], i: usize) -> P3 {
        P3 {
            x: pos[3 * i],
            y: pos[3 * i + 1],
            z: pos[3 * i + 2],
        }
    }
    pub(crate) fn sub(self, o: P3) -> P3 {
        P3 {
            x: self.x - o.x,
            y: self.y - o.y,
            z: self.z - o.z,
        }
    }
    pub(crate) fn div(self, v: f64) -> P3 {
        P3 {
            x: self.x / v,
            y: self.y / v,
            z: self.z / v,
        }
    }
    pub(crate) fn neg(self) -> P3 {
        P3 {
            x: self.x * -1.0,
            y: self.y * -1.0,
            z: self.z * -1.0,
        }
    }
    pub(crate) fn dot(self, o: P3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    pub(crate) fn length(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }
    pub(crate) fn cross(self, o: P3) -> P3 {
        P3 {
            x: self.y * o.z - self.z * o.y,
            y: -self.x * o.z + self.z * o.x,
            z: self.x * o.y - self.y * o.x,
        }
    }
}

/// `ForceField::distance` on a flat position array.
pub(crate) fn distance(pos: &[f64], i: usize, j: usize) -> f64 {
    let (i, j) = if j < i { (j, i) } else { (i, j) };
    let mut res = 0.0;
    for d in 0..3 {
        let t = pos[3 * i + d] - pos[3 * j + d];
        res += t * t;
    }
    res.sqrt()
}

/// Interatomic distances for the most recently evaluated position in one
/// minimization. RDKit's `ForceField` shares this cache between an energy
/// evaluation and the gradient evaluated immediately afterwards.
struct DistanceCache {
    n: usize,
    point: Vec<f64>,
    values: Vec<f64>,
    generations: Vec<u32>,
    generation: u32,
}

impl DistanceCache {
    fn new(n: usize) -> Self {
        Self {
            n,
            point: Vec::new(),
            values: vec![0.0; n * n],
            generations: vec![0; n * n],
            generation: 1,
        }
    }

    fn prepare(&mut self, pos: &[f64]) {
        let same_point = self.point.len() == pos.len()
            && self
                .point
                .iter()
                .zip(pos)
                .all(|(a, b)| a.to_bits() == b.to_bits());
        if !same_point {
            self.point.clear();
            self.point.extend_from_slice(pos);
            self.generation = self.generation.wrapping_add(1);
            if self.generation == 0 {
                self.generations.fill(0);
                self.generation = 1;
            }
        }
    }

    fn get(&mut self, pos: &[f64], i: usize, j: usize) -> f64 {
        let (i, j) = if j < i { (j, i) } else { (i, j) };
        let index = i * self.n + j;
        if self.generations[index] != self.generation {
            let computed = distance(pos, i, j);
            self.values[index] = computed;
            self.generations[index] = self.generation;
            computed
        } else {
            self.values[index]
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct NonbondedTerm {
    i: usize,
    j: usize,
    /// (R_ij_star, epsilon)
    vdw: Option<(f64, f64)>,
    /// (unscreened charge product, screened chargeTerm, is1_4)
    ele: Option<(f64, f64, bool)>,
}

/// RDKit's MMFF94 force field for one molecule (all atoms, RDKit order).
#[derive(Clone, Debug)]
pub struct RdkitMmffField {
    n: usize,
    bonds: Vec<(usize, usize, f64, f64)>,
    angles: Vec<(usize, usize, usize, f64, f64, bool)>,
    stretch_bends: Vec<(usize, usize, usize, f64, f64, f64, f64, f64)>,
    oops: Vec<(usize, usize, usize, usize, f64)>,
    torsions: Vec<(usize, usize, usize, usize, f64, f64, f64)>,
    nonbonded: Vec<NonbondedTerm>,
    /// Which term classes are included (bond, angle, stbn, oop, torsion,
    /// vdW, electrostatic), as `MMFFMolProperties::SetMMFF*Term`.
    terms: [bool; 7],
    dielectric_constant: f64,
}

/// MMFF parameter set (RDKit `mmffVariant`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MmffVariant {
    /// `"MMFF94"`.
    Mmff94,
    /// `"MMFF94s"`: MMFF94 with the MMFF94s out-of-plane and torsion tables.
    Mmff94s,
}

/// Term classes of [`RdkitMmffField`], in RDKit's contribution order.
pub const RDKIT_MMFF_TERMS: [&str; 7] = [
    "bond",
    "angle",
    "stretch_bend",
    "oop",
    "torsion",
    "vdw",
    "ele",
];

/// Diagnostics from one RDKit-compatible MMFF94 minimization.
///
/// [`RdkitMmffField::optimize`] remains the compact compatibility API. Use
/// [`RdkitMmffField::optimize_detailed`] when a caller needs a deterministic
/// convergence signal without evaluating the gradient again.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RdkitMmffOptimization {
    /// RDKit status: 0 converged, 1 iteration limit, negative on failure.
    pub status: i32,
    /// Energy reported with RDKit's post-minimization distance-cache semantics.
    pub energy: f64,
    /// Number of completed BFGS iterations.
    pub iterations: u32,
    /// RDKit's normalized gradient convergence test at the returned point.
    ///
    /// This is `None` only when RDKit declares convergence from a tiny
    /// position step before evaluating the gradient at that new point.
    pub final_gradient_test: Option<f64>,
}

type BondAdjacency = Vec<Vec<(usize, usize)>>;

fn rdkit_adjacency(mol: &Molecule) -> (Vec<chematic_core::BondIdx>, BondAdjacency) {
    let order = rdkit_bond_order_with_added_hs(mol);
    let mut adjacency = vec![Vec::new(); mol.atom_count()];
    for &bond_index in &order {
        let bond = mol.bond(bond_index);
        let (first, second) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        adjacency[first].push((second, bond_index.0 as usize));
        adjacency[second].push((first, bond_index.0 as usize));
    }
    (order, adjacency)
}

fn mmff_bond_terms(
    mol: &Molecule,
    order: &[chematic_core::BondIdx],
    parts: &PreparedParts,
) -> Vec<(usize, usize, f64, f64)> {
    let parameters: FxHashMap<(usize, usize), (f64, f64)> = parts
        .bonds
        .iter()
        .map(|&(i, j, r0, kb)| ((i.min(j), i.max(j)), (r0, kb)))
        .collect();
    order
        .iter()
        .filter_map(|&bond_index| {
            let bond = mol.bond(bond_index);
            let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            parameters
                .get(&(i.min(j), i.max(j)))
                .map(|&(r0, kb)| (i, j, r0, kb))
        })
        .collect()
}

fn mmff_oop_terms(
    adjacency: &BondAdjacency,
    parts: &PreparedParts,
    types: &[u8],
    mmffs: bool,
) -> Vec<(usize, usize, usize, usize, f64)> {
    let parameters: FxHashMap<usize, f64> = parts.oops.iter().copied().collect();
    let mut terms = Vec::new();
    for center in 0..adjacency.len() {
        if adjacency[center].len() != 3 {
            continue;
        }
        let Some(&base_koop) = parameters.get(&center) else {
            continue;
        };
        let koop = if mmffs {
            let Some(value) = crate::mmff94_energy::mmff94s_oop(
                types[center],
                types[adjacency[center][0].0],
                types[adjacency[center][1].0],
                types[adjacency[center][2].0],
            ) else {
                continue;
            };
            value
        } else {
            base_koop
        };
        let idx = [
            adjacency[center][0].0,
            center,
            adjacency[center][1].0,
            adjacency[center][2].0,
        ];
        for permutation in [[0, 1, 2, 3], [0, 1, 3, 2], [2, 1, 3, 0]] {
            terms.push((
                idx[permutation[0]],
                idx[permutation[1]],
                idx[permutation[2]],
                idx[permutation[3]],
                koop,
            ));
        }
    }
    terms
}

fn mmff_torsion_terms(
    mol: &Molecule,
    view: &Molecule,
    adjacency: &BondAdjacency,
    parts: &PreparedParts,
    types: &[u8],
    mmffs: bool,
) -> Vec<(usize, usize, usize, usize, f64, f64, f64)> {
    let parameters: FxHashMap<(usize, usize, usize, usize), (f64, f64, f64)> = parts
        .torsions
        .iter()
        .flat_map(|&(i, j, k, l, v1, v2, v3)| {
            [((i, j, k, l), (v1, v2, v3)), ((l, k, j, i), (v1, v2, v3))]
        })
        .collect();
    let torsion_end = |atom: usize| {
        adjacency[atom].len() > 1
            && !adjacency[atom].iter().any(|&(_, bond)| {
                mol.bond(chematic_core::BondIdx(bond as u32)).order == BondOrder::Triple
            })
    };
    let sp2_or_sp3 = |atom: usize| {
        matches!(
            chematic_smarts::rdkit_hybridization(view, AtomIdx(atom as u32)),
            Some(2) | Some(3)
        )
    };
    let mut terms = Vec::new();
    let mut seen = FxHashSet::default();
    for j in 0..adjacency.len() {
        if !torsion_end(j) {
            continue;
        }
        for &(k, central_bond) in &adjacency[j] {
            if !torsion_end(k)
                || !seen.insert((j.min(k), j.max(k)))
                || !(sp2_or_sp3(j) && sp2_or_sp3(k))
            {
                continue;
            }
            add_mmff_torsions_for_bond(
                &mut terms,
                view,
                adjacency,
                &parameters,
                types,
                mmffs,
                j,
                k,
                central_bond,
            );
        }
    }
    terms
}

#[allow(clippy::too_many_arguments)]
fn add_mmff_torsions_for_bond(
    terms: &mut Vec<(usize, usize, usize, usize, f64, f64, f64)>,
    view: &Molecule,
    adjacency: &BondAdjacency,
    parameters: &FxHashMap<(usize, usize, usize, usize), (f64, f64, f64)>,
    types: &[u8],
    mmffs: bool,
    j: usize,
    k: usize,
    central_bond: usize,
) {
    for &(i, first_bond) in &adjacency[j] {
        if first_bond == central_bond {
            continue;
        }
        for &(l, second_bond) in &adjacency[k] {
            if second_bond == central_bond || second_bond == first_bond || l == i {
                continue;
            }
            if mmffs {
                if let Some(parameter) = crate::mmff94_minimizer::mmff_torsion_term_params_variant(
                    view, i, j, k, l, types, true,
                ) {
                    terms.push((i, j, k, l, parameter.v1, parameter.v2, parameter.v3));
                }
            } else if let Some(&(v1, v2, v3)) = parameters.get(&(i, j, k, l)) {
                terms.push((i, j, k, l, v1, v2, v3));
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn mmff_nonbonded_terms(
    coords: &[[f64; 3]],
    adjacency: &BondAdjacency,
    types: &[u8],
    charges: &[f64],
    threshold: f64,
    ignore_interfrag: bool,
    terms: [bool; 7],
) -> Vec<NonbondedTerm> {
    let n = adjacency.len();
    let topology = topological_distances(n, adjacency);
    let fragments = fragments(n, adjacency);
    let mut output = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if ignore_interfrag && fragments[i] != fragments[j] {
                continue;
            }
            let relation = topology[i * n + j];
            if matches!(relation, 1 | 2) {
                continue;
            }
            let delta = P3 {
                x: coords[i][0] - coords[j][0],
                y: coords[i][1] - coords[j][1],
                z: coords[i][2] - coords[j][2],
            };
            if delta.length() > threshold {
                continue;
            }
            let vdw = terms[5].then(|| vdw_params(types[i], types[j])).flatten();
            let ele = (terms[6] && !is_double_zero(charges[i]) && !is_double_zero(charges[j]))
                .then_some((
                    charges[i] * charges[j],
                    charges[i] * charges[j],
                    relation == 3,
                ));
            if vdw.is_some() || ele.is_some() {
                output.push(NonbondedTerm { i, j, vdw, ele });
            }
        }
    }
    output
}

impl RdkitMmffField {
    /// `MMFFGetMoleculeForceField(mol, MMFFGetMoleculeProperties(mol),
    /// nonBondedThresh, confId, ignoreInterfragInteractions)` with the
    /// conformer `coords` (all atoms of `mol`, which must carry explicit
    /// hydrogens in RDKit `AddHs` order).
    pub fn new(
        mol: &Molecule,
        coords: &[[f64; 3]],
        non_bonded_thresh: f64,
        ignore_interfrag: bool,
    ) -> Result<Self, MinimizerError> {
        Self::with_terms(
            mol,
            coords,
            MmffVariant::Mmff94,
            non_bonded_thresh,
            ignore_interfrag,
            [true; 7],
        )
    }

    /// [`Self::new`] for the given MMFF variant and only the selected term
    /// classes (see [`RDKIT_MMFF_TERMS`]).
    pub fn with_terms(
        mol: &Molecule,
        coords: &[[f64; 3]],
        variant: MmffVariant,
        non_bonded_thresh: f64,
        ignore_interfrag: bool,
        terms: [bool; 7],
    ) -> Result<Self, MinimizerError> {
        let mmffs = variant == MmffVariant::Mmff94s;
        let n = mol.atom_count();
        let (parts, types, view) = rdkit_prepared_parts(mol)?;
        let (order, adj) = rdkit_adjacency(mol);
        let nbrs: Vec<Vec<usize>> = adj
            .iter()
            .map(|v| v.iter().map(|&(a, _)| a).collect())
            .collect();
        let charges =
            mmff94_charges_rdkit_order(mol, &nbrs).map_err(MinimizerError::ChargeCalculation)?;

        let (angles, stretch_bends) = crate::mmff94_minimizer::rdkit_angle_parts(mol, &nbrs)?;
        Ok(RdkitMmffField {
            n,
            bonds: mmff_bond_terms(mol, &order, &parts),
            angles,
            stretch_bends,
            oops: mmff_oop_terms(&adj, &parts, &types, mmffs),
            torsions: mmff_torsion_terms(mol, &view, &adj, &parts, &types, mmffs),
            nonbonded: mmff_nonbonded_terms(
                coords,
                &adj,
                &types,
                &charges,
                non_bonded_thresh,
                ignore_interfrag,
                terms,
            ),
            terms,
            dielectric_constant: 1.0,
        })
    }

    /// Number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.n
    }

    /// `ForceField::calcEnergy(pos)` on a flat `3 * n` position array.
    pub fn energy(&self, pos: &[f64]) -> f64 {
        self.energy_with_distances(pos, pos)
    }

    /// The energy at `pos` with every interatomic distance RDKit takes from
    /// its `ForceField::distance` cache taken at `dpos` instead (see
    /// [`crate::rdkit_bfgs::BfgsOutcome::stale_distance_point`]).
    pub fn energy_with_distances(&self, pos: &[f64], dpos: &[f64]) -> f64 {
        self.energy_with_distance_fn(pos, |i, j| distance(dpos, i, j))
    }

    fn energy_with_distance_fn(
        &self,
        pos: &[f64],
        mut pair_distance: impl FnMut(usize, usize) -> f64,
    ) -> f64 {
        let mut res = 0.0;
        if self.terms[0] && !self.bonds.is_empty() {
            res += self.bond_energy(&mut pair_distance);
        }
        if self.terms[1] && !self.angles.is_empty() {
            res += self.angle_energy(pos, &mut pair_distance);
        }
        if self.terms[2] && !self.stretch_bends.is_empty() {
            res += self.stbn_energy(pos, &mut pair_distance);
        }
        if self.terms[3] && !self.oops.is_empty() {
            res += self.oop_energy(pos);
        }
        if self.terms[4] && !self.torsions.is_empty() {
            res += self.torsion_energy(pos);
        }
        if (self.terms[5] || self.terms[6]) && !self.nonbonded.is_empty() {
            res += self.nonbonded_energy(&mut pair_distance);
        }
        res
    }

    /// `ForceField::calcGrad(pos, grad)`: adds every contribution's gradient
    /// to `grad` (unscaled).
    pub fn gradient(&self, pos: &[f64], grad: &mut [f64]) {
        self.gradient_with_distance_fn(pos, grad, |i, j| distance(pos, i, j));
    }

    fn gradient_with_distance_fn(
        &self,
        pos: &[f64],
        grad: &mut [f64],
        mut pair_distance: impl FnMut(usize, usize) -> f64,
    ) {
        if self.terms[0] {
            self.bond_grad(pos, grad, &mut pair_distance);
        }
        if self.terms[1] {
            self.angle_grad(pos, grad, &mut pair_distance);
        }
        if self.terms[2] {
            self.stbn_grad(pos, grad, &mut pair_distance);
        }
        if self.terms[3] {
            self.oop_grad(pos, grad);
        }
        if self.terms[4] {
            self.torsion_grad(pos, grad);
        }
        if self.terms[5] || self.terms[6] {
            self.nonbonded_grad(pos, grad, &mut pair_distance);
        }
    }

    // ─── Bond stretch ────────────────────────────────────────────────────────

    fn bond_energy(&self, pair_distance: &mut impl FnMut(usize, usize) -> f64) -> f64 {
        let c1 = MDYNE_A_TO_KCAL_MOL;
        let cs = -2.0;
        let c3 = 7.0 / 12.0;
        let mut sum = 0.0;
        for &(i, j, r0, kb) in &self.bonds {
            let dist_term = pair_distance(i, j) - r0;
            let dist_term2 = dist_term * dist_term;
            sum += 0.5 * c1 * kb * dist_term2 * (1.0 + cs * dist_term + c3 * cs * cs * dist_term2);
        }
        sum
    }

    fn bond_grad(
        &self,
        pos: &[f64],
        grad: &mut [f64],
        pair_distance: &mut impl FnMut(usize, usize) -> f64,
    ) {
        let cs = -2.0;
        let c1 = MDYNE_A_TO_KCAL_MOL;
        let c3 = 7.0 / 12.0;
        for &(i, j, r0, kb) in &self.bonds {
            let dist = pair_distance(i, j);
            let dist_term = dist - r0;
            let de_dr = c1
                * kb
                * dist_term
                * (1.0 + 1.5 * cs * dist_term + 2.0 * c3 * cs * cs * dist_term * dist_term);
            for d in 0..3 {
                let dgrad = if dist > 0.0 {
                    de_dr * (pos[3 * i + d] - pos[3 * j + d]) / dist
                } else {
                    kb * 0.01
                };
                grad[3 * i + d] += dgrad;
                grad[3 * j + d] -= dgrad;
            }
        }
    }

    // ─── Angle bend ──────────────────────────────────────────────────────────

    fn angle_energy(
        &self,
        pos: &[f64],
        pair_distance: &mut impl FnMut(usize, usize) -> f64,
    ) -> f64 {
        let mut res = 0.0;
        for &(i, j, k, theta0, ka, linear) in &self.angles {
            let dist1 = pair_distance(i, j);
            let dist2 = pair_distance(j, k);
            let (p1, p2, p3) = (P3::at(pos, i), P3::at(pos, j), P3::at(pos, k));
            let cos_theta = clip_to_one(p1.sub(p2).dot(p3.sub(p2)) / (dist1 * dist2));
            res += angle_bend_energy(theta0, ka, linear, cos_theta);
        }
        res
    }

    fn angle_grad(
        &self,
        pos: &[f64],
        grad: &mut [f64],
        pair_distance: &mut impl FnMut(usize, usize) -> f64,
    ) {
        for &(i, j, k, theta0, ka, linear) in &self.angles {
            let dist = [pair_distance(i, j), pair_distance(j, k)];
            let (p1, p2, p3) = (P3::at(pos, i), P3::at(pos, j), P3::at(pos, k));
            let r = [p1.sub(p2).div(dist[0]), p3.sub(p2).div(dist[1])];
            let cos_theta = clip_to_one(r[0].dot(r[1]));
            let sin_theta_sq = 1.0 - cos_theta * cos_theta;
            let sin_theta = smax(
                if sin_theta_sq > 0.0 {
                    sin_theta_sq.sqrt()
                } else {
                    0.0
                },
                1.0e-8,
            );
            let angle_term = RAD2DEG * cos_theta.acos() - theta0;
            let cb = -0.006981317;
            let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
            let de_dtheta = if linear {
                -MDYNE_A_TO_KCAL_MOL * ka * sin_theta
            } else {
                RAD2DEG * c2 * ka * angle_term * (1.0 + 1.5 * cb * angle_term)
            };
            let dcos = [
                1.0 / dist[0] * (r[1].x - cos_theta * r[0].x),
                1.0 / dist[0] * (r[1].y - cos_theta * r[0].y),
                1.0 / dist[0] * (r[1].z - cos_theta * r[0].z),
                1.0 / dist[1] * (r[0].x - cos_theta * r[1].x),
                1.0 / dist[1] * (r[0].y - cos_theta * r[1].y),
                1.0 / dist[1] * (r[0].z - cos_theta * r[1].z),
            ];
            for d in 0..3 {
                grad[3 * i + d] += de_dtheta * dcos[d] / (-sin_theta);
            }
            for d in 0..3 {
                grad[3 * j + d] += de_dtheta * (-dcos[d] - dcos[d + 3]) / (-sin_theta);
            }
            for d in 0..3 {
                grad[3 * k + d] += de_dtheta * dcos[d + 3] / (-sin_theta);
            }
        }
    }

    // ─── Stretch-bend ────────────────────────────────────────────────────────

    fn stbn_energy(&self, pos: &[f64], pair_distance: &mut impl FnMut(usize, usize) -> f64) -> f64 {
        let mut total = 0.0;
        for &(i, j, k, r1, r2, theta0, k1, k2) in &self.stretch_bends {
            let dist1 = pair_distance(i, j);
            let dist2 = pair_distance(j, k);
            let (p1, p2, p3) = (P3::at(pos, i), P3::at(pos, j), P3::at(pos, k));
            let cos_theta = clip_to_one(p1.sub(p2).dot(p3.sub(p2)) / (dist1 * dist2));
            let delta_theta = RAD2DEG * cos_theta.acos() - theta0;
            let factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * delta_theta;
            let e1 = factor * k1 * (dist1 - r1);
            let e2 = factor * k2 * (dist2 - r2);
            total += e1 + e2;
        }
        total
    }

    fn stbn_grad(
        &self,
        pos: &[f64],
        grad: &mut [f64],
        pair_distance: &mut impl FnMut(usize, usize) -> f64,
    ) {
        for &(i, j, k, rest1, rest2, theta0, fc1, fc2) in &self.stretch_bends {
            let dist1 = pair_distance(i, j);
            let dist2 = pair_distance(j, k);
            let (p1, p2, p3) = (P3::at(pos, i), P3::at(pos, j), P3::at(pos, k));
            let p12 = p1.sub(p2).div(dist1);
            let p32 = p3.sub(p2).div(dist2);
            let c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
            let cos_theta = clip_to_one(p12.dot(p32));
            let sin_theta_sq = 1.0 - cos_theta * cos_theta;
            let sin_theta = smax(sin_theta_sq.sqrt(), 1.0e-8);
            let angle_term = RAD2DEG * cos_theta.acos() - theta0;
            let dist_term = RAD2DEG * (fc1 * (dist1 - rest1) + fc2 * (dist2 - rest2));
            let d1 = 1.0 / dist1 * (p32.x - cos_theta * p12.x);
            let d2 = 1.0 / dist1 * (p32.y - cos_theta * p12.y);
            let d3 = 1.0 / dist1 * (p32.z - cos_theta * p12.z);
            let d4 = 1.0 / dist2 * (p12.x - cos_theta * p32.x);
            let d5 = 1.0 / dist2 * (p12.y - cos_theta * p32.y);
            let d6 = 1.0 / dist2 * (p12.z - cos_theta * p32.z);
            grad[3 * i] += c5 * (p12.x * fc1 * angle_term + d1 / (-sin_theta) * dist_term);
            grad[3 * i + 1] += c5 * (p12.y * fc1 * angle_term + d2 / (-sin_theta) * dist_term);
            grad[3 * i + 2] += c5 * (p12.z * fc1 * angle_term + d3 / (-sin_theta) * dist_term);
            grad[3 * j] += c5
                * ((-p12.x * fc1 - p32.x * fc2) * angle_term
                    + (-d1 - d4) / (-sin_theta) * dist_term);
            grad[3 * j + 1] += c5
                * ((-p12.y * fc1 - p32.y * fc2) * angle_term
                    + (-d2 - d5) / (-sin_theta) * dist_term);
            grad[3 * j + 2] += c5
                * ((-p12.z * fc1 - p32.z * fc2) * angle_term
                    + (-d3 - d6) / (-sin_theta) * dist_term);
            grad[3 * k] += c5 * (p32.x * fc2 * angle_term + d4 / (-sin_theta) * dist_term);
            grad[3 * k + 1] += c5 * (p32.y * fc2 * angle_term + d5 / (-sin_theta) * dist_term);
            grad[3 * k + 2] += c5 * (p32.z * fc2 * angle_term + d6 / (-sin_theta) * dist_term);
        }
    }

    // ─── Out-of-plane bend ───────────────────────────────────────────────────

    fn oop_energy(&self, pos: &[f64]) -> f64 {
        let mut total = 0.0;
        let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        for &(i, j, k, l, koop) in &self.oops {
            let (ip, jp, kp, lp) = (
                P3::at(pos, i),
                P3::at(pos, j),
                P3::at(pos, k),
                P3::at(pos, l),
            );
            let mut rji = ip.sub(jp);
            let mut rjk = kp.sub(jp);
            let mut rjl = lp.sub(jp);
            rji = rji.div(rji.length());
            rjk = rjk.div(rjk.length());
            rjl = rjl.div(rjl.length());
            let mut nv = rji.cross(rjk);
            nv = nv.div(nv.length());
            let sin_chi = clip_to_one(nv.dot(rjl));
            let chi = RAD2DEG * sin_chi.asin();
            total += 0.5 * c2 * koop * chi * chi;
        }
        total
    }

    fn oop_grad(&self, pos: &[f64], grad: &mut [f64]) {
        let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        for &(i, j, k, l, koop) in &self.oops {
            let (ip, jp, kp, lp) = (
                P3::at(pos, i),
                P3::at(pos, j),
                P3::at(pos, k),
                P3::at(pos, l),
            );
            let mut rji = ip.sub(jp);
            let mut rjk = kp.sub(jp);
            let mut rjl = lp.sub(jp);
            let dji = rji.length();
            let djk = rjk.length();
            let djl = rjl.length();
            if is_double_zero(dji) || is_double_zero(djk) || is_double_zero(djl) {
                continue;
            }
            rji = rji.div(dji);
            rjk = rjk.div(djk);
            rjl = rjl.div(djl);
            let mut nv = rji.neg().cross(rjk);
            nv = nv.div(nv.length());
            let sin_chi = clip_to_one(rjl.dot(nv));
            let cos_chi_sq = 1.0 - sin_chi * sin_chi;
            let cos_chi = smax(
                if cos_chi_sq > 0.0 {
                    cos_chi_sq.sqrt()
                } else {
                    0.0
                },
                1.0e-8,
            );
            let chi = RAD2DEG * sin_chi.asin();
            let cos_theta = clip_to_one(rji.dot(rjk));
            let sin_theta_sq = smax(1.0 - cos_theta * cos_theta, 1.0e-8);
            let sin_theta = smax(
                if sin_theta_sq > 0.0 {
                    sin_theta_sq.sqrt()
                } else {
                    0.0
                },
                1.0e-8,
            );
            let de_dchi = RAD2DEG * c2 * koop * chi;
            let t1 = rjl.cross(rjk);
            let t2 = rji.cross(rjl);
            let t3 = rjk.cross(rji);
            let term1 = cos_chi * sin_theta;
            let term2 = sin_chi / (cos_chi * sin_theta_sq);
            let tg1 = [
                (t1.x / term1 - (rji.x - rjk.x * cos_theta) * term2) / dji,
                (t1.y / term1 - (rji.y - rjk.y * cos_theta) * term2) / dji,
                (t1.z / term1 - (rji.z - rjk.z * cos_theta) * term2) / dji,
            ];
            let tg3 = [
                (t2.x / term1 - (rjk.x - rji.x * cos_theta) * term2) / djk,
                (t2.y / term1 - (rjk.y - rji.y * cos_theta) * term2) / djk,
                (t2.z / term1 - (rjk.z - rji.z * cos_theta) * term2) / djk,
            ];
            let tg4 = [
                (t3.x / term1 - rjl.x * sin_chi / cos_chi) / djl,
                (t3.y / term1 - rjl.y * sin_chi / cos_chi) / djl,
                (t3.z / term1 - rjl.z * sin_chi / cos_chi) / djl,
            ];
            for d in 0..3 {
                grad[3 * i + d] += de_dchi * tg1[d];
                grad[3 * j + d] += -de_dchi * (tg1[d] + tg3[d] + tg4[d]);
                grad[3 * k + d] += de_dchi * tg3[d];
                grad[3 * l + d] += de_dchi * tg4[d];
            }
        }
    }

    // ─── Torsion ─────────────────────────────────────────────────────────────

    fn torsion_energy(&self, pos: &[f64]) -> f64 {
        let mut total = 0.0;
        for &(i, j, k, l, v1, v2, v3) in &self.torsions {
            let (ip, jp, kp, lp) = (
                P3::at(pos, i),
                P3::at(pos, j),
                P3::at(pos, k),
                P3::at(pos, l),
            );
            let r1 = ip.sub(jp);
            let r2 = kp.sub(jp);
            let r3 = jp.sub(kp);
            let r4 = lp.sub(kp);
            let t1 = r1.cross(r2);
            let t2 = r3.cross(r4);
            let t1_len = t1.length();
            let t2_len = t2.length();
            let cos_phi = if is_double_zero(t1_len) || is_double_zero(t2_len) {
                0.0
            } else {
                clip_to_one(t1.dot(t2) / (t1_len * t2_len))
            };
            let cos2 = 2.0 * cos_phi * cos_phi - 1.0;
            let cos3 = cos_phi * (2.0 * cos2 - 1.0);
            total += 0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2) + v3 * (1.0 + cos3));
        }
        total
    }

    fn torsion_grad(&self, pos: &[f64], grad: &mut [f64]) {
        for &(i, j, k, l, v1, v2, v3) in &self.torsions {
            let (p1, p2, p3, p4) = (
                P3::at(pos, i),
                P3::at(pos, j),
                P3::at(pos, k),
                P3::at(pos, l),
            );
            // ForceFieldsHelper::computeDihedral
            let r0 = p1.sub(p2);
            let r1 = p3.sub(p2);
            let r2 = r1.neg();
            let r3 = p4.sub(p3);
            let mut t0 = r0.cross(r1);
            let d0 = smax(t0.length(), 1.0e-5);
            t0 = t0.div(d0);
            let mut t1 = r2.cross(r3);
            let d1 = smax(t1.length(), 1.0e-5);
            t1 = t1.div(d1);
            let dp = t0.dot(t1);
            let cos_phi = smax(-1.0, if 1.0 < dp { 1.0 } else { dp });
            let sin_phi_sq = 1.0 - cos_phi * cos_phi;
            let sin_phi = if sin_phi_sq > 0.0 {
                sin_phi_sq.sqrt()
            } else {
                0.0
            };
            let sin2 = 2.0 * sin_phi * cos_phi;
            let sin3 = 3.0 * sin_phi - 4.0 * sin_phi * sin_phi_sq;
            let de_dphi = 0.5 * (-v1 * sin_phi + 2.0 * v2 * sin2 - 3.0 * v3 * sin3);
            let sin_term = -de_dphi
                * (if is_double_zero(sin_phi) {
                    1.0 / cos_phi
                } else {
                    1.0 / sin_phi
                });
            let dc = [
                1.0 / d0 * (t1.x - cos_phi * t0.x),
                1.0 / d0 * (t1.y - cos_phi * t0.y),
                1.0 / d0 * (t1.z - cos_phi * t0.z),
                1.0 / d1 * (t0.x - cos_phi * t1.x),
                1.0 / d1 * (t0.y - cos_phi * t1.y),
                1.0 / d1 * (t0.z - cos_phi * t1.z),
            ];
            let r = [r0, r1, r2, r3];
            grad[3 * i] += sin_term * (dc[2] * r[1].y - dc[1] * r[1].z);
            grad[3 * i + 1] += sin_term * (dc[0] * r[1].z - dc[2] * r[1].x);
            grad[3 * i + 2] += sin_term * (dc[1] * r[1].x - dc[0] * r[1].y);
            grad[3 * j] += sin_term
                * (dc[1] * (r[1].z - r[0].z)
                    + dc[2] * (r[0].y - r[1].y)
                    + dc[4] * (-r[3].z)
                    + dc[5] * (r[3].y));
            grad[3 * j + 1] += sin_term
                * (dc[0] * (r[0].z - r[1].z)
                    + dc[2] * (r[1].x - r[0].x)
                    + dc[3] * (r[3].z)
                    + dc[5] * (-r[3].x));
            grad[3 * j + 2] += sin_term
                * (dc[0] * (r[1].y - r[0].y)
                    + dc[1] * (r[0].x - r[1].x)
                    + dc[3] * (-r[3].y)
                    + dc[4] * (r[3].x));
            grad[3 * k] += sin_term
                * (dc[1] * (r[0].z)
                    + dc[2] * (-r[0].y)
                    + dc[4] * (r[3].z - r[2].z)
                    + dc[5] * (r[2].y - r[3].y));
            grad[3 * k + 1] += sin_term
                * (dc[0] * (-r[0].z)
                    + dc[2] * (r[0].x)
                    + dc[3] * (r[2].z - r[3].z)
                    + dc[5] * (r[3].x - r[2].x));
            grad[3 * k + 2] += sin_term
                * (dc[0] * (r[0].y)
                    + dc[1] * (-r[0].x)
                    + dc[3] * (r[3].y - r[2].y)
                    + dc[4] * (r[2].x - r[3].x));
            grad[3 * l] += sin_term * (dc[4] * r[2].z - dc[5] * r[2].y);
            grad[3 * l + 1] += sin_term * (dc[5] * r[2].x - dc[3] * r[2].z);
            grad[3 * l + 2] += sin_term * (dc[3] * r[2].y - dc[4] * r[2].x);
        }
    }

    // ─── Non-bonded (van der Waals + electrostatics) ─────────────────────────

    fn nonbonded_energy(&self, pair_distance: &mut impl FnMut(usize, usize) -> f64) -> f64 {
        let mut sum = 0.0;
        for t in &self.nonbonded {
            let dist = pair_distance(t.i, t.j);
            if let Some((r_star, well)) = t.vdw {
                sum += vdw_energy(dist, r_star, well);
            }
            if let Some((_, charge_term, is14)) = t.ele {
                let corr_dist = dist + 0.05;
                sum += 332.0716 * charge_term / corr_dist * (if is14 { 0.75 } else { 1.0 });
            }
        }
        sum
    }

    fn nonbonded_grad(
        &self,
        pos: &[f64],
        grad: &mut [f64],
        pair_distance: &mut impl FnMut(usize, usize) -> f64,
    ) {
        let vdw1 = 1.07;
        let vdw1m1 = vdw1 - 1.0;
        let vdw2 = 1.12;
        let vdw2m1 = vdw2 - 1.0;
        let vdw2t7 = vdw2 * 7.0;
        for t in &self.nonbonded {
            let dist = pair_distance(t.i, t.j);
            if dist <= 0.0 {
                // RDKit returns from the whole gradient here.
                if let Some((r_star, _)) = t.vdw {
                    for d in 0..3 {
                        grad[3 * t.i + d] += r_star * 0.01;
                        grad[3 * t.j + d] -= r_star * 0.01;
                    }
                }
                if t.ele.is_some() {
                    for d in 0..3 {
                        grad[3 * t.i + d] += 0.02;
                        grad[3 * t.j + d] -= 0.02;
                    }
                }
                return;
            }
            let mut vdw_grad = 0.0;
            let mut ele_grad = 0.0;
            if let Some((r_star, well)) = t.vdw {
                let q = dist / r_star;
                let q2 = q * q;
                let q6 = q2 * q2 * q2;
                let q7 = q6 * q;
                let q7pvdw2m1 = q7 + vdw2m1;
                let tt = vdw1 / (q + vdw1 - 1.0);
                let t2 = tt * tt;
                let t7 = t2 * t2 * t2 * tt;
                let de_dr = well / r_star
                    * t7
                    * (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1)
                        + ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
                vdw_grad = de_dr / dist;
            }
            if let Some((_, charge_term, is14)) = t.ele {
                let mut corr_dist = dist + 0.05;
                corr_dist *= corr_dist;
                let de_dr =
                    -332.0716 * 1.0 * charge_term / corr_dist * (if is14 { 0.75 } else { 1.0 });
                ele_grad = de_dr / dist;
            }
            let de_dr = vdw_grad + ele_grad;
            for d in 0..3 {
                let dgrad = de_dr * (pos[3 * t.i + d] - pos[3 * t.j + d]);
                grad[3 * t.i + d] += dgrad;
                grad[3 * t.j + d] -= dgrad;
            }
        }
    }

    // ─── Minimizer ───────────────────────────────────────────────────────────

    fn cached_energy(&self, pos: &[f64], cache: &RefCell<DistanceCache>) -> f64 {
        let mut cache = cache.borrow_mut();
        cache.prepare(pos);
        self.energy_with_distance_fn(pos, |i, j| cache.get(pos, i, j))
    }

    fn cached_scaled_gradient(
        &self,
        pos: &[f64],
        grad: &mut [f64],
        cache: &RefCell<DistanceCache>,
    ) -> f64 {
        for g in grad.iter_mut() {
            *g = 0.0;
        }
        let mut cache = cache.borrow_mut();
        cache.prepare(pos);
        self.gradient_with_distance_fn(pos, grad, |i, j| cache.get(pos, i, j));
        crate::rdkit_bfgs::scale_gradient(grad)
    }

    /// `ForceField::minimize(maxIts, forceTol, energyTol)` from the flat
    /// positions `pos` (updated in place): `BFGSOpt::minimize`. Returns
    /// RDKit's status (0 converged, 1 more iterations needed).
    pub fn minimize(&self, pos: &mut [f64], max_its: u32, force_tol: f64) -> i32 {
        if self.is_empty() {
            return 0;
        }
        let cache = RefCell::new(DistanceCache::new(self.n));
        crate::rdkit_bfgs::bfgs_minimize(
            pos,
            force_tol,
            max_its,
            |p| self.cached_energy(p, &cache),
            |p, g| self.cached_scaled_gradient(p, g, &cache),
        )
    }

    fn is_empty(&self) -> bool {
        let any_nb = (self.terms[5] || self.terms[6]) && !self.nonbonded.is_empty();
        !((self.terms[0] && !self.bonds.is_empty())
            || (self.terms[1] && !self.angles.is_empty())
            || (self.terms[2] && !self.stretch_bends.is_empty())
            || (self.terms[3] && !self.oops.is_empty())
            || (self.terms[4] && !self.torsions.is_empty())
            || any_nb)
    }

    /// RDKit `MMFFOptimizeMolecule(mol, maxIters)` (forceTol 1e-4) on the
    /// flat positions: (status, final energy).
    pub fn optimize(&self, pos: &mut [f64], max_iters: u32) -> (i32, f64) {
        let result = self.optimize_detailed(pos, max_iters);
        (result.status, result.energy)
    }

    /// RDKit-compatible MMFF94 optimization plus the already-computed final
    /// convergence diagnostic. This does not add an energy or gradient pass.
    pub fn optimize_detailed(&self, pos: &mut [f64], max_iters: u32) -> RdkitMmffOptimization {
        if self.is_empty() {
            return RdkitMmffOptimization {
                status: 0,
                energy: 0.0,
                iterations: 0,
                final_gradient_test: Some(0.0),
            };
        }
        let cache = RefCell::new(DistanceCache::new(self.n));
        let cached = crate::rdkit_bfgs::bfgs_minimize_detailed_with_cached_energy(
            pos,
            1e-4,
            max_iters,
            |p| self.cached_energy(p, &cache),
            |p, g| self.cached_scaled_gradient(p, g, &cache),
        );
        let energy = match &cached.outcome.stale_distance_point {
            Some(d) => self.energy_with_distances(pos, d),
            None => cached
                .energy_at_returned_position
                .unwrap_or_else(|| self.energy(pos)),
        };
        RdkitMmffOptimization {
            status: cached.outcome.status,
            energy,
            iterations: cached.iterations,
            final_gradient_test: cached.final_gradient_test,
        }
    }
}

fn angle_bend_energy(theta0: f64, ka: f64, linear: bool, cos_theta: f64) -> f64 {
    let angle = RAD2DEG * cos_theta.acos() - theta0;
    let cb = -0.006981317;
    let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    if linear {
        MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cos_theta)
    } else {
        0.5 * c2 * ka * angle * angle * (1.0 + cb * angle)
    }
}

fn vdw_energy(dist: f64, r_star_ij: f64, well_depth: f64) -> f64 {
    let vdw1 = 1.07;
    let vdw1m1 = vdw1 - 1.0;
    let vdw2 = 1.12;
    let vdw2m1 = vdw2 - 1.0;
    let dist2 = dist * dist;
    let dist7 = dist2 * dist2 * dist2 * dist;
    let a_term = vdw1 * r_star_ij / (dist + vdw1m1 * r_star_ij);
    let a_term2 = a_term * a_term;
    let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
    let r2 = r_star_ij * r_star_ij;
    let r7 = r2 * r2 * r2 * r_star_ij;
    let b_term = vdw2 * r7 / (dist7 + vdw2m1 * r7) - 2.0;
    well_depth * a_term7 * b_term
}

/// `MMFFMolProperties::getMMFFVdWParams`: (R_ij_star, epsilon), RDKit's
/// combination rules on the MMFFVDW.PAR per-type values.
fn vdw_params(ti: u8, tj: u8) -> Option<(f64, f64)> {
    const POWER: f64 = 0.25;
    const B: f64 = 0.2;
    const BETA: f64 = 12.0;
    const DARAD: f64 = 0.8;
    const DAEPS: f64 = 0.5;
    let pi = crate::mmff94_energy::mmff94_vdw_energy(ti)?;
    let pj = crate::mmff94_energy::mmff94_vdw_energy(tj)?;
    let ri = pi.a_i * pi.alpha_i.powf(POWER);
    let rj = pj.a_i * pj.alpha_i.powf(POWER);
    // DA codes: 1 = donor, 2 = acceptor.
    let gamma = (ri - rj) / (ri + rj);
    let r_ij = 0.5
        * (ri + rj)
        * (1.0
            + if pi.da == 1 || pj.da == 1 {
                0.0
            } else {
                B * (1.0 - (-BETA * gamma * gamma).exp())
            });
    let r2 = r_ij * r_ij;
    let c4 = 181.16;
    let mut eps = c4 * pi.g_i * pj.g_i * pi.alpha_i * pj.alpha_i
        / (((pi.alpha_i / pi.n_i).sqrt() + (pj.alpha_i / pj.n_i).sqrt()) * r2 * r2 * r2);
    let mut r = r_ij;
    if (pi.da == 1 && pj.da == 2) || (pi.da == 2 && pj.da == 1) {
        r *= DARAD;
        eps *= DAEPS;
    }
    Some((r, eps))
}

/// RDKit's bond numbering for a molecule whose trailing hydrogens were
/// added by `add_hydrogens` (RDKit's `AddHs`): the parent's bonds in RDKit
/// order (SMILES ring closures after the chain bonds), then the added
/// hydrogen bonds in creation order.
pub(crate) fn rdkit_bond_order_with_added_hs(mol: &Molecule) -> Vec<chematic_core::BondIdx> {
    let n = mol.atom_count();
    let mut h0 = n;
    while h0 > 0 {
        let a = mol.atom(AtomIdx((h0 - 1) as u32));
        if a.element.atomic_number() == 1
            && !a.wildcard
            && mol.degree(AtomIdx((h0 - 1) as u32)) == 1
        {
            h0 -= 1;
        } else {
            break;
        }
    }
    let added = |b: &chematic_core::BondIdx| {
        let bond = mol.bond(*b);
        bond.atom1.0 as usize >= h0 || bond.atom2.0 as usize >= h0
    };
    let mut order: Vec<chematic_core::BondIdx> = mol
        .rdkit_bond_order()
        .into_iter()
        .filter(|b| !added(b))
        .collect();
    let mut hs: Vec<chematic_core::BondIdx> = (0..mol.bond_count() as u32)
        .map(chematic_core::BondIdx)
        .filter(added)
        .collect();
    hs.sort_unstable_by_key(|b| b.0);
    order.extend(hs);
    order
}

/// Topological distances (`n * n`, 0 on the diagonal, `u32::MAX` between
/// fragments).
pub(crate) fn topological_distances(n: usize, adj: &[Vec<(usize, usize)>]) -> Vec<u32> {
    let mut d = vec![u32::MAX; n * n];
    for s in 0..n {
        let row = &mut d[s * n..(s + 1) * n];
        row[s] = 0;
        let mut queue = std::collections::VecDeque::from([s]);
        while let Some(a) = queue.pop_front() {
            for &(b, _) in &adj[a] {
                if row[b] == u32::MAX {
                    row[b] = row[a] + 1;
                    queue.push_back(b);
                }
            }
        }
    }
    d
}

pub(crate) fn fragments(n: usize, adj: &[Vec<(usize, usize)>]) -> Vec<usize> {
    let mut id = vec![usize::MAX; n];
    let mut next = 0;
    for s in 0..n {
        if id[s] != usize::MAX {
            continue;
        }
        id[s] = next;
        let mut stack = vec![s];
        while let Some(a) = stack.pop() {
            for &(b, _) in &adj[a] {
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

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// `CC(=O)Nc1ccc(O)cc1` after `AddHs`, RDKit 2026.03.1
    /// `EmbedMolecule(randomSeed=42)` coordinates.
    pub(crate) const PARACETAMOL: &str = "CC(=O)Nc1ccc(O)cc1";
    pub(crate) const COORDS: [[f64; 3]; 20] = [
        [3.762175028044453, 0.30538563337732244, -0.09745815343294277],
        [2.3329369300899403, 0.11610000936619912, -0.5176254421062628],
        [
            2.1306032462731035,
            -0.09593752049575305,
            -1.7310611678300563,
        ],
        [1.3484147775946786, 0.19222362621777822, 0.4794175681721041],
        [
            -0.04651062589597437,
            0.03729534973500948,
            0.2549743150386899,
        ],
        [-0.6601760831775039, -1.210120317536445, 0.3277557282622417],
        [-2.027875532733807, -1.3550238887272026, 0.10667566446336685],
        [-2.822456517029374, -0.2670983319336878, -0.1922446077538659],
        [-4.187714927031642, -0.3854347185261268, -0.416753109003994],
        [-2.215264059602086, 0.990686079372787, -0.2676816459777285],
        [
            -0.8506489653837782,
            1.1179323355827404,
            -0.04458045101770329,
        ],
        [3.8141231841990857, 1.3434050655543417, 0.3024234778452463],
        [4.477057254518884, 0.10425801972676844, -0.9145679992818092],
        [3.909133114355366, -0.42656371299342444, 0.7351586531314641],
        [1.6266288380623064, 0.37726586194308837, 1.4846206751064692],
        [
            -0.03711779601633895,
            -2.0820010342352497,
            0.5648832746587801,
        ],
        [-2.5130630664729114, -2.314409067782831, 0.16030971994515217],
        [-4.834672252546139, -0.3986600801567403, 0.37067475360645646],
        [-2.8406901294426152, 1.8542155191372227, -0.5040194618516627],
        [-0.3648824178056287, 2.096481172374245, -0.10090179197395394],
    ];

    pub(crate) fn paracetamol() -> Molecule {
        chematic_chem::add_hydrogens(&chematic_smiles::parse(PARACETAMOL).unwrap())
    }

    fn check(variant: MmffVariant, energy: f64, opt_energy: f64, x0: f64) {
        let mol = paracetamol();
        let ff =
            RdkitMmffField::with_terms(&mol, &COORDS, variant, 100.0, true, [true; 7]).unwrap();
        let mut pos = COORDS.as_flattened().to_vec();
        let initial = ff.energy(&pos);
        assert!(
            (initial - energy).abs() <= 1e-12,
            "initial energy {initial} differs from RDKit {energy}"
        );
        let (status, e) = ff.optimize(&mut pos, 200);
        assert_eq!(status, 0);
        assert!(
            (e - opt_energy).abs() <= 1e-9,
            "optimized energy {e} differs from RDKit {opt_energy} by {}",
            (e - opt_energy).abs()
        );
        assert!(
            (pos[0] - x0).abs() <= 1e-6,
            "first coordinate {} differs from RDKit {x0}",
            pos[0]
        );

        let mut detailed_pos = COORDS.as_flattened().to_vec();
        let detailed = ff.optimize_detailed(&mut detailed_pos, 200);
        assert_eq!(detailed.status, status);
        assert_eq!(detailed.energy, e);
        assert_eq!(detailed_pos, pos);
        assert!(detailed.iterations > 0);
    }

    /// RDKit 2026.03.1: `MMFFGetMoleculeForceField(...).CalcEnergy()` and
    /// `Minimize(maxIts=200)` (status, energy, first coordinate).
    #[test]
    fn mmff94_matches_rdkit_with_platform_tolerance() {
        check(
            MmffVariant::Mmff94,
            22.565180175565263,
            -12.775974493568183,
            3.6705868960020265,
        );
    }

    #[test]
    fn mmff94s_matches_rdkit_with_platform_tolerance() {
        check(
            MmffVariant::Mmff94s,
            22.959370426002092,
            -11.453744586799015,
            3.66682072353845,
        );
    }

    #[test]
    fn distance_cache_matches_public_energy_and_gradient() {
        let mol = paracetamol();
        let ff = RdkitMmffField::new(&mol, &COORDS, 100.0, true).unwrap();
        let pos = COORDS.as_flattened();
        let cache = RefCell::new(DistanceCache::new(mol.atom_count()));
        assert_eq!(ff.cached_energy(pos, &cache), ff.energy(pos));

        let mut cached_gradient = vec![0.0; pos.len()];
        let mut public_gradient = vec![0.0; pos.len()];
        ff.cached_scaled_gradient(pos, &mut cached_gradient, &cache);
        ff.gradient(pos, &mut public_gradient);
        crate::rdkit_bfgs::scale_gradient(&mut public_gradient);
        assert_eq!(cached_gradient, public_gradient);

        let generation = cache.borrow().generation;
        ff.cached_energy(pos, &cache);
        assert_eq!(cache.borrow().generation, generation);
    }
}

#[cfg(test)]
mod nonbonded_collision_boundary_tests {
    use super::*;

    #[test]
    fn public_default_field_matches_native_ionic_pairs_and_collision_guards() {
        let fixture: serde_json::Value = serde_json::from_str(include_str!(
            "../../../validation/rdkit-2026.03.1-mmff-nonbonded-boundary.json"
        ))
        .unwrap();
        assert_eq!(fixture["rdkit_version"], "2026.03.1");
        for row in fixture["rows"].as_array().unwrap() {
            let source = row["smiles"].as_str().unwrap();
            let mol = chematic_smiles::parse(source).unwrap();
            let types = crate::mmff94_numeric::assign_mmff94_numeric_types(&mol).unwrap();
            let charges = crate::mmff94_numeric::mmff94_charges_numeric(&mol).unwrap();
            for i in 0..2 {
                assert_eq!(
                    u64::from(types[i]),
                    row["atom_types"][i].as_u64().unwrap(),
                    "{source}"
                );
                assert!(
                    (charges[i] - row["charges"][i].as_f64().unwrap()).abs() < 1e-12,
                    "{source}"
                );
            }
            let distance = row["distance"].as_f64().unwrap();
            let ignore = row["ignore_interfrag"].as_bool().unwrap();
            let coords = [[0., 0., 0.], [distance, 0., 0.]];
            let field = RdkitMmffField::new(&mol, &coords, 100., ignore).unwrap();
            assert_eq!(field.num_atoms(), 2);
            let pos = [0., 0., 0., distance, 0., 0.];
            let energy = field.energy(&pos);
            let expected = row["energy"].as_f64().unwrap();
            assert!(
                energy.is_finite() && (energy - expected).abs() < 1e-8 + expected.abs() * 1e-10,
                "{source},d={distance},ignore={ignore}: {energy} != {expected}"
            );
            let mut gradient = [0.; 6];
            field.gradient(&pos, &mut gradient);
            for (axis, value) in row["gradient"].as_array().unwrap().iter().enumerate() {
                let expected = value.as_f64().unwrap();
                assert!(
                    gradient[axis].is_finite()
                        && (gradient[axis] - expected).abs() < 1e-9 + expected.abs() * 1e-10,
                    "{source},d={distance},ignore={ignore},axis={axis}: {} != {expected}",
                    gradient[axis]
                );
            }
            for axis in 0..3 {
                assert!((gradient[axis] + gradient[axis + 3]).abs() < 1e-12);
            }
            if ignore {
                assert_eq!(energy, 0.);
                assert_eq!(gradient, [0.; 6]);
            }
            // At exactly coincident centres RDKit returns a deterministic
            // guard vector, not the derivative of a differentiable geometry.
            if distance >= 1. && !ignore {
                for axis in 0..6 {
                    let mut plus = pos;
                    let mut minus = pos;
                    plus[axis] += 1e-5;
                    minus[axis] -= 1e-5;
                    let finite_difference = (field.energy(&plus) - field.energy(&minus)) / (2e-5);
                    assert!(
                        (gradient[axis] - finite_difference).abs()
                            < 1e-4 + gradient[axis].abs() * 2e-6,
                        "{source},d={distance},axis={axis}"
                    );
                }
            }
        }
        assert_eq!(fixture["rows"].as_array().unwrap().len(), 104);
    }
}
