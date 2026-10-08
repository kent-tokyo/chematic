//! RDKit's UFF force field, operation for operation.
//!
//! [`RdkitUffField`] holds the contributions RDKit 2026.03's
//! `UFF::constructForceField` builds (bond stretch, angle bend including the
//! trigonal-bipyramidal special case, van der Waals, torsion and inversion),
//! in RDKit's order and with RDKit's energy and gradient expressions
//! (`Code/ForceField/UFF/*.cpp`), so `UFFGetMoleculeForceField(...)
//! .CalcEnergy()` / `CalcGrad()` and `UFFOptimizeMolecule` are reproduced bit
//! for bit. Atom typing is RDKit's `UFF::Tools::getAtomLabel` on chematic's
//! RDKit-compatible hybridization, aromaticity and conjugation perception.
//!
//! The molecule must carry explicit hydrogens in RDKit `AddHs` order
//! (chematic's `add_hydrogens`), as RDKit's own UFF expects.

// RDKit's spelling (`x * -1.0`, `-1. * k`) is kept: those are exact
// negations, so rewriting them would not change any result.
#![allow(clippy::neg_multiply, clippy::type_complexity)]

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

use crate::rdkit_mmff::{
    P3, clip_to_one, distance, fragments, is_double_zero, rdkit_bond_order_with_added_hs, smax,
    topological_distances,
};
use crate::rdkit_uff_params::{
    RDKIT_DEFAULT_VALENCE, RDKIT_SYMBOL, RDKIT_VALENCE_LIST, UFF_NEEDS_HYBRID, UFF_PARAMS,
};

/// `UFF::Params::G`, the bond force constant prefactor.
const G: f64 = 332.06;
/// `UFF::Params::lambda`, the bond-order correction scale.
const LAMBDA: f64 = 0.1332;
const ANGLE_CORRECTION_THRESHOLD: f64 = 0.8660;

/// RDKit `UFF::AtomicParams` (`theta0` in radians).
#[derive(Clone, Copy, Debug)]
struct AtomicParams {
    r1: f64,
    theta0: f64,
    x1: f64,
    d1: f64,
    z1: f64,
    v1: f64,
    u1: f64,
    xi: f64,
}

fn params_for(label: &str) -> Option<AtomicParams> {
    let i = UFF_PARAMS.binary_search_by(|(l, _)| l.cmp(&label)).ok()?;
    let v = UFF_PARAMS[i].1;
    Some(AtomicParams {
        r1: v[0],
        theta0: v[1] * std::f64::consts::PI / 180.0,
        x1: v[2],
        d1: v[3],
        z1: v[5],
        v1: v[6],
        u1: v[7],
        xi: v[8],
    })
}

/// RDKit hybridization (`Atom::HybridizationType`) as chematic's SMARTS
/// `^n` code: 0 S, 1 SP, 2 SP2, 3 SP3, 4 SP3D, 5 SP3D2; `None` unspecified.
fn hybridization(mol: &Molecule) -> Vec<Option<u8>> {
    (0..mol.atom_count())
        .map(|i| chematic_smarts::rdkit_hybridization(mol, AtomIdx(i as u32)))
        .collect()
}

fn atomic_number(mol: &Molecule, i: usize) -> usize {
    let a = mol.atom(AtomIdx(i as u32));
    if a.wildcard {
        0
    } else {
        a.element.atomic_number() as usize
    }
}

/// `Bond::getBondTypeAsDouble`.
fn bond_as_double(order: BondOrder) -> f64 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
        _ => 0.0,
    }
}

/// RDKit `getTotalValence` of an atom of an explicit-H molecule (no implicit
/// or bracket hydrogens left): `calculateExplicitValence`.
fn total_valence(mol: &Molecule, i: usize) -> i32 {
    let idx = AtomIdx(i as u32);
    let mut accum = 0.0f64;
    for (_, b) in mol.neighbors(idx) {
        let bond = mol.bond(b);
        accum += match bond.order {
            BondOrder::Dative if bond.atom2 != idx => 0.0,
            o => bond_as_double(o),
        };
    }
    let z = atomic_number(mol, i);
    let atom = mol.atom(idx);
    let mut eff = z;
    let ovalens = RDKIT_VALENCE_LIST[z.min(118)];
    if ovalens[1] != -1 || ovalens[0] != -1 {
        let e = z as i32 - atom.charge as i32;
        if (0..=118).contains(&e) {
            eff = e as usize;
        }
    }
    let dv = RDKIT_DEFAULT_VALENCE[eff] as f64;
    if accum > dv && atom.aromatic {
        let mut pval = dv;
        for &v in &RDKIT_VALENCE_LIST[eff] {
            if v == -1 {
                break;
            }
            if v as f64 > accum {
                break;
            }
            pval = v as f64;
        }
        if accum - pval <= 1.5 {
            accum = pval;
        }
    }
    accum += 0.1;
    accum.round() as i32
}

/// RDKit `UFF::Tools::getAtomLabel` (with `addAtomChargeFlags`'s default
/// `tolerateChargeMismatch = true`).
fn atom_label(mol: &Molecule, i: usize, hyb: Option<u8>, conjugated: bool) -> String {
    let z = atomic_number(mol, i);
    let atom = mol.atom(AtomIdx(i as u32));
    let mut key = RDKIT_SYMBOL[z.min(118)].to_string();
    if key.len() == 1 {
        key.push('_');
    }
    if z != 0 && UFF_NEEDS_HYBRID[z.min(118)] {
        match z {
            12..=15 | 50..=52 | 81..=84 => key.push('3'),
            80 => key.push('1'),
            _ => match hyb {
                Some(1) => key.push('1'),
                Some(2) => {
                    if (atom.aromatic || conjugated) && matches!(z, 6 | 7 | 8 | 16) {
                        key.push('R');
                    } else {
                        key.push('2');
                    }
                }
                Some(3) => key.push('3'),
                Some(4) => key.push('5'),
                Some(5) => key.push('6'),
                _ => {}
            },
        }
    }
    // addAtomChargeFlags
    let suffix = match z {
        29 | 47 => "+1",
        4 | 20 | 25 | 26 | 28 | 46 | 78 => "+2",
        21 | 24 | 27 | 79 | 89 | 96..=103 => "+3",
        2 | 18 | 22 | 36 | 54 | 90..=95 => "+4",
        23 | 41 | 43 | 73 => "+5",
        42 => "+6",
        12 | 30 | 34 | 48 | 52 | 80 | 84 => "+2",
        31 | 33 | 49 | 51 | 81 | 82 | 83 => "+3",
        15 => {
            if total_valence(mol, i) == 3 {
                "+3"
            } else {
                "+5"
            }
        }
        16 if hyb != Some(2) => match total_valence(mol, i) {
            2 => "+2",
            4 => "+4",
            _ => "+6",
        },
        75 => {
            if key == "Re6" {
                key = "Re6+5".to_string();
            } else if key == "Re3" {
                key = "Re3+7".to_string();
            }
            ""
        }
        _ => "",
    };
    key.push_str(suffix);
    if (57..=71).contains(&z) {
        key.push_str("+3");
    }
    key
}

fn conjugated_atoms(mol: &Molecule) -> Vec<bool> {
    let mut conj = vec![false; mol.atom_count()];
    for b in 0..mol.bond_count() {
        if chematic_smarts::rdkit_bond_is_conjugated(mol, BondIdx(b as u32)) {
            let bond = mol.bond(BondIdx(b as u32));
            conj[bond.atom1.0 as usize] = true;
            conj[bond.atom2.0 as usize] = true;
        }
    }
    conj
}

/// RDKit's UFF atom labels (`UFF::Tools::getAtomLabel`) of every atom of an
/// explicit-hydrogen molecule.
pub fn rdkit_uff_atom_labels(mol: &Molecule) -> Vec<String> {
    let view = rdkit_view(mol);
    let mol = &view;
    let hyb = hybridization(mol);
    let conj = conjugated_atoms(mol);
    (0..mol.atom_count())
        .map(|i| atom_label(mol, i, hyb[i], conj[i]))
        .collect()
}

/// RDKit `UFFHasAllMoleculeParams`: every atom's UFF label has parameters.
pub fn rdkit_uff_has_all_params(mol: &Molecule) -> bool {
    rdkit_uff_atom_labels(mol)
        .iter()
        .all(|l| params_for(l).is_some())
}

/// The molecule as RDKit's sanitization perceives it (RDKit aromaticity,
/// cleaned-up functional groups), index for index; `mol` itself when that
/// view cannot be built.
fn rdkit_view(mol: &Molecule) -> Molecule {
    chematic_perception::with_rdkit_parity_view(mol, |v| match v {
        Ok(view)
            if view.atom_count() == mol.atom_count() && view.bond_count() == mol.bond_count() =>
        {
            view.clone()
        }
        _ => mol.clone(),
    })
}

/// RDKit `UFF::Utils::calcBondRestLength(bondOrder, params(label_i),
/// params(label_j))`; `None` when either label has no UFF parameters (RDKit's
/// `getAtomTypes` then yields a null parameter pointer).
pub fn rdkit_uff_bond_rest_length(bond_order: f64, label_i: &str, label_j: &str) -> Option<f64> {
    let pi = params_for(label_i)?;
    let pj = params_for(label_j)?;
    Some(calc_bond_rest_length(bond_order, &pi, &pj))
}

// ─── Parameter formulas (`Code/ForceField/UFF/*.cpp` `Utils`) ───────────────

fn calc_bond_rest_length(bond_order: f64, p1: &AtomicParams, p2: &AtomicParams) -> f64 {
    let ri = p1.r1;
    let rj = p2.r1;
    let r_bo = -LAMBDA * (ri + rj) * bond_order.ln();
    let xi = p1.xi;
    let xj = p2.xi;
    let r_en = ri * rj * (xi.sqrt() - xj.sqrt()) * (xi.sqrt() - xj.sqrt()) / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
}

fn calc_bond_force_constant(rest_length: f64, p1: &AtomicParams, p2: &AtomicParams) -> f64 {
    2.0 * G * p1.z1 * p2.z1 / (rest_length * rest_length * rest_length)
}

/// RDKit `int_pow<n>`: recursive halving.
fn int_pow(x: f64, n: u32) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let half = int_pow(x, n / 2);
            if n.is_multiple_of(2) {
                half * half
            } else {
                half * half * x
            }
        }
    }
}

fn calc_angle_force_constant(
    theta0: f64,
    bo12: f64,
    bo23: f64,
    p1: &AtomicParams,
    p2: &AtomicParams,
    p3: &AtomicParams,
) -> f64 {
    let cos_theta0 = theta0.cos();
    let r12 = calc_bond_rest_length(bo12, p1, p2);
    let r23 = calc_bond_rest_length(bo23, p2, p3);
    let r13 = (r12 * r12 + r23 * r23 - 2. * r12 * r23 * cos_theta0).sqrt();
    let beta = 2. * G / (r12 * r23);
    let pre_factor = beta * p1.z1 * p3.z1 / int_pow(r13, 5);
    let r_term = r12 * r23;
    let inner_bit = 3. * r_term * (1. - cos_theta0 * cos_theta0) - r13 * r13 * cos_theta0;
    pre_factor * r_term * inner_bit
}

fn is_in_group6(num: usize) -> bool {
    matches!(num, 8 | 16 | 34 | 52 | 84)
}

fn equation17(bo23: f64, p2: &AtomicParams, p3: &AtomicParams) -> f64 {
    5. * (p2.u1 * p3.u1).sqrt() * (1. + 4.18 * bo23.ln())
}

/// `calcInversionCoefficientsAndForceConstant`: (K, C0, C1, C2).
fn inversion_coefficients(at2: usize, is_c_bound_to_o: bool) -> (f64, f64, f64, f64) {
    let (mut res, c0, c1, c2);
    if matches!(at2, 6..=8) {
        c0 = 1.0;
        c1 = -1.0;
        c2 = 0.0;
        res = if is_c_bound_to_o { 50.0 } else { 6.0 };
    } else {
        let mut w0 = std::f64::consts::PI / 180.0;
        match at2 {
            15 => w0 *= 84.4339,
            33 => w0 *= 86.9735,
            51 => w0 *= 87.7047,
            83 => w0 *= 90.0,
            _ => {}
        }
        c2 = 1.0;
        c1 = -4.0 * w0.cos();
        c0 = -(c1 * w0.cos() + c2 * (2.0 * w0).cos());
        res = 22.0 / (c0 + c1 + c2);
    }
    res /= 3.0;
    (res, c0, c1, c2)
}

// ─── Contributions ───────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
struct AngleTerm {
    i: usize,
    j: usize,
    k: usize,
    theta0: f64,
    ka: f64,
    order: u32,
    c0: f64,
    c1: f64,
    c2: f64,
}

impl AngleTerm {
    #[allow(clippy::too_many_arguments)]
    fn new(
        i: usize,
        j: usize,
        k: usize,
        bo12: f64,
        bo23: f64,
        p1: &AtomicParams,
        p2: &AtomicParams,
        p3: &AtomicParams,
        order: u32,
    ) -> Self {
        let mut theta0 = p2.theta0;
        let mut order = order;
        if order >= 30 {
            match order {
                30 => theta0 = 150.0 / 180.0 * std::f64::consts::PI,
                35 => theta0 = 60.0 / 180.0 * std::f64::consts::PI,
                40 => theta0 = 135.0 / 180.0 * std::f64::consts::PI,
                45 => theta0 = 90.0 / 180.0 * std::f64::consts::PI,
                _ => {}
            }
            order = 0;
        }
        let ka = calc_angle_force_constant(theta0, bo12, bo23, p1, p2, p3);
        let (mut c0, mut c1, mut c2) = (0.0, 0.0, 0.0);
        if order == 0 {
            let sin_theta0 = theta0.sin();
            let cos_theta0 = theta0.cos();
            c2 = 1. / (4. * smax(sin_theta0 * sin_theta0, 1e-8));
            c1 = -4. * c2 * cos_theta0;
            c0 = c2 * (2. * cos_theta0 * cos_theta0 + 1.);
        }
        AngleTerm {
            i,
            j,
            k,
            theta0,
            ka,
            order,
            c0,
            c1,
            c2,
        }
    }

    fn energy_term(&self, cos_theta: f64, sin_theta_sq: f64) -> f64 {
        let cos2_theta = cos_theta * cos_theta - sin_theta_sq;
        if self.order == 0 {
            self.c0 + self.c1 * cos_theta + self.c2 * cos2_theta
        } else {
            let res = match self.order {
                1 => -cos_theta,
                2 => cos2_theta,
                3 => cos_theta * (cos_theta * cos_theta - 3. * sin_theta_sq),
                4 => {
                    int_pow(cos_theta, 4) - 6. * cos_theta * cos_theta * sin_theta_sq
                        + sin_theta_sq * sin_theta_sq
                }
                _ => 0.0,
            };
            (1. - res) / (self.order * self.order) as f64
        }
    }

    fn theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> f64 {
        let sin2_theta = 2. * sin_theta * cos_theta;
        if self.order == 0 {
            -1. * self.ka * (self.c1 * sin_theta + 2. * self.c2 * sin2_theta)
        } else {
            let d = match self.order {
                1 => -sin_theta,
                2 => sin2_theta,
                3 => sin_theta * (3. - 4. * sin_theta * sin_theta),
                4 => cos_theta * sin_theta * (4. - 8. * sin_theta * sin_theta),
                _ => 0.0,
            };
            d * (self.ka / self.order as f64)
        }
    }

    fn energy(&self, pos: &[f64], dpos: &[f64]) -> f64 {
        let dist1 = distance(dpos, self.i, self.j);
        let dist2 = distance(dpos, self.j, self.k);
        let p1 = P3::at(pos, self.i);
        let p2 = P3::at(pos, self.j);
        let p3 = P3::at(pos, self.k);
        let p12 = p1.sub(p2);
        let p32 = p3.sub(p2);
        let cos_theta = clip_to_one(p12.dot(p32) / (dist1 * dist2));
        let sin_theta_sq = 1. - cos_theta * cos_theta;
        let angle_term = self.energy_term(cos_theta, sin_theta_sq);
        let mut res = self.ka * angle_term;
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            let theta = cos_theta.acos();
            res += (-20.0 * (theta - self.theta0 + 0.25)).exp();
        }
        res
    }

    fn grad(&self, pos: &[f64], grad: &mut [f64]) {
        let dist = [distance(pos, self.i, self.j), distance(pos, self.j, self.k)];
        let p1 = P3::at(pos, self.i);
        let p2 = P3::at(pos, self.j);
        let p3 = P3::at(pos, self.k);
        let r = [p1.sub(p2).div(dist[0]), p3.sub(p2).div(dist[1])];
        let cos_theta = clip_to_one(r[0].dot(r[1]));
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = smax(sin_theta_sq.sqrt(), 1.0e-8);
        let mut de_dtheta = self.theta_deriv(cos_theta, sin_theta);
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            let theta = cos_theta.acos();
            let corr = -20.0 * (-20.0 * (theta - self.theta0 + 0.25)).exp();
            de_dtheta += corr;
        }
        angle_bend_grad(
            &r, &dist, grad, self.i, self.j, self.k, de_dtheta, cos_theta, sin_theta,
        );
    }
}

/// `UFF::Utils::calcAngleBendGrad`.
#[allow(clippy::too_many_arguments)]
fn angle_bend_grad(
    r: &[P3; 2],
    dist: &[f64; 2],
    grad: &mut [f64],
    i: usize,
    j: usize,
    k: usize,
    de_dtheta: f64,
    cos_theta: f64,
    sin_theta: f64,
) {
    let dcos = [
        1.0 / dist[0] * (r[1].x - cos_theta * r[0].x),
        1.0 / dist[0] * (r[1].y - cos_theta * r[0].y),
        1.0 / dist[0] * (r[1].z - cos_theta * r[0].z),
        1.0 / dist[1] * (r[0].x - cos_theta * r[1].x),
        1.0 / dist[1] * (r[0].y - cos_theta * r[1].y),
        1.0 / dist[1] * (r[0].z - cos_theta * r[1].z),
    ];
    let ms = -sin_theta;
    grad[3 * i] += de_dtheta * dcos[0] / ms;
    grad[3 * i + 1] += de_dtheta * dcos[1] / ms;
    grad[3 * i + 2] += de_dtheta * dcos[2] / ms;
    grad[3 * j] += de_dtheta * (-dcos[0] - dcos[3]) / ms;
    grad[3 * j + 1] += de_dtheta * (-dcos[1] - dcos[4]) / ms;
    grad[3 * j + 2] += de_dtheta * (-dcos[2] - dcos[5]) / ms;
    grad[3 * k] += de_dtheta * dcos[3] / ms;
    grad[3 * k + 1] += de_dtheta * dcos[4] / ms;
    grad[3 * k + 2] += de_dtheta * dcos[5] / ms;
}

#[derive(Clone, Copy, Debug)]
struct TorsionTerm {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    force_constant: f64,
    order: u32,
    cos_term: f64,
}

impl TorsionTerm {
    fn energy(&self, pos: &[f64]) -> f64 {
        let p1 = P3::at(pos, self.i);
        let p2 = P3::at(pos, self.j);
        let p3 = P3::at(pos, self.k);
        let p4 = P3::at(pos, self.l);
        // UFF::Utils::calculateCosTorsion
        let r1 = p1.sub(p2);
        let r2 = p3.sub(p2);
        let r3 = p2.sub(p3);
        let r4 = p4.sub(p3);
        let t1 = r1.cross(r2);
        let t2 = r3.cross(r4);
        let d1 = t1.length();
        let d2 = t2.length();
        let cos_phi = if is_double_zero(d1) || is_double_zero(d2) {
            0.0
        } else {
            clip_to_one(t1.dot(t2) / (d1 * d2))
        };
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        let cos_n_phi = match self.order {
            2 => 1.0 - 2.0 * sin_phi_sq,
            3 => cos_phi * (cos_phi * cos_phi - 3. * sin_phi_sq),
            6 => 1.0 + sin_phi_sq * (-32. * sin_phi_sq * sin_phi_sq + 48. * sin_phi_sq - 18.),
            _ => 0.0,
        };
        self.force_constant / 2.0 * (1. - self.cos_term * cos_n_phi)
    }

    fn theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> f64 {
        let sin_theta_sq = sin_theta * sin_theta;
        let res = match self.order {
            2 => 2.0 * sin_theta * cos_theta,
            3 => sin_theta * (3.0 - 4.0 * sin_theta_sq),
            6 => cos_theta * sin_theta * (32.0 * sin_theta_sq * (sin_theta_sq - 1.0) + 6.0),
            _ => 0.0,
        };
        res * (self.force_constant / 2.0 * self.cos_term * -1.0 * self.order as f64)
    }

    fn grad(&self, pos: &[f64], grad: &mut [f64]) {
        let (i, j, k, l) = (self.i, self.j, self.k, self.l);
        let p1 = P3::at(pos, i);
        let p2 = P3::at(pos, j);
        let p3 = P3::at(pos, k);
        let p4 = P3::at(pos, l);
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
        let de_dphi = self.theta_deriv(cos_phi, sin_phi);
        let sin_term = de_dphi
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

#[derive(Clone, Copy, Debug)]
struct InversionTerm {
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    force_constant: f64,
    c0: f64,
    c1: f64,
    c2: f64,
}

impl InversionTerm {
    fn energy(&self, pos: &[f64]) -> f64 {
        let p1 = P3::at(pos, self.i);
        let p2 = P3::at(pos, self.j);
        let p3 = P3::at(pos, self.k);
        let p4 = P3::at(pos, self.l);
        // UFF::Utils::calculateCosY
        const ZERO_TOL: f64 = 1.0e-16;
        let r_ji = p1.sub(p2);
        let r_jk = p3.sub(p2);
        let r_jl = p4.sub(p2);
        let l2_ji = r_ji.dot(r_ji);
        let l2_jk = r_jk.dot(r_jk);
        let l2_jl = r_jl.dot(r_jl);
        let cos_y = if l2_ji < ZERO_TOL || l2_jk < ZERO_TOL || l2_jl < ZERO_TOL {
            0.0
        } else {
            let n = r_ji.cross(r_jk).div(l2_ji.sqrt() * l2_jk.sqrt());
            let l2n = n.dot(n);
            if l2n < ZERO_TOL {
                0.0
            } else {
                n.dot(r_jl) / (l2_jl.sqrt() * l2n.sqrt())
            }
        };
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
        let cos2w = 2.0 * sin_y * sin_y - 1.0;
        self.force_constant * (self.c0 + self.c1 * sin_y + self.c2 * cos2w)
    }

    fn grad(&self, pos: &[f64], grad: &mut [f64]) {
        let p1 = P3::at(pos, self.i);
        let p2 = P3::at(pos, self.j);
        let p3 = P3::at(pos, self.k);
        let p4 = P3::at(pos, self.l);
        let mut r_ji = p1.sub(p2);
        let mut r_jk = p3.sub(p2);
        let mut r_jl = p4.sub(p2);
        let d_ji = r_ji.length();
        let d_jk = r_jk.length();
        let d_jl = r_jl.length();
        if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
            return;
        }
        r_ji = r_ji.div(d_ji);
        r_jk = r_jk.div(d_jk);
        r_jl = r_jl.div(d_jl);
        let mut n = r_ji.neg().cross(r_jk);
        n = n.div(n.length());
        let cos_y = clip_to_one(n.dot(r_jl));
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = smax(sin_y_sq.sqrt(), 1.0e-8);
        let cos_theta = clip_to_one(r_ji.dot(r_jk));
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = smax(sin_theta_sq.sqrt(), 1.0e-8);
        let de_dw = -self.force_constant * (self.c1 * cos_y - 4.0 * self.c2 * cos_y * sin_y);
        let t1 = r_jl.cross(r_jk);
        let t2 = r_ji.cross(r_jl);
        let t3 = r_jk.cross(r_ji);
        let term1 = sin_y * sin_theta;
        let term2 = cos_y / (sin_y * sin_theta_sq);
        let tg1 = [
            (t1.x / term1 - (r_ji.x - r_jk.x * cos_theta) * term2) / d_ji,
            (t1.y / term1 - (r_ji.y - r_jk.y * cos_theta) * term2) / d_ji,
            (t1.z / term1 - (r_ji.z - r_jk.z * cos_theta) * term2) / d_ji,
        ];
        let tg3 = [
            (t2.x / term1 - (r_jk.x - r_ji.x * cos_theta) * term2) / d_jk,
            (t2.y / term1 - (r_jk.y - r_ji.y * cos_theta) * term2) / d_jk,
            (t2.z / term1 - (r_jk.z - r_ji.z * cos_theta) * term2) / d_jk,
        ];
        let tg4 = [
            (t3.x / term1 - r_jl.x * cos_y / sin_y) / d_jl,
            (t3.y / term1 - r_jl.y * cos_y / sin_y) / d_jl,
            (t3.z / term1 - r_jl.z * cos_y / sin_y) / d_jl,
        ];
        for d in 0..3 {
            grad[3 * self.i + d] += de_dw * tg1[d];
            grad[3 * self.j + d] += -de_dw * (tg1[d] + tg3[d] + tg4[d]);
            grad[3 * self.k + d] += de_dw * tg3[d];
            grad[3 * self.l + d] += de_dw * tg4[d];
        }
    }
}

/// RDKit's UFF force field for one molecule (all atoms, RDKit order).
#[derive(Clone, Debug)]
pub struct RdkitUffField {
    n: usize,
    /// (i, j, r0, kb)
    bonds: Vec<(usize, usize, f64, f64)>,
    angles: Vec<AngleTerm>,
    /// (i, j, x_ij, well_depth, thresh)
    vdw: Vec<(usize, usize, f64, f64, f64)>,
    torsions: Vec<TorsionTerm>,
    inversions: Vec<InversionTerm>,
}

impl RdkitUffField {
    /// `UFFGetMoleculeForceField(mol, vdwThresh, confId,
    /// ignoreInterfragInteractions)` with the conformer `coords` (all atoms
    /// of `mol`, explicit hydrogens in RDKit `AddHs` order). Atoms without
    /// UFF parameters are left out of the terms, as RDKit does.
    pub fn new(
        mol: &Molecule,
        coords: &[[f64; 3]],
        vdw_thresh: f64,
        ignore_interfrag: bool,
    ) -> Self {
        let input = mol;
        let view = rdkit_view(input);
        let mol = &view;
        let n = mol.atom_count();
        let hyb = hybridization(mol);
        let conj = conjugated_atoms(mol);
        let params: Vec<Option<AtomicParams>> = (0..n)
            .map(|i| params_for(&atom_label(mol, i, hyb[i], conj[i])))
            .collect();
        let z: Vec<usize> = (0..n).map(|i| atomic_number(mol, i)).collect();

        let order = rdkit_bond_order_with_added_hs(input);
        // RDKit bond index of each chematic bond.
        let mut rdkit_bond_idx = vec![0usize; mol.bond_count()];
        for (pos, b) in order.iter().enumerate() {
            rdkit_bond_idx[b.0 as usize] = pos;
        }
        let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
        for &b in &order {
            let bond = mol.bond(b);
            let (a1, a2) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            adj[a1].push((a2, b.0 as usize));
            adj[a2].push((a1, b.0 as usize));
        }
        let bo = |b: usize| bond_as_double(mol.bond(BondIdx(b as u32)).order);
        let mut field = RdkitUffField {
            n,
            bonds: Vec::new(),
            angles: Vec::new(),
            vdw: Vec::new(),
            torsions: Vec::new(),
            inversions: Vec::new(),
        };

        // addBonds: RDKit's begin atom first (the force constant is not
        // symmetric in floating point). A SMILES ring-closure bond begins at
        // its closing atom in RDKit.
        for &b in &order {
            let bond = mol.bond(b);
            let (i, j) = if input.is_smiles_ring_closure(b) {
                (bond.atom2.0 as usize, bond.atom1.0 as usize)
            } else {
                (bond.atom1.0 as usize, bond.atom2.0 as usize)
            };
            if let (Some(pi), Some(pj)) = (&params[i], &params[j]) {
                let r0 = calc_bond_rest_length(bond_as_double(bond.order), pi, pj);
                let kb = calc_bond_force_constant(r0, pi, pj);
                field.bonds.push((i, j, r0, kb));
            }
        }

        // addAngles
        let rings = chematic_perception::find_symmetrized_sssr(mol);
        let mut ring_sizes: Vec<Vec<usize>> = vec![Vec::new(); n];
        for ring in rings.rings() {
            for a in ring {
                ring_sizes[a.0 as usize].push(ring.len());
            }
        }
        let in_ring = |a: usize, size: usize| ring_sizes[a].contains(&size);
        for j in 0..n {
            let Some(pj) = &params[j] else {
                continue;
            };
            if adj[j].len() == 1 {
                continue;
            }
            for (a, &(i, bi)) in adj[j].iter().enumerate() {
                let Some(pi) = &params[i] else {
                    continue;
                };
                for &(k, bk) in &adj[j][a + 1..] {
                    let Some(pk) = &params[k] else {
                        continue;
                    };
                    if hyb[j] == Some(4) && adj[j].len() == 5 {
                        continue;
                    }
                    let order = match hyb[j] {
                        Some(1) => 1,
                        Some(2) => {
                            let mut o = 3;
                            if in_ring(j, 3) {
                                if in_ring(i, 3) != in_ring(k, 3) {
                                    o = 30;
                                } else if in_ring(i, 3) && in_ring(k, 3) {
                                    o = 35;
                                }
                            } else if in_ring(j, 4) {
                                if in_ring(i, 4) != in_ring(k, 4) {
                                    o = 40;
                                } else if in_ring(i, 4) && in_ring(k, 4) {
                                    o = 45;
                                }
                            }
                            o
                        }
                        Some(5) => 4,
                        _ => 0,
                    };
                    field
                        .angles
                        .push(AngleTerm::new(i, j, k, bo(bi), bo(bk), pi, pj, pk, order));
                }
            }
        }
        // addAngleSpecialCases: trigonal-bipyramidal centres.
        for a in 0..n {
            if hyb[a] == Some(4) && adj[a].len() == 5 {
                field.add_trigonal_bipyramid_angles(a, &adj, &rdkit_bond_idx, coords, &params, &bo);
            }
        }

        // addNonbonded
        let topo = topological_distances(n, &adj);
        let frag = fragments(n, &adj);
        for i in 0..n {
            let Some(pi) = &params[i] else {
                continue;
            };
            for j in (i + 1)..n {
                let Some(pj) = &params[j] else {
                    continue;
                };
                if ignore_interfrag && frag[i] != frag[j] {
                    continue;
                }
                let rel = topo[i * n + j];
                if rel == 1 || rel == 2 {
                    continue;
                }
                let (ci, cj) = (coords[i], coords[j]);
                let dist = P3 {
                    x: ci[0] - cj[0],
                    y: ci[1] - cj[1],
                    z: ci[2] - cj[2],
                }
                .length();
                let x_ij = (pi.x1 * pj.x1).sqrt();
                if dist < vdw_thresh * x_ij {
                    let well = (pi.d1 * pj.d1).sqrt();
                    field.vdw.push((i, j, x_ij, well, 10.0 * x_ij));
                }
            }
        }

        // addTorsions: `[!$(*#*)&!D1]~[!$(*#*)&!D1]` matches in RDKit order.
        let has_triple = |a: usize| {
            adj[a]
                .iter()
                .any(|&(_, b)| mol.bond(BondIdx(b as u32)).order == BondOrder::Triple)
        };
        let torsion_end = |a: usize| adj[a].len() > 1 && !has_triple(a);
        let sp2_or_sp3 = |a: usize| matches!(hyb[a], Some(2) | Some(3));
        let mut seen = std::collections::HashSet::new();
        for j in 0..n {
            if !torsion_end(j) {
                continue;
            }
            for &(k, bond) in &adj[j] {
                if !torsion_end(k) || !seen.insert((j.min(k), j.max(k))) {
                    continue;
                }
                let (Some(pj), Some(pk)) = (&params[j], &params[k]) else {
                    continue;
                };
                if !(sp2_or_sp3(j) && sp2_or_sp3(k)) {
                    continue;
                }
                let start = field.torsions.len();
                let bo23 = bo(bond);
                for &(i, b1) in &adj[j] {
                    if b1 == bond {
                        continue;
                    }
                    for &(l, b2) in &adj[k] {
                        if b2 == bond || b2 == b1 || l == i {
                            continue;
                        }
                        let end_sp2 = hyb[i] == Some(2) || hyb[l] == Some(2);
                        let (fc, ord, cos_term) =
                            torsion_params(bo23, z[j], z[k], hyb[j], hyb[k], pj, pk, end_sp2);
                        field.torsions.push(TorsionTerm {
                            i,
                            j,
                            k,
                            l,
                            force_constant: fc,
                            order: ord,
                            cos_term,
                        });
                    }
                }
                let count = (field.torsions.len() - start) as f64;
                for t in &mut field.torsions[start..] {
                    t.force_constant /= count;
                }
            }
        }

        // addInversions
        for j in 0..n {
            let at2 = z[j];
            if !matches!(at2, 6 | 7 | 8 | 15 | 33 | 51 | 83) || adj[j].len() != 3 {
                continue;
            }
            if matches!(at2, 6..=8) && hyb[j] != Some(2) {
                continue;
            }
            let nb = [adj[j][0].0, adj[j][1].0, adj[j][2].0];
            let bound_to_sp2_o = at2 == 6 && nb.iter().any(|&a| z[a] == 8 && hyb[a] == Some(2));
            let idx = [nb[0], j, nb[1], nb[2]];
            let (k_inv, c0, c1, c2) = inversion_coefficients(at2, bound_to_sp2_o);
            for perm in [[0, 1, 2, 3], [0, 1, 3, 2], [2, 1, 3, 0]] {
                field.inversions.push(InversionTerm {
                    i: idx[perm[0]],
                    j: idx[perm[1]],
                    k: idx[perm[2]],
                    l: idx[perm[3]],
                    force_constant: 1.0 * k_inv,
                    c0,
                    c1,
                    c2,
                });
            }
        }
        field
    }

    /// `Tools::addTrigonalBipyramidAngles` for centre `a`.
    fn add_trigonal_bipyramid_angles(
        &mut self,
        a: usize,
        adj: &[Vec<(usize, usize)>],
        rdkit_bond_idx: &[usize],
        coords: &[[f64; 3]],
        params: &[Option<AtomicParams>],
        bo: &dyn Fn(usize) -> f64,
    ) {
        let Some(pa) = params[a] else {
            return;
        };
        // Bonds of `a` in adjacency order: (bond index, other atom).
        let bonds: Vec<(usize, usize)> = adj[a].iter().map(|&(o, b)| (b, o)).collect();
        let dir = |o: usize| {
            // Point3D::directionVector: (other - self) normalized.
            let p = coords[a];
            let q = coords[o];
            let v = P3 {
                x: q[0] - p[0],
                y: q[1] - p[1],
                z: q[2] - p[2],
            };
            let l = v.length();
            if l == 0.0 { v } else { v.div(l) }
        };
        let mut most_neg = 100.0;
        let mut ax: Option<(usize, usize)> = None;
        for &(b1, o1) in &bonds {
            let v1 = dir(o1);
            for &(b2, o2) in &bonds {
                if rdkit_bond_idx[b2] > rdkit_bond_idx[b1] {
                    let v2 = dir(o2);
                    let dot = v1.dot(v2);
                    if dot < most_neg {
                        most_neg = dot;
                        ax = Some((b1, b2));
                    }
                }
            }
        }
        let Some((ax1, ax2)) = ax else {
            return;
        };
        let eq: Vec<usize> = bonds
            .iter()
            .map(|&(b, _)| b)
            .filter(|&b| b != ax1 && b != ax2)
            .collect();
        if eq.len() < 3 {
            return;
        }
        let other = |b: usize| bonds.iter().find(|&&(x, _)| x == b).unwrap().1;
        let pairs: [(usize, usize, u32); 10] = [
            (ax1, ax2, 2),
            (eq[0], eq[1], 3),
            (eq[0], eq[2], 3),
            (eq[1], eq[2], 3),
            (ax1, eq[0], 0),
            (ax1, eq[1], 0),
            (ax1, eq[2], 0),
            (ax2, eq[0], 0),
            (ax2, eq[1], 0),
            (ax2, eq[2], 0),
        ];
        for (b1, b2, order) in pairs {
            let (i, j) = (other(b1), other(b2));
            if let (Some(pi), Some(pj)) = (&params[i], &params[j]) {
                self.angles
                    .push(AngleTerm::new(i, a, j, bo(b1), bo(b2), pi, &pa, pj, order));
            }
        }
    }

    /// Diagnostic: every term as `(kind, atoms, parameters)` in RDKit's
    /// contribution order (`bond` r0 kb; `angle` theta0 ka order; `vdw` x_ij
    /// D_ij; `torsion` V order cosTerm; `inversion` K C0 C1 C2).
    #[doc(hidden)]
    pub fn debug_terms(&self) -> Vec<(&'static str, Vec<usize>, Vec<f64>)> {
        let mut out = Vec::new();
        for &(i, j, r0, kb) in &self.bonds {
            out.push(("bond", vec![i, j], vec![r0, kb]));
        }
        for a in &self.angles {
            out.push((
                "angle",
                vec![a.i, a.j, a.k],
                vec![a.theta0, a.ka, a.order as f64],
            ));
        }
        for &(i, j, x, d, _) in &self.vdw {
            out.push(("vdw", vec![i, j], vec![x, d]));
        }
        for t in &self.torsions {
            out.push((
                "torsion",
                vec![t.i, t.j, t.k, t.l],
                vec![t.force_constant, t.order as f64, t.cos_term],
            ));
        }
        for t in &self.inversions {
            out.push((
                "inversion",
                vec![t.i, t.j, t.k, t.l],
                vec![t.force_constant, t.c0, t.c1, t.c2],
            ));
        }
        out
    }

    /// Number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.n
    }

    /// `ForceField::calcEnergy(pos)` on flat coordinates.
    pub fn energy(&self, pos: &[f64]) -> f64 {
        self.energy_with_distances(pos, pos)
    }

    /// The energy at `pos` with every interatomic distance RDKit takes from
    /// its `ForceField::distance` cache taken at `dpos` instead (see
    /// [`crate::rdkit_bfgs::BfgsOutcome::stale_distance_point`]).
    pub fn energy_with_distances(&self, pos: &[f64], dpos: &[f64]) -> f64 {
        let mut res = 0.0;
        for &(i, j, r0, kb) in &self.bonds {
            let dist_term = distance(dpos, i, j) - r0;
            res += 0.5 * kb * dist_term * dist_term;
        }
        for a in &self.angles {
            res += a.energy(pos, dpos);
        }
        for &(i, j, x_ij, well, thresh) in &self.vdw {
            let dist = distance(dpos, i, j);
            if dist > thresh || dist <= 0.0 {
                continue;
            }
            let r = x_ij / dist;
            let r6 = int_pow(r, 6);
            let r12 = r6 * r6;
            res += well * (r12 - 2.0 * r6);
        }
        for t in &self.torsions {
            res += t.energy(pos);
        }
        for inv in &self.inversions {
            res += inv.energy(pos);
        }
        res
    }

    /// `ForceField::calcGrad(pos, grad)`: adds the gradient into `grad`.
    pub fn gradient(&self, pos: &[f64], grad: &mut [f64]) {
        for &(i, j, r0, kb) in &self.bonds {
            let dist = distance(pos, i, j);
            let pre_factor = kb * (dist - r0);
            for d in 0..3 {
                let dgrad = if dist > 0.0 {
                    pre_factor * (pos[3 * i + d] - pos[3 * j + d]) / dist
                } else {
                    kb * 0.01
                };
                grad[3 * i + d] += dgrad;
                grad[3 * j + d] -= dgrad;
            }
        }
        for a in &self.angles {
            a.grad(pos, grad);
        }
        for &(i, j, x_ij, well, thresh) in &self.vdw {
            let dist = distance(pos, i, j);
            if dist > thresh {
                continue;
            }
            if dist <= 0.0 {
                for d in 0..3 {
                    grad[3 * i + d] += 100.0;
                    grad[3 * j + d] -= 100.0;
                }
                continue;
            }
            let r = x_ij / dist;
            let r7 = int_pow(r, 7);
            let r13 = int_pow(r, 13);
            let pre_factor = 12. * well / x_ij * (r7 - r13);
            for d in 0..3 {
                let dgrad = pre_factor * (pos[3 * i + d] - pos[3 * j + d]) / dist;
                grad[3 * i + d] += dgrad;
                grad[3 * j + d] -= dgrad;
            }
        }
        for t in &self.torsions {
            t.grad(pos, grad);
        }
        for inv in &self.inversions {
            inv.grad(pos, grad);
        }
    }

    fn is_empty(&self) -> bool {
        self.bonds.is_empty()
            && self.angles.is_empty()
            && self.vdw.is_empty()
            && self.torsions.is_empty()
            && self.inversions.is_empty()
    }

    /// `ForceField::minimize(maxIts, forceTol)` on flat coordinates (updated
    /// in place). Returns RDKit's status (0 converged, 1 more iterations
    /// needed, -1 line-search failure).
    pub fn minimize(&self, pos: &mut [f64], max_its: u32, force_tol: f64) -> i32 {
        if self.is_empty() {
            return 0;
        }
        crate::rdkit_bfgs::bfgs_minimize(
            pos,
            force_tol,
            max_its,
            |p| self.energy(p),
            |p, g| {
                for x in g.iter_mut() {
                    *x = 0.0;
                }
                self.gradient(p, g);
                crate::rdkit_bfgs::scale_gradient(g)
            },
        )
    }

    /// RDKit `UFFOptimizeMolecule(mol, maxIters)` (forceTol 1e-4) on flat
    /// coordinates: (status, final energy).
    pub fn optimize(&self, pos: &mut [f64], max_iters: u32) -> (i32, f64) {
        if self.is_empty() {
            return (0, 0.0);
        }
        let out = crate::rdkit_bfgs::bfgs_minimize_detailed(
            pos,
            1e-4,
            max_iters,
            |p| self.energy(p),
            |p, g| {
                for x in g.iter_mut() {
                    *x = 0.0;
                }
                self.gradient(p, g);
                crate::rdkit_bfgs::scale_gradient(g)
            },
        );
        let energy = match &out.stale_distance_point {
            Some(d) => self.energy_with_distances(pos, d),
            None => self.energy(pos),
        };
        (out.status, energy)
    }
}

/// `TorsionAngleContrib::calcTorsionParams`: (force constant, order,
/// cosTerm).
#[allow(clippy::too_many_arguments)]
fn torsion_params(
    bo23: f64,
    at2: usize,
    at3: usize,
    hyb2: Option<u8>,
    hyb3: Option<u8>,
    p2: &AtomicParams,
    p3: &AtomicParams,
    end_atom_is_sp2: bool,
) -> (f64, u32, f64) {
    if hyb2 == Some(3) && hyb3 == Some(3) {
        let mut fc = (p2.v1 * p3.v1).sqrt();
        let mut order = 3;
        if bo23 == 1.0 && is_in_group6(at2) && is_in_group6(at3) {
            let v2: f64 = if at2 == 8 { 2.0 } else { 6.8 };
            let v3 = if at3 == 8 { 2.0 } else { 6.8 };
            fc = (v2 * v3).sqrt();
            order = 2;
        }
        (fc, order, -1.0)
    } else if hyb2 == Some(2) && hyb3 == Some(2) {
        (equation17(bo23, p2, p3), 2, 1.0)
    } else {
        let mut res = (1.0, 6, 1.0);
        if bo23 == 1.0 {
            if (hyb2 == Some(3) && is_in_group6(at2) && !is_in_group6(at3))
                || (hyb3 == Some(3) && is_in_group6(at3) && !is_in_group6(at2))
            {
                res = (equation17(bo23, p2, p3), 2, -1.0);
            } else if end_atom_is_sp2 {
                res = (2.0, 3, -1.0);
            }
        }
        res
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// RDKit 2026.03.1 `UFFGetMoleculeForceField(vdwThresh=10)`:
    /// `CalcEnergy()` and `Minimize(maxIts=200)` on paracetamol.
    #[test]
    fn uff_matches_rdkit_bitwise() {
        use crate::rdkit_mmff::tests::{COORDS, paracetamol};
        let mol = paracetamol();
        assert!(rdkit_uff_has_all_params(&mol));
        let ff = RdkitUffField::new(&mol, &COORDS, 10.0, true);
        let mut pos = COORDS.as_flattened().to_vec();
        assert_eq!(ff.energy(&pos), 40.395092089572906);
        let (status, e) = ff.optimize(&mut pos, 200);
        assert_eq!((status, e), (0, 19.492409722746174));
        assert_eq!(pos[0], 3.764669178348734);
    }

    #[test]
    fn int_pow_matches_rdkit_recursion() {
        let x = 1.2345678901234567_f64;
        assert_eq!(int_pow(x, 5), (x * x) * (x * x) * x);
        assert_eq!(int_pow(x, 6), ((x * x) * x) * ((x * x) * x));
        let six = ((x * x) * x) * ((x * x) * x);
        assert_eq!(int_pow(x, 13), six * six * x);
    }

    #[test]
    fn param_table_is_sorted_and_has_common_types() {
        assert!(UFF_PARAMS.windows(2).all(|w| w[0].0 < w[1].0));
        for l in [
            "C_3", "C_R", "C_2", "C_1", "N_3", "N_R", "O_3", "O_R", "H_", "S_3+2", "Cl",
        ] {
            assert!(params_for(l).is_some(), "{l}");
        }
    }
}
