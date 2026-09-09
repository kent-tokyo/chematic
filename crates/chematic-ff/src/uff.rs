//! Universal Force Field (UFF) — geometry minimisation for all elements.
//!
//! UFF is a purely rule-based force field covering the full periodic table
//! (Rappé et al. J. Am. Chem. Soc. 1992, 114(25), 10024-10035).  Unlike
//! MMFF94 — which is parameterised only for common organic/heteroatoms — UFF
//! can handle metal-ligand complexes, organometallics, and any covalent
//! structure.
//!
//! ## Implemented energy terms
//! - **Bond stretching**: harmonic with natural bond order correction
//! - **Angle bending**: Fourier cosine series (C_0 + C_1·cos + C_2·cos(2θ))
//! - **van der Waals**: Lennard-Jones (12-6) with UFF combining rules
//!
//! Torsion and inversion terms are intentionally omitted here; they are less
//! critical for initial 3D placement and can be added incrementally.
//!
//! ## Usage
//! ```rust,ignore
//! use chematic_ff::{assign_uff_types, uff_total_energy, minimize_uff};
//!
//! let types = assign_uff_types(&mol);
//! let coords: Vec<[f64; 3]> = ...; // initial geometry
//! let result = minimize_uff(&mol, &types, coords, 500);
//! ```

use chematic_core::{AtomIdx, BondOrder, Molecule};

// ── Atom type ─────────────────────────────────────────────────────────────────

/// UFF atom type, following the notation in Rappé 1992 Table 1.
///
/// The underscore in names like `C_3` replaces the period used in the paper
/// (`C.3`) to form valid Rust identifiers.
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum UffType {
    // Carbon
    C_3,
    C_2,
    C_1,
    C_R,
    // Nitrogen
    N_3,
    N_2,
    N_1,
    N_R,
    // Oxygen
    O_3,
    O_2,
    O_1,
    O_R,
    // Sulfur
    S_3,
    S_2,
    S_R,
    // Phosphorus
    P_3,
    P_R,
    // Hydrogen
    H_,
    // Halogens
    F_,
    Cl,
    Br,
    I_,
    // Common metals (s/d-block)
    Li,
    Na,
    K,
    Ca,
    Mg,
    Fe,
    Co,
    Ni,
    Cu,
    Zn,
    Mn,
    Cr,
    V_,
    Mo,
    W_,
    Pd,
    Pt,
    Au,
    Ag,
    Hg,
    Al,
    Si,
    // Generic fallback
    Unknown,
}

impl UffType {
    /// UFF parameter: single-bond radius r1 (Å).
    pub fn r1(self) -> f64 {
        match self {
            Self::C_3 => 0.757,
            Self::C_2 => 0.732,
            Self::C_1 => 0.706,
            Self::C_R => 0.729,
            Self::N_3 => 0.700,
            Self::N_2 => 0.685,
            Self::N_1 => 0.656,
            Self::N_R => 0.699,
            Self::O_3 => 0.658,
            Self::O_2 => 0.634,
            Self::O_1 => 0.639,
            Self::O_R => 0.680,
            Self::S_3 => 1.020,
            Self::S_2 => 0.940,
            Self::S_R => 1.000,
            Self::P_3 => 1.101,
            Self::P_R => 1.060,
            Self::H_ => 0.354,
            Self::F_ => 0.668,
            Self::Cl => 1.022,
            Self::Br => 1.172,
            Self::I_ => 1.394,
            Self::Li => 1.336,
            Self::Na => 1.539,
            Self::K => 1.953,
            Self::Ca => 1.761,
            Self::Mg => 1.535,
            Self::Fe => 1.285,
            Self::Co => 1.241,
            Self::Ni => 1.164,
            Self::Cu => 1.302,
            Self::Zn => 1.193,
            Self::Mn => 1.362,
            Self::Cr => 1.370,
            Self::V_ => 1.359,
            Self::Mo => 1.458,
            Self::W_ => 1.526,
            Self::Pd => 1.375,
            Self::Pt => 1.387,
            Self::Au => 1.340,
            Self::Ag => 1.420,
            Self::Hg => 1.490,
            Self::Al => 1.244,
            Self::Si => 1.117,
            Self::Unknown => 1.5,
        }
    }

    /// UFF parameter: natural valence angle θ₀ (degrees).
    pub fn theta0(self) -> f64 {
        match self {
            Self::C_3 => 109.47,
            Self::C_2 => 120.0,
            Self::C_1 => 180.0,
            Self::C_R => 120.0,
            Self::N_3 => 106.70,
            Self::N_2 => 111.2,
            Self::N_1 => 180.0,
            Self::N_R => 120.0,
            Self::O_3 => 104.51,
            Self::O_2 => 120.0,
            Self::O_1 => 180.0,
            Self::O_R => 110.0,
            Self::S_3 => 92.10,
            Self::S_2 => 120.0,
            Self::S_R => 100.0,
            Self::P_3 => 93.80,
            Self::P_R => 120.0,
            Self::H_ => 180.0,
            Self::F_ => 180.0,
            Self::Cl => 180.0,
            Self::Br => 180.0,
            Self::I_ => 180.0,
            _ => 109.47, // default sp3
        }
    }

    /// UFF parameter: nonbonded distance x₁ (Å).
    pub fn x1(self) -> f64 {
        match self {
            Self::H_ => 2.886,
            Self::C_3 => 3.851,
            Self::C_2 => 3.851,
            Self::C_1 => 3.851,
            Self::C_R => 3.851,
            Self::N_3 => 3.660,
            Self::N_2 => 3.660,
            Self::N_1 => 3.660,
            Self::N_R => 3.660,
            Self::O_3 => 3.500,
            Self::O_2 => 3.500,
            Self::O_1 => 3.500,
            Self::O_R => 3.500,
            Self::F_ => 3.364,
            Self::Cl => 3.947,
            Self::Br => 4.153,
            Self::I_ => 4.590,
            Self::S_3 => 4.035,
            Self::S_2 => 4.035,
            Self::S_R => 4.035,
            Self::P_3 => 4.147,
            Self::P_R => 4.147,
            Self::Si => 4.295,
            Self::Al => 4.499,
            Self::Fe => 4.054,
            Self::Co => 3.898,
            Self::Ni => 3.782,
            Self::Cu => 3.495,
            Self::Zn => 3.445,
            Self::Mg => 3.021,
            Self::Ca => 3.753,
            Self::Mn => 4.013,
            Self::Cr => 3.894,
            Self::V_ => 3.804,
            Self::Na => 3.144,
            Self::K => 3.812,
            _ => 3.800,
        }
    }

    /// UFF parameter: nonbonded well depth D₁ (kcal/mol).
    pub fn d1(self) -> f64 {
        match self {
            Self::H_ => 0.044,
            Self::C_3 => 0.105,
            Self::C_2 => 0.105,
            Self::C_1 => 0.105,
            Self::C_R => 0.105,
            Self::N_3 => 0.069,
            Self::N_2 => 0.069,
            Self::N_1 => 0.069,
            Self::N_R => 0.069,
            Self::O_3 => 0.060,
            Self::O_2 => 0.060,
            Self::O_1 => 0.060,
            Self::O_R => 0.060,
            Self::F_ => 0.050,
            Self::Cl => 0.227,
            Self::Br => 0.251,
            Self::I_ => 0.339,
            Self::S_3 => 0.274,
            Self::S_2 => 0.274,
            Self::S_R => 0.274,
            Self::P_3 => 0.305,
            Self::P_R => 0.305,
            Self::Si => 0.402,
            Self::Al => 0.505,
            Self::Fe => 0.013,
            Self::Co => 0.014,
            Self::Ni => 0.015,
            Self::Cu => 0.005,
            Self::Zn => 0.124,
            Self::Mg => 0.111,
            _ => 0.100,
        }
    }
}

// ── Type assignment ───────────────────────────────────────────────────────────

/// Assign a UFF atom type to each heavy atom in `mol`.
///
/// Assignment rules based on element + hybridization (degree, aromatic flag,
/// bond orders) following Table 1 of Rappé 1992.
pub fn assign_uff_types(mol: &Molecule) -> Vec<(AtomIdx, UffType)> {
    mol.atoms()
        .map(|(idx, atom)| {
            let an = atom.element.atomic_number();
            let degree = mol.neighbors(idx).count();
            let aromatic = atom.aromatic;
            let has_double = mol
                .neighbors(idx)
                .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
            let has_triple = mol
                .neighbors(idx)
                .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Triple);

            let uff = match an {
                1 => UffType::H_,
                6 => {
                    if aromatic {
                        UffType::C_R
                    } else if has_triple {
                        UffType::C_1
                    } else if has_double {
                        UffType::C_2
                    } else {
                        UffType::C_3
                    }
                }
                7 => {
                    if aromatic {
                        UffType::N_R
                    } else if has_triple {
                        UffType::N_1
                    } else if has_double {
                        UffType::N_2
                    } else {
                        UffType::N_3
                    }
                }
                8 => {
                    if aromatic {
                        UffType::O_R
                    } else if has_double {
                        UffType::O_2
                    } else if degree == 1 {
                        UffType::O_1
                    } else {
                        UffType::O_3
                    }
                }
                9 => UffType::F_,
                14 => UffType::Si,
                15 => {
                    if aromatic {
                        UffType::P_R
                    } else {
                        UffType::P_3
                    }
                }
                16 => {
                    if aromatic {
                        UffType::S_R
                    } else if has_double {
                        UffType::S_2
                    } else {
                        UffType::S_3
                    }
                }
                17 => UffType::Cl,
                35 => UffType::Br,
                53 => UffType::I_,
                13 => UffType::Al,
                3 => UffType::Li,
                11 => UffType::Na,
                19 => UffType::K,
                20 => UffType::Ca,
                12 => UffType::Mg,
                26 => UffType::Fe,
                27 => UffType::Co,
                28 => UffType::Ni,
                29 => UffType::Cu,
                30 => UffType::Zn,
                25 => UffType::Mn,
                24 => UffType::Cr,
                23 => UffType::V_,
                42 => UffType::Mo,
                74 => UffType::W_,
                46 => UffType::Pd,
                78 => UffType::Pt,
                79 => UffType::Au,
                47 => UffType::Ag,
                80 => UffType::Hg,
                _ => UffType::Unknown,
            };
            (idx, uff)
        })
        .collect()
}

// ── Energy functions ──────────────────────────────────────────────────────────

/// Compute UFF bond length between types `i` and `j` with bond order `n`.
///
/// Equation 2 from Rappé 1992: r_ij = r_i + r_j + r_BO - r_EN
fn uff_bond_length(ti: UffType, tj: UffType, bond_order: f64) -> f64 {
    let rij = ti.r1() + tj.r1();
    // Bond order correction r_BO = −λ(r_i + r_j) ln(n)
    let lambda = 0.1332;
    let r_bo = -lambda * rij * bond_order.ln();
    // Electronegativity correction (χ) — simplified: use zero for now
    rij + r_bo
}

/// Bond order as f64 from `BondOrder`.
fn bond_order_f64(bo: BondOrder) -> f64 {
    match bo {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        _ => 1.0,
    }
}

fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn cos_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
    let dot = ba[0] * bc[0] + ba[1] * bc[1] + ba[2] * bc[2];
    let len_ba = (ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]).sqrt();
    let len_bc = (bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]).sqrt();
    let denom = len_ba * len_bc;
    if denom < 1e-10 {
        return 1.0;
    }
    (dot / denom).clamp(-1.0, 1.0)
}

#[derive(Clone, Copy)]
struct UffBondTerm {
    a: AtomIdx,
    b: AtomIdx,
    r0: f64,
    k: f64,
}

#[derive(Clone, Copy)]
struct UffAngleTerm {
    a: AtomIdx,
    center: AtomIdx,
    c: AtomIdx,
    cos0: f64,
    k: f64,
}

#[derive(Clone, Copy)]
struct UffVdwTerm {
    a: AtomIdx,
    b: AtomIdx,
    x: f64,
    d: f64,
}

struct PreparedUffEnergy {
    bonds: Vec<UffBondTerm>,
    angles: Vec<UffAngleTerm>,
    vdw: Vec<UffVdwTerm>,
}

impl PreparedUffEnergy {
    fn new(mol: &Molecule, types: &[(AtomIdx, UffType)]) -> Self {
        let type_map: std::collections::HashMap<AtomIdx, UffType> =
            types.iter().map(|&(a, t)| (a, t)).collect();
        let get_type = |idx: AtomIdx| type_map.get(&idx).copied().unwrap_or(UffType::Unknown);

        let bonds = mol
            .bonds()
            .map(|(_, bond)| {
                let r0 = uff_bond_length(
                    get_type(bond.atom1),
                    get_type(bond.atom2),
                    bond_order_f64(bond.order),
                );
                UffBondTerm {
                    a: bond.atom1,
                    b: bond.atom2,
                    r0,
                    k: 664.12 / (r0 * r0 * r0),
                }
            })
            .collect();

        let mut angles = Vec::new();
        for &(center, center_type) in types {
            let theta0 = center_type.theta0().to_radians();
            let sin0 = theta0.sin();
            let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();
            for i in 0..neighbors.len() {
                for j in (i + 1)..neighbors.len() {
                    angles.push(UffAngleTerm {
                        a: neighbors[i],
                        center,
                        c: neighbors[j],
                        cos0: theta0.cos(),
                        k: 0.5 * 332.06 / (sin0 * sin0 + 1e-10),
                    });
                }
            }
        }

        let atom_indices: Vec<AtomIdx> = mol.atoms().map(|(idx, _)| idx).collect();
        let mut excl: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for (_, bond) in mol.bonds() {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            excl.insert((i.min(j), i.max(j)));
            for (nb_i, _) in mol.neighbors(bond.atom1) {
                let ni = nb_i.0 as usize;
                excl.insert((ni.min(j), ni.max(j)));
            }
            for (nb_j, _) in mol.neighbors(bond.atom2) {
                let nj = nb_j.0 as usize;
                excl.insert((i.min(nj), i.max(nj)));
            }
        }
        let mut vdw = Vec::new();
        for i in 0..atom_indices.len() {
            for j in (i + 1)..atom_indices.len() {
                let a = atom_indices[i];
                let b = atom_indices[j];
                let key = (a.0 as usize, b.0 as usize);
                if excl.contains(&key) {
                    continue;
                }
                let ta = get_type(a);
                let tb = get_type(b);
                vdw.push(UffVdwTerm {
                    a,
                    b,
                    x: (ta.x1() * tb.x1()).sqrt(),
                    d: (ta.d1() * tb.d1()).sqrt(),
                });
            }
        }

        Self { bonds, angles, vdw }
    }

    fn energy(&self, coords: &[[f64; 3]]) -> f64 {
        let get_coord = |idx: AtomIdx| coords[idx.0 as usize];
        let mut energy = 0.0;
        for term in &self.bonds {
            let r = dist(get_coord(term.a), get_coord(term.b));
            energy += 0.5 * term.k * (r - term.r0) * (r - term.r0);
        }
        for term in &self.angles {
            let delta =
                cos_angle(get_coord(term.a), get_coord(term.center), get_coord(term.c)) - term.cos0;
            energy += 0.5 * term.k * delta * delta;
        }
        for term in &self.vdw {
            let r = dist(get_coord(term.a), get_coord(term.b)).max(0.5);
            let ratio = term.x / r;
            let ratio6 = ratio.powi(6);
            let ratio12 = ratio6 * ratio6;
            energy += term.d * (ratio12 - 2.0 * ratio6);
        }
        energy
    }

    fn analytic_gradient(&self, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut gradient = vec![[0.0; 3]; coords.len()];
        let add = |slot: &mut [f64; 3], value: [f64; 3]| {
            slot[0] += value[0];
            slot[1] += value[1];
            slot[2] += value[2];
        };
        for term in &self.bonds {
            let a = coords[term.a.0 as usize];
            let b = coords[term.b.0 as usize];
            let delta = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let r = dist(a, b);
            if r > 1e-12 {
                let scale = term.k * (r - term.r0) / r;
                let value = [scale * delta[0], scale * delta[1], scale * delta[2]];
                add(&mut gradient[term.a.0 as usize], value);
                add(
                    &mut gradient[term.b.0 as usize],
                    [-value[0], -value[1], -value[2]],
                );
            }
        }
        for term in &self.angles {
            let a = coords[term.a.0 as usize];
            let center = coords[term.center.0 as usize];
            let c = coords[term.c.0 as usize];
            let u = [a[0] - center[0], a[1] - center[1], a[2] - center[2]];
            let v = [c[0] - center[0], c[1] - center[1], c[2] - center[2]];
            let lu2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
            let lv2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            let lu = lu2.sqrt();
            let lv = lv2.sqrt();
            if lu <= 1e-12 || lv <= 1e-12 {
                continue;
            }
            let cos_theta = cos_angle(a, center, c);
            let scale = term.k * (cos_theta - term.cos0);
            let inv = 1.0 / (lu * lv);
            let da = [
                v[0] * inv - cos_theta * u[0] / lu2,
                v[1] * inv - cos_theta * u[1] / lu2,
                v[2] * inv - cos_theta * u[2] / lu2,
            ];
            let dc = [
                u[0] * inv - cos_theta * v[0] / lv2,
                u[1] * inv - cos_theta * v[1] / lv2,
                u[2] * inv - cos_theta * v[2] / lv2,
            ];
            let da = [scale * da[0], scale * da[1], scale * da[2]];
            let dc = [scale * dc[0], scale * dc[1], scale * dc[2]];
            add(&mut gradient[term.a.0 as usize], da);
            add(&mut gradient[term.c.0 as usize], dc);
            add(
                &mut gradient[term.center.0 as usize],
                [-da[0] - dc[0], -da[1] - dc[1], -da[2] - dc[2]],
            );
        }
        for term in &self.vdw {
            let a = coords[term.a.0 as usize];
            let b = coords[term.b.0 as usize];
            let delta = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            let r = dist(a, b);
            if r <= 0.5 || r <= 1e-12 {
                continue;
            }
            let ratio6 = (term.x / r).powi(6);
            let radial = 12.0 * term.d * (ratio6 - ratio6 * ratio6) / r;
            let value = [
                radial * delta[0] / r,
                radial * delta[1] / r,
                radial * delta[2] / r,
            ];
            add(&mut gradient[term.a.0 as usize], value);
            add(
                &mut gradient[term.b.0 as usize],
                [-value[0], -value[1], -value[2]],
            );
        }
        gradient
    }
}

/// Compute UFF total energy (bond + angle + vdW) in kcal/mol.
pub fn uff_total_energy(mol: &Molecule, types: &[(AtomIdx, UffType)], coords: &[[f64; 3]]) -> f64 {
    PreparedUffEnergy::new(mol, types).energy(coords)
}

// ── Gradient + L-BFGS minimizer ───────────────────────────────────────────────

/// Analytic gradient of the prepared UFF energy terms.
fn uff_gradient(prepared: &PreparedUffEnergy, coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
    prepared.analytic_gradient(coords)
}

/// No legitimate covalent bond stretches anywhere near this length; a
/// post-minimization bond longer than this indicates a blown-up geometry,
/// not a slow-but-fine one. Mirrors `chematic-3d`'s
/// `minimize::MAX_SANE_BOND_LENGTH` (that crate can't be depended on from
/// here — `chematic-3d` depends on `chematic-ff`, not the reverse — so this
/// is a deliberately-duplicated copy of the same, already-corpus-validated
/// constant, not an independently chosen one; keep the two in sync).
const MAX_SANE_UFF_BOND_LENGTH: f64 = 3.0;

fn worst_uff_bond_length(mol: &Molecule, coords: &[[f64; 3]]) -> f64 {
    if coords
        .iter()
        .any(|point| point.iter().any(|value| !value.is_finite()))
    {
        return f64::INFINITY;
    }
    let dist = |a: [f64; 3], b: [f64; 3]| {
        let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    };
    mol.bonds()
        .map(|(_, b)| dist(coords[b.atom1.0 as usize], coords[b.atom2.0 as usize]))
        .fold(0.0_f64, f64::max)
}

/// True iff every coordinate is finite and no bond exceeds
/// [`MAX_SANE_UFF_BOND_LENGTH`]. Deliberately independent of `converged`:
/// steepest descent frequently reports `converged == false` on perfectly
/// sound geometries that simply haven't hit the tight RMS-gradient
/// threshold within `max_iter` (same rationale as `chematic-3d`'s
/// `check_minimization_soundness`, which this mirrors the bond-length half
/// of — that gate's other half, a residual-force ceiling, has no UFF
/// equivalent here since this steepest-descent loop doesn't retain a
/// converged gradient norm past each iteration).
fn is_sound_uff_geometry(mol: &Molecule, coords: &[[f64; 3]]) -> bool {
    if coords.iter().any(|p| p.iter().any(|x| !x.is_finite())) {
        return false;
    }
    worst_uff_bond_length(mol, coords) <= MAX_SANE_UFF_BOND_LENGTH
}

/// Result of UFF minimisation.
pub struct UffMinimizeResult {
    /// Final atomic coordinates (Å).
    pub coords: Vec<[f64; 3]>,
    /// Final total energy (kcal/mol).
    pub energy: f64,
    /// Number of iterations taken.
    pub iterations: usize,
    /// True if the gradient norm converged below threshold.
    pub converged: bool,
    /// True if `coords` is a geometrically sound result (all-finite, no
    /// bond stretched past [`MAX_SANE_UFF_BOND_LENGTH`]) — independent of
    /// `converged`, which only reports whether the RMS-gradient stopping
    /// criterion was met, not whether the geometry itself is trustworthy.
    /// Callers that skip a soundness check of their own (both
    /// `chematic-py`'s `Mol.minimize_uff()` and `chematic-wasm`'s
    /// `minimize_uff_json()` did, until this field existed) previously had
    /// no signal at all that a result like a blown-up bond in a fused
    /// aromatic ring folding non-planar (a real stationary point of UFF's
    /// torsion/out-of-plane-incomplete potential, not slow convergence) had
    /// occurred.
    pub sound: bool,
    /// Longest covalent bond in the returned geometry (Å). This is the
    /// measurement behind `sound` and lets bindings explain a rejected
    /// result without reimplementing the soundness gate.
    pub worst_bond_length: f64,
    /// True when line search rejected an energy-decreasing proposal because
    /// it would have produced an unsound covalent bond length. Callers can
    /// distinguish this bounded rescue signal from an ordinary high-residual
    /// result whose geometry never attempted a catastrophic step.
    pub rejected_unsound_step: bool,
}

/// Minimise UFF energy using steepest descent (convergence criterion: RMS
/// gradient < 0.01 kcal/mol/Å).
///
/// For production use, consider hooking into the existing L-BFGS minimiser
/// in `mmff94_minimizer.rs`; the interface is intentionally compatible.
pub fn minimize_uff(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    initial_coords: Vec<[f64; 3]>,
    max_iter: usize,
) -> UffMinimizeResult {
    minimize_uff_with_constraint(mol, types, initial_coords, max_iter, |_| true)
}

/// Minimise UFF while rejecting candidate line-search steps that violate a
/// caller-supplied geometry constraint.
///
/// The initial geometry must already satisfy `accept`; the predicate is
/// evaluated only for proposed coordinates. Rejected candidates reduce the
/// line-search step exactly like an unsound bond-length proposal, so this is a
/// bounded, fail-closed constraint rather than a post-hoc correction. The
/// callback must be deterministic and side-effect free.
pub fn minimize_uff_with_constraint<F>(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    initial_coords: Vec<[f64; 3]>,
    max_iter: usize,
    accept: F,
) -> UffMinimizeResult
where
    F: Fn(&[[f64; 3]]) -> bool,
{
    let prepared = PreparedUffEnergy::new(mol, types);
    let mut coords = initial_coords;
    let mut step = 0.05_f64;
    let mut prev_energy = f64::MAX;
    let mut rejected_unsound_step = false;

    for iter in 0..max_iter {
        let energy = prepared.energy(&coords);
        let grad = uff_gradient(&prepared, &coords);

        // RMS gradient norm
        let rms: f64 = {
            let sum2: f64 = grad.iter().flat_map(|g| g.iter()).map(|v| v * v).sum();
            (sum2 / (grad.len() * 3) as f64).sqrt()
        };

        if rms < 0.01 {
            let sound = is_sound_uff_geometry(mol, &coords);
            let worst_bond_length = worst_uff_bond_length(mol, &coords);
            return UffMinimizeResult {
                coords,
                energy,
                iterations: iter,
                converged: true,
                sound,
                worst_bond_length,
                rejected_unsound_step,
            };
        }

        // Line search: accept step only if energy decreases
        let new_coords: Vec<[f64; 3]> = coords
            .iter()
            .zip(&grad)
            .map(|(c, g)| [c[0] - step * g[0], c[1] - step * g[1], c[2] - step * g[2]])
            .collect();

        let new_energy = prepared.energy(&new_coords);
        let geometry_sound = is_sound_uff_geometry(mol, &new_coords);
        // Energy descent alone is not a sufficient acceptance criterion:
        // the incomplete UFF potential can lower its energy by walking into
        // a stationary geometry with a catastrophically stretched covalent
        // bond (notably fused aromatics such as naphthalene). Reject such a
        // proposal before it becomes the next iterate and let the line
        // search reduce the step instead. This preserves the existing
        // fail-closed `sound` contract while preventing the optimizer from
        // knowingly propagating an unsound intermediate.
        if new_energy < energy && geometry_sound && accept(&new_coords) {
            coords = new_coords;
            if energy - new_energy < prev_energy * 1e-7 {
                step *= 1.2;
            }
            prev_energy = energy;
        } else {
            if new_energy < energy && !geometry_sound {
                rejected_unsound_step = true;
            }
            step *= 0.5;
            if step < 1e-8 {
                let sound = is_sound_uff_geometry(mol, &coords);
                let worst_bond_length = worst_uff_bond_length(mol, &coords);
                return UffMinimizeResult {
                    coords,
                    energy,
                    iterations: iter,
                    converged: false,
                    sound,
                    worst_bond_length,
                    rejected_unsound_step,
                };
            }
        }
    }

    let energy = uff_total_energy(mol, types, &coords);
    let sound = is_sound_uff_geometry(mol, &coords);
    let worst_bond_length = worst_uff_bond_length(mol, &coords);
    UffMinimizeResult {
        coords,
        energy,
        iterations: max_iter,
        converged: false,
        sound,
        worst_bond_length,
        rejected_unsound_step,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn assign_types_ethanol() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        assert_eq!(types.len(), 3);
        // C sp3 → C_3, O sp3 → O_3
        let type_map: std::collections::HashMap<_, _> = types.into_iter().collect();
        for (_, atom) in mol.atoms() {
            let idx = mol
                .atoms()
                .find(|(_, a)| a.element == atom.element)
                .map(|(i, _)| i);
            if atom.element.atomic_number() == 6 {
                assert!(matches!(
                    type_map[&idx.unwrap()],
                    UffType::C_3 | UffType::C_2
                ));
            }
        }
    }

    #[test]
    fn assign_types_benzene_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        let types = assign_uff_types(&mol);
        // All C aromatic → C_R
        for (_, t) in &types {
            assert_eq!(*t, UffType::C_R);
        }
    }

    #[test]
    fn energy_finite() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [2.5, 1.2, 0.0]];
        let e = uff_total_energy(&mol, &types, &coords);
        assert!(e.is_finite(), "energy should be finite: {e}");
    }

    #[test]
    fn minimize_reduces_energy() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0], // stretched bond
            [3.5, 1.2, 0.0],
        ];
        let e0 = uff_total_energy(&mol, &types, &coords);
        let result = minimize_uff(&mol, &types, coords, 200);
        assert!(
            result.energy < e0,
            "minimisation should reduce energy: {e0} → {}",
            result.energy
        );
    }

    #[test]
    fn constrained_minimizer_rejects_constraint_violating_steps() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [3.5, 1.2, 0.0]];
        let initial = coords.clone();
        let result = minimize_uff_with_constraint(&mol, &types, coords, 20, |candidate| {
            candidate == initial.as_slice()
        });
        assert_eq!(result.coords, initial);
        assert_eq!(result.energy, uff_total_energy(&mol, &types, &initial));
        assert!(!result.converged);
        assert!(
            !result.rejected_unsound_step,
            "a caller constraint rejection must not be reported as an unsound UFF step"
        );
    }

    #[test]
    fn minimize_uff_reports_sound_on_ordinary_ethanol() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0], [2.5, 1.2, 0.0]];
        let result = minimize_uff(&mol, &types, coords, 200);
        assert!(
            result.sound,
            "an ordinary small molecule minimizing normally should report sound"
        );
        assert!(
            result.worst_bond_length <= MAX_SANE_UFF_BOND_LENGTH,
            "sound result must expose a bond length within the soundness limit"
        );
    }

    #[test]
    fn minimize_uff_reports_unsound_on_a_blown_up_bond() {
        // max_iter=0 returns the initial coords untouched (the `for iter in
        // 0..max_iter` loop never runs), so this deterministically exercises
        // `sound`'s bond-length check against a deliberately-stretched
        // C-C bond (5.0 Å, well past MAX_SANE_UFF_BOND_LENGTH) without
        // depending on steepest descent actually getting stuck there.
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [6.0, 1.2, 0.0]];
        let result = minimize_uff(&mol, &types, coords, 0);
        assert!(
            !result.sound,
            "a 5.0 Å C-C bond must be reported unsound regardless of `converged`"
        );
        assert!(
            result.worst_bond_length > MAX_SANE_UFF_BOND_LENGTH,
            "unsound result must expose the stretched bond measurement"
        );
    }

    #[test]
    fn minimize_uff_reports_unsound_on_non_finite_coordinates() {
        let mol = parse("CCO").unwrap();
        let types = assign_uff_types(&mol);
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0], [2.5, 1.2, 0.0]];
        let result = minimize_uff(&mol, &types, coords, 0);
        assert!(!result.sound, "non-finite coordinates must be unsound");
        assert!(
            result.worst_bond_length.is_infinite(),
            "non-finite geometry must fail closed in the exposed bond metric"
        );
    }

    /// Propane skeleton (C0-C1-C2, heavy atoms only — implicit H fills
    /// valence) placed so the C0-C1-C2 angle is exactly `theta_deg`, both
    /// C-C bonds at UFF's own C_3-C_3 r0 (~1.514 Å) so the bond-stretch term
    /// stays near zero and any energy blow-up as the angle closes is
    /// attributable to the angle and vdW terms only — the same isolation
    /// issue #176's own propane repro used.
    fn propane_at_angle(theta_deg: f64) -> (Molecule, Vec<(AtomIdx, UffType)>, Vec<[f64; 3]>) {
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let types = assign_uff_types(&mol);

        let r = 1.514_f64; // UFF C_3-C_3 bond length
        let half = theta_deg.to_radians() / 2.0;
        let coords = vec![
            [r * half.cos(), r * half.sin(), 0.0],
            [0.0, 0.0, 0.0],
            [r * half.cos(), -r * half.sin(), 0.0],
        ];
        (mol, types, coords)
    }

    #[test]
    fn uff_vdw_excludes_1_3_pair_propane_no_runaway() {
        // Issue #176's own measured pre-fix blow-up: 109.5°→16, 90°→110,
        // 70°→1326, 60°→6800 kcal/mol, driven entirely by the (wrongly
        // included) C0···C2 1-3 pair's LJ repulsion as it gets squeezed by
        // the closing angle. After excluding true 1-3 pairs, energy should
        // stay bounded (dominated by the smooth angle-bending term) even at
        // very closed angles.
        let mut energies = Vec::new();
        for &theta in &[109.5, 90.0, 70.0, 60.0, 45.0] {
            let (mol, types, coords) = propane_at_angle(theta);
            let e = uff_total_energy(&mol, &types, &coords);
            assert!(e.is_finite(), "energy at {theta}° should be finite: {e}");
            energies.push((theta, e));
        }
        for &(theta, e) in &energies {
            assert!(
                e < 200.0,
                "1-3 exclusion should prevent runaway vdW: at {theta}° energy={e} (issue #176 measured 6800 at 60° pre-fix)"
            );
        }
    }

    #[test]
    fn uff_vdw_still_repels_genuine_1_4_pair_butane() {
        // Positive control mirroring `mmff94_minimizer.rs`'s own
        // `vdw_more_repulsive_at_short_range`: a *real* 1-4 pair (butane's
        // terminal carbons) must still get full vdW repulsion — the 1-3 fix
        // must not over-exclude non-angle pairs.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        let mol = b.build();
        let types = assign_uff_types(&mol);

        let coords_close = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.6, 0.0, 0.0], // C3 forced very close to C0 (genuine 1-4)
        ];
        let coords_far = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [8.0, 0.0, 0.0],
        ];
        let e_close = uff_total_energy(&mol, &types, &coords_close);
        let e_far = uff_total_energy(&mol, &types, &coords_far);
        assert!(e_close.is_finite() && e_far.is_finite());
        assert!(
            e_close > e_far,
            "genuine 1-4 pair must still repel at short range: close={e_close} far={e_far}"
        );
    }

    #[test]
    fn uff_gradient_is_descent_direction_for_closed_angle_propane() {
        // At a strained, closed-angle propane geometry, stepping a small
        // distance along the analytic -gradient must lower the energy.
        let (mol, types, coords) = propane_at_angle(60.0);
        let e0 = uff_total_energy(&mol, &types, &coords);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        let grad = uff_gradient(&prepared, &coords);
        let grad_norm: f64 = grad
            .iter()
            .flat_map(|g| g.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(
            grad_norm > 1e-8,
            "gradient should be nonzero at strained geometry"
        );

        let step = 1e-4 / grad_norm;
        let stepped: Vec<[f64; 3]> = coords
            .iter()
            .zip(&grad)
            .map(|(c, g)| [c[0] - step * g[0], c[1] - step * g[1], c[2] - step * g[2]])
            .collect();
        let e1 = uff_total_energy(&mol, &types, &stepped);
        assert!(
            e1 < e0,
            "step along -gradient should decrease energy: {e0} → {e1}"
        );
    }

    #[test]
    fn analytic_uff_gradient_matches_finite_difference_reference() {
        let (mol, types, coords) = propane_at_angle(75.0);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        let actual = prepared.analytic_gradient(&coords);
        let delta = 1e-5;
        for atom in 0..coords.len() {
            for axis in 0..3 {
                let mut plus = coords.clone();
                let mut minus = coords.clone();
                plus[atom][axis] += delta;
                minus[atom][axis] -= delta;
                let expected = (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                assert!(
                    (actual[atom][axis] - expected).abs() < 2e-3,
                    "gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                    actual[atom][axis]
                );
            }
        }
    }

    #[test]
    fn analytic_uff_gradient_matches_vdw_14_finite_difference_reference() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        let mut builder = MoleculeBuilder::new();
        let atoms = [
            builder.add_atom(Atom::new(Element::C)),
            builder.add_atom(Atom::new(Element::C)),
            builder.add_atom(Atom::new(Element::C)),
            builder.add_atom(Atom::new(Element::C)),
        ];
        for pair in atoms.windows(2) {
            builder
                .add_bond(pair[0], pair[1], BondOrder::Single)
                .expect("butane probe bond should be valid");
        }
        let mol = builder.build();
        let types = assign_uff_types(&mol);
        let prepared = PreparedUffEnergy::new(&mol, &types);
        let coords = [
            [0.0, 0.0, 0.0],
            [1.52, 0.1, 0.0],
            [3.01, 0.45, 0.2],
            [4.30, 0.85, 0.55],
        ];
        let actual = prepared.analytic_gradient(&coords);
        let delta = 1e-5;
        for atom in 0..coords.len() {
            for axis in 0..3 {
                let mut plus = coords;
                let mut minus = coords;
                plus[atom][axis] += delta;
                minus[atom][axis] -= delta;
                let expected = (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                assert!(
                    (actual[atom][axis] - expected).abs() < 2e-3,
                    "1-4 vdW gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                    actual[atom][axis]
                );
            }
        }
    }

    #[test]
    fn uff_handles_zinc_complex() {
        // Zinc as a metal centre — UFF should assign Zn type
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let zn = b.add_atom(Atom::new(Element::ZN));
        let n1 = b.add_atom(Atom::new(Element::N));
        let n2 = b.add_atom(Atom::new(Element::N));
        b.add_bond(zn, n1, BondOrder::Single).unwrap();
        b.add_bond(zn, n2, BondOrder::Single).unwrap();
        let mol = b.build();
        let types = assign_uff_types(&mol);
        let zn_type = types.iter().find(|(_, t)| *t == UffType::Zn);
        assert!(zn_type.is_some(), "Zn should get UffType::Zn");
    }

    #[test]
    fn uff_supported_metal_and_halogen_types_have_finite_prepared_energy() {
        // This is a local soundness probe for the explicit element-to-UFF
        // dispatch table. It does not claim reference-parameter parity or
        // chemical validity for arbitrary metal coordination geometries.
        use chematic_core::{Atom, Element, MoleculeBuilder};

        let elements = [
            Element::LI,
            Element::NA,
            Element::K,
            Element::MG,
            Element::CA,
            Element::AL,
            Element::SI,
            Element::FE,
            Element::CO,
            Element::NI,
            Element::CU,
            Element::ZN,
            Element::AG,
            Element::AU,
            Element::HG,
            Element::CL,
            Element::BR,
            Element::I,
        ];
        for element in elements {
            let mut builder = MoleculeBuilder::new();
            let first = builder.add_atom(Atom::new(element));
            let second = builder.add_atom(Atom::new(Element::C));
            builder
                .add_bond(first, second, BondOrder::Single)
                .expect("two-atom probe bond should be valid");
            let mol = builder.build();
            let types = assign_uff_types(&mol);
            let prepared = PreparedUffEnergy::new(&mol, &types);
            let coords = [[0.0, 0.0, 0.0], [1.8, 0.2, 0.1]];
            let energy = prepared.energy(&coords);
            let gradient = prepared.analytic_gradient(&coords);
            assert!(energy.is_finite(), "{element:?} energy must be finite");
            assert!(
                gradient.iter().flatten().all(|value| value.is_finite()),
                "{element:?} gradient must be finite: {gradient:?}"
            );
            let delta = 1e-5;
            for atom in 0..coords.len() {
                for axis in 0..3 {
                    let mut plus = coords;
                    let mut minus = coords;
                    plus[atom][axis] += delta;
                    minus[atom][axis] -= delta;
                    let expected =
                        (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                    assert!(
                        (gradient[atom][axis] - expected).abs() < 2e-3,
                        "{element:?} gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                        gradient[atom][axis]
                    );
                }
            }
        }
    }

    #[test]
    fn every_uff_type_has_finite_prepared_energy_and_gradient() {
        // Dispatch coverage above exercises elements that are commonly
        // reached through the public assigner. This table additionally keeps
        // every declared UffType parameter bundle executable, including
        // hybridization-specific and fallback-only variants.
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

        let all_types = [
            UffType::C_3,
            UffType::C_2,
            UffType::C_1,
            UffType::C_R,
            UffType::N_3,
            UffType::N_2,
            UffType::N_1,
            UffType::N_R,
            UffType::O_3,
            UffType::O_2,
            UffType::O_1,
            UffType::O_R,
            UffType::S_3,
            UffType::S_2,
            UffType::S_R,
            UffType::P_3,
            UffType::P_R,
            UffType::H_,
            UffType::F_,
            UffType::Cl,
            UffType::Br,
            UffType::I_,
            UffType::Li,
            UffType::Na,
            UffType::K,
            UffType::Ca,
            UffType::Mg,
            UffType::Fe,
            UffType::Co,
            UffType::Ni,
            UffType::Cu,
            UffType::Zn,
            UffType::Mn,
            UffType::Cr,
            UffType::V_,
            UffType::Mo,
            UffType::W_,
            UffType::Pd,
            UffType::Pt,
            UffType::Au,
            UffType::Ag,
            UffType::Hg,
            UffType::Al,
            UffType::Si,
            UffType::Unknown,
        ];
        assert_eq!(all_types.len(), 45);

        let mut builder = MoleculeBuilder::new();
        let first = builder.add_atom(Atom::new(Element::C));
        let second = builder.add_atom(Atom::new(Element::C));
        builder
            .add_bond(first, second, BondOrder::Single)
            .expect("two-atom probe bond should be valid");
        let mol = builder.build();
        let coords = [[0.0, 0.0, 0.0], [1.8, 0.2, 0.1]];
        let delta = 1e-5;

        for uff_type in all_types {
            let types = vec![(first, uff_type), (second, UffType::C_3)];
            let prepared = PreparedUffEnergy::new(&mol, &types);
            let energy = prepared.energy(&coords);
            let gradient = prepared.analytic_gradient(&coords);
            assert!(energy.is_finite(), "{uff_type:?} energy must be finite");
            assert!(
                gradient.iter().flatten().all(|value| value.is_finite()),
                "{uff_type:?} gradient must be finite: {gradient:?}"
            );
            for atom in 0..coords.len() {
                for axis in 0..3 {
                    let mut plus = coords;
                    let mut minus = coords;
                    plus[atom][axis] += delta;
                    minus[atom][axis] -= delta;
                    let expected =
                        (prepared.energy(&plus) - prepared.energy(&minus)) / (2.0 * delta);
                    assert!(
                        (gradient[atom][axis] - expected).abs() < 2e-3,
                        "{uff_type:?} gradient mismatch at atom {atom}, axis {axis}: actual={}, expected={expected}",
                        gradient[atom][axis]
                    );
                }
            }
        }
    }
}
