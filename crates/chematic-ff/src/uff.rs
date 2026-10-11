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
//! - **Aromatic inversion**: a bounded planar-restraint term for trigonal
//!   aromatic centers, preventing fused aromatic systems from settling into
//!   non-planar stationary points of the incomplete potential
//! - **Torsion**: bounded periodic terms for the common sp3/sp2 central-bond
//!   classes, with a fail-closed zero contribution for unsupported classes
//!
//! Higher-coordinate metal torsion parameter parity remains intentionally
//! deferred; the implemented organic classes are an incremental UFF slice,
//! not a claim of complete RDKit energy parity.
//!
//! ## Usage
//! ```rust,ignore
//! use chematic_ff::{assign_uff_types, uff_total_energy, minimize_uff};
//!
//! let types = assign_uff_types(&mol);
//! let coords: Vec<[f64; 3]> = ...; // initial geometry
//! let result = minimize_uff(&mol, &types, coords, 500);
//! ```

#[allow(unused_imports)]
use crate::fmath::DetMath;
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

    /// UFF torsional barrier parameter V1 (kcal/mol).
    ///
    /// These are the atom-local values used by the UFF sp3--sp3 rule. A zero
    /// value is intentional for atom classes whose torsion is determined by
    /// the U1 sp2 rule instead.
    fn v1(self) -> f64 {
        match self {
            Self::C_3 => 2.119,
            Self::N_3 => 0.450,
            Self::O_3 => 0.018,
            Self::S_3 => 0.484,
            _ => 0.0,
        }
    }

    /// UFF torsional parameter U1 (kcal/mol), used by the sp2--sp2 rule.
    fn u1(self) -> f64 {
        match self {
            Self::C_3
            | Self::C_2
            | Self::C_1
            | Self::C_R
            | Self::N_3
            | Self::N_2
            | Self::N_1
            | Self::N_R
            | Self::O_3
            | Self::O_2
            | Self::O_1
            | Self::O_R
            | Self::F_
            | Self::Cl
            | Self::Br
            | Self::I_ => 2.0,
            Self::S_3 | Self::S_2 | Self::S_R | Self::P_3 | Self::P_R => 1.25,
            _ => 0.0,
        }
    }

    fn is_sp3(self) -> bool {
        matches!(
            self,
            Self::C_3 | Self::N_3 | Self::O_3 | Self::S_3 | Self::P_3
        )
    }

    fn is_sp2(self) -> bool {
        matches!(
            self,
            Self::C_2
                | Self::C_R
                | Self::N_2
                | Self::N_R
                | Self::O_2
                | Self::O_R
                | Self::S_2
                | Self::S_R
                | Self::P_R
        )
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
    let r_bo = -lambda * rij * bond_order.dln();
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

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn scale_vec(v: [f64; 3], scale: f64) -> [f64; 3] {
    [v[0] * scale, v[1] * scale, v[2] * scale]
}

fn add_vec(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
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

#[derive(Clone, Copy)]
struct UffInversionTerm {
    center: AtomIdx,
    a: AtomIdx,
    b: AtomIdx,
    c: AtomIdx,
    /// Force constant for the squared signed volume (kcal/mol/Å⁶).
    k: f64,
}

#[derive(Clone, Copy)]
struct UffTorsionTerm {
    i: AtomIdx,
    j: AtomIdx,
    k: AtomIdx,
    l: AtomIdx,
    barrier: f64,
    periodicity: u8,
    phase_cosine: f64,
}

struct PreparedUffEnergy {
    bonds: Vec<UffBondTerm>,
    angles: Vec<UffAngleTerm>,
    vdw: Vec<UffVdwTerm>,
    inversions: Vec<UffInversionTerm>,
    torsions: Vec<UffTorsionTerm>,
}

impl PreparedUffEnergy {
    fn new(mol: &Molecule, types: &[(AtomIdx, UffType)]) -> Self {
        let type_map: std::collections::HashMap<AtomIdx, UffType> =
            types.iter().map(|&(a, t)| (a, t)).collect();
        let bonds = prepare_uff_bonds(mol, &type_map);
        let (angles, inversions) = prepare_uff_angles(mol, types);
        let torsions = prepare_uff_torsions(mol, &type_map);
        let vdw = prepare_uff_vdw(mol, &type_map);
        Self {
            bonds,
            angles,
            vdw,
            inversions,
            torsions,
        }
    }

    fn inversion_volume(term: UffInversionTerm, coords: &[[f64; 3]]) -> f64 {
        let center = coords[term.center.0 as usize];
        let u = sub(coords[term.a.0 as usize], center);
        let v = sub(coords[term.b.0 as usize], center);
        let w = sub(coords[term.c.0 as usize], center);
        dot(u, cross(v, w))
    }

    fn torsion_cosine(term: UffTorsionTerm, coords: &[[f64; 3]]) -> Option<f64> {
        let b1 = sub(coords[term.j.0 as usize], coords[term.i.0 as usize]);
        let b2 = sub(coords[term.k.0 as usize], coords[term.j.0 as usize]);
        let b3 = sub(coords[term.l.0 as usize], coords[term.k.0 as usize]);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let n1_len = dot(n1, n1).sqrt();
        let n2_len = dot(n2, n2).sqrt();
        if n1_len <= 1e-12 || n2_len <= 1e-12 {
            return None;
        }
        Some((dot(n1, n2) / (n1_len * n2_len)).clamp(-1.0, 1.0))
    }

    fn torsion_cosine_gradient(
        term: UffTorsionTerm,
        coords: &[[f64; 3]],
    ) -> Option<(f64, [[f64; 3]; 4])> {
        let b1 = sub(coords[term.j.0 as usize], coords[term.i.0 as usize]);
        let b2 = sub(coords[term.k.0 as usize], coords[term.j.0 as usize]);
        let b3 = sub(coords[term.l.0 as usize], coords[term.k.0 as usize]);
        let n1 = cross(b1, b2);
        let n2 = cross(b2, b3);
        let n1_len = dot(n1, n1).sqrt();
        let n2_len = dot(n2, n2).sqrt();
        if n1_len <= 1e-12 || n2_len <= 1e-12 {
            return None;
        }
        let c = (dot(n1, n2) / (n1_len * n2_len)).clamp(-1.0, 1.0);
        let nh1 = scale_vec(n1, 1.0 / n1_len);
        let nh2 = scale_vec(n2, 1.0 / n2_len);
        let ga = scale_vec(sub(nh2, scale_vec(nh1, c)), 1.0 / n1_len);
        let gb = scale_vec(sub(nh1, scale_vec(nh2, c)), 1.0 / n2_len);
        let db1 = cross(b2, ga);
        let db2 = add_vec(cross(ga, b1), cross(b3, gb));
        let db3 = cross(gb, b2);
        Some((c, [scale_vec(db1, -1.0), sub(db1, db2), sub(db2, db3), db3]))
    }

    fn torsion_energy(term: UffTorsionTerm, coords: &[[f64; 3]]) -> f64 {
        let Some(cos_phi) = Self::torsion_cosine(term, coords) else {
            return 0.0;
        };
        let cos_n_phi = match term.periodicity {
            2 => 2.0 * cos_phi * cos_phi - 1.0,
            3 => 4.0 * cos_phi.powi(3) - 3.0 * cos_phi,
            6 => {
                let cos_3 = 4.0 * cos_phi.powi(3) - 3.0 * cos_phi;
                2.0 * cos_3 * cos_3 - 1.0
            }
            _ => return 0.0,
        };
        0.5 * term.barrier * (1.0 - term.phase_cosine * cos_n_phi)
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
        for &term in &self.inversions {
            let volume = Self::inversion_volume(term, coords);
            energy += term.k * volume * volume;
        }
        for &term in &self.torsions {
            energy += Self::torsion_energy(term, coords);
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
        for &term in &self.inversions {
            let center = coords[term.center.0 as usize];
            let u = sub(coords[term.a.0 as usize], center);
            let v = sub(coords[term.b.0 as usize], center);
            let w = sub(coords[term.c.0 as usize], center);
            let volume = dot(u, cross(v, w));
            let scale = 2.0 * term.k * volume;
            let da = scale_vec(cross(v, w), scale);
            let db = scale_vec(cross(w, u), scale);
            let dc = scale_vec(cross(u, v), scale);
            add(&mut gradient[term.a.0 as usize], da);
            add(&mut gradient[term.b.0 as usize], db);
            add(&mut gradient[term.c.0 as usize], dc);
            add(
                &mut gradient[term.center.0 as usize],
                [
                    -da[0] - db[0] - dc[0],
                    -da[1] - db[1] - dc[1],
                    -da[2] - db[2] - dc[2],
                ],
            );
        }
        for &term in &self.torsions {
            let Some((cos_phi, basis)) = Self::torsion_cosine_gradient(term, coords) else {
                continue;
            };
            let d_cos_n = match term.periodicity {
                2 => 4.0 * cos_phi,
                3 => 12.0 * cos_phi * cos_phi - 3.0,
                6 => 192.0 * cos_phi.powi(5) - 192.0 * cos_phi.powi(3) + 36.0 * cos_phi,
                _ => continue,
            };
            let scale = -0.5 * term.barrier * term.phase_cosine * d_cos_n;
            for (atom, direction) in [term.i, term.j, term.k, term.l].into_iter().zip(basis) {
                let slot = &mut gradient[atom.0 as usize];
                slot[0] += scale * direction[0];
                slot[1] += scale * direction[1];
                slot[2] += scale * direction[2];
            }
        }
        gradient
    }
}

fn prepared_uff_type(
    types: &std::collections::HashMap<AtomIdx, UffType>,
    atom: AtomIdx,
) -> UffType {
    types.get(&atom).copied().unwrap_or(UffType::Unknown)
}

fn prepare_uff_bonds(
    mol: &Molecule,
    types: &std::collections::HashMap<AtomIdx, UffType>,
) -> Vec<UffBondTerm> {
    mol.bonds()
        .map(|(_, bond)| {
            let r0 = uff_bond_length(
                prepared_uff_type(types, bond.atom1),
                prepared_uff_type(types, bond.atom2),
                bond_order_f64(bond.order),
            );
            UffBondTerm {
                a: bond.atom1,
                b: bond.atom2,
                r0,
                k: 664.12 / (r0 * r0 * r0),
            }
        })
        .collect()
}

fn prepare_uff_angles(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
) -> (Vec<UffAngleTerm>, Vec<UffInversionTerm>) {
    let mut angles = Vec::new();
    let mut inversions = Vec::new();
    for &(center, center_type) in types {
        let theta0 = center_type.theta0().to_radians();
        let sin0 = theta0.dsin();
        let neighbors = mol
            .neighbors(center)
            .map(|(atom, _)| atom)
            .collect::<Vec<_>>();
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                angles.push(UffAngleTerm {
                    a: neighbors[i],
                    center,
                    c: neighbors[j],
                    cos0: theta0.dcos(),
                    k: 0.5 * 332.06 / (sin0 * sin0 + 1e-10),
                });
            }
        }
        if matches!(center_type, UffType::C_R | UffType::N_R | UffType::S_R) && neighbors.len() == 3
        {
            inversions.push(UffInversionTerm {
                center,
                a: neighbors[0],
                b: neighbors[1],
                c: neighbors[2],
                k: 0.5,
            });
        }
    }
    (angles, inversions)
}

fn prepare_uff_torsions(
    mol: &Molecule,
    types: &std::collections::HashMap<AtomIdx, UffType>,
) -> Vec<UffTorsionTerm> {
    let mut torsions = Vec::new();
    for (_, bond) in mol.bonds() {
        let (j, k) = (bond.atom1, bond.atom2);
        let (tj, tk) = (prepared_uff_type(types, j), prepared_uff_type(types, k));
        let Some((barrier, periodicity, phase_cosine)) =
            uff_torsion_barrier(tj, tk, bond_order_f64(bond.order))
        else {
            continue;
        };
        let left = mol
            .neighbors(j)
            .map(|(atom, _)| atom)
            .filter(|&atom| atom != k);
        let right = mol
            .neighbors(k)
            .map(|(atom, _)| atom)
            .filter(|&atom| atom != j)
            .collect::<Vec<_>>();
        for i in left {
            for &l in &right {
                torsions.push(UffTorsionTerm {
                    i,
                    j,
                    k,
                    l,
                    barrier,
                    periodicity,
                    phase_cosine,
                });
            }
        }
    }
    torsions
}

fn uff_torsion_barrier(tj: UffType, tk: UffType, order: f64) -> Option<(f64, u8, f64)> {
    let params = if tj.is_sp3() && tk.is_sp3() {
        ((tj.v1() * tk.v1()).sqrt(), 3, -1.0)
    } else if tj.is_sp2() && tk.is_sp2() {
        (
            5.0 * (tj.u1() * tk.u1()).sqrt() * (1.0 + 4.18 * order.dln()),
            2,
            1.0,
        )
    } else if (tj.is_sp2() && tk.is_sp3()) || (tj.is_sp3() && tk.is_sp2()) {
        (1.0, 6, 1.0)
    } else {
        return None;
    };
    (params.0 > 0.0 && params.0.is_finite()).then_some(params)
}

fn prepare_uff_vdw(
    mol: &Molecule,
    types: &std::collections::HashMap<AtomIdx, UffType>,
) -> Vec<UffVdwTerm> {
    let atoms = mol.atoms().map(|(atom, _)| atom).collect::<Vec<_>>();
    let excluded = uff_nonbonded_exclusions(mol);
    let mut terms = Vec::new();
    for (position, &a) in atoms.iter().enumerate() {
        for &b in &atoms[position + 1..] {
            let key = (a.0 as usize, b.0 as usize);
            if excluded.contains(&key) {
                continue;
            }
            let (ta, tb) = (prepared_uff_type(types, a), prepared_uff_type(types, b));
            terms.push(UffVdwTerm {
                a,
                b,
                x: (ta.x1() * tb.x1()).sqrt(),
                d: (ta.d1() * tb.d1()).sqrt(),
            });
        }
    }
    terms
}

fn uff_nonbonded_exclusions(mol: &Molecule) -> std::collections::HashSet<(usize, usize)> {
    let mut excluded = std::collections::HashSet::new();
    for (_, bond) in mol.bonds() {
        let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        excluded.insert((i.min(j), i.max(j)));
        for (neighbor, _) in mol.neighbors(bond.atom1) {
            let neighbor = neighbor.0 as usize;
            excluded.insert((neighbor.min(j), neighbor.max(j)));
        }
        for (neighbor, _) in mol.neighbors(bond.atom2) {
            let neighbor = neighbor.0 as usize;
            excluded.insert((i.min(neighbor), i.max(neighbor)));
        }
    }
    excluded
}

/// Compute UFF total energy (bond + angle + torsion + inversion + vdW) in
/// kcal/mol.
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

/// Trust radius for one UFF steepest-descent proposal (Å per atom).
///
/// A close non-bonded contact can produce a perfectly finite UFF vdW
/// gradient above 1e12 kcal/mol/Å. Multiplying that directly by the nominal
/// line-search step sends the first proposal many orders of magnitude beyond
/// molecular scale; repeatedly halving the nominal step then reaches the
/// `1e-8` floor before a physically local proposal is ever tried. Limiting
/// the largest atomic displacement keeps the proposal local without clipping
/// the energy or gradient, and the existing energy-decrease, geometry, and
/// caller-constraint gates still decide whether it is accepted.
const MAX_UFF_PROPOSAL_DISPLACEMENT: f64 = 0.1;

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
/// [`MAX_SANE_UFF_BOND_LENGTH`]. This is checked independently while the
/// optimizer runs, and also gates the public `converged` result:
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
    /// True if the gradient norm converged below threshold on a geometrically
    /// sound result. A low RMS gradient on an unsound stationary point is not
    /// reported as convergence; inspect [`Self::sound`] and
    /// [`Self::worst_bond_length`] for the diagnostic outcome.
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
    let mut trust_radius = MAX_UFF_PROPOSAL_DISPLACEMENT;
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
                // Even with the bounded organic torsion and aromatic
                // inversion terms, unsupported metal-class terms and other
                // incomplete UFF coverage can leave a low-gradient but
                // geometrically absurd stationary point. Do not expose that
                // state as successful convergence to direct callers.
                converged: sound,
                sound,
                worst_bond_length,
                rejected_unsound_step,
            };
        }

        // Line search: accept step only if energy decreases
        let max_gradient_norm = grad
            .iter()
            .map(|g| g[0].dhypot(g[1]).dhypot(g[2]))
            .fold(0.0_f64, f64::max);
        let proposal_step = if max_gradient_norm.is_finite() && max_gradient_norm > 0.0 {
            step.min(trust_radius / max_gradient_norm)
        } else {
            step
        };
        let new_coords: Vec<[f64; 3]> = coords
            .iter()
            .zip(&grad)
            .map(|(c, g)| {
                [
                    c[0] - proposal_step * g[0],
                    c[1] - proposal_step * g[1],
                    c[2] - proposal_step * g[2],
                ]
            })
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
            trust_radius = (trust_radius * 1.2).min(MAX_UFF_PROPOSAL_DISPLACEMENT);
            prev_energy = energy;
        } else {
            if new_energy < energy && !geometry_sound {
                rejected_unsound_step = true;
            }
            step *= 0.5;
            trust_radius *= 0.5;
            if step < 1e-8 || trust_radius < 1e-8 {
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
#[path = "uff_tests.rs"]
mod embedded_tests;
