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

/// Per-term UFF energy, for diagnosing which term dominates a given geometry
/// (issues #185/#188: `generate_coords`' starting geometries have enormous
/// total energy, and knowing whether that's bond/angle/vdW-dominated is
/// needed to pick a fix). `total` is a separately-accumulated value, not
/// required to be bit-identical to [`uff_total_energy`]'s single running
/// sum -- see that function's own doc comment for why the two are computed
/// independently rather than one delegating to the other.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UffEnergyBreakdown {
    pub bond: f64,
    pub angle: f64,
    pub vdw: f64,
    pub total: f64,
}

/// Same three energy terms as [`uff_total_energy`], accumulated into
/// separate per-term running sums instead of one combined sum, so a caller
/// can see which term dominates. Deliberately a **separate, independent**
/// implementation of the same three loops rather than a refactor of
/// `uff_total_energy` into "compute breakdown, then sum" -- floating-point
/// addition is not associative, so accumulating bond/angle/vdW into three
/// sums and adding those at the end is not guaranteed bit-identical to one
/// single running total accumulated bond-then-angle-then-vdW in sequence.
/// `uff_total_energy` is a published, load-bearing function already
/// exercised by production soundness gates in `chematic-3d`; this function
/// exists purely for diagnostics and must not risk perturbing that one by
/// even the last bit, so it is not refactored to share an accumulator with
/// it. `total` here is expected to numerically agree with
/// `uff_total_energy`'s result to well within any threshold this crate or
/// its consumers use (see `uff_energy_breakdown_total_matches_uff_total_energy`
/// in the test module) -- it is a genuinely separate value, verified close,
/// not assumed bit-equal.
pub fn uff_energy_breakdown(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    coords: &[[f64; 3]],
) -> UffEnergyBreakdown {
    let type_map: std::collections::HashMap<AtomIdx, UffType> =
        types.iter().map(|&(a, t)| (a, t)).collect();
    let get_type = |idx: AtomIdx| type_map.get(&idx).copied().unwrap_or(UffType::Unknown);
    let get_coord = |idx: AtomIdx| coords[idx.0 as usize];

    let mut bond_energy = 0.0;
    for (_, bond) in mol.bonds() {
        let ti = get_type(bond.atom1);
        let tj = get_type(bond.atom2);
        let n = bond_order_f64(bond.order);
        let r0 = uff_bond_length(ti, tj, n);
        let r = dist(get_coord(bond.atom1), get_coord(bond.atom2));
        let k = 664.12 / (r0 * r0 * r0);
        bond_energy += 0.5 * k * (r - r0) * (r - r0);
    }

    let mut angle_energy = 0.0;
    for (center_idx, center_type) in types {
        let theta0_deg = center_type.theta0();
        let theta0 = theta0_deg.to_radians();
        let cos0 = theta0.cos();
        let sin0 = theta0.sin();

        let neighbors: Vec<AtomIdx> = mol.neighbors(*center_idx).map(|(nb, _)| nb).collect();
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let cos_theta = cos_angle(
                    get_coord(neighbors[i]),
                    get_coord(*center_idx),
                    get_coord(neighbors[j]),
                );
                let delta = cos_theta - cos0;
                let k_angle = 0.5 * 332.06 / (sin0 * sin0 + 1e-10);
                angle_energy += 0.5 * k_angle * delta * delta;
            }
        }
    }

    let mut vdw_energy = 0.0;
    let atom_indices: Vec<AtomIdx> = mol.atoms().map(|(idx, _)| idx).collect();
    let n = atom_indices.len();
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
    for i in 0..n {
        for j in (i + 1)..n {
            let ai = atom_indices[i];
            let aj = atom_indices[j];
            let (ai_u, aj_u) = (ai.0 as usize, aj.0 as usize);
            if excl.contains(&(ai_u.min(aj_u), ai_u.max(aj_u))) {
                continue;
            }
            let ti = get_type(ai);
            let tj = get_type(aj);
            let x_ij = (ti.x1() * tj.x1()).sqrt();
            let d_ij = (ti.d1() * tj.d1()).sqrt();
            let r = dist(get_coord(ai), get_coord(aj)).max(0.5);
            let ratio = x_ij / r;
            let ratio6 = ratio.powi(6);
            let ratio12 = ratio6 * ratio6;
            vdw_energy += d_ij * (ratio12 - 2.0 * ratio6);
        }
    }

    UffEnergyBreakdown {
        bond: bond_energy,
        angle: angle_energy,
        vdw: vdw_energy,
        total: bond_energy + angle_energy + vdw_energy,
    }
}

/// Worst (longest) covalent bond length in `coords`, mirroring
/// `chematic-3d`'s own `worst_bond_length_vec`/`MAX_SANE_BOND_LENGTH`
/// convention exactly (this crate has no dependency on `chematic-3d`, so
/// this is a from-scratch but semantically-identical helper, not a shared
/// call) -- diagnostics that report this alongside `chematic-3d`'s own
/// soundness gate stay directly comparable.
pub fn uff_worst_bond_length(mol: &Molecule, coords: &[[f64; 3]]) -> f64 {
    mol.bonds()
        .map(|(_, b)| dist(coords[b.atom1.0 as usize], coords[b.atom2.0 as usize]))
        .fold(0.0_f64, f64::max)
}

/// Compute UFF total energy (bond + angle + vdW) in kcal/mol.
pub fn uff_total_energy(mol: &Molecule, types: &[(AtomIdx, UffType)], coords: &[[f64; 3]]) -> f64 {
    let type_map: std::collections::HashMap<AtomIdx, UffType> =
        types.iter().map(|&(a, t)| (a, t)).collect();
    let get_type = |idx: AtomIdx| type_map.get(&idx).copied().unwrap_or(UffType::Unknown);
    let get_coord = |idx: AtomIdx| coords[idx.0 as usize];

    let mut energy = 0.0;

    // ── Bond stretching ───────────────────────────────────────────────────
    // E_bond = k_ij/2 * (r - r0)^2   with k_ij = 664.12 * Z*_i * Z*_j / r0^3
    for (_, bond) in mol.bonds() {
        let ti = get_type(bond.atom1);
        let tj = get_type(bond.atom2);
        let n = bond_order_f64(bond.order);
        let r0 = uff_bond_length(ti, tj, n);
        let r = dist(get_coord(bond.atom1), get_coord(bond.atom2));
        // Force constant: simplified Badger's rule
        let k = 664.12 / (r0 * r0 * r0);
        energy += 0.5 * k * (r - r0) * (r - r0);
    }

    // ── Angle bending ─────────────────────────────────────────────────────
    // For sp3 / sp2 / sp centres use different Fourier expansion
    for (center_idx, center_type) in types {
        let theta0_deg = center_type.theta0();
        let theta0 = theta0_deg.to_radians();
        let cos0 = theta0.cos();
        let sin0 = theta0.sin();

        let neighbors: Vec<AtomIdx> = mol.neighbors(*center_idx).map(|(nb, _)| nb).collect();
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let cos_theta = cos_angle(
                    get_coord(neighbors[i]),
                    get_coord(*center_idx),
                    get_coord(neighbors[j]),
                );
                // Fourier: E = k/n^2 * C0 + C1*cos + C2*cos(2θ)
                // Simplified harmonic in cos space:
                let delta = cos_theta - cos0;
                let k_angle = 0.5 * 332.06 / (sin0 * sin0 + 1e-10);
                energy += 0.5 * k_angle * delta * delta;
            }
        }
    }

    // ── van der Waals (Lennard-Jones 12-6) ────────────────────────────────
    // Only 1-4+ pairs get vdW (1-2 bonded and 1-3 angle-bonded pairs are
    // excluded). `(i+2)..n` is an atom-*index* heuristic and only matches a
    // graph-based 1-3 exclusion for an unbranched chain visited in bond
    // order; build the real exclusion set from the bond graph instead
    // (mirrors `mmff94_minimizer.rs`'s `vdw_energy`, which already does this
    // correctly): every bonded pair, plus every pair that shares a common
    // bonded neighbor (i.e. is an angle apex away from each other).
    let atom_indices: Vec<AtomIdx> = mol.atoms().map(|(idx, _)| idx).collect();
    let n = atom_indices.len();
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
    for i in 0..n {
        for j in (i + 1)..n {
            let ai = atom_indices[i];
            let aj = atom_indices[j];
            let (ai_u, aj_u) = (ai.0 as usize, aj.0 as usize);
            if excl.contains(&(ai_u.min(aj_u), ai_u.max(aj_u))) {
                continue;
            }

            let ti = get_type(ai);
            let tj = get_type(aj);

            // UFF combining rules: x_ij = sqrt(x_i * x_j), D_ij = sqrt(D_i * D_j)
            let x_ij = (ti.x1() * tj.x1()).sqrt();
            let d_ij = (ti.d1() * tj.d1()).sqrt();

            let r = dist(get_coord(ai), get_coord(aj)).max(0.5);
            let ratio = x_ij / r;
            let ratio6 = ratio.powi(6);
            let ratio12 = ratio6 * ratio6;
            energy += d_ij * (ratio12 - 2.0 * ratio6);
        }
    }

    energy
}

// ── Gradient + L-BFGS minimizer ───────────────────────────────────────────────

/// Numerical gradient of UFF total energy with step δ = 1e-4 Å. `pub` so
/// diagnostic/experimental minimizer variants (issues #185/#188) can build
/// on the exact same gradient `minimize_uff` itself uses, rather than
/// re-deriving a numerically-similar-but-not-identical one externally.
pub fn uff_gradient(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    coords: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    const DELTA: f64 = 1e-4;
    let n = coords.len();
    let mut grad = vec![[0.0_f64; 3]; n];
    let mut perturbed = coords.to_vec();
    for i in 0..n {
        for k in 0..3 {
            perturbed[i][k] += DELTA;
            let ep = uff_total_energy(mol, types, &perturbed);
            perturbed[i][k] -= 2.0 * DELTA;
            let em = uff_total_energy(mol, types, &perturbed);
            perturbed[i][k] += DELTA;
            grad[i][k] = (ep - em) / (2.0 * DELTA);
        }
    }
    grad
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
    let mut coords = initial_coords;
    let mut step = 0.05_f64;
    let mut prev_energy = f64::MAX;

    for iter in 0..max_iter {
        let energy = uff_total_energy(mol, types, &coords);
        let grad = uff_gradient(mol, types, &coords);

        // RMS gradient norm
        let rms: f64 = {
            let sum2: f64 = grad.iter().flat_map(|g| g.iter()).map(|v| v * v).sum();
            (sum2 / (grad.len() * 3) as f64).sqrt()
        };

        if rms < 0.01 {
            return UffMinimizeResult {
                coords,
                energy,
                iterations: iter,
                converged: true,
            };
        }

        // Line search: accept step only if energy decreases
        let new_coords: Vec<[f64; 3]> = coords
            .iter()
            .zip(&grad)
            .map(|(c, g)| [c[0] - step * g[0], c[1] - step * g[1], c[2] - step * g[2]])
            .collect();

        let new_energy = uff_total_energy(mol, types, &new_coords);
        if new_energy < energy {
            coords = new_coords;
            if energy - new_energy < prev_energy * 1e-7 {
                step *= 1.2;
            }
            prev_energy = energy;
        } else {
            step *= 0.5;
            if step < 1e-8 {
                return UffMinimizeResult {
                    coords,
                    energy,
                    iterations: iter,
                    converged: false,
                };
            }
        }
    }

    let energy = uff_total_energy(mol, types, &coords);
    UffMinimizeResult {
        coords,
        energy,
        iterations: max_iter,
        converged: false,
    }
}

/// One outer-loop iteration's diagnostic snapshot from
/// [`minimize_uff_with_trace`] -- issues #185/#188 (catastrophic bond-length
/// blowup on some starting geometries, sound convergence on others) need
/// per-iteration energy/gradient/step-size/geometry history to diagnose,
/// which [`minimize_uff`]'s plain [`UffMinimizeResult`] summary can't show.
#[derive(Debug, Clone, Copy)]
pub struct UffTraceStep {
    pub iter: usize,
    /// Total energy at the start of this iteration (before this iteration's
    /// step is attempted) -- the same value `minimize_uff`'s own loop uses
    /// for its accept/reject decision.
    pub energy: f64,
    pub bond_energy: f64,
    pub angle_energy: f64,
    pub vdw_energy: f64,
    /// RMS gradient norm -- the same quantity `minimize_uff` checks against
    /// its 0.01 convergence threshold.
    pub grad_rms: f64,
    /// Largest single gradient component (kcal/mol/Å) -- comparable to
    /// `chematic-3d`'s own `max_residual_force` finite-difference metric.
    pub grad_max: f64,
    /// The step size used to compute this iteration's attempted move,
    /// before any post-hoc 1.2x growth or 0.5x shrink is applied.
    pub step_size: f64,
    /// Whether this iteration's attempted step was accepted (energy
    /// decreased) or rejected (step will be halved next iteration).
    pub accepted: bool,
    pub worst_bond_length: f64,
    /// Largest single-atom displacement (Å) this iteration's attempted step
    /// would produce, whether or not it was ultimately accepted.
    pub max_atom_displacement: f64,
}

/// Identical algorithm to [`minimize_uff`] (same steepest-descent loop, same
/// accept/reject/step-growth/shrink rules, same convergence threshold),
/// additionally recording one [`UffTraceStep`] per outer-loop iteration.
///
/// Deliberately a **separate, parallel implementation** rather than
/// `minimize_uff` calling into this one internally (or vice versa) --
/// `minimize_uff` is a published, load-bearing function already relied on by
/// `chematic-3d`'s force-field bridge; changing its control flow to thread a
/// trace-collection callback through it (even one that's a no-op by default)
/// is exactly the kind of restructuring this diagnostic work must not risk
/// perturbing. The regression test
/// `traced_minimizer_matches_untraced_minimizer_byte_for_byte` proves the two
/// loops haven't drifted apart by asserting identical
/// `iterations`/`converged`/`energy`/`coords` on the same inputs, rather than
/// relying on the two implementations merely *looking* parallel.
pub fn minimize_uff_with_trace(
    mol: &Molecule,
    types: &[(AtomIdx, UffType)],
    initial_coords: Vec<[f64; 3]>,
    max_iter: usize,
) -> (UffMinimizeResult, Vec<UffTraceStep>) {
    let mut coords = initial_coords;
    let mut step = 0.05_f64;
    let mut prev_energy = f64::MAX;
    let mut trace = Vec::new();

    for iter in 0..max_iter {
        let energy = uff_total_energy(mol, types, &coords);
        let breakdown = uff_energy_breakdown(mol, types, &coords);
        let grad = uff_gradient(mol, types, &coords);

        let rms: f64 = {
            let sum2: f64 = grad.iter().flat_map(|g| g.iter()).map(|v| v * v).sum();
            (sum2 / (grad.len() * 3) as f64).sqrt()
        };
        let grad_max = grad
            .iter()
            .flat_map(|g| g.iter())
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        let worst_bond = uff_worst_bond_length(mol, &coords);

        if rms < 0.01 {
            trace.push(UffTraceStep {
                iter,
                energy,
                bond_energy: breakdown.bond,
                angle_energy: breakdown.angle,
                vdw_energy: breakdown.vdw,
                grad_rms: rms,
                grad_max,
                step_size: step,
                accepted: true,
                worst_bond_length: worst_bond,
                max_atom_displacement: 0.0,
            });
            return (
                UffMinimizeResult {
                    coords,
                    energy,
                    iterations: iter,
                    converged: true,
                },
                trace,
            );
        }

        let new_coords: Vec<[f64; 3]> = coords
            .iter()
            .zip(&grad)
            .map(|(c, g)| [c[0] - step * g[0], c[1] - step * g[1], c[2] - step * g[2]])
            .collect();
        let max_disp = coords
            .iter()
            .zip(&new_coords)
            .map(|(c, nc)| dist(*c, *nc))
            .fold(0.0_f64, f64::max);

        let new_energy = uff_total_energy(mol, types, &new_coords);
        let accepted = new_energy < energy;

        trace.push(UffTraceStep {
            iter,
            energy,
            bond_energy: breakdown.bond,
            angle_energy: breakdown.angle,
            vdw_energy: breakdown.vdw,
            grad_rms: rms,
            grad_max,
            step_size: step,
            accepted,
            worst_bond_length: worst_bond,
            max_atom_displacement: max_disp,
        });

        if accepted {
            coords = new_coords;
            if energy - new_energy < prev_energy * 1e-7 {
                step *= 1.2;
            }
            prev_energy = energy;
        } else {
            step *= 0.5;
            if step < 1e-8 {
                return (
                    UffMinimizeResult {
                        coords,
                        energy,
                        iterations: iter,
                        converged: false,
                    },
                    trace,
                );
            }
        }
    }

    let energy = uff_total_energy(mol, types, &coords);
    (
        UffMinimizeResult {
            coords,
            energy,
            iterations: max_iter,
            converged: false,
        },
        trace,
    )
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
        // Gradient check following this crate's existing pattern (finite-
        // difference gradient, no analytic gradient exists in chematic-ff —
        // see PR #169's own finding on this): at a strained, closed-angle
        // propane geometry, stepping a small distance along -gradient must
        // lower the energy.
        let (mol, types, coords) = propane_at_angle(60.0);
        let e0 = uff_total_energy(&mol, &types, &coords);
        let grad = uff_gradient(&mol, &types, &coords);
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

    /// `uff_energy_breakdown` is a deliberately separate accumulator from
    /// `uff_total_energy` (see that function's own doc comment for why) --
    /// this proves the two still agree numerically on a real, non-trivial
    /// strained geometry, not just that both compile.
    #[test]
    fn uff_energy_breakdown_total_matches_uff_total_energy() {
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let types = assign_uff_types(&mol);
        // Deliberately distorted coordinates (not a real embedding) so bond,
        // angle, and vdW terms are all simultaneously nonzero and large.
        let coords: Vec<[f64; 3]> = (0..mol.atom_count())
            .map(|i| {
                let t = i as f64;
                [t.sin() * 2.0, t.cos() * 1.3, (t * 0.7).sin()]
            })
            .collect();
        let total = uff_total_energy(&mol, &types, &coords);
        let breakdown = uff_energy_breakdown(&mol, &types, &coords);
        assert!(
            (breakdown.total - total).abs() < 1e-9 * total.abs().max(1.0),
            "breakdown.total ({}) should numerically agree with uff_total_energy ({})",
            breakdown.total,
            total
        );
        assert_eq!(
            breakdown.bond + breakdown.angle + breakdown.vdw,
            breakdown.total
        );
    }

    /// The traced and untraced minimizers must be genuinely parallel
    /// implementations, not just similar-looking ones -- assert byte-for-byte
    /// identical `UffMinimizeResult`s across several starting geometries and
    /// budgets, including ones that terminate via each of the three return
    /// paths (converged, step-collapsed, budget-exhausted).
    #[test]
    fn traced_minimizer_matches_untraced_minimizer_byte_for_byte() {
        let cases: Vec<(&str, usize)> = vec![
            ("CCO", 200),
            ("c1ccc2ccccc2c1", 50),
            ("c1ccc2ccccc2c1", 5000),
            ("CCCCCC", 30),
        ];
        for (smi, max_iter) in cases {
            let mol = parse(smi).expect("parse");
            let types = assign_uff_types(&mol);
            let coords: Vec<[f64; 3]> = (0..mol.atom_count())
                .map(|i| [1.5 * i as f64, 0.0, 0.0])
                .collect();

            let untraced = minimize_uff(&mol, &types, coords.clone(), max_iter);
            let (traced, steps) = minimize_uff_with_trace(&mol, &types, coords, max_iter);

            assert_eq!(untraced.iterations, traced.iterations, "smi={smi}");
            assert_eq!(untraced.converged, traced.converged, "smi={smi}");
            assert_eq!(untraced.energy, traced.energy, "smi={smi}");
            assert_eq!(untraced.coords, traced.coords, "smi={smi}");
            assert!(
                !steps.is_empty(),
                "smi={smi} should record at least one step"
            );
            // The loop body pushes exactly one trace entry per iteration it runs,
            // then either returns early (converged, or step collapsed below
            // 1e-8 -- both report `iterations` as the 0-indexed loop variable at
            // return, so one more entry than that index) or falls through after
            // exhausting `max_iter` iterations (0-indexed loop ran exactly
            // `max_iter` times, and `iterations` is set to `max_iter` itself
            // outside the loop -- so entries == iterations there, not +1).
            let budget_exhausted = !traced.converged && traced.iterations == max_iter;
            let expected_len = if budget_exhausted {
                traced.iterations
            } else {
                traced.iterations + 1
            };
            assert_eq!(
                steps.len(),
                expected_len,
                "smi={smi} max_iter={max_iter} converged={} iterations={}",
                traced.converged,
                traced.iterations
            );
        }
    }
}
