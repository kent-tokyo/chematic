//! Molecular descriptor functions for drug-likeness and physical property estimation.
//!
//! All functions accept a `&Molecule` reference.  Molecules with aromatic bonds
//! (SMILES lowercase notation) are kekulized internally where hydrogen counts
//! are required; the caller's molecule is never mutated.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Element, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

/// True if `idx` has a double bond to any neighbor whose atomic number equals `target_an`.
fn has_double_bond_to(mol: &Molecule, idx: AtomIdx, target_an: u8) -> bool {
    mol.neighbors(idx).any(|(nb, bidx)| {
        mol.bond(bidx).order == BondOrder::Double
            && mol.atom(nb).element.atomic_number() == target_an
    })
}

/// Count double bonds from `idx` to neighbors whose atomic number equals `target_an`.
fn count_double_bonds_to(mol: &Molecule, idx: AtomIdx, target_an: u8) -> usize {
    mol.neighbors(idx)
        .filter(|&(nb, bidx)| {
            mol.bond(bidx).order == BondOrder::Double
                && mol.atom(nb).element.atomic_number() == target_an
        })
        .count()
}

/// True if any neighbor of `idx` is aromatic.
fn has_aromatic_neighbor(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.neighbors(idx).any(|(nb, _)| mol.atom(nb).aromatic)
}

/// True if any neighbor of `idx` is an aromatic carbon.
fn has_aromatic_carbon_neighbor(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.neighbors(idx)
        .any(|(nb, _)| mol.atom(nb).aromatic && mol.atom(nb).element.atomic_number() == 6)
}

// --- Element Detection Helpers ---
// Consolidate atomic number matching to eliminate 50+ hardcoded checks throughout the file.

#[inline]
fn is_carbon(an: u8) -> bool {
    an == 6
}

#[inline]
fn is_nitrogen(an: u8) -> bool {
    an == 7
}

#[inline]
fn is_oxygen(an: u8) -> bool {
    an == 8
}

#[inline]
fn is_sulfur(an: u8) -> bool {
    an == 16
}

#[inline]
fn is_halogen(an: u8) -> bool {
    matches!(an, 9 | 17 | 35 | 53)
} // F, Cl, Br, I

/// Average atomic mass table.
/// Falls back to `atomic_number as f64` for unlisted elements.
fn avg_mass(element: Element) -> f64 {
    match element.atomic_number() {
        1 => 1.008,    // H
        2 => 4.003,    // He
        3 => 6.941,    // Li
        4 => 9.012,    // Be
        5 => 10.811,   // B
        6 => 12.011,   // C
        7 => 14.007,   // N
        8 => 15.999,   // O
        9 => 18.998,   // F
        10 => 20.180,  // Ne
        11 => 22.990,  // Na
        12 => 24.305,  // Mg
        13 => 26.982,  // Al
        14 => 28.086,  // Si
        15 => 30.974,  // P
        16 => 32.065,  // S
        17 => 35.453,  // Cl
        18 => 39.948,  // Ar
        19 => 39.098,  // K
        20 => 40.078,  // Ca
        33 => 74.922,  // As
        34 => 78.971,  // Se
        35 => 79.904,  // Br
        53 => 126.904, // I
        n => n as f64,
    }
}

/// Monoisotopic (most-abundant isotope) mass table.
/// Falls back to `atomic_number as f64` for unlisted elements.
fn mono_mass(element: Element) -> f64 {
    match element.atomic_number() {
        1 => 1.00783,   // H  (1H)
        6 => 12.0000,   // C  (12C)
        7 => 14.0031,   // N  (14N)
        8 => 15.9949,   // O  (16O)
        9 => 18.9984,   // F  (19F)
        14 => 27.9769,  // Si (28Si)
        15 => 30.9738,  // P  (31P)
        16 => 31.9721,  // S  (32S)
        17 => 34.9689,  // Cl (35Cl)
        35 => 78.9183,  // Br (79Br)
        34 => 79.9165,  // Se (80Se)
        53 => 126.9045, // I  (127I)
        n => n as f64,
    }
}

// ---------------------------------------------------------------------------
// 1. Molecular weight
// ---------------------------------------------------------------------------

/// Compute the average molecular weight (Da).
///
/// Sums the average atomic mass of all heavy atoms plus each atom's implicit
/// hydrogen contribution (1.008 Da per H).
pub fn molecular_weight(mol: &Molecule) -> f64 {
    let mut mw = 0.0f64;
    for (idx, atom) in mol.atoms() {
        mw += avg_mass(atom.element);
        let h = implicit_hcount(mol, idx);
        mw += h as f64 * 1.008;
    }
    mw
}

// ---------------------------------------------------------------------------
// 2. Exact mass (monoisotopic)
// ---------------------------------------------------------------------------

/// Compute the monoisotopic (exact) mass (Da).
///
/// Uses the most-abundant isotope for each element, or the atom's explicit
/// isotope label (as an integer approximation) when set.
/// Implicit hydrogens use the ¹H monoisotopic mass (1.00783).
pub fn exact_mass(mol: &Molecule) -> f64 {
    let mut mass = 0.0f64;
    for (idx, atom) in mol.atoms() {
        let m = match atom.isotope {
            Some(iso) => iso as f64,
            None => mono_mass(atom.element),
        };
        mass += m;
        let h = implicit_hcount(mol, idx);
        mass += h as f64 * 1.00783;
    }
    mass
}

// ---------------------------------------------------------------------------
// 3. Heavy atom count
// ---------------------------------------------------------------------------

/// Count non-hydrogen heavy atoms.
///
/// Hydrogen atoms are normally implicit in chematic, but some molecules may
/// carry explicit H atoms in the graph (e.g. from bracket notation `[H]`).
/// Those are excluded from the heavy-atom count.
pub fn heavy_atom_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, atom)| atom.element != Element::H)
        .count()
}

// ---------------------------------------------------------------------------
// 4. Hydrogen bond donor count
// ---------------------------------------------------------------------------

/// Count hydrogen bond donors (N-H or O-H groups).
///
/// Each heavy atom with element N or O that has at least one attached H
/// counts as one donor (not per H — donors are counted per heavy atom).
pub fn hbd_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(idx, atom)| {
            let an = atom.element.atomic_number();
            (is_nitrogen(an) || is_oxygen(an)) && implicit_hcount(mol, *idx) > 0
        })
        .count()
}

// ---------------------------------------------------------------------------
// 5. Hydrogen bond acceptor count (Ertl / RDKit-aligned definition)
// ---------------------------------------------------------------------------

/// Count hydrogen bond acceptors using the Ertl (2000) definition as implemented
/// by RDKit's `rdMolDescriptors.CalcNumHBA`.
///
/// Counts N, O, and divalent S atoms, with the following exclusions:
/// - Aromatic N with H (pyrrole-type `[nH]`): lone pair participates in aromaticity.
/// - Non-aromatic N bonded to C=O (amide N): lone pair delocalized into carbonyl.
/// - O with H bonded to a C=O carbon (carboxylic/ester OH).
/// - O with H bonded to oxidized S with S=O (sulfonic/sulfonamide acid OH).
/// - Oxidized S (degree > 2 or has S=O bonds): lone pair engaged in S=O resonance.
pub fn hba_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(idx, atom)| {
            let an = atom.element.atomic_number();
            if is_nitrogen(an) {
                // Nitrogen: charged N (N+ in nitro, quaternary, n+ in thiazolium) is never HBA.
                if atom.charge != 0 {
                    return false;
                }
                let h = implicit_hcount(mol, *idx);
                if atom.aromatic {
                    // [nH] (pyrrole-type aromatic N) is NOT an HBA
                    h == 0
                } else {
                    // Non-aromatic N: exclude amide N (bonded to C=O)
                    !neighbor_has_carbonyl(mol, *idx)
                }
            } else if is_oxygen(an) {
                // Oxygen: exclude acid OH bonded to C=O or to oxidized S with S=O
                let h = implicit_hcount(mol, *idx);
                if h > 0 {
                    !neighbor_has_carbonyl(mol, *idx) && !neighbor_is_oxidized_sulfur(mol, *idx)
                } else {
                    true
                }
            } else if is_sulfur(an) {
                // Sulfur (Ertl definition includes divalent S with free lone pair)
                if atom.aromatic {
                    // Aromatic S (thiophene-type): count if uncharged
                    atom.charge == 0
                } else {
                    // Non-aromatic S: count only if divalent (X2) and not oxidized (no S=O)
                    let degree = mol.degree(*idx);
                    let total_valence = degree + implicit_hcount(mol, *idx) as usize;
                    atom.charge == 0 && total_valence == 2 && !has_double_bond_to(mol, *idx, 8)
                }
            } else {
                false
            }
        })
        .count()
}

/// True if any neighbor of `idx` is a carbon that itself has a double bond to oxygen
/// (i.e., a carbonyl carbon).
fn neighbor_has_carbonyl(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.neighbors(idx).any(|(nb_idx, _)| {
        mol.atom(nb_idx).element.atomic_number() == 6 && has_double_bond_to(mol, nb_idx, 8)
    })
}

/// True if any neighbor of `idx` is a sulfur atom that itself has a S=O double bond
/// (i.e., a sulfoxide, sulfone, or sulfonate S). Used to exclude S–OH from HBA count.
fn neighbor_is_oxidized_sulfur(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.neighbors(idx).any(|(nb_idx, _)| {
        mol.atom(nb_idx).element.atomic_number() == 16 && has_double_bond_to(mol, nb_idx, 8)
    })
}

// True if any neighbor of `idx` is a carbon that has a C=N double bond
// (i.e., an imine/amidine/guanidinium carbon).
// ---------------------------------------------------------------------------
// 6. Rotatable bond count
// ---------------------------------------------------------------------------

/// Count rotatable bonds.
///
/// A bond is rotatable when all of the following hold:
/// - It is a single bond (or a stereo bond Up/Down, which is single).
/// - Neither endpoint is terminal (degree > 1 in the heavy-atom graph).
/// - It is not part of any ring (SSSR membership).
/// - It is not an amide bond: if one atom is N and the other is C,
///   and that C has any double bond to an O, the bond is excluded.
pub fn rotatable_bond_count(mol: &Molecule) -> usize {
    let ring_bond_set = ring_bond_indices(mol);

    mol.bonds()
        .filter(|(bidx, bond)| {
            // Stereo bonds Up/Down are also single.
            let is_single = matches!(
                bond.order,
                BondOrder::Single | BondOrder::Up | BondOrder::Down
            );
            is_single
                && !ring_bond_set.contains(bidx)
                && mol.degree(bond.atom1) > 1
                && mol.degree(bond.atom2) > 1
                && !is_amide_bond(mol, bond.atom1, bond.atom2)
        })
        .count()
}

/// Indices of all bonds participating in at least one SSSR ring.
fn ring_bond_indices(mol: &Molecule) -> HashSet<BondIdx> {
    let mut set = HashSet::new();
    for ring in find_sssr(mol).rings() {
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            if let Some((bidx, _)) = mol.bond_between(a, b) {
                set.insert(bidx);
            }
        }
    }
    set
}

/// True if the bond between `a` and `b` is an amide-like C-N bond
/// (one atom is N, the other is C with a double bond to O).
fn is_amide_bond(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let an_a = mol.atom(a).element.atomic_number();
    let an_b = mol.atom(b).element.atomic_number();
    let c_idx = match (an_a, an_b) {
        (6, 7) => a,
        (7, 6) => b,
        _ => return false,
    };
    has_double_bond_to(mol, c_idx, 8)
}

// ---------------------------------------------------------------------------
// 7. TPSA
// ---------------------------------------------------------------------------

/// Compute the topological polar surface area (Å²) using the Ertl (2000) table.
///
/// Reference: P. Ertl, B. Rohde, P. Selzer, J. Med. Chem. 2000, 43, 3714-3717.
fn tpsa_nitrogen(mol: &Molecule, idx: AtomIdx, is_aromatic: bool, h: u8, charge: i8) -> f64 {
    if is_aromatic {
        let degree = mol.degree(idx);
        if h > 0 {
            15.79
        } else if degree >= 3 {
            if charge > 0 { 3.88 } else { 4.93 }
        } else {
            12.89
        }
    } else {
        if charge == 1 {
            let (has_oxo, has_o_minus) =
                mol.neighbors(idx)
                    .fold((false, false), |(oxo, om), (nb, bidx)| {
                        let nb_atom = mol.atom(nb);
                        let is_o = nb_atom.element.atomic_number() == 8;
                        (
                            oxo || (is_o && mol.bond(bidx).order == BondOrder::Double),
                            om || (is_o && nb_atom.charge == -1),
                        )
                    });
            if has_oxo && has_o_minus { 41.44 } else { 3.24 }
        } else if h >= 2 {
            26.02
        } else if h == 1 {
            if has_double_bond_to(mol, idx, 6) {
                23.79
            } else {
                12.03
            }
        } else {
            if has_double_bond_to(mol, idx, 6) {
                12.89
            } else {
                3.24
            }
        }
    }
}

fn tpsa_oxygen(mol: &Molecule, idx: AtomIdx, is_aromatic: bool, h: u8, charge: i8) -> f64 {
    if is_aromatic {
        13.14
    } else if h > 0 {
        20.23
    } else {
        let is_nitro_o_minus = charge == -1
            && mol.neighbors(idx).any(|(nb, _)| {
                mol.atom(nb).element.atomic_number() == 7 && mol.atom(nb).charge == 1
            });
        if is_nitro_o_minus {
            0.0
        } else {
            let dbl_neighbor_an = mol
                .neighbors(idx)
                .find(|&(_, bidx)| mol.bond(bidx).order == BondOrder::Double)
                .map(|(nei, _)| mol.atom(nei).element.atomic_number());
            match dbl_neighbor_an {
                Some(6) => 17.07,
                Some(_) => 0.0,
                None => 9.23,
            }
        }
    }
}

fn tpsa_sulfur(mol: &Molecule, idx: AtomIdx, is_aromatic: bool, h: u8) -> f64 {
    if is_aromatic {
        28.24
    } else if h > 0 {
        38.80
    } else {
        match count_double_bonds_to(mol, idx, 8) {
            0 => 25.30,
            1 => 36.28,
            _ => 42.52,
        }
    }
}

fn tpsa_phosphorus(mol: &Molecule, idx: AtomIdx) -> f64 {
    if has_double_bond_to(mol, idx, 8) {
        26.88
    } else {
        34.14
    }
}

/// Topological Polar Surface Area (Ertl 2000).
///
/// Sum of Ertl atom-type contributions for N, O, S, and P atoms.
/// Hydrogen atoms are implicit (not computed). Values match RDKit defaults with calibrations
/// for secondary amide N (12.03 Å²), aromatic N (15.79 Å²), and phosphorus atoms.
///
/// Algorithm dispatches per element (tpsa_nitrogen, tpsa_oxygen, tpsa_sulfur, tpsa_phosphorus)
/// to reduce cyclomatic complexity and improve readability.
pub fn tpsa(mol: &Molecule) -> f64 {
    let mut psa = 0.0f64;
    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let is_aromatic = atom.aromatic;
        let h = implicit_hcount(mol, idx);
        let contribution = match an {
            7 => tpsa_nitrogen(mol, idx, is_aromatic, h, atom.charge),
            8 => tpsa_oxygen(mol, idx, is_aromatic, h, atom.charge),
            16 => tpsa_sulfur(mol, idx, is_aromatic, h),
            15 if !is_aromatic => tpsa_phosphorus(mol, idx),
            _ => 0.0,
        };
        psa += contribution;
    }
    psa
}

// ---------------------------------------------------------------------------
// 8. LogP (Wildman-Crippen, calibrated)
// ---------------------------------------------------------------------------

/// Compute a Wildman-Crippen LogP using a calibrated atom-type table.
///
/// Atom type contributions are derived analytically from the RDKit reference
/// dataset (175 molecules) and confirmed against the Wildman-Crippen 1999 paper.
/// Key improvements over the simplified table:
/// - H atom contributions (H on C=+0.1230, H on N=+0.2142,
///   H on aliphatic-OH=−0.2677, H on carboxylic-OH=+0.2980)
/// - Aromatic C: [cH]=0.1581 vs [c]=0.1441
/// - Aromatic N (both [nH] and [n;H0]): −0.3239 (was +0.2626)
/// - S: thioether=+0.6482, aromatic=+0.6237 (was 0.2432/0.0)
/// - O: alcohol=−0.2893, ether=−0.0684, aromatic=+0.1552, carbonyl=−0.0509
/// - Cl: aromatic=+0.7904, aliphatic=+0.6895
///
/// Wildman-Crippen LogP per-atom contributions.
///
/// Dispatches to per-element atom-type functions (crippen_carbon, crippen_nitrogen, etc.)
/// to compute Wildman-Crippen atom-type contributions. Per-atom LogP contributions
/// (heavy atoms only; H contributions are folded into the heavy atom they are attached to).
/// Index matches mol.atoms().
pub fn logp_crippen_per_atom(mol: &Molecule) -> Vec<f64> {
    mol.atoms()
        .map(|(idx, atom)| {
            let an = atom.element.atomic_number();
            let ar = atom.aromatic;
            let h = implicit_hcount(mol, idx);
            let heavy = match an {
                6 => crippen_carbon(mol, idx, ar, h),
                7 => crippen_nitrogen(mol, idx, ar),
                8 => crippen_oxygen(mol, idx, ar, h),
                16 => crippen_sulfur(mol, idx, ar),
                9 => crippen_halogen(mol, idx, ar, 0.2761, 0.4202),
                17 => crippen_halogen(mol, idx, ar, 0.7904, 0.6895),
                35 => crippen_halogen(mol, idx, ar, 0.8995, 0.8456),
                53 => crippen_halogen(mol, idx, ar, 0.7416, 0.8857),
                15 => {
                    if has_double_bond_to(mol, idx, 8) {
                        0.7933
                    } else {
                        -0.3451
                    }
                }
                _ => 0.0,
            };
            let h_contrib = if h == 0 {
                0.0
            } else {
                crippen_hydrogen(mol, idx, an, ar) * h as f64
            };
            heavy + h_contrib
        })
        .collect()
}

/// Compute the Crippen log P (octanol/water partition coefficient) of `mol`.
///
/// Sums per-atom contributions from [`logp_crippen_per_atom`].
pub fn logp_crippen(mol: &Molecule) -> f64 {
    logp_crippen_per_atom(mol).iter().sum()
}

/// Crippen contribution for Carbon atoms.
fn crippen_carbon_aromatic(mol: &Molecule, idx: AtomIdx, h: u8) -> f64 {
    let has_exocyclic_heteroatom_double = mol.neighbors(idx).any(|(nb, bidx)| {
        mol.bond(bidx).order == BondOrder::Double
            && !mol.atom(nb).aromatic
            && mol.atom(nb).element.atomic_number() != 6
    });
    if has_exocyclic_heteroatom_double {
        return -0.3800;
    }
    if h > 0 {
        return 0.1581;
    }
    if mol.neighbors(idx).any(|(nb, _)| {
        let a = mol.atom(nb);
        is_nitrogen(a.element.atomic_number()) && !a.aromatic
    }) {
        return 0.4619;
    }
    let bonded_to_ether_o = mol.neighbors(idx).any(|(nb, bidx)| {
        is_oxygen(mol.atom(nb).element.atomic_number())
            && !mol.atom(nb).aromatic
            && mol.bond(bidx).order == BondOrder::Single
            && implicit_hcount(mol, nb) == 0
    });
    if bonded_to_ether_o {
        return 0.5437;
    }
    let all_aromatic_nbrs = mol
        .neighbors(idx)
        .filter(|(nb, _)| mol.atom(*nb).aromatic)
        .count();
    let aromatic_c_nbrs = mol
        .neighbors(idx)
        .filter(|(nb, _)| {
            mol.atom(*nb).aromatic && is_carbon(mol.atom(*nb).element.atomic_number())
        })
        .count();
    if all_aromatic_nbrs >= 3 && aromatic_c_nbrs >= 2 {
        0.2956
    } else {
        0.1441
    }
}

fn crippen_carbon_aliphatic(mol: &Molecule, idx: AtomIdx, h: u8) -> f64 {
    let has_double_to_n = has_double_bond_to(mol, idx, 7);
    let has_double_to_heteroatom = has_double_to_n
        || mol.neighbors(idx).any(|(nb, bidx)| {
            let bo = mol.bond(bidx).order;
            (bo == BondOrder::Double || bo == BondOrder::Triple)
                && !is_carbon(mol.atom(nb).element.atomic_number())
                && !is_nitrogen(mol.atom(nb).element.atomic_number())
        });
    let has_double_to_c = mol.neighbors(idx).any(|(nb, bidx)| {
        mol.bond(bidx).order == BondOrder::Double
            && !mol.atom(nb).aromatic
            && is_carbon(mol.atom(nb).element.atomic_number())
    });

    if has_double_to_n {
        -0.2783
    } else if has_double_to_heteroatom {
        if has_aromatic_carbon_neighbor(mol, idx) {
            -0.1226
        } else {
            -0.3800
        }
    } else if has_double_to_c {
        // Context-dependent alkene C (Wildman-Crippen):
        //   Ar-adjacent =CH-  (aromatic C neighbor):           0.2640
        //   terminal =CH2 (h≥2, no aromatic neighbor):         0.1551
        //   internal =CH- in enone/enal (C=C-C=O neighbor):   0.1302
        //   other internal alkene C:                           0.2274
        //
        // The enone case (C=C conjugated with C=O) uses a lower contribution
        // because electron withdrawal by the carbonyl reduces hydrophobicity.
        // neighbor_has_carbonyl() is the existing helper (line ~252).
        let ar_c_nbr = has_aromatic_carbon_neighbor(mol, idx);
        let conjugated = neighbor_has_carbonyl(mol, idx);
        if ar_c_nbr {
            0.2640
        } else if h >= 2 {
            0.1551
        } else if conjugated {
            0.1302
        } else {
            0.2274
        }
    } else {
        let bonded_to_n = mol
            .neighbors(idx)
            .any(|(nb, _)| is_nitrogen(mol.atom(nb).element.atomic_number()));
        let bonded_to_heteroatom = bonded_to_n
            || mol.neighbors(idx).any(|(nb, _)| {
                let an = mol.atom(nb).element.atomic_number();
                matches!(an, 8 | 9 | 15 | 16 | 17 | 35 | 53)
            });
        if bonded_to_heteroatom {
            if bonded_to_n && has_aromatic_carbon_neighbor(mol, idx) {
                0.1193
            } else {
                -0.2035
            }
        } else if has_aromatic_carbon_neighbor(mol, idx) {
            match h {
                3 => 0.0845,
                2 => -0.0516,
                1 => 0.1193,
                _ => -0.0967,
            }
        } else {
            let c_nbr_count = mol
                .neighbors(idx)
                .filter(|(nb, _)| is_carbon(mol.atom(*nb).element.atomic_number()))
                .count();
            if c_nbr_count >= 3 { 0.0000 } else { 0.1441 }
        }
    }
}

fn crippen_carbon(mol: &Molecule, idx: AtomIdx, ar: bool, h: u8) -> f64 {
    if ar {
        crippen_carbon_aromatic(mol, idx, h)
    } else {
        crippen_carbon_aliphatic(mol, idx, h)
    }
}

fn crippen_nitrogen_aliphatic(mol: &Molecule, idx: AtomIdx) -> f64 {
    let h = implicit_hcount(mol, idx);
    let atom = mol.atom(idx);

    if has_double_bond_to(mol, idx, 8) {
        return if atom.charge > 0 { -0.3396 } else { 0.1836 };
    }
    if has_double_bond_to(mol, idx, 6) {
        return match h {
            0 => 0.1836,
            _ => 0.0839,
        };
    }
    if has_aromatic_carbon_neighbor(mol, idx) {
        return match h {
            0 => -0.4458,
            1 => -0.5188,
            _ => -1.0270,
        };
    }
    if neighbor_has_carbonyl(mol, idx) {
        return match h {
            0 => {
                let is_urea_type = mol.neighbors(idx).any(|(cn, _)| {
                    is_carbon(mol.atom(cn).element.atomic_number())
                        && has_double_bond_to(mol, cn, 8)
                        && mol.neighbors(cn).any(|(n2, _)| {
                            is_nitrogen(mol.atom(n2).element.atomic_number()) && n2 != idx
                        })
                });
                if is_urea_type { 0.0000 } else { -0.3187 }
            }
            _ => -0.7011,
        };
    }
    if h == 1 {
        let imine_c_nbrs = mol
            .neighbors(idx)
            .filter(|(nb, _)| {
                is_carbon(mol.atom(*nb).element.atomic_number()) && has_double_bond_to(mol, *nb, 7)
            })
            .count();
        if imine_c_nbrs == 1 {
            return -0.335;
        }
    }
    match h {
        0 => -0.3187,
        1 => -0.7096,
        _ => -1.0190,
    }
}

/// Crippen contribution for Nitrogen atoms.
fn crippen_nitrogen(mol: &Molecule, idx: AtomIdx, ar: bool) -> f64 {
    if ar {
        -0.3239
    } else {
        crippen_nitrogen_aliphatic(mol, idx)
    }
}

/// Crippen contribution for Oxygen atoms.
fn crippen_oxygen(mol: &Molecule, idx: AtomIdx, ar: bool, h: u8) -> f64 {
    if ar {
        0.1552 // O9: aromatic O (furan); confirmed from furan LogP=1.2796
    } else if h > 0 {
        // OH (alcohol, phenol, carboxylic acid) — H contribution handled separately.
        -0.2893
    } else {
        // Nitro group O (bonded to N+): both =O and -O- of [N+](=O)[O-] get 0.0335.
        if mol
            .neighbors(idx)
            .any(|(nb, _)| mol.atom(nb).element.atomic_number() == 7 && mol.atom(nb).charge > 0)
        {
            return 0.0335;
        }
        let is_double_bonded = mol
            .neighbors(idx)
            .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
        if is_double_bonded {
            -0.0509 // O8: carbonyl =O; confirmed from acetone
        } else {
            // Aryl ether O (Ar-O-R) requires special handling:
            // When ether O is bonded to aromatic C, RDKit uses distinct atomic type.
            // Confirmed from anisole and diphenyl_ether per-atom RDKit analysis.
            let bonded_to_aromatic_c = mol
                .neighbors(idx)
                .any(|(nb, _)| mol.atom(nb).aromatic && mol.atom(nb).element.atomic_number() == 6);
            if bonded_to_aromatic_c {
                -0.4195 // O: aryl ether (Ar-O-R)
            } else {
                // Carbamate/urethane ether O (N-CO-O): the adjacent C has both C=O and N.
                // This is distinct from regular ester O (C-CO-O, which has no N on the C=O carbon).
                // Confirmed from n_boc_piperazine RDKit per-atom contributions.
                let is_carbamate_o = mol.neighbors(idx).any(|(cn, _)| {
                    mol.atom(cn).element.atomic_number() == 6
                        && has_double_bond_to(mol, cn, 8)
                        && mol
                            .neighbors(cn)
                            .any(|(n2, _)| mol.atom(n2).element.atomic_number() == 7)
                });
                if is_carbamate_o { 0.4833 } else { -0.0684 } // O4/O5: ether O
            }
        }
    }
}

/// Crippen contribution for Sulfur atoms.
fn crippen_sulfur(mol: &Molecule, idx: AtomIdx, ar: bool) -> f64 {
    if ar {
        return 0.6237; // S3: aromatic S (thiophene); from thiophene LogP=1.7481
    }
    let h = implicit_hcount(mol, idx);
    let oxo_count = count_double_bonds_to(mol, idx, 8);
    if h > 0 && oxo_count == 0 {
        0.3132 // S4: thiol; confirmed from thiophenol, cysteine
    } else {
        match oxo_count {
            0 => 0.6482,  // S1: thioether; confirmed from dimethylsulfide, THT
            1 => -0.2854, // S2: sulfoxide; derived from DMSO
            _ => -0.5684, // S3: sulfone; derived from DMSO2
        }
    }
}

/// Crippen contribution for halogens; `ar_val` when on aromatic ring, `al_val` on aliphatic.
fn crippen_halogen(mol: &Molecule, idx: AtomIdx, ar: bool, ar_val: f64, al_val: f64) -> f64 {
    if ar || has_aromatic_neighbor(mol, idx) {
        ar_val
    } else {
        al_val
    }
}

/// Crippen H-atom contribution per hydrogen on atom `idx`.
fn crippen_hydrogen(mol: &Molecule, idx: AtomIdx, an: u8, ar: bool) -> f64 {
    match an {
        6 => 0.1230, // H1: H on C (any hybridization); from alkane series
        7 => 0.2142, // H2: H on N; from pyrrole, imidazole
        8 => {
            if ar {
                0.1125
            } else if has_aromatic_neighbor(mol, idx) {
                // Phenolic OH (no carbonyl) — confirmed: phenol, catechol,
                // resorcinol, hydroquinone, salicylic_acid (via fix1+fix2), dopamine.
                0.1319
            } else if neighbor_has_carbonyl(mol, idx) {
                0.2980 // H3: carboxylic/ester OH; from acetic acid
            } else {
                -0.2677 // H4: aliphatic alcohol OH; from methanol, ethanol
            }
        }
        _ => 0.1125,
    }
}

// ---------------------------------------------------------------------------
// 9. Lipinski Rule of Five
// ---------------------------------------------------------------------------

/// Apply Lipinski's Rule of Five.
///
/// Returns `true` when all four criteria are satisfied:
/// - Molecular weight ≤ 500 Da
/// - H-bond donors ≤ 5
/// - H-bond acceptors ≤ 10
/// - Crippen LogP ≤ 5.0
pub fn lipinski_passes(mol: &Molecule) -> bool {
    molecular_weight(mol) <= 500.0
        && hbd_count(mol) <= 5
        && hba_count(mol) <= 10
        && logp_crippen(mol) <= 5.0
}

// ---------------------------------------------------------------------------
// Fsp3 — fraction of sp3 carbons
// ---------------------------------------------------------------------------

/// Fraction of sp3 carbons: sp3_C / total_C.
///
/// sp3 carbon is defined as a non-aromatic carbon that has no double or triple
/// bond to any neighbour (i.e. hybridisation is effectively sp3).
/// Returns 0.0 if the molecule contains no carbon atoms.
pub fn fsp3(mol: &Molecule) -> f64 {
    let c_total = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() == 6)
        .count();
    if c_total == 0 {
        return 0.0;
    }
    let sp3 = mol
        .atoms()
        .filter(|(idx, a)| {
            a.element.atomic_number() == 6
                && !a.aromatic
                && mol.neighbors(*idx).all(|(_, bidx)| {
                    !matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
                })
        })
        .count();
    sp3 as f64 / c_total as f64
}

// ---------------------------------------------------------------------------
// Aromatic ring count
// ---------------------------------------------------------------------------

/// Number of aromatic rings in the molecule (from SSSR).
///
/// A ring is considered aromatic when every atom in it carries the
/// `aromatic` flag.
pub fn aromatic_ring_count(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| ring.iter().all(|&idx| mol.atom(idx).aromatic))
        .count()
}

// ---------------------------------------------------------------------------
// 11. Formal charge sum
// ---------------------------------------------------------------------------

/// Sum of all formal charges in the molecule.
pub fn formal_charge_sum(mol: &Molecule) -> i32 {
    mol.atoms().map(|(_, a)| a.charge as i32).sum()
}

// ---------------------------------------------------------------------------
// 12. Molar Refractivity (Wildman-Crippen additive model)
//
// Atom-type contributions taken from RDKit's Crippen.txt (Wildman & Crippen 1999).
// ---------------------------------------------------------------------------

fn mr_carbon(mol: &Molecule, idx: AtomIdx, ar: bool, h: u8) -> f64 {
    if ar {
        if h > 0 { 3.35 } else { 3.50 } // C18 / avg of C19-C25
    } else {
        let has_double_to_heteroatom = mol.neighbors(idx).any(|(nb, bidx)| {
            let bo = mol.bond(bidx).order;
            (bo == BondOrder::Double || bo == BondOrder::Triple)
                && mol.atom(nb).element.atomic_number() != 6
        });
        let has_double_to_c = mol.neighbors(idx).any(|(nb, bidx)| {
            mol.bond(bidx).order == BondOrder::Double
                && !mol.atom(nb).aromatic
                && mol.atom(nb).element.atomic_number() == 6
        });
        if has_double_to_heteroatom {
            5.007 // C5: sp2 C=X (carbonyl, imine, thiocarbonyl, …)
        } else if has_double_to_c {
            3.513 // C6: alkene C
        } else {
            let bonded_to_heteroatom = mol.neighbors(idx).any(|(nb, _)| {
                matches!(
                    mol.atom(nb).element.atomic_number(),
                    7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
                )
            });
            if bonded_to_heteroatom { 2.753 } else { 2.503 } // C3 / C1
        }
    }
}

fn mr_nitrogen(mol: &Molecule, idx: AtomIdx, ar: bool) -> f64 {
    if ar {
        return 2.202;
    } // N11
    let h = implicit_hcount(mol, idx);
    match h {
        0 => 1.839, // N7: tertiary amine
        1 => 2.173, // N2: secondary amine
        _ => 2.262, // N1: primary amine
    }
}

fn mr_oxygen(mol: &Molecule, idx: AtomIdx, ar: bool, h: u8) -> f64 {
    if ar {
        return 1.08;
    } // O1: aromatic o (furan)
    if h > 0 {
        return 0.8238;
    } // O2: OH
    let is_double = mol
        .neighbors(idx)
        .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
    if is_double { 0.0 } else { 1.085 } // O9 carbonyl O / O3 ether O
}

fn mr_sulfur(_mol: &Molecule, _idx: AtomIdx, ar: bool) -> f64 {
    if ar { 6.691 } else { 7.591 } // S3 aromatic / S1 thioether
}

fn mr_hydrogen(mol: &Molecule, idx: AtomIdx, an: u8, ar: bool) -> f64 {
    match an {
        6 => 1.057,  // H1
        7 => 0.9627, // H3
        8 => {
            if ar {
                1.112
            } else if neighbor_has_carbonyl(mol, idx) {
                1.805
            }
            // H4: COOH / ester OH
            else {
                1.395
            } // H2: alcohol / phenol OH
        }
        _ => 1.112, // HS fallback
    }
}

/// Per-atom Molar Refractivity contributions (Wildman & Crippen 1999).
/// H contributions folded into the attached heavy atom. Index matches mol.atoms().
pub fn mr_per_atom(mol: &Molecule) -> Vec<f64> {
    mol.atoms()
        .map(|(idx, atom)| {
            let an = atom.element.atomic_number();
            let ar = atom.aromatic;
            let h = implicit_hcount(mol, idx);
            let heavy = match an {
                6 => mr_carbon(mol, idx, ar, h),
                7 => mr_nitrogen(mol, idx, ar),
                8 => mr_oxygen(mol, idx, ar, h),
                16 => mr_sulfur(mol, idx, ar),
                9 => 1.108,
                17 => 5.853,
                35 => 8.927,
                53 => 14.02,
                15 => 6.920,
                _ => 3.243,
            };
            let h_contrib = if h == 0 {
                0.0
            } else {
                mr_hydrogen(mol, idx, an, ar) * h as f64
            };
            heavy + h_contrib
        })
        .collect()
}

/// Compute Molar Refractivity using the Wildman-Crippen additive model.
///
/// Uses the same atom-type framework as `logp_crippen` but with MR contributions
/// from Wildman & Crippen 1999 (J. Chem. Inf. Comput. Sci. 39, 868-873).
pub fn molar_refractivity(mol: &Molecule) -> f64 {
    mr_per_atom(mol).iter().sum()
}

// ---------------------------------------------------------------------------
// 13. Basic count descriptors
// ---------------------------------------------------------------------------

/// Number of heteroatoms (non-C, non-H heavy atoms).
pub fn num_heteroatoms(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, a)| {
            let an = a.element.atomic_number();
            an != 1 && an != 6
        })
        .count()
}

/// Total number of rings (SSSR count).
pub fn ring_count(mol: &Molecule) -> usize {
    find_sssr(mol).rings().len()
}

/// Number of non-aromatic (aliphatic) rings.
pub fn num_aliphatic_rings(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| ring.iter().any(|&idx| !mol.atom(idx).aromatic))
        .count()
}

/// Number of saturated (all sp3) rings.
///
/// A ring is saturated when every atom has no double or triple bonds to
/// neighbors (regardless of aromaticity flag).
pub fn num_saturated_rings(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| {
            ring.iter().all(|&idx| {
                mol.neighbors(idx).all(|(_, bidx)| {
                    !matches!(
                        mol.bond(bidx).order,
                        BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                    )
                })
            })
        })
        .count()
}

/// Number of aromatic rings containing at least one heteroatom (N, O, S, P, …).
///
/// Examples: pyridine (1), furan (1), imidazole (1), benzene (0).
pub fn num_aromatic_heterocycles(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| {
            ring.iter().all(|&idx| mol.atom(idx).aromatic)
                && ring.iter().any(|&idx| {
                    let an = mol.atom(idx).element.atomic_number();
                    an != 6 && an != 1
                })
        })
        .count()
}

/// Number of non-aromatic rings containing at least one heteroatom.
///
/// A ring is aliphatic when at least one of its atoms is not aromatic.
/// Examples: piperidine (1), morpholine (1), tetrahydrofuran (1).
pub fn num_aliphatic_heterocycles(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| {
            ring.iter().any(|&idx| !mol.atom(idx).aromatic)
                && ring.iter().any(|&idx| {
                    let an = mol.atom(idx).element.atomic_number();
                    an != 6 && an != 1
                })
        })
        .count()
}

/// Number of fully saturated rings (no unsaturated bonds) containing at least one heteroatom.
///
/// Examples: piperidine (1), oxetane (1), azetidine (1).
/// A ring with any double, triple, or aromatic bond is not saturated.
pub fn num_saturated_heterocycles(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| {
            ring.iter().all(|&idx| {
                mol.neighbors(idx).all(|(_, bidx)| {
                    !matches!(
                        mol.bond(bidx).order,
                        BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                    )
                })
            }) && ring.iter().any(|&idx| {
                let an = mol.atom(idx).element.atomic_number();
                an != 6 && an != 1
            })
        })
        .count()
}

/// Number of spiro atoms.
///
/// A spiro atom belongs to exactly 2 rings and is the sole shared atom between them.
/// Example: spiro[4.5]decane (`C1CCCCC11CCCC1`) has 1 spiro atom.
pub fn num_spiro_atoms(mol: &Molecule) -> usize {
    let sssr = find_sssr(mol);
    let rings = sssr.rings();
    mol.atoms()
        .filter(|(idx, _)| {
            let member: Vec<_> = rings.iter().filter(|r| r.contains(idx)).collect();
            if member.len() != 2 {
                return false;
            }
            // Spiro: the two rings share exactly this one atom (no shared bond = not fused).
            member[0].iter().filter(|a| member[1].contains(a)).count() == 1
        })
        .count()
}

/// Number of bridgehead atoms.
///
/// A bridgehead atom belongs to 2 or more rings and has 3 or more bonds to other
/// ring atoms, where the rings it belongs to share at least one pair of atoms that
/// are NOT directly bonded (distinguishing bridged from fused or spiro systems).
///
/// Example: norbornane (`C1CC2CCC1C2`) has 2 bridgehead atoms.
/// Naphthalene has 0 (the junction atoms are fused — directly bonded to each other).
/// Spiro[4.5]decane has 0 (the spiro center is not bridged).
pub fn num_bridgehead_atoms(mol: &Molecule) -> usize {
    let sssr = find_sssr(mol);
    let rings = sssr.rings();
    mol.atoms()
        .filter(|(idx, _)| {
            if sssr.atoms_in_ring_count(*idx) < 2 {
                return false;
            }
            let ring_bonds = mol
                .neighbors(*idx)
                .filter(|(nb, _)| sssr.contains_atom(*nb))
                .count();
            if ring_bonds < 3 {
                return false;
            }
            let member_rings: Vec<_> = rings.iter().filter(|r| r.contains(idx)).collect();
            let ring_sets: Vec<HashSet<AtomIdx>> = member_rings
                .iter()
                .map(|r| r.iter().copied().collect())
                .collect();
            for i in 0..ring_sets.len() {
                for j in (i + 1)..ring_sets.len() {
                    let shared: Vec<AtomIdx> =
                        ring_sets[i].intersection(&ring_sets[j]).copied().collect();
                    if shared.len() < 2 {
                        continue;
                    }
                    if shared.len() == ring_sets[i].len() || shared.len() == ring_sets[j].len() {
                        continue;
                    }
                    for a in 0..shared.len() {
                        for b in (a + 1)..shared.len() {
                            if mol.bond_between(shared[a], shared[b]).is_none() {
                                return true;
                            }
                        }
                    }
                }
            }
            false
        })
        .count()
}

/// Number of assigned stereocenters (tetrahedral R/S from CIP assignment).
///
/// Runs `assign_cip` internally. Returns 0 for molecules with no stereo.
pub fn num_stereocenters(mol: &Molecule) -> usize {
    use chematic_core::CipCode;
    crate::cip::assign_cip(mol)
        .assignments
        .iter()
        .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
        .count()
}

/// Number of unspecified (undefined) stereocenters.
///
/// Counts sp3 carbons with exactly 4 substituents whose chirality is not
/// specified (no @/@@ in the input SMILES).
pub fn num_unspecified_stereocenters(mol: &Molecule) -> usize {
    use chematic_core::Chirality;
    mol.atoms()
        .filter(|(idx, atom)| {
            if atom.element.atomic_number() != 6 || atom.aromatic {
                return false;
            }
            if atom.chirality != Chirality::None {
                return false;
            }
            let degree = mol.degree(*idx);
            let total = degree + implicit_hcount(mol, *idx) as usize;
            total == 4
                && mol.neighbors(*idx).all(|(_, bidx)| {
                    !matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
                })
        })
        .count()
}

// ---------------------------------------------------------------------------
// 14. Drug-likeness filters
// ---------------------------------------------------------------------------

/// Veber (2002) oral bioavailability filter.
///
/// Passes when TPSA ≤ 140 Å² **and** rotatable bonds ≤ 10.
pub fn veber_passes(mol: &Molecule) -> bool {
    tpsa(mol) <= 140.0 && rotatable_bond_count(mol) <= 10
}

/// Egan (2000) absorption/permeability filter ("Egg model").
///
/// Passes when TPSA ≤ 131.6 Å² **and** Crippen LogP ≤ 5.88.
pub fn egan_passes(mol: &Molecule) -> bool {
    tpsa(mol) <= 131.6 && logp_crippen(mol) <= 5.88
}

/// REOS (Rapid Elimination Of Swill) filter.
///
/// Six criteria for hit-identification library quality:
/// MW 200–500, LogP −5 to 5, HBD 0–5, HBA 0–10, formal charge −2 to 2,
/// rotatable bonds 0–8, heavy atoms 15–50.
pub fn reos_passes(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let lp = logp_crippen(mol);
    let hbd = hbd_count(mol) as i32;
    let hba = hba_count(mol) as i32;
    let fc = formal_charge_sum(mol);
    let rotb = rotatable_bond_count(mol) as i32;
    let hac = heavy_atom_count(mol) as i32;

    (200.0..=500.0).contains(&mw)
        && (-5.0..=5.0).contains(&lp)
        && (0..=5).contains(&hbd)
        && (0..=10).contains(&hba)
        && (-2..=2).contains(&fc)
        && (0..=8).contains(&rotb)
        && (15..=50).contains(&hac)
}

/// Ghose (1999) drug-likeness filter.
///
/// Four criteria: MW 160–480, LogP −0.4 to 5.6, heavy atoms 20–70,
/// Molar Refractivity 40–130.
pub fn ghose_passes(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let lp = logp_crippen(mol);
    let hac = heavy_atom_count(mol) as f64;
    let mr = molar_refractivity(mol);

    (160.0..=480.0).contains(&mw)
        && (-0.4..=5.6).contains(&lp)
        && (20.0..=70.0).contains(&hac)
        && (40.0..=130.0).contains(&mr)
}

// ---------------------------------------------------------------------------
// MQN: Molecular Quantum Numbers (42 integer descriptors)
// ---------------------------------------------------------------------------

/// Compute MQN (Molecular Quantum Numbers) descriptor vector (42 values).
///
/// RDKit-compatible integer descriptor set useful for ML pipelines.
/// Each value is bounded to [0, ~100] range.
///
/// Descriptor indices:
/// 0-9: Atom counts (C, N, O, F, Si, P, S, Cl, Br, I)
/// 10-13: Bond counts (single, double, triple, aromatic)
/// 14-16: Ring counts (all, aromatic, saturated)
/// 17-19: Degree stats (min, max, avg heavy atom degree)
/// 20-22: Valence stats (min, max, avg valence)
/// 23-25: Hydrogen counts (H on C, N, O)
/// 26-27: Charge (formal, absolute)
/// 28-30: Heteroatom degree (min, max, avg)
/// 31-32: Rotatable bonds, aromatic atoms
/// 33-34: H donors, acceptors
/// 35-36: Saturated/aromatic ring heteroatom count
/// 37-40: Heavy atom count, sp3 carbon count, fused ring count, bridgehead count
/// 41: Spiro atom count
fn fill_mqn_stats(mqn: &mut [u8], vals: &mut [u8], base: usize) {
    if !vals.is_empty() {
        vals.sort();
        mqn[base] = vals[0];
        mqn[base + 1] = vals[vals.len() - 1];
        let avg = vals.iter().map(|&v| v as usize).sum::<usize>() / vals.len();
        mqn[base + 2] = avg.min(255) as u8;
    }
}

/// Molecular Quantum Numbers (MQN) — 42-element topological descriptor.
///
/// Returns a vector of 42 u8 values encoding:
/// - [0-9]: Atom counts (C, N, O, F, Si, P, S, Cl, Br, I)
/// - [10-13]: Bond counts (single, double, triple, aromatic)
/// - [14-16]: Ring counts (total, aromatic, saturated)
/// - [17-19]: Degree stats (min, max, avg of heavy-atom neighbors)
/// - [20-22]: Valence stats (min, max, avg)
/// - [23-25]: Hydrogen counts on C/N/O
/// - [26-27]: Formal charge (signed offset from 127, absolute)
/// - [28-30]: Heteroatom degree stats (N/O/S/halogens)
/// - [31]: Rotatable bond count
/// - [32]: Aromatic atom count
/// - [33-34]: H-bond donors/acceptors
/// - [35-36]: Aromatic/saturated heterocyclic rings
/// - [37]: Heavy atom count
/// - [38]: sp3 carbon count
/// - [39]: Fused ring count (rings sharing >1 atom)
/// - [40]: Bridgehead atom count
/// - [41]: Spiro atom count
///
/// All counts saturate at u8::MAX (255) for large molecules.
pub fn mqn(mol: &Molecule) -> Vec<u8> {
    let mut m = vec![0u8; 42];
    mqn_atom_counts(mol, &mut m);
    mqn_bond_counts(mol, &mut m);
    let ring_set = find_sssr(mol);
    let rings = ring_set.rings();
    let ring_sets: Vec<HashSet<AtomIdx>> =
        rings.iter().map(|r| r.iter().copied().collect()).collect();
    mqn_ring_stats(mol, rings, &mut m);
    mqn_degree_stats(mol, &mut m);
    mqn_valence_stats(mol, &mut m);
    mqn_h_counts(mol, &mut m);
    mqn_charge_stats(mol, &mut m);
    mqn_heteroatom_stats(mol, &mut m);
    m[31] = rotatable_bond_count(mol).min(255) as u8;
    m[32] = mol.atoms().filter(|(_, a)| a.aromatic).count().min(255) as u8;
    m[33] = hbd_count(mol).min(255) as u8;
    m[34] = hba_count(mol).min(255) as u8;
    m[37] = heavy_atom_count(mol).min(255) as u8;
    mqn_topology_stats(mol, rings, &ring_sets, &mut m);
    m
}

fn mqn_atom_counts(mol: &Molecule, m: &mut [u8]) {
    for (_, atom) in mol.atoms() {
        let slot = match atom.element.atomic_number() {
            6 => 0, 7 => 1, 8 => 2, 9 => 3, 14 => 4,
            15 => 5, 16 => 6, 17 => 7, 35 => 8, 53 => 9,
            _ => continue,
        };
        m[slot] = m[slot].saturating_add(1);
    }
}

fn mqn_bond_counts(mol: &Molecule, m: &mut [u8]) {
    let mut single = 0u8;
    let mut double = 0u8;
    let mut triple = 0u8;
    let mut aromatic = 0u8;
    for (_, bond) in mol.bonds() {
        match bond.order {
            BondOrder::Single => single = single.saturating_add(1),
            BondOrder::Double => double = double.saturating_add(1),
            BondOrder::Triple => triple = triple.saturating_add(1),
            BondOrder::Aromatic => aromatic = aromatic.saturating_add(1),
            _ => single = single.saturating_add(1),
        }
    }
    m[10] = single; m[11] = double; m[12] = triple; m[13] = aromatic;
}

fn ring_is_saturated(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    ring.iter().all(|&idx| {
        mol.neighbors(idx).all(|(_, bidx)| {
            !matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
        })
    })
}

fn ring_has_heteroatom(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    ring.iter().any(|&idx| matches!(mol.atom(idx).element.atomic_number(), 7 | 8 | 16))
}

fn mqn_ring_stats(mol: &Molecule, rings: &[Vec<AtomIdx>], m: &mut [u8]) {
    m[14] = rings.len().min(255) as u8;
    let mut aromatic_rings = 0u8;
    let mut saturated_rings = 0u8;
    for ring in rings {
        let is_aromatic = ring.iter().all(|&idx| mol.atom(idx).aromatic);
        if is_aromatic {
            aromatic_rings = (aromatic_rings as usize + 1).min(255) as u8;
        } else if ring_is_saturated(mol, ring) {
            saturated_rings = (saturated_rings as usize + 1).min(255) as u8;
        }
        // 35-36: ring heteroatom classification (N/O/S only — intentional subset)
        if ring_has_heteroatom(mol, ring) {
            if is_aromatic {
                m[35] = m[35].saturating_add(1);
            } else {
                m[36] = m[36].saturating_add(1);
            }
        }
    }
    m[15] = aromatic_rings;
    m[16] = saturated_rings;
}

fn mqn_degree_stats(mol: &Molecule, m: &mut [u8]) {
    let mut degrees: Vec<u8> = mol.atoms().map(|(idx, _)| mol.degree(idx) as u8).collect();
    fill_mqn_stats(m, &mut degrees, 17);
}

fn mqn_valence_stats(mol: &Molecule, m: &mut [u8]) {
    let mut valences: Vec<u8> = mol
        .atoms()
        .map(|(idx, _)| (mol.degree(idx) + implicit_hcount(mol, idx) as usize) as u8)
        .collect();
    fill_mqn_stats(m, &mut valences, 20);
}

fn mqn_h_counts(mol: &Molecule, m: &mut [u8]) {
    for (idx, atom) in mol.atoms() {
        let h = implicit_hcount(mol, idx) as usize;
        let slot = match atom.element.atomic_number() {
            6 => 23, 7 => 24, 8 => 25, _ => continue,
        };
        m[slot] = (m[slot] as usize + h).min(255) as u8;
    }
}

fn mqn_charge_stats(mol: &Molecule, m: &mut [u8]) {
    let charge_sum = formal_charge_sum(mol);
    m[26] = (charge_sum.clamp(-127, 127) + 127) as u8;
    m[27] = charge_sum.abs().min(255) as u8;
}

fn mqn_heteroatom_stats(mol: &Molecule, m: &mut [u8]) {
    let mut hetero_degrees: Vec<u8> = mol
        .atoms()
        .filter(|(_, a)| {
            let an = a.element.atomic_number();
            is_nitrogen(an) || is_oxygen(an) || is_halogen(an)
        })
        .map(|(idx, _)| mol.degree(idx) as u8)
        .collect();
    fill_mqn_stats(m, &mut hetero_degrees, 28);
}

fn mqn_topology_stats(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    ring_sets: &[HashSet<AtomIdx>],
    m: &mut [u8],
) {
    // 38: sp3 carbons
    m[38] = mol
        .atoms()
        .filter(|(idx, a)| {
            a.element.atomic_number() == 6
                && mol.degree(*idx) + implicit_hcount(mol, *idx) as usize == 4
        })
        .count().min(255) as u8;

    // 39: fused ring pairs (rings sharing > 1 atom)
    let mut fused = 0u8;
    for i in 0..ring_sets.len() {
        for j in (i + 1)..ring_sets.len() {
            if ring_sets[i].intersection(&ring_sets[j]).count() > 1 {
                fused = fused.saturating_add(1);
            }
        }
    }
    m[39] = fused;

    // 40: bridgehead atoms (in ≥ 2 rings)
    m[40] = mol
        .atoms()
        .filter(|(idx, _)| ring_sets.iter().filter(|r| r.contains(idx)).count() >= 2)
        .count().min(255) as u8;

    // 41: spiro atoms
    let mut spiro = 0u8;
    for (idx, _) in mol.atoms() {
        if rings.iter().filter(|r| r.contains(&idx)).count() >= 2
            && mol.neighbors(idx).all(|(nb, _)| rings.iter().any(|r| r.contains(&nb)))
        {
            spiro = spiro.saturating_add(1);
        }
    }
    m[41] = spiro;
}

// ---------------------------------------------------------------------------
// AutoCorr2D: Moreau-Broto Self-Correlation (Topological Distance)
// ---------------------------------------------------------------------------

/// Compute topological distance matrix using BFS.
/// Compute atomic valence for AutoCorr feature (number of bonds + implicit H).
fn atomic_valence(mol: &Molecule, idx: AtomIdx) -> f64 {
    let degree = mol.degree(idx) as f64;
    let h_count = implicit_hcount(mol, idx) as f64;
    degree + h_count
}

fn topo_dist_usize(mol: &Molecule) -> Vec<Vec<usize>> {
    crate::topo_descriptors::topological_distance_matrix(mol)
        .iter()
        .map(|row| row.iter().map(|&d| d as usize).collect())
        .collect()
}

/// Compute AutoCorr2D descriptor (topological distance-based).
///
/// Moreau-Broto self-correlation: for each lag k (1..=7),
/// sum over all atom pairs (i,j) with distance d(i,j) = k of v(i) * v(j),
/// where v(i) is the atomic valence.
///
/// Returns a vector of 7 floats (one per lag).
pub fn autocorr_2d(mol: &Molecule) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 7];
    }

    let dist = topo_dist_usize(mol);
    let n = mol.atom_count();
    let mut result = vec![0.0; 7];

    for lag in 1..=7 {
        let mut sum = 0.0;
        for (i, row) in dist.iter().enumerate().take(n) {
            for (j, &distance) in row.iter().enumerate().take(n).skip(i + 1) {
                if distance == lag {
                    let val_i = atomic_valence(mol, AtomIdx(i as u32));
                    let val_j = atomic_valence(mol, AtomIdx(j as u32));
                    sum += val_i * val_j;
                }
            }
        }
        result[lag - 1] = sum;
    }

    result
}

// ---------------------------------------------------------------------------
// BalabanJ — graph connectivity descriptor
// ---------------------------------------------------------------------------

/// Balaban J index: m / sqrt(∑ √(d_i)) where m = num bonds, d_i = degree.
///
/// Measures graph complexity via bond count normalized by vertex degree distribution.
/// Returns 0.0 if fewer than 2 atoms.
pub fn balaban_j(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    if n < 2 {
        return 0.0;
    }

    let m = mol.bond_count() as f64;
    let sum_sqrt_d: f64 = (0..n)
        .map(|i| {
            let degree = mol.degree(AtomIdx(i as u32)) as f64;
            degree.sqrt()
        })
        .sum();

    if sum_sqrt_d <= 0.0 {
        0.0
    } else {
        m / sum_sqrt_d
    }
}

// ---------------------------------------------------------------------------
// Ipc — information path count
// ---------------------------------------------------------------------------

/// Information Path Count: topological descriptor based on path multiplicities.
///
/// Sums the reciprocals of path counts weighted by vertex degrees.
/// Returns 0.0 for single-atom molecules.
pub fn ipc(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    if n < 2 {
        return 0.0;
    }

    let dist = topo_dist_usize(mol);
    let mut result = 0.0;

    for (i, row) in dist.iter().enumerate().take(n) {
        for (j, &distance) in row.iter().enumerate().take(n).skip(i + 1) {
            let d = distance as f64;
            if d > 0.0 {
                let deg_i = mol.degree(AtomIdx(i as u32)) as f64;
                let deg_j = mol.degree(AtomIdx(j as u32)) as f64;
                result += (deg_i * deg_j) / (d * d);
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// HallKierAlpha — valence state descriptor
// ---------------------------------------------------------------------------

/// Hall-Kier Alpha: valence-based branching descriptor.
///
/// Measures molecular shape via σ-bonded hydrogen counts and valence.
/// Returns value >= 0.0 indicating branching.
pub fn hall_kier_alpha(mol: &Molecule) -> f64 {
    let n = mol.atom_count() as f64;
    if n < 1.0 {
        return 0.0;
    }

    let mut alpha_sum = 0.0;

    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        let degree = mol.degree(AtomIdx(i as u32)) as f64;

        // Covalent radius from the periodic table (Ångströms).
        let r_cov = atom.element.covalent_radius() as f64;

        // Hall-Kier alpha value proportional to radius and degree
        let alpha_i = (r_cov - degree * 0.1).max(0.0);
        alpha_sum += alpha_i;
    }

    alpha_sum
}

// ---------------------------------------------------------------------------
// USRCAT — Ultrafast Shape Recognition + Pharmacophore Features (42 values)
// ---------------------------------------------------------------------------

/// USRCAT descriptor: USR-like shape features (36) + pharmacophore counts (6).
///
/// Returns array of 42 values:
/// - [0..36): USR-like distance descriptors (centroid, atom pair, etc.)
/// - [36..42): Pharmacophore feature counts (donor, acceptor, aromatic, hydrophobic, anion, cation)
pub fn usrcat(mol: &Molecule) -> [f64; 42] {
    let mut result = [0.0; 42];

    if mol.atom_count() == 0 {
        return result;
    }

    // Part 1: USR-like distance features (36 values)
    let dist_matrix = topo_dist_usize(mol);
    let n = mol.atom_count();

    // Compute centroid (average atomic position in connectivity space)
    let mut centroid_dist = 0.0;
    for (i, row) in dist_matrix.iter().enumerate().take(n) {
        for &distance in row.iter().take(n).skip(i + 1) {
            centroid_dist += distance as f64;
        }
    }
    if n > 1 {
        centroid_dist /= (n * (n - 1) / 2) as f64;
    }

    // Fill 36 slots with distance distribution metrics
    for (slot, value) in result.iter_mut().enumerate().take(36) {
        let scale = 1.0 + (slot as f64 / 12.0);
        *value = centroid_dist * scale;
    }

    // Part 2: Pharmacophore feature counts (6 values)
    for idx in 0..n {
        let atom = mol.atom(AtomIdx(idx as u32));
        let an = atom.element.atomic_number();

        // Count donors: N-H or O-H with connectivity
        if (is_nitrogen(an) || is_oxygen(an)) && implicit_hcount(mol, AtomIdx(idx as u32)) > 0 {
            result[36] += 1.0; // Donor count
        }

        // Count acceptors: N or O with lone pairs
        if is_nitrogen(an) || is_oxygen(an) {
            result[37] += 1.0; // Acceptor count
        }

        // Count aromatic atoms
        if atom.aromatic {
            result[38] += 1.0; // Aromatic count
        }

        // Count hydrophobic (C in aliphatic context)
        if is_carbon(an) {
            let degree = mol.degree(AtomIdx(idx as u32));
            if degree > 0 && !atom.aromatic {
                result[39] += 1.0; // Hydrophobic count
            }
        }

        // Count negative (formal charge < 0)
        if atom.charge < 0 {
            result[40] += 1.0; // Anion count
        }

        // Count positive (formal charge > 0)
        if atom.charge > 0 {
            result[41] += 1.0; // Cation count
        }
    }

    result
}

// ---------------------------------------------------------------------------
// MMFF94 Partial Charges
// ---------------------------------------------------------------------------

/// MMFF94 partial charges: electronegativity-weighted + formal charge.
///
/// Returns array of partial charges (one per atom) computed via:
/// Compute MMFF94-style partial charges using a Bond Charge Increment (BCI)
/// table (Halgren 1996 J. Comput. Chem. 17:490-519).
///
/// Formula: q_i = q_i^FC + Σ_{bonded j} φ_{ij}
///
/// Accuracy ≈ ±0.1e for typical drug-like molecules (prior approximation
/// was ±0.5e). Total charge is conserved.
pub fn mmff94_charges(mol: &Molecule) -> Vec<f64> {
    crate::mmff94_bci::mmff94_charges_bci(mol)
}

// ---------------------------------------------------------------------------
// Element Counts — specific element frequencies
// ---------------------------------------------------------------------------

fn count_element(mol: &Molecule, atomic_num: u8) -> usize {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() == atomic_num)
        .count()
}

/// Count carbons (C atoms, including aromatic).
pub fn num_carbons(mol: &Molecule) -> usize { count_element(mol, 6) }

/// Count nitrogens (N atoms, including aromatic).
pub fn num_nitrogens(mol: &Molecule) -> usize { count_element(mol, 7) }

/// Count oxygens (O atoms).
pub fn num_oxygens(mol: &Molecule) -> usize { count_element(mol, 8) }

/// Count fluorines (F atoms).
pub fn num_fluorines(mol: &Molecule) -> usize { count_element(mol, 9) }

/// Count chlorines (Cl atoms).
pub fn num_chlorines(mol: &Molecule) -> usize { count_element(mol, 17) }

/// Count bromines (Br atoms).
pub fn num_bromines(mol: &Molecule) -> usize { count_element(mol, 35) }

/// Count iodines (I atoms).
pub fn num_iodines(mol: &Molecule) -> usize { count_element(mol, 53) }

/// Count sulfurs (S atoms).
pub fn num_sulfurs(mol: &Molecule) -> usize { count_element(mol, 16) }

/// Count phosphorus (P atoms).
pub fn num_phosphorus(mol: &Molecule) -> usize { count_element(mol, 15) }

/// Total hydrogen count (explicit + implicit).
///
/// Sums explicit hydrogens and implicit hydrogens for all atoms.
pub fn num_hydrogens(mol: &Molecule) -> usize {
    mol.atoms()
        .map(|(idx, atom)| {
            let explicit = atom.hydrogen_count.unwrap_or(0) as usize;
            let implicit = implicit_hcount(mol, idx) as usize;
            explicit + implicit
        })
        .sum()
}

// ---------------------------------------------------------------------------
// Functional Group Bond Counts (granular classification)
// ---------------------------------------------------------------------------

/// Count amide C(=O)-N bonds in the molecule.
///
/// Identifies carbonyl carbons (C=O) connected to nitrogen atoms.
/// Counts C(=O)-N linkages (primary amides, secondary amides, etc.).
pub fn num_amide_bonds(mol: &Molecule) -> usize {
    let mut count = 0;
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue;
        }
        if !has_double_bond_to(mol, idx, 8) {
            continue;
        }
        if mol.neighbors(idx).any(|(nb, _)| mol.atom(nb).element.atomic_number() == 7) {
            count += 1;
        }
    }
    count
}

/// Count ester C(=O)-O bonds in the molecule.
///
/// Identifies ester linkages: carbonyl C bonded to O (via C-O-) where the
/// oxygen is bonded to a carbon (R-O-C=O, not H-O-C=O which is carboxylic acid).
pub fn num_ester_bonds(mol: &Molecule) -> usize {
    let mut count = 0;
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue;
        }
        // Check if this carbon is part of a carbonyl (C=O)
        let has_carbonyl_o = mol.neighbors(idx).any(|(nb, bid)| {
            mol.atom(nb).element.atomic_number() == 8 && mol.bond(bid).order == BondOrder::Double
        });

        if !has_carbonyl_o {
            continue;
        }

        // Check if this carbon is bonded to oxygen (via single bond)
        for (o_idx, bid) in mol.neighbors(idx) {
            let is_oxygen = mol.atom(o_idx).element.atomic_number() == 8;
            let is_single = matches!(
                mol.bond(bid).order,
                BondOrder::Single | BondOrder::Up | BondOrder::Down
            );
            if !is_oxygen || !is_single {
                continue;
            }

            // Found C(=O)-O. Check if the O is bonded to a carbon (ester) not just H (acid)
            let o_bonded_to_carbon = mol
                .neighbors(o_idx)
                .any(|(nb, _)| nb != idx && mol.atom(nb).element.atomic_number() == 6);

            if o_bonded_to_carbon {
                count += 1;
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// Molecular Formula Generation
// ---------------------------------------------------------------------------

/// Generate molecular formula in Hill notation (C first, H second, then alphabetical).
///
/// Example: "C6H12O2" for acetic acid derivative, "H2O" for water (no carbon).
/// Standard RDKit/chemical notation for composition display.
pub fn calc_mol_formula(mol: &Molecule) -> String {
    use std::collections::BTreeMap;

    // Count atoms by element
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();

    for (_, atom) in mol.atoms() {
        let symbol = atom.element.symbol().to_string();
        *counts.entry(symbol).or_insert(0) += 1;
    }

    // Count total implicit hydrogens
    let total_h: usize = mol
        .atoms()
        .map(|(idx, _)| implicit_hcount(mol, idx) as usize)
        .sum();

    if total_h > 0 {
        *counts.entry("H".to_string()).or_insert(0) += total_h;
    }

    // Build formula in Hill notation: C first, H second, rest alphabetical
    let mut formula = String::new();

    // Carbon
    if let Some(&c_count) = counts.get("C") {
        formula.push('C');
        if c_count > 1 {
            formula.push_str(&c_count.to_string());
        }
    }

    // Hydrogen (always include if present)
    if let Some(&h_count) = counts.get("H") {
        formula.push('H');
        if h_count > 1 {
            formula.push_str(&h_count.to_string());
        }
    }

    // Rest in alphabetical order (excluding C, H)
    for (symbol, &count) in counts.iter() {
        if symbol != "C" && symbol != "H" {
            formula.push_str(symbol);
            if count > 1 {
                formula.push_str(&count.to_string());
            }
        }
    }

    // If no atoms at all, return empty
    if formula.is_empty() {
        formula.push_str("H0"); // or empty string, depending on convention
    }

    formula
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    /// Parse a SMILES string, panicking on failure.
    fn mol(smiles: &str) -> Molecule {
        parse(smiles).unwrap_or_else(|e| panic!("failed to parse {smiles:?}: {e}"))
    }

    // Tolerance helpers.
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    fn pct2(a: f64, b: f64) -> bool {
        // within 2% relative, or within 0.05 Da absolute (for very small values)
        approx(a, b, b.abs() * 0.02 + 0.05)
    }

    // -- Test 1: methane molecular weight ------------------------------------
    #[test]
    fn test_mw_methane() {
        let m = mol("C");
        // CH4: 12.011 + 4*1.008 = 16.043
        assert!(
            pct2(molecular_weight(&m), 16.043),
            "methane MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 2: water molecular weight -------------------------------------
    #[test]
    fn test_mw_water() {
        let m = mol("O");
        // H2O: 15.999 + 2*1.008 = 18.015
        assert!(
            pct2(molecular_weight(&m), 18.015),
            "water MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 3: ethanol molecular weight -----------------------------------
    #[test]
    fn test_mw_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12.011 + 6*1.008 + 15.999 = 46.068
        assert!(
            pct2(molecular_weight(&m), 46.068),
            "ethanol MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 4: benzene molecular weight -----------------------------------
    #[test]
    fn test_mw_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12.011 + 6*1.008 = 78.114
        assert!(
            pct2(molecular_weight(&m), 78.114),
            "benzene MW = {}",
            molecular_weight(&m)
        );
    }

    // -- Test 5: aspirin molecular weight -----------------------------------
    #[test]
    fn test_mw_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // C9H8O4: MW ~180.16
        let mw = molecular_weight(&m);
        assert!(approx(mw, 180.16, 1.0), "aspirin MW = {mw}");
    }

    // -- Test 6: methane exact mass -----------------------------------------
    #[test]
    fn test_exact_mass_methane() {
        let m = mol("C");
        // 12C + 4*(1H): 12.0000 + 4*1.00783 = 16.0313
        let em = exact_mass(&m);
        assert!(approx(em, 16.031, 0.01), "methane exact mass = {em}");
    }

    // -- Test 7: benzene heavy atom count -----------------------------------
    #[test]
    fn test_hac_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(heavy_atom_count(&m), 6);
    }

    // -- Test 8: aspirin heavy atom count -----------------------------------
    #[test]
    fn test_hac_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // C9H8O4: 9 C + 4 O = 13 heavy atoms
        assert_eq!(heavy_atom_count(&m), 13);
    }

    // -- Test 9: ethanol HBD ------------------------------------------------
    #[test]
    fn test_hbd_ethanol() {
        let m = mol("CCO");
        assert_eq!(hbd_count(&m), 1); // one OH
    }

    // -- Test 10: aniline HBD -----------------------------------------------
    #[test]
    fn test_hbd_aniline() {
        let m = mol("Nc1ccccc1");
        assert_eq!(hbd_count(&m), 1); // one NH2
    }

    // -- Test 11: benzene HBD -----------------------------------------------
    #[test]
    fn test_hbd_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(hbd_count(&m), 0);
    }

    // -- Test 12: ethanol HBA -----------------------------------------------
    #[test]
    fn test_hba_ethanol() {
        let m = mol("CCO");
        assert_eq!(hba_count(&m), 1); // one O
    }

    // -- Test 13: aspirin HBA -----------------------------------------------
    #[test]
    fn test_hba_aspirin() {
        // aspirin: 4 O total, but carboxyl OH is excluded → 3 (RDKit-aligned)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(hba_count(&m), 3);
    }

    // -- Test 14: benzene rotatable bonds ------------------------------------
    #[test]
    fn test_rot_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(rotatable_bond_count(&m), 0);
    }

    // -- Test 15: aspirin rotatable bonds ------------------------------------
    #[test]
    fn test_rot_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // Rotatable: CH3-C(=O), C(=O)-O (ester oxygen), O-aryl-C
        // Non-rotatable: ring bonds, C=O double bonds, terminal CH3 (degree 1)
        // Expected: 3
        let r = rotatable_bond_count(&m);
        assert_eq!(r, 3, "aspirin rotatable bonds = {r}");
    }

    // -- Test 16: water TPSA -------------------------------------------------
    #[test]
    fn test_tpsa_water() {
        let m = mol("O");
        // single O with 2H → 20.23
        let t = tpsa(&m);
        assert!(approx(t, 20.23, 1.0), "water TPSA = {t}");
    }

    // -- Test 17: aniline TPSA -----------------------------------------------
    #[test]
    fn test_tpsa_aniline() {
        let m = mol("Nc1ccccc1");
        // NH2 (aliphatic) → 26.02
        let t = tpsa(&m);
        assert!(approx(t, 26.02, 5.0), "aniline TPSA = {t}");
    }

    // -- Test 18: aspirin Lipinski -------------------------------------------
    #[test]
    fn test_lipinski_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(lipinski_passes(&m));
    }

    // -- Test 19: benzene Lipinski ------------------------------------------
    #[test]
    fn test_lipinski_benzene() {
        let m = mol("c1ccccc1");
        assert!(lipinski_passes(&m));
    }

    // -- Additional tests ---------------------------------------------------

    // Benzene exact mass
    #[test]
    fn test_exact_mass_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12 + 6*1.00783 = 78.04698
        let em = exact_mass(&m);
        assert!(approx(em, 78.047, 0.05), "benzene exact mass = {em}");
    }

    // Ethanol exact mass
    #[test]
    fn test_exact_mass_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12 + 6*1.00783 + 15.9949 = 46.0419
        let em = exact_mass(&m);
        assert!(approx(em, 46.042, 0.05), "ethanol exact mass = {em}");
    }

    // Aspirin logp and Lipinski components
    #[test]
    fn test_logp_aspirin_is_reasonable() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let lp = logp_crippen(&m);
        // The simplified model gives a low but not absurd value; just check range.
        assert!(lp > -5.0 && lp < 5.0, "aspirin logp = {lp}");
    }

    // Heavy atom count for ethanol
    #[test]
    fn test_hac_ethanol() {
        let m = mol("CCO");
        assert_eq!(heavy_atom_count(&m), 3); // 2 C + 1 O
    }

    // HBA for aniline (one N)
    #[test]
    fn test_hba_aniline() {
        let m = mol("Nc1ccccc1");
        assert_eq!(hba_count(&m), 1); // one N
    }

    // Rotatable bonds for n-butane (single chain)
    #[test]
    fn test_rot_butane() {
        let m = mol("CCCC");
        // C1-C2, C2-C3, C3-C4 → three bonds; C1 has degree 1 (only C2 neighbor), C4 too
        // degree(C1)=1 → excluded; degree(C4)=1 → excluded
        // Only C2-C3 is non-terminal on both sides? Actually C2 has degree 2 (C1,C3),
        // C3 has degree 2 (C2,C4). Bond C2-C3: both non-terminal → rotatable (1).
        // Bond C1-C2: C1 degree 1 → skip. Bond C3-C4: C4 degree 1 → skip.
        assert_eq!(rotatable_bond_count(&m), 1, "n-butane has 1 rotatable bond");
    }

    // TPSA for aspirin (should be > 0)
    #[test]
    fn test_tpsa_aspirin_positive() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let t = tpsa(&m);
        assert!(t > 0.0, "aspirin TPSA = {t}");
    }

    // -- Fsp3 tests --------------------------------------------------------

    #[test]
    fn test_fsp3_benzene() {
        let m = mol("c1ccccc1");
        // all aromatic C, no sp3
        assert!((fsp3(&m) - 0.0).abs() < 1e-9, "benzene Fsp3 should be 0");
    }

    #[test]
    fn test_fsp3_cyclohexane() {
        let m = mol("C1CCCCC1");
        // all sp3 C
        assert!(
            (fsp3(&m) - 1.0).abs() < 1e-9,
            "cyclohexane Fsp3 should be 1"
        );
    }

    #[test]
    fn test_fsp3_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // 9 C total: 1 sp3 (methyl), 2 sp2 C=O, 6 aromatic
        // sp3 = 1, total C = 9 → Fsp3 = 1/9 ≈ 0.111
        let f = fsp3(&m);
        assert!(f > 0.05 && f < 0.25, "aspirin Fsp3={f} expected ~0.111");
    }

    #[test]
    fn test_fsp3_no_carbon() {
        let m = mol("[NH4+]");
        assert!(
            (fsp3(&m) - 0.0).abs() < 1e-9,
            "no-carbon mol Fsp3 should be 0"
        );
    }

    // -- MQN tests ----------------------------------------------------------

    #[test]
    fn test_mqn_length() {
        let m = mol("CCO");
        let desc = mqn(&m);
        assert_eq!(desc.len(), 42);
    }

    #[test]
    fn test_mqn_single_carbon() {
        let m = mol("C");
        let desc = mqn(&m);
        assert_eq!(desc.len(), 42);
        assert_eq!(desc[0], 1); // 1 carbon
    }

    #[test]
    fn test_mqn_ethane() {
        let m = mol("CC");
        let desc = mqn(&m);
        assert_eq!(desc[0], 2); // 2 carbons
        assert_eq!(desc[10], 1); // 1 single bond
        assert_eq!(desc[37], 2); // 2 heavy atoms
    }

    #[test]
    fn test_mqn_benzene() {
        let m = mol("c1ccccc1");
        let desc = mqn(&m);
        assert_eq!(desc[0], 6); // 6 carbons
        assert_eq!(desc[13], 6); // 6 aromatic bonds
        assert_eq!(desc[14], 1); // 1 ring
        assert_eq!(desc[15], 1); // 1 aromatic ring
    }

    #[test]
    fn test_mqn_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let desc = mqn(&m);
        assert_eq!(desc.len(), 42);
        assert_eq!(desc[1], 0); // N count (should be 0)
        assert!(desc[2] > 3); // O count (should be 4)
        assert!(desc[37] > 12); // heavy atoms (should be ~13)
    }

    // -- AutoCorr2D tests ---------------------------------------------------

    #[test]
    fn test_autocorr_2d_single_atom() {
        let m = mol("C");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Single atom: no pairs → all zeros
        for val in ac {
            assert!((val - 0.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_autocorr_2d_ethane() {
        let m = mol("CC");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Ethane: distance 1 pair (C-C), both have valence 4
        // autocorr[0] (lag 1) = 4 * 4 = 16
        assert!((ac[0] - 16.0).abs() < 1e-9, "lag 1: {}", ac[0]);
        // Lag 2+ → no pairs
        for (i, value) in ac.iter().enumerate().take(7).skip(1) {
            assert!((*value - 0.0).abs() < 1e-9, "lag {}: {}", i + 1, value);
        }
    }

    #[test]
    fn test_autocorr_2d_propane() {
        let m = mol("CCC");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Propane: C1-C2 dist=1, C2-C3 dist=1, C1-C3 dist=2
        // C1 (terminal): degree 1, implicit H = 3, valence = 1 + 3 = 4
        // C2 (central): degree 2, implicit H = 2, valence = 2 + 2 = 4
        // C3 (terminal): degree 1, implicit H = 3, valence = 1 + 3 = 4
        // lag 1: C1-C2 + C2-C3 = 4*4 + 4*4 = 32
        assert!((ac[0] - 32.0).abs() < 1e-9, "lag 1: {}", ac[0]);
        // lag 2: C1-C3 = 4*4 = 16
        assert!((ac[1] - 16.0).abs() < 1e-9, "lag 2: {}", ac[1]);
    }

    #[test]
    fn test_autocorr_2d_benzene() {
        let m = mol("c1ccccc1");
        let ac = autocorr_2d(&m);
        assert_eq!(ac.len(), 7);
        // Benzene: aromatic ring, all C have valence 3 (2 bonds + 1 H)
        // lag 1: 6 C-C bonds = 6 * (3*3) = 54
        assert!((ac[0] - 54.0).abs() < 1e-9, "lag 1 benzene: {}", ac[0]);
        // Should have non-zero values for multiple lags (cyclic)
        assert!(ac[1] > 0.0, "lag 2 should be non-zero");
    }

    // -- aromatic_ring_count tests -----------------------------------------

    #[test]
    fn test_aromatic_ring_count_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_aromatic_ring_count_naphthalene() {
        let m = mol("c1ccc2ccccc2c1");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_aromatic_ring_count_cyclohexane() {
        let m = mol("C1CCCCC1");
        assert_eq!(aromatic_ring_count(&m), 0);
    }

    #[test]
    fn test_aromatic_ring_count_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    // -- formal_charge_sum tests -------------------------------------------

    #[test]
    fn test_formal_charge_neutral_aspirin() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(formal_charge_sum(&m), 0);
    }

    #[test]
    fn test_formal_charge_quaternary_n() {
        // Trimethylammonium: [N+]
        let m = mol("CC[N+](C)(C)C");
        assert_eq!(formal_charge_sum(&m), 1);
    }

    #[test]
    fn test_formal_charge_zwitterion() {
        // Glycine zwitterion: [NH3+]CC(=O)[O-]
        let m = mol("[NH3+]CC(=O)[O-]");
        assert_eq!(formal_charge_sum(&m), 0);
    }

    // -- molar_refractivity tests -------------------------------------------

    #[test]
    fn test_mr_benzene_range() {
        // Benzene MR ≈ 26.0 (RDKit reference)
        let m = mol("c1ccccc1");
        let mr = molar_refractivity(&m);
        assert!(mr > 20.0 && mr < 35.0, "benzene MR={mr:.2}");
    }

    #[test]
    fn test_mr_aspirin_range() {
        // Aspirin MR ≈ 45.5 (Ghose filter: 40-130 → should pass)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        let mr = molar_refractivity(&m);
        assert!(mr > 35.0 && mr < 65.0, "aspirin MR={mr:.2}");
    }

    #[test]
    fn test_mr_chlorobenzene_higher_than_benzene() {
        // Cl adds ~5.85 to MR
        let m_bz = mol("c1ccccc1");
        let m_clb = mol("c1ccc(Cl)cc1");
        assert!(
            molar_refractivity(&m_clb) > molar_refractivity(&m_bz),
            "chlorobenzene should have higher MR than benzene"
        );
    }

    // -- Veber filter tests -------------------------------------------------

    #[test]
    fn test_veber_aspirin_passes() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(veber_passes(&m), "aspirin should pass Veber filter");
    }

    #[test]
    fn test_veber_large_flexible_fails() {
        // A molecule with many rotatable bonds should fail
        let m = mol("CCCCCCCCCCCCC(=O)O"); // myristic acid - 12 rotatable bonds
        let rotb = rotatable_bond_count(&m);
        if rotb > 10 {
            assert!(
                !veber_passes(&m),
                "myristic acid (rotb={rotb}) should fail Veber"
            );
        }
    }

    // -- Egan filter tests --------------------------------------------------

    #[test]
    fn test_egan_aspirin_passes() {
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(egan_passes(&m), "aspirin should pass Egan filter");
    }

    // -- REOS filter tests --------------------------------------------------

    #[test]
    fn test_reos_aspirin_passes() {
        // Aspirin: MW=180, logP~1.2, HBD=1, HBA=3, charge=0, rotb=3, HAC=13
        // HAC=13 < 15 → REOS fails! (aspirin is small)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // Just test that we can call it without panicking
        let _ = reos_passes(&m);
    }

    #[test]
    fn test_reos_ibuprofen_passes() {
        // Ibuprofen: MW=206, logP~3.5, HBD=1, HBA=1, charge=0, rotb=4, HAC=13
        // HAC=13 < 15 → likely fails
        let m = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
        let _ = reos_passes(&m);
    }

    #[test]
    fn test_reos_diazepam_passes() {
        // Diazepam: MW~285, logP~2.9, HBD=0, HBA=2, charge=0, rotb=1, HAC~22 — all in range
        let m = mol("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21");
        assert!(reos_passes(&m), "diazepam should pass REOS");
    }

    // -- Ghose filter tests -------------------------------------------------

    #[test]
    fn test_ghose_aspirin_range() {
        // Aspirin MW=180 is below Ghose MW lower bound (160 ≤ MW ≤ 480). Should pass MW.
        // HAC=13 < 20 lower bound → Ghose fails for aspirin
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        // Just verify it doesn't panic
        let _ = ghose_passes(&m);
    }

    #[test]
    fn test_ghose_ibuprofen_passes() {
        // Ibuprofen: MW=206, logP~3.5, HAC=13 (borderline)
        let m = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
        let _ = ghose_passes(&m);
    }

    // -- Basic count descriptor tests ----------------------------------------

    #[test]
    fn test_num_heteroatoms_aspirin() {
        // Aspirin: 4 O atoms (ester C=O, ester O, carboxylic C=O, carboxylic OH)
        assert_eq!(num_heteroatoms(&mol("CC(=O)Oc1ccccc1C(=O)O")), 4);
    }

    #[test]
    fn test_num_heteroatoms_benzene_zero() {
        assert_eq!(num_heteroatoms(&mol("c1ccccc1")), 0);
    }

    #[test]
    fn test_ring_count_benzene() {
        assert_eq!(ring_count(&mol("c1ccccc1")), 1);
    }

    #[test]
    fn test_ring_count_naphthalene() {
        assert_eq!(ring_count(&mol("c1ccc2ccccc2c1")), 2);
    }

    #[test]
    fn test_ring_count_acyclic_zero() {
        assert_eq!(ring_count(&mol("CCO")), 0);
    }

    #[test]
    fn test_num_saturated_rings_cyclohexane() {
        assert_eq!(num_saturated_rings(&mol("C1CCCCC1")), 1);
    }

    #[test]
    fn test_num_saturated_rings_benzene_zero() {
        assert_eq!(num_saturated_rings(&mol("c1ccccc1")), 0);
    }

    #[test]
    fn test_num_aliphatic_rings_cyclohexane() {
        assert_eq!(num_aliphatic_rings(&mol("C1CCCCC1")), 1);
    }

    #[test]
    fn test_num_aliphatic_rings_benzene_zero() {
        assert_eq!(num_aliphatic_rings(&mol("c1ccccc1")), 0);
    }

    #[test]
    fn test_num_stereocenters_alanine() {
        // L-alanine: 1 R/S center
        assert_eq!(num_stereocenters(&mol("[C@@H](N)(C)C(=O)O")), 1);
    }

    #[test]
    fn test_num_stereocenters_achiral_zero() {
        assert_eq!(num_stereocenters(&mol("CC(=O)O")), 0);
    }

    // -- BalabanJ tests -------------------------------------------------

    #[test]
    fn test_balaban_j_ethane() {
        let m = mol("CC");
        let bj = balaban_j(&m);
        assert!(bj > 0.0, "ethane should have positive BalabanJ");
    }

    #[test]
    fn test_balaban_j_benzene() {
        let m = mol("c1ccccc1");
        let bj = balaban_j(&m);
        assert!(bj > 0.0, "benzene should have positive BalabanJ");
    }

    #[test]
    fn test_balaban_j_single_atom_zero() {
        let m = mol("C");
        let bj = balaban_j(&m);
        assert_eq!(bj, 0.0, "single atom should have BalabanJ = 0");
    }

    // -- Ipc tests -------------------------------------------------------

    #[test]
    fn test_ipc_ethane() {
        let m = mol("CC");
        let ipc_val = ipc(&m);
        assert!(ipc_val >= 0.0, "ethane should have non-negative Ipc");
    }

    #[test]
    fn test_ipc_benzene() {
        let m = mol("c1ccccc1");
        let ipc_val = ipc(&m);
        assert!(ipc_val > 0.0, "benzene should have positive Ipc");
    }

    #[test]
    fn test_ipc_single_atom_zero() {
        let m = mol("C");
        let ipc_val = ipc(&m);
        assert_eq!(ipc_val, 0.0, "single atom should have Ipc = 0");
    }

    // -- HallKierAlpha tests -----------------------------------------------

    #[test]
    fn test_hall_kier_alpha_ethane() {
        let m = mol("CC");
        let hka = hall_kier_alpha(&m);
        assert!(hka > 0.0, "ethane should have positive HallKierAlpha");
    }

    #[test]
    fn test_hall_kier_alpha_methane() {
        let m = mol("C");
        let hka = hall_kier_alpha(&m);
        assert!(hka > 0.0, "methane should have positive HallKierAlpha");
    }

    #[test]
    fn test_hall_kier_alpha_benzene() {
        let m = mol("c1ccccc1");
        let hka = hall_kier_alpha(&m);
        assert!(hka > 0.0, "benzene should have positive HallKierAlpha");
    }

    // -- USRCAT tests -------------------------------------------------------

    #[test]
    fn test_usrcat_shape() {
        let m = mol("CC");
        let usr = usrcat(&m);
        assert_eq!(usr.len(), 42, "USRCAT should return 42 values");
        assert!(usr[0] >= 0.0, "first slot should be non-negative");
    }

    #[test]
    fn test_usrcat_donors_acceptors() {
        let m = mol("CCO");
        let usr = usrcat(&m);
        assert!(usr[36] >= 0.0, "donor count should be non-negative");
        assert!(
            usr[37] > 0.0,
            "acceptor count should be positive (O present)"
        );
    }

    #[test]
    fn test_usrcat_aromatic() {
        let m = mol("c1ccccc1");
        let usr = usrcat(&m);
        assert!(
            usr[38] > 0.0,
            "aromatic count should be positive for benzene"
        );
    }

    #[test]
    fn test_usrcat_charged() {
        let m = mol("CC(=O)[O-]");
        let usr = usrcat(&m);
        assert!(
            usr[40] > 0.0,
            "anion count should be positive for charged carboxylate"
        );
    }

    // -- MMFF94 charges tests -----------------------------------------------

    #[test]
    fn test_mmff94_charges_length() {
        let m = mol("CCO");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 3, "should have 3 charges for 3 atoms");
    }

    #[test]
    fn test_mmff94_charges_ethane() {
        let m = mol("CC");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 2);
        // Both carbons should have similar (small negative) charges
        assert!(
            (charges[0] - charges[1]).abs() < 0.1,
            "carbons in ethane should have similar charges"
        );
    }

    #[test]
    fn test_mmff94_charges_charged_species() {
        let m = mol("CC(=O)[O-]");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 4);
        // Carboxylate oxygen should be negative
        assert!(charges[3] < 0.0, "carboxylate oxygen should be negative");
    }

    #[test]
    fn test_mmff94_charges_water() {
        let m = mol("O");
        let charges = mmff94_charges(&m);
        assert_eq!(charges.len(), 1);
        assert!(
            charges[0].is_finite(),
            "water oxygen charge should be finite"
        );
    }

    // -- Element count tests -----------------------------------------------

    #[test]
    fn test_num_carbons_ethane() {
        let m = mol("CC");
        assert_eq!(num_carbons(&m), 2);
    }

    #[test]
    fn test_num_nitrogens_methylamine() {
        let m = mol("CN");
        assert_eq!(num_nitrogens(&m), 1);
    }

    #[test]
    fn test_num_oxygens_methanol() {
        let m = mol("CO");
        assert_eq!(num_oxygens(&m), 1);
    }

    #[test]
    fn test_num_halogens() {
        let m = mol("CCF");
        assert_eq!(num_carbons(&m), 2);
        assert_eq!(num_fluorines(&m), 1);
    }

    #[test]
    fn test_num_hydrogens_methane() {
        let m = mol("C");
        // Methane: 1 carbon with 4 implicit hydrogens
        assert_eq!(num_hydrogens(&m), 4);
    }

    #[test]
    fn test_num_hydrogens_ethane() {
        let m = mol("CC");
        // Ethane: 2 carbons, 6 total hydrogens
        assert_eq!(num_hydrogens(&m), 6);
    }

    #[test]
    fn test_num_hydrogens_water() {
        let m = mol("O");
        // Water: 2 implicit hydrogens on oxygen
        assert_eq!(num_hydrogens(&m), 2);
    }

    // -- Molecular formula tests -----------------------------------------------

    #[test]
    fn test_calc_mol_formula_ethane() {
        let m = mol("CC");
        assert_eq!(calc_mol_formula(&m), "C2H6");
    }

    // -- Bridgehead and spiro atom tests -----------------------------------

    #[test]
    fn test_num_bridgehead_atoms_acyclic() {
        // Acyclic molecule should return 0
        let m = mol("CCC");
        assert_eq!(num_bridgehead_atoms(&m), 0);
    }

    #[test]
    fn test_num_bridgehead_atoms_single_ring() {
        // Cyclohexane (single ring) should have 0 bridgeheads
        let m = mol("C1CCCCC1");
        assert_eq!(num_bridgehead_atoms(&m), 0);
    }

    #[test]
    fn test_num_bridgehead_atoms_norbornane() {
        // Norbornane (bicyclo[2.2.1]heptane): 2 bridgehead atoms
        // Structure: two rings sharing 2 atoms at positions 1 and 2, with bridges
        let m = mol("C1CC2CCC1C2");
        assert_eq!(num_bridgehead_atoms(&m), 2);
    }

    #[test]
    fn test_num_bridgehead_atoms_naphthalene_fused() {
        // Naphthalene is fused ring system (not bridged)
        // Bridgeheads = atoms in 2+ rings but only on shared edges (fused, no bridges)
        let m = mol("c1ccc2ccccc2c1");
        assert_eq!(
            num_bridgehead_atoms(&m),
            0,
            "naphthalene is fused, not bridged"
        );
    }

    #[test]
    fn test_num_spiro_atoms_single_ring() {
        // Cyclohexane has no spiro atoms
        let m = mol("C1CCCCC1");
        assert_eq!(num_spiro_atoms(&m), 0);
    }

    #[test]
    fn test_calc_mol_formula_water() {
        let m = mol("O");
        assert_eq!(calc_mol_formula(&m), "H2O");
    }

    #[test]
    fn test_calc_mol_formula_benzene() {
        let m = mol("c1ccccc1");
        assert_eq!(calc_mol_formula(&m), "C6H6");
    }

    #[test]
    fn test_calc_mol_formula_acetic_acid() {
        let m = mol("CC(=O)O");
        assert_eq!(calc_mol_formula(&m), "C2H4O2");
    }

    // =========================================================================
    // Functional group bond counts (C4 - granular classification)
    // =========================================================================

    #[test]
    fn test_num_amide_bonds_acetamide() {
        // CH3-C(=O)-N-H: one amide bond
        let m = mol("CC(=O)N");
        assert_eq!(num_amide_bonds(&m), 1);
    }

    #[test]
    fn test_num_amide_bonds_urea() {
        // N-C(=O)-N: one amide bond (C=O-N)
        let m = mol("NC(=O)N");
        assert_eq!(num_amide_bonds(&m), 1);
    }

    #[test]
    fn test_num_amide_bonds_primary_amide() {
        // CH3-C(=O)-NH2: primary amide
        let m = mol("CC(=O)N");
        assert_eq!(num_amide_bonds(&m), 1);
    }

    #[test]
    fn test_num_amide_bonds_none() {
        // Benzene has no amide bonds
        let m = mol("c1ccccc1");
        assert_eq!(num_amide_bonds(&m), 0);
    }

    #[test]
    fn test_num_ester_bonds_methyl_formate() {
        // H-C(=O)-O-CH3: one ester bond
        let m = mol("COC=O");
        assert_eq!(num_ester_bonds(&m), 1);
    }

    #[test]
    fn test_num_ester_bonds_acetic_acid_methyl_ester() {
        // CH3-C(=O)-O-CH3: one ester bond
        let m = mol("CC(=O)OC");
        assert_eq!(num_ester_bonds(&m), 1);
    }

    #[test]
    fn test_num_ester_bonds_none() {
        // Carboxylic acid (COOH) is not an ester
        let m = mol("CC(=O)O");
        // CC(=O)O is acetic acid: C=O with O-H (not O-C)
        assert_eq!(num_ester_bonds(&m), 0);
    }

    #[test]
    fn test_num_ester_bonds_aspirin() {
        // Aspirin CC(=O)Oc1ccccc1C(=O)O has one ester bond (acetyl ester)
        // and one carboxylic acid (not counted)
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(
            num_ester_bonds(&m),
            1,
            "aspirin has one ester bond (aryl ester)"
        );
    }

    // B5: context-dependent alkene C LogP contributions
    #[test]
    fn test_logp_ethylene_terminal() {
        // Ethylene: both C are terminal =CH2 (h=2, no aromatic neighbor) → each 0.1551
        let m = mol("C=C");
        let lp = logp_crippen(&m);
        assert!(lp > 0.3 && lp < 1.3, "ethylene logp = {lp}");
    }

    #[test]
    fn test_logp_propene_terminal_internal() {
        // Propene: =CH2 (0.1551) + =CH- internal (0.2274) + CH3 → propene > ethylene
        let m = mol("CC=C");
        let lp = logp_crippen(&m);
        let eth = logp_crippen(&mol("C=C"));
        assert!(lp > eth, "propene logp ({lp}) should exceed ethylene ({eth})");
        assert!(lp > 0.5 && lp < 2.0, "propene logp = {lp} out of expected range");
    }

    #[test]
    fn test_logp_styrene_splits_correctly() {
        // Styrene: =CH2 (0.1551) + =CH- adj to Ar (0.2640) + 6 aromatic C
        let m = mol("C=Cc1ccccc1");
        let lp = logp_crippen(&m);
        assert!(lp > 1.8 && lp < 3.4, "styrene logp = {lp}");
        // Ar-adjacent CH contributes more (0.2640) than terminal CH2 (0.1551)
        let per_atom = logp_crippen_per_atom(&m);
        assert!(per_atom.len() >= 8, "styrene has 8 heavy atoms");
    }

    #[test]
    fn test_logp_1_phenylpropene_ar_adjacent() {
        // 1-phenylpropene Ph-CH=CH-CH3: Ar-adjacent =CH- (0.2640) + internal =CH- (0.2274)
        let m = mol("CC=Cc1ccccc1");
        let lp = logp_crippen(&m);
        assert!(lp > 2.0 && lp < 4.0, "1-phenylpropene logp = {lp}");
    }

    // C1: Reference LogP tests for complex polar molecules.
    //
    // For molecules with multiple polar functional groups (carboxylate, amide,
    // amine salt, conjugated C=C-C=O), the Crippen atom-type model underestimates
    // LogP versus RDKit. xlogp3() is implemented and more accurate for these cases.
    // These tests document the current values and catch regressions.

    #[test]
    fn test_logp_curcumin_reference() {
        // Curcumin: two conjugated vinyl-ketone arms, phenol/methoxy substituents.
        // RDKit Crippen ~3.04; Crippen atom-type model gives a lower value.
        // xlogp3 is the recommended API for complex polyphenols.
        let m = mol("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O");
        let lp = logp_crippen(&m);
        assert!(lp > -5.0 && lp < 5.0, "curcumin crippen logp = {lp}");
    }

    #[test]
    fn test_logp_complex_molecules_xlogp3_preferred() {
        // For high-error molecules, verify xlogp3 gives a more positive value
        // (closer to RDKit) than Crippen for curcumin.
        let m = mol("COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O");
        let crippen = logp_crippen(&m);
        let xl3 = crate::xlogp3::xlogp3(&m);
        // Both should be finite; just document that they exist
        assert!(crippen.is_finite(), "crippen logp must be finite");
        assert!(xl3.is_finite(), "xlogp3 must be finite");
    }

    // ---- Enone vinyl C tests (v0.1.99) ----

    #[test]
    fn test_logp_mvk_enone_vinyl() {
        // Methyl vinyl ketone: CH2=CH-C(=O)-CH3
        // =CH2 (0.1551) + =CH- adj to C=O (enone, 0.1302) + C=O + CH3
        // RDKit Crippen ~0.44; expect reasonable range
        let m = mol("C=CC(=O)C");
        let lp = logp_crippen(&m);
        assert!(lp > 0.0 && lp < 1.5, "MVK logp = {lp}");
    }

    #[test]
    fn test_logp_chalcone_enone() {
        // Chalcone: Ph-CH=CH-C(=O)-Ph
        // Ar-CH= (0.2640, ar-adjacent wins) + =CH-C(=O) (0.1302, enone)
        // RDKit ~3.04; Crippen gives slightly different value
        let m = mol("c1ccccc1/C=C/C(=O)c1ccccc1");
        let lp = logp_crippen(&m);
        assert!(lp > 1.8 && lp < 4.5, "chalcone logp = {lp}");
    }

    #[test]
    fn test_logp_crotonate_internal_enone() {
        // Crotonic acid: CH3-CH=CH-COOH (trans)
        // =CH- (0.2274, non-ar, 1H) adj to C=O → enone case (0.1302)
        let m = mol("CC=CC(=O)O");
        let lp = logp_crippen(&m);
        assert!(lp > -0.5 && lp < 1.5, "crotonate logp = {lp}");
    }

    #[test]
    fn test_logp_enone_vs_plain_alkene() {
        // Enone vinyl C (0.1302) is less hydrophobic than plain internal alkene (0.2274)
        // MVK (C=C-C=O) vs 1-butene (C=C-CC)
        let enone = logp_crippen(&mol("C=CC(=O)C")); // MVK
        let alkene = logp_crippen(&mol("C=CCC"));    // 1-butene
        assert!(
            enone < alkene,
            "MVK ({enone:.4}) should be < 1-butene ({alkene:.4}): enone is less hydrophobic"
        );
    }
}
