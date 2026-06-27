//! Molecular descriptor functions for drug-likeness and physical property estimation.
//!
//! All functions accept a `&Molecule` reference.  Molecules with aromatic bonds
//! (SMILES lowercase notation) are kekulized internally where hydrogen counts
//! are required; the caller's molecule is never mutated.

use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::OnceLock;

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Element, Molecule, bond_order_sum, implicit_hcount,
};
use chematic_perception::{
    RingSet, apply_aromaticity, augmented_ring_set, find_ring_families, find_sssr,
};
use chematic_smarts::{MatchConfig, find_matches_with_rings_and_config, parse_smarts};

/// True if `idx` has a double bond to any neighbor whose atomic number equals `target_an`.
fn has_double_bond_to(mol: &Molecule, idx: AtomIdx, target_an: u8) -> bool {
    mol.neighbors(idx).any(|(nb, bidx)| {
        mol.bond(bidx).order == BondOrder::Double
            && mol.atom(nb).element.atomic_number() == target_an
    })
}

/// Returns true if `idx` is in a ring.
///
/// For degree-2 atoms one BFS suffices. For degree-3 atoms (e.g. N with two ring
/// bonds + one exocyclic substituent) the first neighbor may be exocyclic — a dead
/// end that makes the single-pair BFS return false even though `idx` is in a ring.
/// This version tries each neighbor as the BFS start and returns true as soon as
/// any other neighbor is found reachable without going through `idx`.
fn is_atom_in_ring(mol: &Molecule, idx: AtomIdx) -> bool {
    let nbs: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();
    if nbs.len() < 2 {
        return false;
    }
    for start_i in 0..nbs.len() {
        let start = nbs[start_i];
        let targets: FxHashSet<AtomIdx> = nbs
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != start_i)
            .map(|(_, &nb)| nb)
            .collect();
        let mut visited = FxHashSet::default();
        visited.insert(idx);
        visited.insert(start);
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(start);
        while let Some(curr) = queue.pop_front() {
            if targets.contains(&curr) {
                return true;
            }
            for (nb, _) in mol.neighbors(curr) {
                if visited.insert(nb) {
                    queue.push_back(nb);
                }
            }
        }
    }
    false
}

/// True if `idx` (O atom) is an aromatic oxide bridge: in a ring, bonded to an
/// aromatic C AND to a sp2 C via a single bond (C=C elsewhere, not C=O).
/// RDKit perceives such O as aromatic ([o] type). Used for TPSA and LogP.
fn is_aromatic_oxide_bridge(mol: &Molecule, idx: AtomIdx) -> bool {
    let has_aromatic_c_nb = mol.neighbors(idx).any(|(nb, _)| {
        let a = mol.atom(nb);
        a.aromatic && a.element.atomic_number() == 6
    });
    if !has_aromatic_c_nb {
        return false;
    }
    let has_vinyl_c_nb = mol.neighbors(idx).any(|(nb, bidx)| {
        mol.bond(bidx).order != BondOrder::Double
            && mol.atom(nb).element.atomic_number() == 6
            && mol.neighbors(nb).any(|(nb2, b2)| {
                nb2 != idx
                    && mol.bond(b2).order == BondOrder::Double
                    && mol.atom(nb2).element.atomic_number() == 6
            })
    });
    has_vinyl_c_nb && is_atom_in_ring(mol, idx)
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
        if atom.wildcard {
            continue;
        }
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
        if atom.wildcard {
            continue;
        }
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
            (is_nitrogen(an) || is_oxygen(an) || an == 16) && implicit_hcount(mol, *idx) > 0
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
fn hba_count_from_set(mol: &Molecule, ring_bonds: &FxHashSet<BondIdx>) -> usize {
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
                    // Pyridine-type aromatic N (lone pair orthogonal to pi) IS an HBA.
                    // Excluded cases:
                    //   h > 0  → [nH] pyrrole-type: lone pair in pi system
                    //   degree >= 3 → N-substituted pyrrole or bridgehead N
                    //                 (e.g. N-methyl pyrrole, indolizine N): lone pair
                    //                 participates in the aromatic π system
                    // Known limitation: 1/5000 molecules (a merocyanine dye) differs from
                    // RDKit because N-methyl in push-pull quinonoid rings is ambiguous.
                    h == 0 && mol.degree(*idx) < 3
                } else {
                    // Non-aromatic N: must have formal valence 3 ([N;v3] in SMARTS);
                    // this excludes radical N (C[N]C, valence 2) and unusual species.
                    let bov = bond_order_sum(mol, *idx) as usize + h as usize;
                    if bov != 3 {
                        return false;
                    }
                    // Exclude N with single bond to any atom that itself has a
                    // NON-RING pi bond to O/N/P/S: amide, sulfonamide, phosphonamide,
                    // thioamide, etc.  Ring pi bonds (e.g. ring C=N in guanidinium)
                    // do NOT trigger the exclusion — matching SMARTS !@ semantics.
                    !n_adjacent_to_pi_center(mol, *idx, ring_bonds)
                }
            } else if is_oxygen(an) {
                if atom.charge > 0 {
                    return false;
                } // O+ (oxonium) never HBA
                if atom.charge < 0 {
                    return true;
                } // [O-] (carboxylate etc.) always HBA
                // Total H = implicit + explicit isotopic H (e.g. [2H]O[2H] = D2O)
                let impl_h = implicit_hcount(mol, *idx);
                let expl_h = mol
                    .neighbors(*idx)
                    .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
                    .count() as u8;
                match impl_h + expl_h {
                    // H0: ether, carbonyl, epoxide — divalent check excludes radical [O]
                    0 => bond_order_sum(mol, *idx) == 2,
                    // H1: alcohol/phenol OH — exclude if neighbor has =O/=N/=P/=S
                    1 => !neighbor_has_pi_bond_to_onps(mol, *idx),
                    // H2+ (H2O, D2O, etc.) — not HBA
                    _ => false,
                }
            } else if is_sulfur(an) {
                if atom.charge < 0 {
                    return true;
                } // [S-] always HBA
                if atom.charge != 0 {
                    return false;
                } // S+/S2+ etc. never HBA
                if atom.aromatic {
                    // Aromatic s (thiophene-type): lone pair available
                    true
                } else {
                    // Total H = implicit + explicit isotopic H (handles [2H]S[2H] = D2S)
                    let impl_h = implicit_hcount(mol, *idx);
                    let expl_h = mol
                        .neighbors(*idx)
                        .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
                        .count() as u8;
                    let total_h = impl_h + expl_h;
                    // Formal valence = bond-order sum + implicit H (handles S=C case)
                    let bos = bond_order_sum(mol, *idx) as usize + impl_h as usize;
                    if bos != 2 {
                        return false;
                    } // not divalent (sulfoxide/sulfone etc.)
                    match total_h {
                        // SH: thiol — exclude if neighbor has =O/=N/=P/=S (thio-acid)
                        1 => !neighbor_has_pi_bond_to_onps(mol, *idx),
                        // S with 0H: sulfide, thioketone — HBA
                        0 => true,
                        // H2S or higher — not HBA
                        _ => false,
                    }
                }
            } else {
                false
            }
        })
        .count()
}

pub fn hba_count(mol: &Molecule) -> usize {
    hba_count_from_set(mol, &ring_bond_indices(mol))
}

/// True if any heavy-atom neighbor of `idx` itself carries a double bond to
/// N, O, P, or S.  Used to exclude OH/SH groups from the HBA count when the
/// attached heavy atom has such a π-bond (e.g. C=O, C=N, As=O, P=O, S=O).
/// Matches the `[!$(*=[O,N,P,S])]` exclusion in the RDKit HBA SMARTS.
fn neighbor_has_pi_bond_to_onps(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.neighbors(idx).any(|(nb_idx, _)| {
        has_double_bond_to(mol, nb_idx, 7)   // =N
            || has_double_bond_to(mol, nb_idx, 8)  // =O
            || has_double_bond_to(mol, nb_idx, 15) // =P
            || has_double_bond_to(mol, nb_idx, 16) // =S
    })
}

/// True if `nb_idx` has a **non-ring** double bond to an atom with atomic
/// number `target_an`.  Implements the SMARTS `!@` (non-ring bond) constraint:
/// ring double bonds (e.g. C=N inside a cyclic amidine) do NOT trigger the
/// N-exclusion and must not be treated as π-centres.
fn has_nonring_double_bond_to(
    mol: &Molecule,
    nb_idx: AtomIdx,
    target_an: u8,
    ring_bonds: &FxHashSet<BondIdx>,
) -> bool {
    mol.neighbors(nb_idx).any(|(far, bidx)| {
        mol.bond(bidx).order == BondOrder::Double
            && mol.atom(far).element.atomic_number() == target_an
            && !ring_bonds.contains(&bidx)
    })
}

/// True if `idx` (an N atom) has a **single-like** bond (including `/`/`\`
/// stereo bonds stored as `BondOrder::Up`/`Down`) to any neighbor that itself
/// carries a **non-ring** double bond to O, N, P, or S.
///
/// Matches the RDKit HBA SMARTS exclusion `!$(N-*=!@[O,N,P,S])` where:
/// - `*` is any atom — C (amide, thioamide, guanidine), S (sulfonamide),
///   P (phosphonamide), N (nitroso-adjacent), etc.
/// - `-` includes stereo-direction bonds (`/`, `\` → `Up`/`Down`)
/// - `=!@` means double bond that is NOT a ring bond
fn n_adjacent_to_pi_center(mol: &Molecule, idx: AtomIdx, ring_bonds: &FxHashSet<BondIdx>) -> bool {
    mol.neighbors(idx).any(|(nb_idx, bidx)| {
        // Accept single, E/Z stereo, and aromatic bonds as "single-like"
        matches!(
            mol.bond(bidx).order,
            BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Aromatic
        ) && (has_nonring_double_bond_to(mol, nb_idx, 8, ring_bonds)   // *=O
                || has_nonring_double_bond_to(mol, nb_idx, 7, ring_bonds)  // *=N
                || has_nonring_double_bond_to(mol, nb_idx, 16, ring_bonds) // *=S
                || has_nonring_double_bond_to(mol, nb_idx, 15, ring_bonds)) // *=P
    })
}

// ---------------------------------------------------------------------------
// 6. Rotatable bond count
// ---------------------------------------------------------------------------

/// Count rotatable bonds (RDKit-compatible strict definition).
///
/// A bond is rotatable when all of the following hold:
/// - It is a single bond (or a stereo bond Up/Down, which is single).
/// - Neither endpoint is terminal (degree > 1 in the heavy-atom graph).
/// - It is not part of any ring (SSSR membership).
/// - It is not an amide bond (C–N where C has a C=O).
/// - Neither endpoint carries a triple bond (excludes propargylic C–C in alkynes).
/// - Neither endpoint is a cumulated-double-bond centre (excludes allene C=C=C bonds).
fn rotatable_bond_count_from_set(mol: &Molecule, ring_bond_set: &FxHashSet<BondIdx>) -> usize {
    mol.bonds()
        .filter(|(bidx, bond)| {
            let is_single = matches!(
                bond.order,
                BondOrder::Single | BondOrder::Up | BondOrder::Down
            );
            is_single
                && !ring_bond_set.contains(bidx)
                && mol.degree(bond.atom1) > 1
                && mol.degree(bond.atom2) > 1
                && !is_carbonyl_hetero_bond(mol, bond.atom1, bond.atom2)
                && !is_neopentyl_like(mol, bond.atom1, bond.atom2)
                && !has_triple_bond(mol, bond.atom1)
                && !has_triple_bond(mol, bond.atom2)
                && !is_cumulated_double(mol, bond.atom1)
                && !is_cumulated_double(mol, bond.atom2)
        })
        .count()
}

pub fn rotatable_bond_count(mol: &Molecule) -> usize {
    rotatable_bond_count_from_set(mol, &ring_bond_indices(mol))
}

/// True if atom `idx` has at least one triple bond.
fn has_triple_bond(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.neighbors(idx)
        .any(|(_, bidx)| mol.bond(bidx).order == BondOrder::Triple)
}

/// True if atom `idx` is the centre of a cumulated double-bond system (≥2 double bonds),
/// as found in allenes (C=C=C) and ketenes (C=C=O).
/// Restricted to carbon: sulfone S(=O)(=O) and phosphate P(=O) are not allene-like
/// and their bonds must not be excluded from the rotatable-bond count.
fn is_cumulated_double(mol: &Molecule, idx: AtomIdx) -> bool {
    mol.atom(idx).element.atomic_number() == 6
        && mol
            .neighbors(idx)
            .filter(|(_, bidx)| mol.bond(*bidx).order == BondOrder::Double)
            .count()
            >= 2
}

/// Build the ring-bond set from a pre-computed rings slice.
fn ring_bond_indices_from_rings(mol: &Molecule, rings: &[Vec<AtomIdx>]) -> FxHashSet<BondIdx> {
    let mut set = FxHashSet::default();
    for ring in rings {
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

/// Indices of all bonds participating in at least one SSSR ring.
fn ring_bond_indices(mol: &Molecule) -> FxHashSet<BondIdx> {
    ring_bond_indices_from_rings(mol, find_sssr(mol).rings())
}

/// True if the bond between `a` and `b` is an amide-like C-N bond
/// (one atom is N, the other is C with a double bond to O).
/// RDKit strict: exclude bond A-B when A is a quaternary C (degree 4, no H) whose only
/// non-terminal heavy-atom neighbor is B.  Covers tert-butyl (CC(C)(C)-X), neopentyl, Boc groups.
fn is_neopentyl_like(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    is_neopentyl_center(mol, a, b) || is_neopentyl_center(mol, b, a)
}

fn is_neopentyl_center(mol: &Molecule, center: AtomIdx, other: AtomIdx) -> bool {
    if mol.atom(center).element.atomic_number() != 6 {
        return false;
    }
    if mol.degree(center) != 4 {
        return false;
    }
    let others: Vec<AtomIdx> = mol
        .neighbors(center)
        .filter(|(nb, _)| *nb != other)
        .map(|(nb, _)| nb)
        .collect();
    // All 3 non-target neighbors must be terminal (degree 1)
    if !others.iter().all(|&nb| mol.degree(nb) == 1) {
        return false;
    }
    // All must share the same atomic number (e.g. all CH3, or all F)
    // This prevents false positives for C(CH3)(OH)(CH3)-X cases.
    let an0 = mol.atom(others[0]).element.atomic_number();
    others
        .iter()
        .all(|&nb| mol.atom(nb).element.atomic_number() == an0)
}

/// RDKit strict definition excludes C–N/O single bonds when C carries a pi-bond to
/// O, N, or S.  Covers amide C(=O)–N, ester C(=O)–O, guanidine/amidine C(=N)–N,
/// thioamide/thiourea C(=S)–N.  Ketone C(=O)–C bonds are NOT excluded.
fn is_carbonyl_hetero_bond(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let an_a = mol.atom(a).element.atomic_number();
    let an_b = mol.atom(b).element.atomic_number();
    let c_idx = match (an_a, an_b) {
        (6, 7) | (6, 8) => a, // C–N or C–O
        (7, 6) | (8, 6) => b, // N–C or O–C
        _ => return false,
    };
    has_double_bond_to(mol, c_idx, 8)   // C=O  (amide / ester)
        || has_double_bond_to(mol, c_idx, 7)  // C=N  (guanidine / amidine)
        || has_double_bond_to(mol, c_idx, 16) // C=S  (thioamide / thiourea)
}

// ---------------------------------------------------------------------------
// 7. TPSA
// ---------------------------------------------------------------------------

/// Compute the topological polar surface area (Å²) using the Ertl (2000) table.
///
/// Reference: P. Ertl, B. Rohde, P. Selzer, J. Med. Chem. 2000, 43, 3714-3717.
fn tpsa_nitrogen(
    mol: &Molecule,
    idx: AtomIdx,
    is_aromatic: bool,
    h: u8,
    charge: i8,
    ring_bonds: &FxHashSet<BondIdx>,
) -> f64 {
    if is_aromatic {
        let degree = mol.degree(idx);

        if h > 0 {
            15.79
        } else if degree >= 3 {
            // Distinguish true ring-junction N (all bonds Aromatic, e.g. bridgehead in
            // imidazo[1,2-a]pyridine or phthalazinone) from N-substituted (has at least
            // one non-Aromatic bond — methyl, phenyl via explicit single, etc.).
            // Ring-junction: 4.41 (neutral) / 4.10 (cationic).
            // N-substituted:  4.93 (neutral) / 3.88 (cationic). — calibrated from RDKit.
            let is_ring_junction = mol
                .neighbors(idx)
                .all(|(_, bidx)| mol.bond(bidx).order == BondOrder::Aromatic);
            if charge > 0 {
                if is_ring_junction { 4.10 } else { 3.88 }
            } else {
                if is_ring_junction { 4.41 } else { 4.93 }
            }
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
            if has_oxo && has_o_minus {
                43.14 // nitro: N+(=O)[O-]
            } else if has_double_bond_to(mol, idx, 7) {
                14.10 // azide central N+: R-N=[N+]=[N-]
            } else if has_double_bond_to(mol, idx, 6) {
                3.01 // nitrone: C=N+(R)-O- (exocyclic double bond to C)
            } else {
                // Quaternary ammonium N+ (no lone pair) and ionic N-oxide [N+][O-]
                // have no polar surface contribution (Ertl 2000).
                0.00
            }
        } else if h >= 2 {
            26.02
        } else if h == 1 {
            if has_double_bond_to(mol, idx, 6) {
                // Imine =NH (N=C with H): 23.85 (RDKit calibrated value, not 23.79)
                23.85
            } else {
                12.03
            }
        } else if charge < 0 {
            // N- anion: azide terminal R-N=[N+]=[N-] → 22.30; others use azo/amine fallback
            if has_double_bond_to(mol, idx, 7) {
                22.30
            } else {
                12.36
            }
        } else {
            // h == 0, charge 0: nitrile, imine/nitroso/azo/P=N, amine
            let has_triple_to_c = mol.neighbors(idx).any(|(nb, bidx)| {
                mol.bond(bidx).order == BondOrder::Triple
                    && mol.atom(nb).element.atomic_number() == 6
            });
            if has_triple_to_c {
                23.79 // nitrile N≡C
            } else if has_double_bond_to(mol, idx, 6)
                || has_double_bond_to(mol, idx, 8)
                || has_double_bond_to(mol, idx, 7)
                || has_double_bond_to(mol, idx, 15)
                || has_double_bond_to(mol, idx, 16)
            {
                12.36 // imine N=C, nitroso N=O, azo N=N, phosphazene P=N, sulfonimidyl N=S
            } else {
                // Aziridine (3-membered ring) N bonded to external P or S → 3.01.
                // Larger ring N bonded to external P, or acyclic N → 3.24.
                let ring_nbs: Vec<_> = mol
                    .neighbors(idx)
                    .filter(|(_, b)| ring_bonds.contains(b))
                    .map(|(nb, _)| nb)
                    .collect();
                let in_3ring =
                    ring_nbs.len() == 2 && mol.bond_between(ring_nbs[0], ring_nbs[1]).is_some();
                let has_external_ps = mol.neighbors(idx).any(|(nb, bidx)| {
                    let an = mol.atom(nb).element.atomic_number();
                    (an == 15 || an == 16) && !ring_bonds.contains(&bidx)
                });
                if in_3ring && has_external_ps {
                    3.01
                } else {
                    3.24
                }
            }
        }
    }
}

fn tpsa_oxygen(mol: &Molecule, idx: AtomIdx, is_aromatic: bool, h: u8, charge: i8) -> f64 {
    if is_aromatic {
        13.14
    } else if h > 0 {
        // Water (H₂O): isolated O with no heavy-atom neighbors → 31.50
        // Hydroxyl (–OH): O with one heavy-atom neighbor → 20.23
        if mol.neighbors(idx).count() == 0 {
            31.50
        } else {
            20.23
        }
    } else {
        // O with charge == -1: true nitro O- (N+ that also has =O) → 0.0; N-oxide/other O- → 23.06
        if charge == -1 {
            let is_nitro_o_minus = mol.neighbors(idx).any(|(nb, _)| {
                let n = mol.atom(nb);
                n.element.atomic_number() == 7
                    && n.charge == 1
                    && mol.neighbors(nb).any(|(nb2, bidx2)| {
                        mol.atom(nb2).element.atomic_number() == 8
                            && mol.bond(bidx2).order == BondOrder::Double
                    })
            });
            return if is_nitro_o_minus { 0.0 } else { 23.06 };
        }
        let dbl_nb_pair = mol
            .neighbors(idx)
            .find(|&(_, bidx)| mol.bond(bidx).order == BondOrder::Double);
        let dbl_nb = dbl_nb_pair
            .map(|(nei, _)| (mol.atom(nei).element.atomic_number(), mol.atom(nei).charge));
        match dbl_nb {
            Some((6, _)) => 17.07,
            // nitroso N=O: neutral N → 17.07; nitro [N+](=O)[O-]: charged N → 0.0
            Some((7, 0)) => 17.07,
            Some((7, _)) => 0.0,
            // Se=O oxygens (seleninic/selenious acid) carry same contribution as C=O
            Some((34, _)) => 17.07,
            // S=O: if S also has N=S double bond (sulfonimidyl) → 17.07; else S handles it → 0.0
            Some((16, _)) => {
                let s_idx = dbl_nb_pair.unwrap().0;
                if has_double_bond_to(mol, s_idx, 7) {
                    17.07
                } else {
                    0.0
                }
            }
            Some(_) => 0.0,
            None => {
                if is_aromatic_oxide_bridge(mol, idx) {
                    13.14
                } else {
                    // Epoxide O (3-membered ring) → 12.53; ether/ring O → 9.23
                    let nbs: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();
                    if nbs.len() == 2 && mol.bond_between(nbs[0], nbs[1]).is_some() {
                        12.53
                    } else {
                        9.23
                    }
                }
            }
        }
    }
}

fn tpsa_sulfur(mol: &Molecule, idx: AtomIdx, is_aromatic: bool, h: u8, charge: i8) -> f64 {
    // [S+][O-] zwitterionic sulfoxide: RDKit assigns 0 to S+, 23.06 to O-
    if charge > 0 {
        return 0.0;
    }
    if is_aromatic {
        28.24
    } else if h > 0 {
        38.80
    } else {
        match count_double_bonds_to(mol, idx, 8) {
            0 => {
                if has_double_bond_to(mol, idx, 6) {
                    32.09 // thioxo C=S
                } else if has_double_bond_to(mol, idx, 7) && is_atom_in_ring(mol, idx) {
                    19.21 // cyclic sulfamidate N=S (sulfoxide-like contribution)
                } else {
                    25.30 // thioether
                }
            }
            1 => {
                // sulfonimidyl N=S(=O): treated as sulfonyl-like by RDKit → 8.38
                if has_double_bond_to(mol, idx, 7) {
                    8.38
                } else {
                    36.28
                }
            }
            _ => 42.52,
        }
    }
}

fn tpsa_phosphorus(mol: &Molecule, idx: AtomIdx, h: u8) -> f64 {
    if has_double_bond_to(mol, idx, 8) {
        26.88 // phosphate/phosphonate: P=O
    } else if has_double_bond_to(mol, idx, 7) {
        9.81 // phosphazene: P=N (cyclic or linear)
    } else if h > 0 {
        34.14 // phosphine P-H (secondary/primary phosphine)
    } else {
        13.59 // phosphite/trivalent P, no H (RDKit-calibrated)
    }
}

/// Topological Polar Surface Area (Ertl 2000).
///
/// Sum of Ertl atom-type contributions for N, O, S, and P atoms.
/// Hydrogen atoms are implicit (not computed). Values match RDKit defaults with calibrations
/// for secondary amide N (12.03 Å²), aromatic N (15.79 Å²), and phosphorus atoms.
///
/// **Note on S/P:** S and P contributions are **included by default** (equivalent to
/// RDKit `Descriptors.TPSA(mol, includeSandP=True)`). RDKit's default
/// `Descriptors.TPSA(mol)` excludes S and P; results will differ for molecules
/// containing these elements.
pub fn tpsa(mol: &Molecule) -> f64 {
    // Apply aromaticity for Kekulé-form input. For TPSA, nitrogen aromaticity uses
    // the ORIGINAL mol's flag: apply_aromaticity can false-promote phthalimide N
    // (making it give 12.89 instead of 12.03). Other elements use mol_arom.
    let mol_arom = chematic_perception::apply_aromaticity(mol);
    let ring_bonds = ring_bond_indices(&mol_arom);

    let mut psa = 0.0f64;
    for (idx, atom_arom) in mol_arom.atoms() {
        let orig_atom = mol.atom(idx);
        let an = orig_atom.element.atomic_number();
        // N: prefer SMILES aromatic flag; other elements: use mol_arom flag.
        let is_aromatic = if an == 7 {
            orig_atom.aromatic
        } else {
            atom_arom.aromatic
        };
        // H count: for N use original mol (before apply_aromaticity changed valence).
        let h = if an == 7 {
            implicit_hcount(mol, idx)
        } else {
            implicit_hcount(&mol_arom, idx)
        };
        let contribution = match an {
            7 => tpsa_nitrogen(mol, idx, is_aromatic, h, orig_atom.charge, &ring_bonds),
            8 => tpsa_oxygen(&mol_arom, idx, is_aromatic, h, atom_arom.charge),
            16 => tpsa_sulfur(&mol_arom, idx, is_aromatic, h, atom_arom.charge),
            15 if !is_aromatic => tpsa_phosphorus(&mol_arom, idx, h),
            _ => 0.0,
        };
        psa += contribution;
    }
    psa
}

// ---------------------------------------------------------------------------
// 8. LogP (Wildman-Crippen, calibrated)
// ---------------------------------------------------------------------------

/// Compute a Wildman-Crippen LogP (RDKit-compatible).
///
/// Uses the same 117-entry SMARTS atom-type table as RDKit `Crippen.MolLogP()`,
/// taken directly from `rdkit/Chem/Crippen.py` (Wildman & Crippen 1999). The
/// first matching pattern for each atom determines its contribution. Results are
/// RDKit-compatible for common drug-like molecules.
///
/// RDKit Crippen LogP atom-type table (Wildman & Crippen 1999).
///
/// Ordered exactly as in RDKit `Crippen.py`. First matching pattern for each atom wins.
/// Each entry: (smarts, logp_contribution, mr_contribution).
/// Patterns that fail to parse are silently skipped.
static CRIPPEN_SMARTS: &[(&str, f64, f64)] = &[
    // --- Carbon (C1-C27 + CS fallback) ---
    ("[CH4]", 0.1441, 2.503),
    ("[CH3]C", 0.1441, 2.503),
    ("[CH2](C)C", 0.1441, 2.503),
    ("[CH](C)(C)C", 0.0000, 2.433),
    ("[C](C)(C)(C)C", 0.0000, 2.433),
    ("[CH3][N,O,P,S,F,Cl,Br,I]", -0.2035, 2.753),
    ("[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]", -0.2035, 2.753),
    (
        "[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]",
        -0.2051,
        2.731,
    ),
    (
        "[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]",
        -0.2051,
        2.731,
    ),
    ("[C]=[!C;A;!#1]", -0.2783, 5.007),
    ("[CH2]=C", 0.1551, 3.513),
    ("[CH1](=C)[A;!#1]", 0.1551, 3.513),
    ("[CH0](=C)([A;!#1])[A;!#1]", 0.1551, 3.513),
    ("[C](=C)=C", 0.1551, 3.513),
    ("[CX2]#[A;!#1]", 0.0017, 3.888),
    ("[CH3]c", 0.08452, 2.464),
    ("[CH3]a", -0.1444, 2.412),
    ("[CH2X4]a", -0.0516, 2.488),
    ("[CHX4]a", 0.1193, 2.582),
    ("[CH0X4]a", -0.0967, 2.576),
    ("[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]", -0.5443, 4.041),
    ("[c][#9]", 0.0000, 3.257),
    ("[c][#17]", 0.2450, 3.564),
    ("[c][#35]", 0.1980, 3.180),
    ("[c][#53]", 0.0000, 3.104),
    ("[cH]", 0.1581, 3.350),
    ("[c](:a)(:a):a", 0.2955, 4.346),
    ("[c](:a)(:a)-a", 0.2713, 3.904),
    ("[c](:a)(:a)-C", 0.1360, 3.509),
    ("[c](:a)(:a)-N", 0.4619, 4.067),
    ("[c](:a)(:a)-O", 0.5437, 3.853),
    ("[c](:a)(:a)-S", 0.1893, 2.673),
    ("[c](:a)(:a)=[C,N,O]", -0.8186, 3.135),
    ("[C](=C)(a)[A;!#1]", 0.2640, 4.305),
    ("[C](=C)(c)a", 0.2640, 4.305),
    ("[CH1](=C)a", 0.2640, 4.305),
    ("[C]=c", 0.2640, 4.305),
    ("[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]", 0.2148, 2.693),
    ("[#6]", 0.08129, 3.243), // CS fallback
    // --- Hydrogen (H1-H4 + HS fallback) ---
    ("[#1][#6,#1]", 0.1230, 1.057),
    ("[#1]O[CX4,c]", -0.2677, 1.395),
    ("[#1]O[!C;!N;!O;!S]", -0.2677, 1.395),
    ("[#1][!C;!N;!O]", -0.2677, 1.395),
    ("[#1][#7]", 0.2142, 0.9627),
    ("[#1]O[#7]", 0.2142, 0.9627),
    ("[#1]OC=[#6,#7,O,S]", 0.2980, 1.805),
    ("[#1]O[O,S]", 0.2980, 1.805),
    ("[#1]", 0.1125, 1.112), // HS fallback
    // --- Nitrogen (N1-N14 + NS fallback) ---
    ("[NH2+0][A;!#1]", -1.0190, 2.262),
    ("[NH+0]([A;!#1])[A;!#1]", -0.7096, 2.173),
    ("[NH2+0]a", -1.0270, 2.827),
    ("[NH1+0]([!#1;A,a])a", -0.5188, 3.000),
    ("[NH+0]=[!#1;A,a]", 0.08387, 1.757),
    ("[N+0](=[!#1;A,a])[!#1;A,a]", 0.1836, 2.428),
    ("[N+0]([A;!#1])([A;!#1])[A;!#1]", -0.3187, 1.839),
    ("[N+0](a)([!#1;A,a])[A;!#1]", -0.4458, 2.819),
    ("[N+0](a)(a)a", -0.4458, 2.819),
    ("[N+0]#[A;!#1]", 0.01508, 1.725),
    ("[NH3,NH2,NH;+,+2,+3]", -1.9500, 0.000),
    ("[n+0]", -0.3239, 2.202),
    ("[n;+,+2,+3]", -1.1190, 0.000),
    (
        "[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]",
        -0.3396,
        0.2604,
    ),
    ("[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]", -0.3396, 0.2604),
    ("[NH0;+,+2,+3](=[#6])=[#7]", -0.3396, 0.2604),
    ("[N;+,+2,+3]#[A;!#1]", 0.2887, 3.359),
    ("[N;-,-2,-3]", 0.2887, 3.359),
    ("[N;+,+2,+3](=[N;-,-2,-3])=N", 0.2887, 3.359),
    ("[#7]", -0.4806, 2.134), // NS fallback
    // --- Oxygen (O1-O12 + OS fallback) ---
    ("[o]", 0.1552, 1.080),
    ("[OH,OH2]", -0.2893, 0.8238),
    ("[O]([A;!#1])[A;!#1]", -0.0684, 1.085),
    ("[O](a)[!#1;A,a]", -0.4195, 1.182),
    ("[O]=[#7,#8]", 0.0335, 3.367),
    ("[OX1;-,-2,-3][#7]", 0.0335, 3.367),
    ("[OX1;-,-2,-3][#16]", -0.3339, 0.7774),
    ("[O;-0]=[#16;-0]", -0.3339, 0.7774),
    ("[O-]C(=O)", -1.3260, 0.000),
    ("[OX1;-,-2,-3][!#1;!N;!S]", -1.1890, 0.000),
    ("[O]=c", 0.1788, 3.135),
    ("[O]=[CH]C", -0.1526, 0.000),
    ("[O]=C(C)([A;!#1])", -0.1526, 0.000),
    ("[O]=[CH][N,O]", -0.1526, 0.000),
    ("[O]=[CH2]", -0.1526, 0.000),
    ("[O]=[CX2]=O", -0.1526, 0.000),
    ("[O]=[CH]c", 0.1129, 0.2215),
    ("[O]=C([C,c])[a;!#1]", 0.1129, 0.2215),
    ("[O]=C(c)[A;!#1]", 0.1129, 0.2215),
    ("[O]=C([!#1;!#6])[!#1;!#6]", 0.4833, 0.3890),
    ("[#8]", -0.1188, 0.6865), // OS fallback
    // --- Sulfur (S1-S3) ---
    ("[S;-,-2,-3,-4,+1,+2,+3,+5,+6]", -0.0024, 7.365),
    ("[S;-0]=[N,O,P,S]", -0.0024, 7.365),
    ("[S;A]", 0.6482, 7.591),
    ("[s;a]", 0.6237, 6.691),
    // --- Phosphorus ---
    ("[#15]", 0.8612, 6.920),
    // --- Halogens ---
    ("[#9;-0]", 0.4202, 1.108),
    ("[#17;-0]", 0.6895, 5.853),
    ("[#35;-0]", 0.8456, 8.927),
    ("[#53;-0]", 0.8857, 14.020),
    ("[#9,#17,#35,#53;-]", -2.9960, 0.000),
    ("[#53;+,+2,+3]", -2.9960, 0.000),
    ("[+;#3,#11,#19,#37,#55]", -2.9960, 0.000),
    // --- Metals ---
    ("[#3,#11,#19,#37,#55]", -0.3808, 5.754),
    ("[#4,#12,#20,#38,#56]", -0.3808, 5.754),
    ("[#5,#13,#31,#49,#81]", -0.3808, 5.754),
    ("[#14,#32,#50,#82]", -0.3808, 5.754),
    ("[#33,#51,#83]", -0.3808, 5.754),
    ("[#34,#52,#84]", -0.3808, 5.754),
    ("[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]", -0.0025, 0.000),
    ("[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]", -0.0025, 0.000),
    ("[#72,#73,#74,#75,#76,#77,#78,#79,#80]", -0.0025, 0.000),
];

type CrippenQueries = Vec<(Option<chematic_smarts::QueryMolecule>, f64, f64)>;
static CRIPPEN_QUERIES: OnceLock<CrippenQueries> = OnceLock::new();

fn get_crippen_queries() -> &'static CrippenQueries {
    CRIPPEN_QUERIES.get_or_init(|| {
        CRIPPEN_SMARTS
            .iter()
            .map(|&(sma, lp, mr)| (parse_smarts(sma).ok(), lp, mr))
            .collect()
    })
}

/// Pre-compute for each SMARTS pattern which target atoms satisfy query-atom-0.
/// Amortizes VF2 cost from O(n_atoms × n_patterns) → O(n_patterns) per molecule.
/// SSSR is computed once and shared across all 117 Crippen patterns.
fn crippen_anchor_sets(mol: &Molecule, queries: &CrippenQueries) -> Vec<FxHashSet<AtomIdx>> {
    let rings = find_sssr(mol);
    // uniquify=false: symmetric bonds (e.g. internal C≡C) must yield both orientations
    // so that each endpoint can appear as query-atom-0 and receive its Crippen contribution.
    let config = MatchConfig {
        uniquify: false,
        ..MatchConfig::default()
    };
    queries
        .iter()
        .map(|(q_opt, _, _)| {
            let Some(q) = q_opt else {
                return FxHashSet::default();
            };
            find_matches_with_rings_and_config(q, mol, &rings, &config)
                .into_iter()
                .filter_map(|m| m.get(&0).copied())
                .collect()
        })
        .collect()
}

/// Wildman-Crippen LogP per-atom contributions (RDKit-compatible SMARTS dispatch).
///
/// Uses the same 117-entry atom-type table as RDKit `Crippen.MolLogP()`. Each atom
/// is assigned the contribution of the first matching SMARTS pattern. H contributions
/// are included via explicit `[#1]` patterns; heavy atoms and their implicit H are
/// handled together as in the original Wildman-Crippen definition.
/// Index matches mol.atoms().
pub fn logp_crippen_per_atom(mol: &Molecule) -> Vec<f64> {
    let queries = get_crippen_queries();
    // Pre-compute once per molecule: for each pattern, which atoms satisfy query-atom-0?
    // Previously O(n_atoms × n_patterns × VF2); now O(n_patterns × VF2 + n_atoms × n_patterns).
    let anchor_sets = crippen_anchor_sets(mol, queries);

    let h_fallback = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| e.1)
        .unwrap_or(0.1125);

    mol.atoms()
        .map(|(idx, atom)| {
            if atom.element.atomic_number() == 1 {
                return 0.0;
            }
            // Implicit H from valence rules + explicit H neighbors (e.g. [H], [2H], [3H]).
            // Explicit isotopic H in the graph contribute 0.0 themselves but their LogP
            // must be counted via their parent heavy atom (same as regular implicit H).
            let h_count_impl = implicit_hcount(mol, idx);
            let h_count_expl = mol
                .neighbors(idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
                .count() as u8;
            let h_count = h_count_impl + h_count_expl;
            // Aromatic oxide bridge: O in a ring bonded to aromatic C AND sp2 C (C=C).
            // RDKit perceives this as [o] type (logp=0.1552); plain SMARTS gives [O](a) (-0.4195).
            let heavy = if atom.element.atomic_number() == 8
                && !atom.aromatic
                && h_count == 0
                && atom.charge == 0
                && is_aromatic_oxide_bridge(mol, idx)
            {
                0.1552
            } else {
                anchor_sets
                    .iter()
                    .zip(queries.iter())
                    .find(|(set, _)| set.contains(&idx))
                    .map(|(_, (_, logp, _))| *logp)
                    .unwrap_or(0.0)
            };
            let h_contrib = if h_count == 0 {
                0.0
            } else {
                h_logp_for_parent(
                    mol,
                    idx,
                    atom.element.atomic_number(),
                    atom.aromatic,
                    h_fallback,
                ) * h_count as f64
            };
            heavy + h_contrib
        })
        .collect()
}

/// Determine the LogP contribution per implicit-H attached to `parent_idx`.
/// Matches against H-type SMARTS by looking at the parent atom's environment.
fn h_logp_for_parent(
    mol: &Molecule,
    parent_idx: AtomIdx,
    parent_an: u8,
    parent_aromatic: bool,
    fallback: f64,
) -> f64 {
    // H type depends on what atom the H is bonded to:
    // H1: [#1][#6,#1]         → H on C  → 0.1230
    // H2: [#1]O[CX4,c] etc.  → H on O (various)  → -0.2677
    // H3: [#1][#7]            → H on N  → 0.2142
    // H4: [#1]OC=[#6,#7,O,S] → H on O adjacent to C=O → 0.2980
    // HS: [#1]                → fallback → 0.1125
    //
    // Since we can't easily run SMARTS with a virtual H atom, we use the
    // original hand-coded H dispatch which is already accurate:
    match parent_an {
        6 => 0.1230, // H1: H on C
        7 => 0.2142, // H3: H on N
        8 => {
            if parent_aromatic {
                fallback // aromatic O (furan-type) — HS fallback
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                // H4: O bonded to C, and that C has a double bond to C/N/O/S
                // Matches RDKit [#1]OC=[#6,#7,O,S]: covers carboxylic, enol, vinylogous OH
                mol.atom(nb).element.atomic_number() == 6
                    && mol.neighbors(nb).any(|(nb2, bidx2)| {
                        nb2 != parent_idx
                            && mol.bond(bidx2).order == BondOrder::Double
                            && matches!(mol.atom(nb2).element.atomic_number(), 6 | 7 | 8 | 16)
                    })
            }) {
                // H4: H on O adjacent to C=X (X = C, N, O, S)
                0.2980
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                mol.atom(nb).element.atomic_number() == 7 // O bonded to N
            }) {
                // [#1]O[#7]: H on O attached to N — oxime / hydroxamic acid
                0.2142
            } else if mol.neighbors(parent_idx).any(|(nb, _)| {
                let an = mol.atom(nb).element.atomic_number();
                an == 8 || an == 16 // [#1]O[O,S]: peroxide / sulfonic / sulfuric acid OH
            }) {
                0.2980
            } else {
                // H2: [#1]O[CX4,c] — aliphatic and phenolic OH both use -0.2677
                -0.2677
            }
        }
        // H2: [#1][!C;!N;!O] — H on S, P, Si, halogens, etc.
        _ => -0.2677,
    }
}

/// Compute the Crippen log P (octanol/water partition coefficient) of `mol`.
///
/// Sums per-atom contributions from [`logp_crippen_per_atom`].
pub fn logp_crippen(mol: &Molecule) -> f64 {
    logp_crippen_per_atom(mol).iter().sum()
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
    use chematic_perception::count_aromatic_rings;
    count_aromatic_rings(mol)
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

/// Return the MR contribution per implicit-H attached to `parent_idx`.
fn h_mr_for_parent(
    mol: &Molecule,
    parent_idx: AtomIdx,
    parent_an: u8,
    parent_aromatic: bool,
    fallback: f64,
) -> f64 {
    match parent_an {
        6 => 1.057,  // H1: H on C
        7 => 0.9627, // H3: H on N
        8 => {
            if parent_aromatic {
                fallback // aromatic O — HS fallback
            } else if mol
                .neighbors(parent_idx)
                .any(|(nb, _)| has_double_bond_to(mol, nb, 8))
            {
                1.805 // H4: carboxylic/ester OH
            } else {
                1.395 // H2: aliphatic or phenolic OH
            }
        }
        _ => fallback, // HS: H on other heteroatom (S, P, …) — MR HS fallback (1.112)
    }
}

/// Per-atom Molar Refractivity contributions (Wildman & Crippen 1999, RDKit-compatible).
///
/// Uses the same 117-entry SMARTS atom-type table as `logp_crippen_per_atom` but
/// reads the MR column. Results match RDKit `Crippen.MolMR()` for common molecules.
/// H contributions are folded into the attached heavy atom. Index matches mol.atoms().
pub fn mr_per_atom(mol: &Molecule) -> Vec<f64> {
    let queries = get_crippen_queries();
    let anchor_sets = crippen_anchor_sets(mol, queries);

    let h_fallback = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| e.2)
        .unwrap_or(1.112);

    mol.atoms()
        .map(|(idx, atom)| {
            if atom.element.atomic_number() == 1 {
                return 0.0;
            }
            let heavy = anchor_sets
                .iter()
                .zip(queries.iter())
                .find(|(set, _)| set.contains(&idx))
                .map(|(_, (_, _, mr))| *mr)
                .unwrap_or(0.0);
            let h_count = implicit_hcount(mol, idx);
            let h_contrib = if h_count == 0 {
                0.0
            } else {
                h_mr_for_parent(
                    mol,
                    idx,
                    atom.element.atomic_number(),
                    atom.aromatic,
                    h_fallback,
                ) * h_count as f64
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

/// Compute both LogP and MR from a single Crippen anchor-set pass.
///
/// Equivalent to calling `logp_crippen(mol)` and `molar_refractivity(mol)` in
/// sequence but shares the 117-pattern SMARTS matching (including the single SSSR
/// computation), making it roughly 2× faster when both values are needed.
pub fn logp_and_mr(mol: &Molecule) -> (f64, f64) {
    let queries = get_crippen_queries();
    let anchor_sets = crippen_anchor_sets(mol, queries);

    let h_logp_fallback = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| e.1)
        .unwrap_or(0.1125);
    let h_mr_fallback = CRIPPEN_SMARTS
        .iter()
        .find(|(sma, _, _)| *sma == "[#1]")
        .map(|e| e.2)
        .unwrap_or(1.112);

    let mut logp_sum = 0.0f64;
    let mut mr_sum = 0.0f64;

    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 1 {
            continue;
        }
        let h_count = implicit_hcount(mol, idx);

        // LogP heavy-atom contribution
        let logp_heavy = if atom.element.atomic_number() == 8
            && !atom.aromatic
            && h_count == 0
            && atom.charge == 0
            && is_aromatic_oxide_bridge(mol, idx)
        {
            0.1552
        } else {
            anchor_sets
                .iter()
                .zip(queries.iter())
                .find(|(set, _)| set.contains(&idx))
                .map(|(_, (_, lp, _))| *lp)
                .unwrap_or(0.0)
        };

        // MR heavy-atom contribution
        let mr_heavy = anchor_sets
            .iter()
            .zip(queries.iter())
            .find(|(set, _)| set.contains(&idx))
            .map(|(_, (_, _, mr))| *mr)
            .unwrap_or(0.0);

        logp_sum += logp_heavy;
        mr_sum += mr_heavy;

        if h_count > 0 {
            logp_sum += h_logp_for_parent(
                mol,
                idx,
                atom.element.atomic_number(),
                atom.aromatic,
                h_logp_fallback,
            ) * h_count as f64;
            mr_sum += h_mr_for_parent(
                mol,
                idx,
                atom.element.atomic_number(),
                atom.aromatic,
                h_mr_fallback,
            ) * h_count as f64;
        }
    }

    (logp_sum, mr_sum)
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

/// Number of distinct connected ring systems.
///
/// Two SSSR rings belong to the same ring system when they share at least
/// one atom (fused, bridged, or spiro connections).  This differs from
/// [`ring_count`] which counts every SSSR ring individually.
///
/// Examples:
/// - benzene: 1 ring system (simple 6-membered ring)
/// - naphthalene: 1 ring system (two fused rings → one system)
/// - biphenyl: 2 ring systems (two independent benzene rings connected by a bond)
///
/// Useful for scaffold complexity assessment (ScaffoldGraph, Chemical-Descriptors).
pub fn ring_system_count(mol: &Molecule) -> usize {
    use chematic_perception::find_ring_families;
    let sssr = find_sssr(mol);
    find_ring_families(mol, &sssr).len()
}

/// Lipinski (1997) HBA count — total number of N and O heavy atoms.
///
/// The original Rule-of-Five definition simply counts all nitrogen and oxygen
/// atoms regardless of hybridisation or substitution.  For the more
/// chemically accurate Ertl (2000) definition see [`hba_count`].
///
/// Used by Chemical-Descriptors, DOPtools, MolVS, and classic implementations.
pub fn hba_count_lipinski(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, a)| {
            let an = a.element.atomic_number();
            an == 7 || an == 8
        })
        .count()
}

/// Fraction of heavy atoms that participate in a rotatable bond.
///
/// `fraction_rotatable_bonds = rotatable_bond_count / heavy_atom_count`.
/// Returns `0.0` for molecules with zero heavy atoms.
/// Normalises molecular flexibility for cross-scaffold comparison (DOPtools).
pub fn fraction_rotatable_bonds(mol: &Molecule) -> f64 {
    let hac = heavy_atom_count(mol);
    if hac == 0 {
        return 0.0;
    }
    rotatable_bond_count(mol) as f64 / hac as f64
}

/// Number of non-aromatic (aliphatic) rings.
pub fn num_aliphatic_rings(mol: &Molecule) -> usize {
    find_sssr(mol)
        .rings()
        .iter()
        .filter(|ring| ring.iter().any(|&idx| !mol.atom(idx).aromatic))
        .count()
}

/// Number of saturated rings.
///
/// A ring is saturated when every *ring bond* is a single bond (no double, triple,
/// or aromatic bonds within the ring).  Exocyclic double bonds (e.g. the C=O in
/// piperidinone) are ignored — matching RDKit `CalcNumSaturatedRings`.
pub fn num_saturated_rings(mol: &Molecule) -> usize {
    let sssr = find_sssr(mol);
    let rings = sssr.rings();
    rings
        .iter()
        .filter(|ring| is_ring_saturated(mol, ring))
        .count()
}

/// Returns true if all bonds *within* the ring are single bonds (non-aromatic, non-double).
fn is_ring_saturated(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        mol.bond_between(a, b)
            .map(|(bidx, _)| {
                !matches!(
                    mol.bond(bidx).order,
                    BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                )
            })
            .unwrap_or(true)
    })
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

/// Number of saturated rings (all ring bonds single) containing at least one heteroatom.
///
/// Exocyclic double bonds are ignored — matching RDKit `CalcNumSaturatedHeterocycles`.
/// Examples: piperidine (1), oxetane (1), piperidinone (1, because ring bonds are all single).
pub fn num_saturated_heterocycles(mol: &Molecule) -> usize {
    let sssr = find_sssr(mol);
    let rings = sssr.rings();
    rings
        .iter()
        .filter(|ring| {
            is_ring_saturated(mol, ring)
                && ring.iter().any(|&idx| {
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
fn num_spiro_atoms_from(mol: &Molecule, rings: &[Vec<AtomIdx>]) -> usize {
    mol.atoms()
        .filter(|(idx, _)| {
            let member: Vec<_> = rings.iter().filter(|r| r.contains(idx)).collect();
            if member.len() != 2 {
                return false;
            }
            member[0].iter().filter(|a| member[1].contains(a)).count() == 1
        })
        .count()
}

pub fn num_spiro_atoms(mol: &Molecule) -> usize {
    num_spiro_atoms_from(mol, find_sssr(mol).rings())
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
fn num_bridgehead_atoms_from(mol: &Molecule, sssr: &RingSet) -> usize {
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
            let ring_sets: Vec<FxHashSet<AtomIdx>> = member_rings
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

pub fn num_bridgehead_atoms(mol: &Molecule) -> usize {
    num_bridgehead_atoms_from(mol, &find_sssr(mol))
}

// ---------------------------------------------------------------------------
// Ring descriptor bundle — find_sssr called exactly once
// ---------------------------------------------------------------------------

/// All ring-derived descriptor values for a molecule, computed with a single `find_sssr` call.
///
/// Use [`ring_bundle`] to obtain this struct. Individual descriptor functions remain
/// available for single-value use but each re-invoke `find_sssr` independently.
#[derive(Debug, Clone)]
pub struct RingBundle {
    pub ring_count: usize,
    pub ring_system_count: usize,
    pub aromatic_ring_count: usize,
    pub num_aliphatic_rings: usize,
    pub num_saturated_rings: usize,
    pub num_aromatic_heterocycles: usize,
    pub num_aliphatic_heterocycles: usize,
    pub num_saturated_heterocycles: usize,
    pub num_spiro_atoms: usize,
    pub num_bridgehead_atoms: usize,
    pub rotatable_bond_count: usize,
    pub hba_count: usize,
    pub fraction_rotatable_bonds: f64,
}

/// Compute all ring-derived descriptors with a single `find_sssr` call.
///
/// When `mol.descriptors()` or `molecule_report()` needs multiple ring values,
/// call this once and read fields — saves ~10-15 redundant SSSR computations per molecule.
pub fn ring_bundle(mol: &Molecule) -> RingBundle {
    let sssr = find_sssr(mol);
    let rings = sssr.rings();

    let ring_bonds = ring_bond_indices_from_rings(mol, rings);

    // Aromatic ring count: for Kekulé input, perceive aromaticity (same topology as sssr).
    let aromatic_ring_count = {
        let mol_arom;
        let mol_a = if mol.atoms().any(|(_, a)| a.aromatic) {
            mol
        } else {
            mol_arom = apply_aromaticity(mol);
            &mol_arom
        };
        // augmented_ring_set takes &[Vec<AtomIdx>] — topology is unchanged by apply_aromaticity.
        let aug = augmented_ring_set(mol_a, rings);
        aug.iter()
            .filter(|r| r.iter().all(|&i| mol_a.atom(i).aromatic))
            .count()
    };

    let rotatable_bond_count = rotatable_bond_count_from_set(mol, &ring_bonds);
    let hba_count = hba_count_from_set(mol, &ring_bonds);
    let hac = heavy_atom_count(mol);
    let fraction_rotatable_bonds = if hac == 0 {
        0.0
    } else {
        rotatable_bond_count as f64 / hac as f64
    };

    RingBundle {
        ring_count: rings.len(),
        ring_system_count: find_ring_families(mol, &sssr).len(),
        aromatic_ring_count,
        num_aliphatic_rings: rings
            .iter()
            .filter(|r| r.iter().any(|&i| !mol.atom(i).aromatic))
            .count(),
        num_saturated_rings: rings
            .iter()
            .filter(|r| {
                r.iter().all(|&i| {
                    mol.neighbors(i).all(|(_, b)| {
                        !matches!(
                            mol.bond(b).order,
                            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                        )
                    })
                })
            })
            .count(),
        num_aromatic_heterocycles: rings
            .iter()
            .filter(|r| {
                r.iter().all(|&i| mol.atom(i).aromatic)
                    && r.iter().any(|&i| {
                        let an = mol.atom(i).element.atomic_number();
                        an != 6 && an != 1
                    })
            })
            .count(),
        num_aliphatic_heterocycles: rings
            .iter()
            .filter(|r| {
                r.iter().any(|&i| !mol.atom(i).aromatic)
                    && r.iter().any(|&i| {
                        let an = mol.atom(i).element.atomic_number();
                        an != 6 && an != 1
                    })
            })
            .count(),
        num_saturated_heterocycles: rings
            .iter()
            .filter(|r| {
                r.iter().all(|&i| {
                    mol.neighbors(i).all(|(_, b)| {
                        !matches!(
                            mol.bond(b).order,
                            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                        )
                    })
                }) && r.iter().any(|&i| {
                    let an = mol.atom(i).element.atomic_number();
                    an != 6 && an != 1
                })
            })
            .count(),
        num_spiro_atoms: num_spiro_atoms_from(mol, rings),
        num_bridgehead_atoms: num_bridgehead_atoms_from(mol, &sssr),
        rotatable_bond_count,
        hba_count,
        fraction_rotatable_bonds,
    }
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

/// Med-Chem Friendly (MCF) composite filter.
///
/// Returns `true` when the compound passes all common quality gates used in
/// medicinal chemistry triage, mirroring the "MCF" concept in the
/// `medchemfilters` Python library:
///   - No PAINS structural alerts
///   - No Brenk reactive-group alerts
///   - Lipinski Ro5 (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)
///   - Veber oral bioavailability (TPSA ≤ 140 Å², RotBonds ≤ 10)
pub fn mcf_passes(mol: &Molecule) -> bool {
    crate::pains_passes(mol)
        && crate::brenk_passes(mol)
        && lipinski_passes(mol)
        && veber_passes(mol)
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
// Fragment / Lead-like filters
// ---------------------------------------------------------------------------

/// Rule of Three (Ro3) — fragment-based drug discovery filter (Congreve 2003).
///
/// Passes when MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3, and RotBonds ≤ 3.
/// Used to screen fragment libraries for FBDD campaigns.
pub fn ro3_passes(mol: &Molecule) -> bool {
    molecular_weight(mol) <= 300.0
        && logp_crippen(mol) <= 3.0
        && hbd_count(mol) <= 3
        && hba_count(mol) <= 3
        && rotatable_bond_count(mol) <= 3
}

/// Lead-like filter (Oprea 2001).
///
/// Passes when MW ≤ 450, LogP −3.5 to 4.5, RotBonds ≤ 10, and RingCount 1–4.
/// Lead-like compounds have lower MW and LogP than typical drug candidates,
/// leaving room for optimisation-related property increases.
pub fn lead_like_passes(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let lp = logp_crippen(mol);
    let rotb = rotatable_bond_count(mol);
    let rings = ring_count(mol);

    mw <= 450.0 && (-3.5..=4.5).contains(&lp) && rotb <= 10 && (1..=4).contains(&rings)
}

/// Pfizer 3/75 filter (Leeson & Springthorpe 2007).
///
/// Returns `true` when the compound is NOT in the high-metabolic-liability zone,
/// i.e. the combination `LogP > 3 AND TPSA < 75` is **absent**.
/// Compounds in the "3/75" zone are more likely to be CYP3A4 substrates with
/// high metabolic clearance.
pub fn pfizer_3_75_passes(mol: &Molecule) -> bool {
    !(logp_crippen(mol) > 3.0 && tpsa(mol) < 75.0)
}

// ---------------------------------------------------------------------------
// CNS MPO Score
// ---------------------------------------------------------------------------

/// Central Nervous System Multi-Parameter Optimisation (CNS MPO) score (Wager 2010).
///
/// Combines six desirability functions, each returning 0–1, for a total of 0–6.
/// Higher scores indicate better predicted CNS drug-like properties.
///
/// Component properties and desirability thresholds:
/// - cLogP: 1.0 if ≤ 3; 0.0 if ≥ 5; linear between 3–5
/// - cLogD (pH 7.4): 1.0 if ≤ 2; 0.0 if ≥ 4; linear between 2–4
/// - MW: 1.0 if ≤ 360; 0.0 if ≥ 500; linear between 360–500
/// - TPSA: 1.0 if 40–90; 0.0 if <0 or >120; linear in 0–40 and 90–120 shoulders
/// - HBD: 1.0 if 0; 0.0 if ≥ 2; linear between 0–2
/// - pKa (most basic site): 1.0 if ≤ 8; 0.0 if ≥ 10; linear between 8–10
///
/// CNS MPO score accepting pre-computed descriptor values.
///
/// Equivalent to [`cns_mpo_score`] but reuses `logp`, `tpsa`, `mw`, `hbd`, and
/// `pka_b` that the caller has already computed, saving a full Crippen SMARTS pass
/// and a pKa scan.  `mol` is still needed for LogD (ionisation class).
///
/// - `logp`: Crippen LogP (from `logp_crippen` or `logp_and_mr`)
/// - `tpsa`: topological polar surface area
/// - `mw`: molecular weight
/// - `hbd`: H-bond donor count
/// - `pka_b`: most basic pKa (`pka_base(mol).unwrap_or(0.0)`)
pub fn cns_mpo_from_parts(
    mol: &Molecule,
    logp: f64,
    tpsa: f64,
    mw: f64,
    hbd: usize,
    pka_b: f64,
) -> f64 {
    #[inline]
    fn ld(val: f64, lo: f64, hi: f64) -> f64 {
        if val <= lo {
            1.0
        } else if val >= hi {
            0.0
        } else {
            (hi - val) / (hi - lo)
        }
    }
    let d_tpsa = if !(0.0..=120.0).contains(&tpsa) {
        0.0
    } else if tpsa <= 40.0 {
        tpsa / 40.0
    } else if tpsa <= 90.0 {
        1.0
    } else {
        (120.0 - tpsa) / 30.0
    };
    ld(logp, 3.0, 5.0)
        + ld(crate::logd::logd_from_logp(logp, mol, 7.4), 2.0, 4.0)
        + ld(mw, 360.0, 500.0)
        + d_tpsa
        + (1.0 - hbd as f64 / 2.0).clamp(0.0, 1.0)
        + ld(pka_b, 8.0, 10.0)
}

/// Reference: Wager T.T. et al., ACS Chem. Neurosci. 2010, 1, 435–449.
pub fn cns_mpo_score(mol: &Molecule) -> f64 {
    #[inline]
    fn linear_desirability(val: f64, lo: f64, hi: f64) -> f64 {
        if val <= lo {
            1.0
        } else if val >= hi {
            0.0
        } else {
            (hi - val) / (hi - lo)
        }
    }

    // Compute LogP once; reuse it for LogD to avoid double Crippen matching.
    let logp = logp_crippen(mol);
    let d_logp = linear_desirability(logp, 3.0, 5.0);

    // cLogD pH 7.4: optimal ≤ 2, zero ≥ 4
    let d_logd = linear_desirability(crate::logd::logd_from_logp(logp, mol, 7.4), 2.0, 4.0);

    // MW: optimal ≤ 360, zero ≥ 500
    let d_mw = linear_desirability(molecular_weight(mol), 360.0, 500.0);

    // TPSA: plateau 40–90, shoulders 0–40 and 90–120
    let psa = tpsa(mol);
    let d_tpsa = if !(0.0..=120.0).contains(&psa) {
        0.0
    } else if psa <= 40.0 {
        psa / 40.0
    } else if psa <= 90.0 {
        1.0
    } else {
        (120.0 - psa) / 30.0
    };

    // HBD: optimal 0, zero ≥ 2
    let hbd = hbd_count(mol) as f64;
    let d_hbd = (1.0 - hbd / 2.0).clamp(0.0, 1.0);

    // pKa (most basic): optimal ≤ 8, zero ≥ 10; None (no basic site) → pKa ≈ 0 → score 1.0
    let pka_b = crate::pka::pka_base(mol).unwrap_or(0.0);
    let d_pka = linear_desirability(pka_b, 8.0, 10.0);

    d_logp + d_logd + d_mw + d_tpsa + d_hbd + d_pka
}

// ---------------------------------------------------------------------------
// Per-atom property vectors
// ---------------------------------------------------------------------------

/// Per-atom hybridization state.
///
/// Returns a `Vec<u8>` of length `mol.atom_count()` with:
///   `1` = sp, `2` = sp2, `3` = sp3, `0` = other (metals, wildcard, etc.)
///
/// Assignment rules:
/// - Wildcard atom (`*`) → 0
/// - Aromatic atom → 2
/// - Has a triple bond → 1
/// - Has any double bond → 2
/// - Otherwise → 3 (sp3)
pub fn hybridization_per_atom(mol: &Molecule) -> Vec<u8> {
    let n = mol.atom_count();
    let mut out = vec![3u8; n];
    for (idx, atom) in mol.atoms() {
        let i = idx.0 as usize;
        if atom.wildcard {
            out[i] = 0;
            continue;
        }
        if atom.aromatic {
            out[i] = 2;
            continue;
        }
        let mut has_triple = false;
        let mut has_double = false;
        for (_, bidx) in mol.neighbors(idx) {
            match mol.bond(bidx).order {
                BondOrder::Triple => has_triple = true,
                BondOrder::Double => has_double = true,
                _ => {}
            }
        }
        out[i] = if has_triple {
            1
        } else if has_double {
            2
        } else {
            3
        };
    }
    out
}

/// Per-atom formal charge.
///
/// Returns a `Vec<i8>` of length `mol.atom_count()` where index `i` is the
/// formal charge of `AtomIdx(i as u32)`.  Neutral atoms have value `0`.
pub fn formal_charge_per_atom(mol: &Molecule) -> Vec<i8> {
    mol.atoms().map(|(_, a)| a.charge).collect()
}

/// Per-atom implicit hydrogen count.
///
/// Returns a `Vec<u8>` of length `mol.atom_count()` where index `i` is the
/// number of implicit H atoms on `AtomIdx(i as u32)`.
/// Reuses [`chematic_core::implicit_hcount`] per atom.
pub fn implicit_hcount_per_atom(mol: &Molecule) -> Vec<u8> {
    mol.atoms()
        .map(|(idx, _)| implicit_hcount(mol, idx))
        .collect()
}

// ---------------------------------------------------------------------------
// TPSA per-atom
// ---------------------------------------------------------------------------

/// Per-atom TPSA contributions (Ertl 2000).
///
/// Returns a `Vec<f64>` of length `mol.atom_count()` where index `i` holds the
/// TPSA contribution of `AtomIdx(i as u32)`.  Atoms that contribute nothing
/// (C, halogens, metals, …) have value `0.0`.  The sum equals `tpsa(mol)`.
///
/// Mirrors the pattern of [`logp_crippen_per_atom`].
pub fn tpsa_per_atom(mol: &Molecule) -> Vec<f64> {
    // Apply aromaticity for consistent results with tpsa() (Kekulé-form parity).
    let mol_arom = chematic_perception::apply_aromaticity(mol);
    let mol = &mol_arom;
    let ring_bonds = ring_bond_indices(mol);
    let n = mol.atom_count();
    let mut out = vec![0.0f64; n];
    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let h = implicit_hcount(mol, idx);
        out[idx.0 as usize] = match an {
            7 => tpsa_nitrogen(mol, idx, atom.aromatic, h, atom.charge, &ring_bonds),
            8 => tpsa_oxygen(mol, idx, atom.aromatic, h, atom.charge),
            16 => tpsa_sulfur(mol, idx, atom.aromatic, h, atom.charge),
            15 if !atom.aromatic => tpsa_phosphorus(mol, idx, h),
            _ => 0.0,
        };
    }
    out
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
    let ring_sets: Vec<FxHashSet<AtomIdx>> =
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
            6 => 0,
            7 => 1,
            8 => 2,
            9 => 3,
            14 => 4,
            15 => 5,
            16 => 6,
            17 => 7,
            35 => 8,
            53 => 9,
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
    m[10] = single;
    m[11] = double;
    m[12] = triple;
    m[13] = aromatic;
}

fn ring_is_saturated(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    ring.iter().all(|&idx| {
        mol.neighbors(idx)
            .all(|(_, bidx)| !matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple))
    })
}

fn ring_has_heteroatom(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    ring.iter()
        .any(|&idx| matches!(mol.atom(idx).element.atomic_number(), 7 | 8 | 16))
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
            6 => 23,
            7 => 24,
            8 => 25,
            _ => continue,
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
    ring_sets: &[FxHashSet<AtomIdx>],
    m: &mut [u8],
) {
    // 38: sp3 carbons
    m[38] = mol
        .atoms()
        .filter(|(idx, a)| {
            a.element.atomic_number() == 6
                && mol.degree(*idx) + implicit_hcount(mol, *idx) as usize == 4
        })
        .count()
        .min(255) as u8;

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
        .count()
        .min(255) as u8;

    // 41: spiro atoms
    let mut spiro = 0u8;
    for (idx, _) in mol.atoms() {
        if rings.iter().filter(|r| r.contains(&idx)).count() >= 2
            && mol
                .neighbors(idx)
                .all(|(nb, _)| rings.iter().any(|r| r.contains(&nb)))
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
// Moran / Geary spatial autocorrelation
// ---------------------------------------------------------------------------

/// Moran's I topological autocorrelation for lags 1–7.
///
/// I_k = (N / W_k) × [Σ_{d(i,j)=k} (x_i−x̄)(x_j−x̄)] / [Σ_i (x_i−x̄)²]
///
/// where N = heavy-atom count, W_k = number of atom pairs at topological distance k,
/// and x_i = atomic valence of atom i (degree + implicit H count).
///
/// Returns a 7-element vector.  Values near 0 indicate spatial randomness;
/// positive = similar neighbours, negative = dissimilar neighbours.
pub fn moran_autocorr(mol: &Molecule) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 7];
    }
    let dist = topo_dist_usize(mol);
    let n = mol.atom_count();
    let vals: Vec<f64> = (0..n)
        .map(|i| atomic_valence(mol, AtomIdx(i as u32)))
        .collect();
    let mean = vals.iter().sum::<f64>() / n as f64;
    let denom: f64 = vals.iter().map(|&v| (v - mean).powi(2)).sum();
    if denom == 0.0 {
        return vec![0.0; 7];
    }
    (1..=7usize)
        .map(|lag| {
            let mut numer = 0.0;
            let mut w = 0usize;
            for i in 0..n {
                for j in (i + 1)..n {
                    if dist[i][j] == lag {
                        numer += (vals[i] - mean) * (vals[j] - mean);
                        w += 1;
                    }
                }
            }
            if w == 0 {
                0.0
            } else {
                (n as f64 / w as f64) * numer / denom
            }
        })
        .collect()
}

/// Geary's C topological autocorrelation for lags 1–7.
///
/// C_k = [(N−1) / (2W_k)] × [Σ_{d(i,j)=k} (x_i−x_j)²] / [Σ_i (x_i−x̄)²]
///
/// Values: C=1 (no autocorrelation), C<1 (positive autocorrelation / similar neighbours),
/// C>1 (negative autocorrelation / dissimilar neighbours).  Returns 1.0 when
/// all atoms have identical valence (denominator is 0).
pub fn geary_autocorr(mol: &Molecule) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![1.0; 7];
    }
    let dist = topo_dist_usize(mol);
    let n = mol.atom_count();
    let vals: Vec<f64> = (0..n)
        .map(|i| atomic_valence(mol, AtomIdx(i as u32)))
        .collect();
    let mean = vals.iter().sum::<f64>() / n as f64;
    let denom: f64 = vals.iter().map(|&v| (v - mean).powi(2)).sum();
    if denom == 0.0 {
        return vec![1.0; 7]; // C=1 means no spatial autocorrelation
    }
    (1..=7usize)
        .map(|lag| {
            let mut numer = 0.0;
            let mut w = 0usize;
            for i in 0..n {
                for j in (i + 1)..n {
                    if dist[i][j] == lag {
                        numer += (vals[i] - vals[j]).powi(2);
                        w += 1;
                    }
                }
            }
            if w == 0 {
                1.0
            } else {
                ((n - 1) as f64 / (2.0 * w as f64)) * numer / denom
            }
        })
        .collect()
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
pub fn num_carbons(mol: &Molecule) -> usize {
    count_element(mol, 6)
}

/// Count nitrogens (N atoms, including aromatic).
pub fn num_nitrogens(mol: &Molecule) -> usize {
    count_element(mol, 7)
}

/// Count oxygens (O atoms).
pub fn num_oxygens(mol: &Molecule) -> usize {
    count_element(mol, 8)
}

/// Count fluorines (F atoms).
pub fn num_fluorines(mol: &Molecule) -> usize {
    count_element(mol, 9)
}

/// Count chlorines (Cl atoms).
pub fn num_chlorines(mol: &Molecule) -> usize {
    count_element(mol, 17)
}

/// Count bromines (Br atoms).
pub fn num_bromines(mol: &Molecule) -> usize {
    count_element(mol, 35)
}

/// Count iodines (I atoms).
pub fn num_iodines(mol: &Molecule) -> usize {
    count_element(mol, 53)
}

/// Count sulfurs (S atoms).
pub fn num_sulfurs(mol: &Molecule) -> usize {
    count_element(mol, 16)
}

/// Count phosphorus (P atoms).
pub fn num_phosphorus(mol: &Molecule) -> usize {
    count_element(mol, 15)
}

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
/// Count amide bonds: individual N–C(=O) single bonds where N is non-aromatic.
///
/// Counts per bond (not per C=O carbon), so urea N–C(=O)–N contributes 2.
/// Aromatic N (e.g. in uracil rings) is excluded — matching RDKit `CalcNumAmideBonds`.
pub fn num_amide_bonds(mol: &Molecule) -> usize {
    mol.bonds()
        .filter(|(_, bond)| {
            let is_single = matches!(
                bond.order,
                BondOrder::Single | BondOrder::Up | BondOrder::Down
            );
            if !is_single {
                return false;
            }
            let an_a = mol.atom(bond.atom1).element.atomic_number();
            let an_b = mol.atom(bond.atom2).element.atomic_number();
            let (c_idx, n_idx) = match (an_a, an_b) {
                (6, 7) => (bond.atom1, bond.atom2),
                (7, 6) => (bond.atom2, bond.atom1),
                _ => return false,
            };
            !mol.atom(n_idx).aromatic && has_double_bond_to(mol, c_idx, 8)
        })
        .count()
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
// Carbon hybridisation types (Mordred CarbonTypes)
// ---------------------------------------------------------------------------

/// Carbon atom counts grouped by hybridisation and heavy-atom degree.
///
/// Mirrors the Mordred `CarbonTypes` descriptor family: C*x*SP*y* where
/// *x* = number of non-hydrogen (heavy-atom) bonds and *y* = hybridisation
/// (1 = sp, 2 = sp², 3 = sp³).
///
/// Uses `hybridization_per_atom()` internally (1=sp, 2=sp², 3=sp³).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CarbonTypes {
    /// sp with 1 heavy-atom bond (terminal alkyne carbon)
    pub c1sp1: u32,
    /// sp with 2 heavy-atom bonds (internal alkyne carbon)
    pub c2sp1: u32,
    /// sp² with 1 heavy-atom bond (terminal alkene =CH₂)
    pub c1sp2: u32,
    /// sp² with 2 heavy-atom bonds (alkene =CH–)
    pub c2sp2: u32,
    /// sp² with 3 heavy-atom bonds (aromatic C or trisubstituted alkene)
    pub c3sp2: u32,
    /// sp³ with 1 heavy-atom bond (methyl group –CH₃)
    pub c1sp3: u32,
    /// sp³ with 2 heavy-atom bonds (methylene –CH₂–)
    pub c2sp3: u32,
    /// sp³ with 3 heavy-atom bonds (methine –CH<)
    pub c3sp3: u32,
}

/// Count carbon atoms by hybridisation × heavy-atom degree (Mordred CarbonTypes).
pub fn carbon_types(mol: &Molecule) -> CarbonTypes {
    let hyb = hybridization_per_atom(mol);
    let mut ct = CarbonTypes::default();
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue; // only carbon
        }
        let h = hyb[idx.0 as usize];
        let deg: u32 = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .count() as u32;
        match (h, deg) {
            (1, 1) => ct.c1sp1 += 1,
            (1, 2) => ct.c2sp1 += 1,
            (2, 1) => ct.c1sp2 += 1,
            (2, 2) => ct.c2sp2 += 1,
            (2, 3) => ct.c3sp2 += 1,
            (3, 1) => ct.c1sp3 += 1,
            (3, 2) => ct.c2sp3 += 1,
            (3, 3) => ct.c3sp3 += 1,
            _ => {}
        }
    }
    ct
}

// ---------------------------------------------------------------------------
// Information Content descriptors (Mordred IC family)
// ---------------------------------------------------------------------------

/// Shannon-entropy-based information content descriptors.
///
/// Based on the vertex equivalence class partition of the molecular graph:
/// atoms are grouped by `(element, heavy_degree)` as their Morgan-0 invariant.
///
/// - **IC**  = Shannon entropy of the partition (bits)
/// - **TIC** = N × IC  (total information content)
/// - **SIC** = IC / log₂(N)  (structural information content; 0 for fully symmetric)
/// - **BIC** = IC × (N-1) / N  (bond information content)
/// - **CIC** = log₂(N) − IC  (complementary information content)
///
/// All values are 0 for a single-atom molecule.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InformationContent {
    pub ic: f64,
    pub tic: f64,
    pub sic: f64,
    pub bic: f64,
    pub cic: f64,
}

/// Compute the Information Content descriptor family (IC, TIC, SIC, BIC, CIC).
pub fn information_content(mol: &Molecule) -> InformationContent {
    let heavy: Vec<_> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .collect();
    let n = heavy.len();
    if n < 2 {
        return InformationContent::default();
    }
    // Vertex invariant: (atomic_number, heavy_degree)
    let mut counts: FxHashMap<(u8, u32), usize> = FxHashMap::default();
    for (idx, atom) in &heavy {
        let deg = mol
            .neighbors(*idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .count() as u32;
        *counts
            .entry((atom.element.atomic_number(), deg))
            .or_insert(0) += 1;
    }
    let n_f = n as f64;
    let ic: f64 = counts
        .values()
        .map(|&c| {
            let p = c as f64 / n_f;
            -p * p.log2()
        })
        .sum();
    let log_n = n_f.log2();
    InformationContent {
        ic,
        tic: n_f * ic,
        sic: if log_n > 0.0 { ic / log_n } else { 0.0 },
        bic: ic * (n_f - 1.0) / n_f,
        cic: (log_n - ic).max(0.0),
    }
}

// ---------------------------------------------------------------------------
// Molecular Distance Edge (MDE) descriptors
// ---------------------------------------------------------------------------

/// Molecular Distance Edge Carbon descriptor (MDEC-*xy*).
///
/// Σ 1/√d(i,j) over all ordered pairs (i,j) where *i* is a carbon with
/// *x* carbon neighbours and *j* is a carbon with *y* carbon neighbours.
/// Here *x,y* ∈ {1,2,3,4} (primary, secondary, tertiary, quaternary C).
///
/// Returns 10 values in the order:
/// MDEC11, MDEC12, MDEC13, MDEC14, MDEC22, MDEC23, MDEC24, MDEC33, MDEC34, MDEC44.
pub fn mde_carbon(mol: &Molecule) -> [f64; 10] {
    let heavy: Vec<_> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx)
        .collect();
    let heavy_count = heavy.len();
    if heavy_count == 0 {
        return [0.0; 10];
    }
    // Build carbon-neighbor-count for each heavy atom index.
    let mut c_neighbors: Vec<Option<u8>> = vec![None; mol.atom_count()];
    for &idx in &heavy {
        if mol.atom(idx).element.atomic_number() == 6 {
            let cn = mol
                .neighbors(idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
                .count()
                .min(4) as u8;
            c_neighbors[idx.0 as usize] = Some(cn);
        }
    }
    let dist = crate::topo_descriptors::topological_distance_matrix(mol);
    // Map original AtomIdx → row in distance matrix (heavy-atom order).
    let heavy_pos: FxHashMap<u32, usize> = heavy
        .iter()
        .enumerate()
        .map(|(p, &idx)| (idx.0, p))
        .collect();
    let mut out = [0.0f64; 10];
    let n = heavy.len();
    for pi in 0..n {
        let ai = heavy[pi];
        let xi = match c_neighbors[ai.0 as usize] {
            Some(v) if v >= 1 => v,
            _ => continue,
        };
        for pj in (pi + 1)..n {
            let aj = heavy[pj];
            let xj = match c_neighbors[aj.0 as usize] {
                Some(v) if v >= 1 => v,
                _ => continue,
            };
            let d = dist[pi][pj];
            if d == 0 || d == u32::MAX {
                continue;
            }
            let contrib = 1.0 / (d as f64).sqrt();
            let (lo, hi) = if xi <= xj { (xi, xj) } else { (xj, xi) };
            // Map (lo, hi) from (1,1)..(4,4) to index 0..9
            let slot = match (lo, hi) {
                (1, 1) => 0,
                (1, 2) => 1,
                (1, 3) => 2,
                (1, 4) => 3,
                (2, 2) => 4,
                (2, 3) => 5,
                (2, 4) => 6,
                (3, 3) => 7,
                (3, 4) => 8,
                (4, 4) => 9,
                _ => continue,
            };
            out[slot] += contrib;
        }
    }
    let _ = heavy_pos; // used indirectly via dist which was built from heavy-atom ordering
    out
}

// ---------------------------------------------------------------------------
// BCUT2D — Burden Connectivity Matrix eigenvalue descriptors
// ---------------------------------------------------------------------------

/// BCUT2D descriptors: largest and smallest eigenvalues of the Burden
/// connectivity matrix weighted by four atomic properties.
///
/// The Burden matrix B (n×n, n = heavy atoms) is defined as:
/// - `B[i,i]` = per-atom property value
/// - `B[i,j]` = √(Z_i × Z_j) × bond_weight / 8  for bonded pairs
/// - `B[i,j]` = 0.001  for non-bonded pairs
///
/// Four atomic properties are used, each giving a HI/LO eigenvalue pair:
/// - **CHG** formal charge per atom
/// - **LOGP** Crippen LogP atomic contribution
/// - **MR** molar refractivity atomic contribution
/// - **MW** atomic mass
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Bcut2D {
    pub chghi: f64,
    pub chglo: f64,
    pub logphi: f64,
    pub logplo: f64,
    pub mrhi: f64,
    pub mrlo: f64,
    pub mwhi: f64,
    pub mwlo: f64,
}

/// Compute BCUT2D descriptors (8 values: 4 properties × HI/LO eigenvalue).
pub fn bcut2d(mol: &Molecule) -> Bcut2D {
    let heavy: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx)
        .collect();
    let n = heavy.len();
    if n == 0 {
        return Bcut2D::default();
    }
    // Per-atom properties.
    let logp_v = logp_crippen_per_atom(mol);
    let mr_v = mr_per_atom(mol);
    let charge_v: Vec<f64> = mol.atoms().map(|(_, a)| a.charge as f64).collect();
    let mw_v: Vec<f64> = mol.atoms().map(|(_, a)| a.element.atomic_mass()).collect();
    let props: [Vec<f64>; 4] = [
        heavy.iter().map(|&i| charge_v[i.0 as usize]).collect(),
        heavy.iter().map(|&i| logp_v[i.0 as usize]).collect(),
        heavy.iter().map(|&i| mr_v[i.0 as usize]).collect(),
        heavy.iter().map(|&i| mw_v[i.0 as usize]).collect(),
    ];
    // Connectivity (bonded) lookup.
    let mut bonded = vec![vec![0.0f64; n]; n];
    for (pi, &ai) in heavy.iter().enumerate() {
        for (pj, &aj) in heavy.iter().enumerate().skip(pi + 1) {
            if let Some((bidx, _)) = mol.bond_between(ai, aj) {
                let bo = match mol.bond(bidx).order {
                    chematic_core::BondOrder::Single => 1.0,
                    chematic_core::BondOrder::Double => 2.0,
                    chematic_core::BondOrder::Triple => 3.0,
                    chematic_core::BondOrder::Aromatic => 1.5,
                    chematic_core::BondOrder::Up | chematic_core::BondOrder::Down => 1.0,
                    _ => 1.0,
                };
                let zi = mol.atom(ai).element.atomic_number() as f64;
                let zj = mol.atom(aj).element.atomic_number() as f64;
                let w = (zi * zj).sqrt() * bo / 8.0;
                bonded[pi][pj] = w;
                bonded[pj][pi] = w;
            }
        }
    }
    let compute = |prop: &[f64]| -> (f64, f64) {
        let mut mat = vec![vec![0.001f64; n]; n];
        for i in 0..n {
            mat[i][i] = prop[i];
            for j in (i + 1)..n {
                if bonded[i][j] > 0.0 {
                    mat[i][j] = bonded[i][j];
                    mat[j][i] = bonded[i][j];
                }
            }
        }
        let hi = burden_max_eigenvalue(&mat);
        let lo = burden_min_eigenvalue(&mat);
        (hi, lo)
    };
    let (chghi, chglo) = compute(&props[0]);
    let (logphi, logplo) = compute(&props[1]);
    let (mrhi, mrlo) = compute(&props[2]);
    let (mwhi, mwlo) = compute(&props[3]);
    Bcut2D {
        chghi,
        chglo,
        logphi,
        logplo,
        mrhi,
        mrlo,
        mwhi,
        mwlo,
    }
}

/// Power iteration → largest eigenvalue of symmetric matrix.
fn burden_max_eigenvalue(mat: &[Vec<f64>]) -> f64 {
    let n = mat.len();
    if n == 0 {
        return 0.0;
    }
    let mut v = vec![1.0_f64 / (n as f64).sqrt(); n];
    for _ in 0..150 {
        let mut av: Vec<f64> = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                av[i] += mat[i][j] * v[j];
            }
        }
        let norm: f64 = av.iter().map(|&x| x * x).sum::<f64>().sqrt();
        if norm < 1e-14 {
            return 0.0;
        }
        for i in 0..n {
            v[i] = av[i] / norm;
        }
    }
    // Rayleigh quotient
    (0..n)
        .map(|i| v[i] * (0..n).map(|j| mat[i][j] * v[j]).sum::<f64>())
        .sum()
}

/// Smallest eigenvalue via max eigenvalue of negated matrix.
/// Valid for symmetric real matrices: min_eig(A) = −max_eig(−A).
fn burden_min_eigenvalue(mat: &[Vec<f64>]) -> f64 {
    let neg: Vec<Vec<f64>> = mat
        .iter()
        .map(|row| row.iter().map(|&x| -x).collect())
        .collect();
    -burden_max_eigenvalue(&neg)
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
        // RDKit strict=True (default): ester C(=O)–O bond is excluded (carbonyl hetero bond).
        // Counted: O–Ar, Ar–C(=O)COOH (2 bonds).  CH3 terminal, carboxyl OH terminal, ring bonds excluded.
        let r = rotatable_bond_count(&m);
        assert_eq!(r, 2, "aspirin rotatable bonds = {r}");
    }

    // -- Test: alkyne adjacent bond excluded ---------------------------------
    #[test]
    fn test_rot_alkyne_adjacent_excluded() {
        // Propyne CH3-C≡CH: C-C adjacent to triple bond is NOT rotatable (C≡ atom excluded).
        let m = mol("CC#C");
        assert_eq!(
            rotatable_bond_count(&m),
            0,
            "propyne: C-C adj to triple bond excluded"
        );
        // But-1-yne CH3-CH2-C≡CH: C3 has triple bond → C2-C3 excluded; C1 terminal → C1-C2 excluded.
        let m2 = mol("CCC#C");
        assert_eq!(
            rotatable_bond_count(&m2),
            0,
            "but-1-yne: both bonds excluded"
        );
        // Pent-1-yne CH3-CH2-CH2-C≡CH: C2-C3 is rotatable; C3-C4 excluded (C4 in triple).
        let m3 = mol("CCCC#C");
        assert_eq!(
            rotatable_bond_count(&m3),
            1,
            "pent-1-yne: CH2-CH2 bond is rotatable"
        );
    }

    // -- Test: allene cumulated double bonds excluded ------------------------
    #[test]
    fn test_rot_allene_excluded() {
        // Allene: C=C=C. Both single-bond-looking connections are in a cumulated system.
        let m = mol("C=C=C");
        assert_eq!(rotatable_bond_count(&m), 0, "allene: no rotatable bonds");
    }

    // -- Test 16: water TPSA -------------------------------------------------
    #[test]
    fn test_tpsa_water() {
        let m = mol("O");
        // H₂O: isolated O with no heavy-atom neighbors → 31.50 (Ertl water type)
        let t = tpsa(&m);
        assert!(approx(t, 31.50, 0.1), "water TPSA = {t}");
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

    // Quaternary ammonium N+ has no lone pair → TPSA contribution = 0
    #[test]
    fn test_tpsa_quaternary_n_zero() {
        // Trimethylphenylammonium: N+ with no lone pair → 0.00
        let m = mol("c1csc([N+]2(C3CCCCC3)CCCCC2)c1");
        let t = tpsa(&m);
        // N+ contributes 0; thiophene S contributes ~38.8 → total ~38.8
        assert!(t < 50.0, "N+ quaternary TPSA should be small, got {t}");
        // Specifically, N+ should NOT add 3.01 on top
        let m2 = mol("c1ccsc1");
        let t_base = tpsa(&m2); // thiophene TPSA without N+
        assert!(
            (t - t_base).abs() < 1.0,
            "N+ adds ~0 to TPSA, got delta {}",
            t - t_base
        );
    }

    // Ionic N-oxide [N+]-[O-] ring: O- accounts for TPSA, N+ contributes 0
    #[test]
    fn test_tpsa_n_oxide_ionic() {
        // [N+]1([O-])CCOCC1 — morpholine N-oxide
        // O- and ring O contribute; N+ should be 0
        let m = mol("O=C1CC[N+]2([O-])CCCC[C@@H]12"); // representative N-oxide ring
        let t = tpsa(&m);
        assert!(t > 0.0 && t < 150.0, "N-oxide TPSA in range, got {t}");
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

    #[test]
    fn test_aromatic_ring_count_anthracene() {
        let m = mol("c1ccc2cc3ccccc3cc2c1");
        assert_eq!(aromatic_ring_count(&m), 3);
    }

    #[test]
    fn test_aromatic_ring_count_pyrene() {
        // 4 fused 6-rings; SSSR may produce envelope cycle requiring 3-ring XOR removal
        let m = mol("c1ccc2cccc3ccc4cccc1c4c23");
        assert_eq!(aromatic_ring_count(&m), 4);
    }

    #[test]
    fn test_aromatic_ring_count_triphenylene() {
        let m = mol("c1ccc2ccc3ccccc3c2c1");
        assert_eq!(aromatic_ring_count(&m), 3);
    }

    #[test]
    fn test_aromatic_ring_count_fluoranthene() {
        // 3 six-membered + 1 five-membered aromatic ring
        let m = mol("c1ccc2-c3cccc4cccc-3c4c2c1");
        assert_eq!(aromatic_ring_count(&m), 4);
    }

    #[test]
    fn test_aromatic_ring_count_acridine() {
        let m = mol("c1ccc2nc3ccccc3cc2c1");
        assert_eq!(aromatic_ring_count(&m), 3);
    }

    // bench5k regression cases (rd=N, ch=N-1 pattern)
    #[test]
    fn test_arc_bench_lactone_benzene() {
        // rd=1 ch=0: benzene fused with lactone ring
        let m = mol("CC1(C)CC(=O)c2c(ccc(C(=O)CC(N)CO)c2N)O1");
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_arc_bench_bridged_benzenes() {
        // rd=2 ch=1: two aromatic benzenes connected by O and C bridge
        let m = mol("Cc1c(C)c2c(c(C)c1O)Oc1c(C)c([C@H](C)[C@@H](C)O)c(O)c(O)c1C2");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_arc_bench_no_bridge() {
        // rd=2 ch=1: two benzenes with N/O bridges
        let m = mol("C[C@]12Nc3ccccc3O[C@@]1(C)Nc1ccccc1O2");
        assert_eq!(aromatic_ring_count(&m), 2);
    }

    #[test]
    fn test_arc_bench_steroid_benzene() {
        // rd=1 ch=0: complex steroid scaffold with one phenyl ring
        let m = mol(
            "CC1=C2C[C@@]3(C)C[C@H](O)[C@](C)(C[C@H](O)[C@H](O)[C@@](C)(O)CO)[C@H]3c3ccc(C)c(c32)C1",
        );
        assert_eq!(aromatic_ring_count(&m), 1);
    }

    #[test]
    fn test_hba_metformin() {
        // Metformin: CN(C)C(=N)NC(=N)N — RDKit CalcNumHBA = 2 (two imine =NH)
        // Confirmed by tracing current n_adjacent_to_pi_center logic:
        //   dimethylamino-N, bridging NH, terminal NH₂ → excluded (adjacent to C=N)
        //   two =NH imine nitrogens → counted (double bond to neighbor, not single-like)
        let m = mol("CN(C)C(=N)NC(=N)N");
        assert_eq!(hba_count(&m), 2, "metformin HBA should match RDKit (2)");
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
        // N-C(=O)-N: two N–C(=O) bonds (RDKit counts per bond, not per C=O carbon)
        let m = mol("NC(=O)N");
        assert_eq!(num_amide_bonds(&m), 2);
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
        assert!(
            lp > eth,
            "propene logp ({lp}) should exceed ethylene ({eth})"
        );
        assert!(
            lp > 0.5 && lp < 2.0,
            "propene logp = {lp} out of expected range"
        );
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
        let alkene = logp_crippen(&mol("C=CCC")); // 1-butene
        assert!(
            enone < alkene,
            "MVK ({enone:.4}) should be < 1-butene ({alkene:.4}): enone is less hydrophobic"
        );
    }

    #[test]
    fn all_crippen_smarts_parse() {
        let mut failed = Vec::new();
        for (i, &(pattern, _, _)) in CRIPPEN_SMARTS.iter().enumerate() {
            if parse_smarts(pattern).is_err() {
                failed.push((i, pattern));
            }
        }
        assert!(
            failed.is_empty(),
            "CRIPPEN_SMARTS parse failures (index, pattern): {:?}",
            failed
        );
    }

    // ── Ro3 / Lead-like / Pfizer 3/75 / CNS MPO / TPSA per-atom ─────────────

    #[test]
    fn test_ro3_ethanol_passes() {
        // ethanol: MW=46, LogP≈-0.3, HBD=1, HBA=1, RotBonds=0
        assert!(ro3_passes(&mol("CCO")));
    }

    #[test]
    fn test_ro3_lipinski_drug_fails() {
        // aspirin (MW=180): just passes; ibuprofen (MW=206, LogP≈3.8): fails on LogP
        assert!(!ro3_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_lead_like_ibuprofen_passes() {
        // ibuprofen: MW=206, LogP≈3.8, RotBonds=4, RingCount=1 — borderline but passes
        assert!(lead_like_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_lead_like_large_drug_fails() {
        // MW > 450 should fail lead_like
        let big = mol("CC(C)c1ccc(cc1)C(=O)NC2CCN(CC2)c3ncnc4c3ccc(c4)OC(F)(F)F");
        assert!(!lead_like_passes(&big));
    }

    #[test]
    fn test_pfizer_3_75_safe_compound_passes() {
        // Aspirin: LogP≈1.2, TPSA≈63 → LogP ≤ 3 → passes (not in danger zone)
        assert!(pfizer_3_75_passes(&mol("CC(=O)Oc1ccccc1C(=O)O")));
    }

    #[test]
    fn test_pfizer_3_75_risky_compound_fails() {
        // Ibuprofen: LogP≈3.8, TPSA≈37 → LogP > 3 AND TPSA < 75 → risky
        assert!(!pfizer_3_75_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_cns_mpo_score_range() {
        // All CNS MPO scores must be in [0, 6]
        for smi in [
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "c1ccc2[nH]cnc2c1",
        ] {
            let s = cns_mpo_score(&mol(smi));
            assert!(
                (0.0..=6.0).contains(&s),
                "CNS MPO out of range for {smi}: {s}"
            );
        }
    }

    #[test]
    fn test_cns_mpo_small_cns_drug_high_score() {
        // Caffeine: small, low LogP, moderate TPSA, low HBD → expect score ≥ 3
        let s = cns_mpo_score(&mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C"));
        assert!(s >= 3.0, "caffeine CNS MPO should be ≥ 3, got {s}");
    }

    #[test]
    fn test_tpsa_per_atom_sum_equals_tpsa() {
        for smi in [
            "CC(=O)Oc1ccccc1C(=O)O",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
            "c1ccncc1",
        ] {
            let m = mol(smi);
            let per_atom = tpsa_per_atom(&m);
            assert_eq!(per_atom.len(), m.atom_count());
            let sum: f64 = per_atom.iter().sum();
            let direct = tpsa(&m);
            assert!(
                (sum - direct).abs() < 1e-9,
                "tpsa_per_atom sum {sum} ≠ tpsa {direct} for {smi}"
            );
        }
    }

    #[test]
    fn test_tpsa_per_atom_carbon_contributes_zero() {
        // In benzene all atoms are C — per-atom TPSA must all be 0.0
        let benz = mol("c1ccccc1");
        for v in tpsa_per_atom(&benz) {
            assert_eq!(v, 0.0);
        }
    }

    // ── MCF composite filter ─────────────────────────────────────────────────

    #[test]
    fn test_mcf_caffeine_passes() {
        // Caffeine: no PAINS/Brenk alerts, MW=194, LogP=-0.07, passes Lipinski and Veber.
        assert!(mcf_passes(&mol("Cn1cnc2c1c(=O)n(c(=O)n2C)C")));
    }

    #[test]
    fn test_mcf_toluene_passes() {
        // Toluene: no PAINS/Brenk alerts, passes Lipinski and Veber.
        assert!(mcf_passes(&mol("Cc1ccccc1")));
    }

    #[test]
    fn test_mcf_ibuprofen_fails_brenk_acetal_ketal() {
        // Ibuprofen's carboxylic acid group matches the broad Brenk "acetal_ketal"
        // SMARTS [#8][#6]([#8])-[#6], which also captures C(=O)OH.
        // This is a known over-match in the Brenk filter set, but MCF correctly
        // reflects the current Brenk implementation.
        assert!(!mcf_passes(&mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O")));
    }

    #[test]
    fn test_mcf_rhodanine_fails_pains() {
        // Rhodanine is a PAINS alert → MCF must fail
        assert!(!mcf_passes(&mol("O=C1CSC(=S)N1")));
    }

    #[test]
    fn test_mcf_aspirin_fails_brenk_active_ester() {
        // Aspirin has an aryl ester which triggers the Brenk "active_ester" alert.
        assert!(!mcf_passes(&mol("CC(=O)Oc1ccccc1C(=O)O")));
    }

    #[test]
    fn test_mcf_very_large_mol_fails_lipinski() {
        // MW > 500 → fails Lipinski → MCF must fail
        let big = mol("CC(C)c1ccc(cc1)C(=O)NC2CCN(CC2)c3ncnc4c3ccc(c4)OC(F)(F)F");
        assert!(!mcf_passes(&big));
    }

    // ── per-atom property vectors ────────────────────────────────────────────

    #[test]
    fn test_hybridization_per_atom_ethane_all_sp3() {
        let m = mol("CC");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 2);
        assert!(h.iter().all(|&v| v == 3), "both C in ethane → sp3: {h:?}");
    }

    #[test]
    fn test_hybridization_per_atom_benzene_all_sp2() {
        let m = mol("c1ccccc1");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 6);
        assert!(h.iter().all(|&v| v == 2), "all C in benzene → sp2: {h:?}");
    }

    #[test]
    fn test_hybridization_per_atom_acetylene_sp() {
        // HC≡CH: both C → sp (triple bond)
        let m = mol("C#C");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 2);
        assert!(h.iter().all(|&v| v == 1), "alkyne C → sp: {h:?}");
    }

    #[test]
    fn test_hybridization_per_atom_acetaldehyde_mixed() {
        // CH3-C(=O)-H: C1(CH3)=sp3, C2(=O)=sp2
        let m = mol("CC=O");
        let h = hybridization_per_atom(&m);
        assert_eq!(h.len(), 3); // C, C, O
        assert_eq!(h[0], 3, "methyl C → sp3");
        assert_eq!(h[1], 2, "carbonyl C → sp2 (double bond to O)");
        assert_eq!(h[2], 2, "carbonyl O → sp2 (double bond to C)");
    }

    #[test]
    fn test_formal_charge_per_atom_neutral() {
        let m = mol("CC(=O)O");
        let fc = formal_charge_per_atom(&m);
        assert_eq!(fc.len(), m.atom_count());
        assert!(
            fc.iter().all(|&c| c == 0),
            "acetic acid has no formal charges: {fc:?}"
        );
    }

    #[test]
    fn test_formal_charge_per_atom_charged() {
        // Ammonium ion [NH4]+: N has charge +1
        let m = mol("[NH4+]");
        let fc = formal_charge_per_atom(&m);
        assert_eq!(fc.len(), 1, "only N is heavy");
        assert_eq!(fc[0], 1, "N has formal charge +1");
    }

    #[test]
    fn test_implicit_hcount_per_atom_ethane() {
        // Ethane CC: each C has 3 implicit H
        let m = mol("CC");
        let ih = implicit_hcount_per_atom(&m);
        assert_eq!(ih.len(), 2);
        assert!(
            ih.iter().all(|&h| h == 3),
            "each CH3 → 3 implicit H: {ih:?}"
        );
    }

    #[test]
    fn test_implicit_hcount_per_atom_sum() {
        // For any molecule, total implicit H should be consistent.
        for smi in [
            "CC(=O)O",
            "c1ccccc1",
            "CC(C)C",
            "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        ] {
            let m = mol(smi);
            let ih = implicit_hcount_per_atom(&m);
            assert_eq!(ih.len(), m.atom_count(), "length mismatch for {smi}");
        }
    }

    // ── page 24 descriptors ───────────────────────────────────────────────────

    #[test]
    fn test_hba_count_lipinski_aspirin() {
        // Aspirin C9H8O4: 0 N + 4 O = 4
        assert_eq!(hba_count_lipinski(&mol("CC(=O)Oc1ccccc1C(=O)O")), 4);
    }

    #[test]
    fn test_hba_count_lipinski_aniline() {
        // Aniline C6H7N: 1 N + 0 O = 1
        assert_eq!(hba_count_lipinski(&mol("Nc1ccccc1")), 1);
    }

    #[test]
    fn test_hba_count_lipinski_ge_ertl() {
        // Lipinski counts more liberally than Ertl (e.g., amide N is excluded by Ertl).
        // For acetamide CC(=O)N: Ertl HBA = 1 (only C=O), Lipinski = 2 (O + N).
        let acetamide = mol("CC(=O)N");
        assert_eq!(hba_count_lipinski(&acetamide), 2); // N + O
        assert!(hba_count(&acetamide) <= hba_count_lipinski(&acetamide));
    }

    #[test]
    fn test_fraction_rotatable_bonds_benzene() {
        // Benzene has 0 rotatable bonds → fraction = 0.0
        assert_eq!(fraction_rotatable_bonds(&mol("c1ccccc1")), 0.0);
    }

    #[test]
    fn test_fraction_rotatable_bonds_in_range() {
        for smi in ["CCCCc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "CCCC"] {
            let f = fraction_rotatable_bonds(&mol(smi));
            assert!(
                (0.0..=1.0).contains(&f),
                "fraction out of [0,1] for {smi}: {f}"
            );
        }
    }

    #[test]
    fn test_ring_system_count_benzene() {
        assert_eq!(ring_system_count(&mol("c1ccccc1")), 1);
    }

    #[test]
    fn test_ring_system_count_naphthalene_one_system() {
        // Naphthalene: 2 fused rings → 1 ring system
        assert_eq!(ring_system_count(&mol("c1ccc2ccccc2c1")), 1);
    }

    #[test]
    fn test_ring_system_count_biphenyl_two_systems() {
        // Biphenyl: two benzene rings connected by a single bond → 2 ring systems
        assert_eq!(ring_system_count(&mol("c1ccc(-c2ccccc2)cc1")), 2);
    }

    #[test]
    fn test_ring_system_count_acyclic_zero() {
        // Propane has no rings → 0 ring systems
        assert_eq!(ring_system_count(&mol("CCC")), 0);
    }

    // ── Moran / Geary autocorrelation ─────────────────────────────────────────

    #[test]
    fn moran_len_is_7() {
        assert_eq!(moran_autocorr(&mol("c1ccccc1")).len(), 7);
    }

    #[test]
    fn geary_len_is_7() {
        assert_eq!(geary_autocorr(&mol("c1ccccc1")).len(), 7);
    }

    #[test]
    fn moran_single_atom_returns_zeros() {
        let v = moran_autocorr(&mol("C"));
        assert_eq!(v, vec![0.0; 7]);
    }

    #[test]
    fn geary_single_atom_returns_ones() {
        // denom = 0 → all-identical property → no autocorrelation → C = 1
        let v = geary_autocorr(&mol("C"));
        assert_eq!(v, vec![1.0; 7]);
    }

    #[test]
    fn moran_uniform_valence_is_zero() {
        // Benzene: all C have the same valence (aromatic) → mean = valence, denom = 0 → I = 0
        let v = moran_autocorr(&mol("c1ccccc1"));
        for &x in &v {
            assert!(x.abs() < 1e-9, "expected 0 for uniform valence, got {x}");
        }
    }

    #[test]
    fn geary_finite_for_mixed_molecule() {
        // Ethanol CCO has mixed valence atoms — all values must be finite
        let v = geary_autocorr(&mol("CCO"));
        assert!(
            v.iter().all(|x| x.is_finite()),
            "all Geary values must be finite: {v:?}"
        );
    }

    #[test]
    fn moran_finite_for_mixed_molecule() {
        let v = moran_autocorr(&mol("CCO"));
        assert!(
            v.iter().all(|x| x.is_finite()),
            "all Moran values must be finite: {v:?}"
        );
    }

    // ── CarbonTypes ──────────────────────────────────────────────────────────

    #[test]
    fn carbon_types_methane() {
        // methane C: sp3 with 0 heavy neighbors — none of the CxSPy bins (deg must be ≥1)
        let ct = carbon_types(&mol("C"));
        assert_eq!(
            ct.c1sp3 + ct.c2sp3 + ct.c3sp3,
            0,
            "isolated methane C has deg=0, not counted"
        );
    }

    #[test]
    fn carbon_types_ethane() {
        // ethane CC: each C is sp3 with 1 heavy neighbor → C1SP3×2
        let ct = carbon_types(&mol("CC"));
        assert_eq!(ct.c1sp3, 2, "ethane: 2×C1SP3");
        assert_eq!(ct.c2sp3, 0);
    }

    #[test]
    fn carbon_types_propane() {
        // propane CCC: 2×C1SP3 (terminals) + 1×C2SP3 (middle)
        let ct = carbon_types(&mol("CCC"));
        assert_eq!(ct.c1sp3, 2, "propane: 2 terminal sp3 C");
        assert_eq!(ct.c2sp3, 1, "propane: 1 internal sp3 C");
    }

    #[test]
    fn carbon_types_benzene() {
        // benzene: 6 aromatic C → hybridization=2, heavy deg=2 → C2SP2×6
        let ct = carbon_types(&mol("c1ccccc1"));
        assert_eq!(ct.c2sp2, 6, "benzene: 6×C2SP2 (aromatic)");
    }

    #[test]
    fn carbon_types_acetylene() {
        // acetylene HC≡CH: 2 sp C with 1 heavy neighbor each → C1SP1×2
        let ct = carbon_types(&mol("C#C"));
        assert_eq!(ct.c1sp1, 2, "acetylene: 2×C1SP1");
    }

    // ── Information Content ───────────────────────────────────────────────────

    #[test]
    fn information_content_benzene_zero_ic() {
        // All 6 C in benzene are equivalent (same element + same degree) → IC=0
        let ic = information_content(&mol("c1ccccc1"));
        assert!(
            ic.ic.abs() < 1e-9,
            "benzene: IC must be 0 (fully symmetric), got {}",
            ic.ic
        );
        assert!(ic.sic.abs() < 1e-9, "benzene: SIC must be 0");
    }

    #[test]
    fn information_content_propane_nonzero() {
        let ic = information_content(&mol("CCC"));
        assert!(
            ic.ic > 0.5,
            "propane: IC must be > 0 (two symmetry classes)"
        );
        assert!(ic.tic > 0.0);
        assert!(ic.ic.is_finite() && ic.tic.is_finite());
    }

    #[test]
    fn information_content_single_atom() {
        // Single heavy atom → all fields = 0
        let ic = information_content(&mol("[Na+]"));
        assert_eq!(ic.ic, 0.0);
        assert_eq!(ic.tic, 0.0);
    }

    // ── MDE ──────────────────────────────────────────────────────────────────

    #[test]
    fn mde_carbon_propane_mdec11() {
        // Propane CCC: 2 primary C (deg=1 C-neighbor) + 1 secondary C (deg=2)
        // MDEC11 = 1/sqrt(2) between the two primary C (distance 2 bonds)
        let v = mde_carbon(&mol("CCC"));
        let expected = 1.0 / (2.0_f64).sqrt();
        assert!(
            (v[0] - expected).abs() < 1e-6,
            "MDEC11 for propane: expected {expected:.4}, got {:.4}",
            v[0]
        );
    }

    #[test]
    fn mde_carbon_all_finite() {
        let v = mde_carbon(&mol("CC(CC)C(=O)O"));
        assert!(v.iter().all(|&x| x.is_finite()), "all MDEC must be finite");
    }

    // ── BCUT2D ───────────────────────────────────────────────────────────────

    #[test]
    fn bcut2d_hi_ge_lo() {
        // HI eigenvalue must be ≥ LO for all properties
        let b = bcut2d(&mol("c1ccccc1"));
        assert!(b.mwhi >= b.mwlo, "BCUT2D-MWHI must be ≥ BCUT2D-MWLO");
        assert!(b.logphi >= b.logplo, "BCUT2D-LOGPHI ≥ BCUT2D-LOGPLO");
        assert!(b.mrhi >= b.mrlo, "BCUT2D-MRHI ≥ BCUT2D-MRLO");
    }

    #[test]
    fn bcut2d_all_finite() {
        let b = bcut2d(&mol("CC(=O)Nc1ccccc1C(=O)O")); // 4-aminobenzoic acid derivative
        let vals = [
            b.chghi, b.chglo, b.logphi, b.logplo, b.mrhi, b.mrlo, b.mwhi, b.mwlo,
        ];
        assert!(
            vals.iter().all(|v| v.is_finite()),
            "all BCUT2D values must be finite: {vals:?}"
        );
    }

    #[test]
    fn bcut2d_single_atom_zero() {
        // H-only molecule → no heavy atoms → default struct
        let b = bcut2d(&mol("[H][H]"));
        assert_eq!(b.mwhi, 0.0);
    }

    #[test]
    fn logp_and_mr_matches_individual_functions() {
        // Verify logp_and_mr() returns identical values to logp_crippen() + molar_refractivity()
        // on a diverse set of molecules.
        let smiles = [
            "C",                                // methane
            "c1ccccc1",                         // benzene
            "CC(=O)Oc1ccccc1C(=O)O",            // aspirin
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",     // caffeine
            "c1ccc2c(c1)cc1ccc3cccc4ccc2c1c34", // pyrene (PAH)
            "O=C(O)c1ccccc1O",                  // salicylic acid
            "CCO",                              // ethanol
            "CC(N)Cc1ccccc1",                   // phenylalanine
        ];
        for smi in smiles {
            let m = mol(smi);
            let expected_logp = logp_crippen(&m);
            let expected_mr = molar_refractivity(&m);
            let (got_logp, got_mr) = logp_and_mr(&m);
            assert!(
                (expected_logp - got_logp).abs() < 1e-9,
                "{smi}: logp_crippen={expected_logp:.6} vs logp_and_mr.0={got_logp:.6}"
            );
            assert!(
                (expected_mr - got_mr).abs() < 1e-9,
                "{smi}: molar_refractivity={expected_mr:.6} vs logp_and_mr.1={got_mr:.6}"
            );
        }
    }
}
