//! Molecular descriptor functions for drug-likeness and physical property estimation.
//!
//! All functions accept a `&Molecule` reference.  Molecules with aromatic bonds
//! (SMILES lowercase notation) are kekulized internally where hydrogen counts
//! are required; the caller's molecule is never mutated.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, BondIdx, Element, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

/// Average atomic mass table.
/// Falls back to `atomic_number as f64` for unlisted elements.
fn avg_mass(element: Element) -> f64 {
    match element.atomic_number() {
        1  => 1.008,   // H
        2  => 4.003,   // He
        3  => 6.941,   // Li
        4  => 9.012,   // Be
        5  => 10.811,  // B
        6  => 12.011,  // C
        7  => 14.007,  // N
        8  => 15.999,  // O
        9  => 18.998,  // F
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
        n  => n as f64,
    }
}

/// Monoisotopic (most-abundant isotope) mass table.
/// Falls back to `atomic_number as f64` for unlisted elements.
fn mono_mass(element: Element) -> f64 {
    match element.atomic_number() {
        1  => 1.00783,   // H  (1H)
        6  => 12.0000,   // C  (12C)
        7  => 14.0031,   // N  (14N)
        8  => 15.9949,   // O  (16O)
        9  => 18.9984,   // F  (19F)
        14 => 27.9769,   // Si (28Si)
        15 => 30.9738,   // P  (31P)
        16 => 31.9721,   // S  (32S)
        17 => 34.9689,   // Cl (35Cl)
        35 => 78.9183,   // Br (79Br)
        34 => 79.9165,   // Se (80Se)
        53 => 126.9045,  // I  (127I)
        n  => n as f64,
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
            (an == 7 || an == 8) && implicit_hcount(mol, *idx) > 0
        })
        .count()
}

// ---------------------------------------------------------------------------
// 5. Hydrogen bond acceptor count (Lipinski style)
// ---------------------------------------------------------------------------

/// Count hydrogen bond acceptors (Lipinski definition).
///
/// Simple N + O count: every nitrogen or oxygen atom (regardless of H count).
pub fn hba_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, atom)| {
            let an = atom.element.atomic_number();
            an == 7 || an == 8
        })
        .count()
}

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
    let rings = find_sssr(mol);

    // Build the set of bond indices that belong to at least one ring.
    let mut ring_bond_set: HashSet<BondIdx> = HashSet::new();
    for ring in rings.rings() {
        for i in 0..ring.len() {
            let a = ring[i];
            let b = ring[(i + 1) % ring.len()];
            if let Some((bidx, _)) = mol.bond_between(a, b) {
                ring_bond_set.insert(bidx);
            }
        }
    }

    let mut count = 0usize;
    for (bidx, bond) in mol.bonds() {
        // Must be a single bond (stereo bonds Up/Down are also single).
        let is_single = matches!(bond.order, BondOrder::Single | BondOrder::Up | BondOrder::Down);
        if !is_single {
            continue;
        }

        // Not in a ring.
        if ring_bond_set.contains(&bidx) {
            continue;
        }

        let a1 = bond.atom1;
        let a2 = bond.atom2;

        // Both endpoints must be non-terminal.
        if mol.degree(a1) <= 1 || mol.degree(a2) <= 1 {
            continue;
        }

        // Exclude amide bonds: C-N bond where the C has a double bond to O.
        if is_amide_bond(mol, a1, a2) {
            continue;
        }

        count += 1;
    }

    count
}

/// Return true if the bond between `a` and `b` is an amide-like C-N bond.
///
/// Condition: one atom is N, the other is C, and that C has at least one
/// double bond to an oxygen neighbor.
fn is_amide_bond(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let atom_a = mol.atom(a);
    let atom_b = mol.atom(b);

    let c_idx = if atom_a.element.atomic_number() == 6
        && atom_b.element.atomic_number() == 7
    {
        a
    } else if atom_b.element.atomic_number() == 6
        && atom_a.element.atomic_number() == 7
    {
        b
    } else {
        return false;
    };

    // Check whether the C has any double bond to an O neighbor.
    mol.neighbors(c_idx).any(|(nb, nbidx)| {
        let bond = mol.bond(nbidx);
        mol.atom(nb).element.atomic_number() == 8
            && bond.order == BondOrder::Double
    })
}

// ---------------------------------------------------------------------------
// 7. TPSA
// ---------------------------------------------------------------------------

/// Compute the topological polar surface area (Å²) using the Ertl (2000) table.
///
/// Reference: P. Ertl, B. Rohde, P. Selzer, J. Med. Chem. 2000, 43, 3714-3717.
pub fn tpsa(mol: &Molecule) -> f64 {
    let mut psa = 0.0f64;

    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let is_aromatic = atom.aromatic;
        let h = implicit_hcount(mol, idx);

        let contribution = match an {
            // Nitrogen
            7 => {
                if is_aromatic {
                    if h > 0 {
                        13.97 // aromatic NH (pyrrole-type)
                    } else {
                        12.89 // aromatic N= (pyridine-type)
                    }
                } else {
                    // aliphatic N
                    if h >= 2 {
                        26.02 // NH2
                    } else if h == 1 {
                        12.03 // NH
                    } else {
                        3.24  // tertiary N
                    }
                }
            }
            // Oxygen
            8 => {
                if is_aromatic {
                    13.14
                } else if h > 0 {
                    20.23 // OH
                } else {
                    // Distinguish carbonyl O (C=O, 17.07 Å²) from ether O (C-O-C, 9.23 Å²).
                    // A carbonyl O has a double bond to its neighbor; an ether O has only single bonds.
                    let is_carbonyl = mol.neighbors(idx).any(|(_, bidx)| {
                        mol.bond(bidx).order == BondOrder::Double
                    });
                    if is_carbonyl { 17.07 } else { 9.23 }
                }
            }
            // Sulfur
            16 => {
                if is_aromatic {
                    0.0
                } else if h > 0 {
                    38.80 // SH
                } else {
                    25.30 // S
                }
            }
            // Phosphorus
            15 => {
                if !is_aromatic {
                    34.14
                } else {
                    0.0
                }
            }
            _ => 0.0,
        };

        psa += contribution;
    }

    psa
}

// ---------------------------------------------------------------------------
// 8. LogP (Wildman-Crippen, simplified)
// ---------------------------------------------------------------------------

/// Compute an approximate Wildman-Crippen LogP.
///
/// Uses a simplified atom-type table.  Each atom is classified by element,
/// aromaticity, and implicit hydrogen count; unlisted types contribute 0.0.
pub fn logp_crippen(mol: &Molecule) -> f64 {
    let mut logp = 0.0f64;

    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let is_aromatic = atom.aromatic;
        let h = implicit_hcount(mol, idx);

        // Determine whether this C is sp2 (has any double or triple bond, or aromatic).
        let is_sp2_c = an == 6
            && !is_aromatic
            && mol.neighbors(idx).any(|(_, bidx)| {
                matches!(
                    mol.bond(bidx).order,
                    BondOrder::Double | BondOrder::Triple
                )
            });

        let contrib = match an {
            // Carbon
            6 => {
                if is_aromatic {
                    0.1441
                } else if is_sp2_c {
                    // sp2 C: match by H count
                    match h {
                        0 => -0.2150, // carbonyl C, no H
                        1 => -0.1477, // =CH-
                        _ => -0.1477, // fallback sp2
                    }
                } else {
                    // sp3 C
                    match h {
                        0 => -0.2035, // quaternary
                        1 => -0.2051, // CH
                        2 => -0.1321, // CH2
                        _ => -0.0880, // CH3 (3 or more H)
                    }
                }
            }
            // Nitrogen
            7 => {
                if is_aromatic {
                    0.2626
                } else {
                    match h {
                        0 => -0.4806, // tertiary N
                        1 => -0.5188, // secondary N
                        _ => -0.7323, // primary NH2
                    }
                }
            }
            // Oxygen
            8 => {
                if is_aromatic {
                    -0.1188
                } else {
                    0.1552 // OH or ether O (simplified)
                }
            }
            // Sulfur
            16 => {
                if is_aromatic {
                    0.0000
                } else {
                    0.2432
                }
            }
            // Halogens
            9  => 0.4202, // F
            17 => 0.6895, // Cl
            35 => 0.8456, // Br
            53 => 0.8857, // I
            // Phosphorus
            15 => 0.0000,
            _  => 0.0,
        };

        logp += contrib;
    }

    logp
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
        assert!(pct2(molecular_weight(&m), 16.043), "methane MW = {}", molecular_weight(&m));
    }

    // -- Test 2: water molecular weight -------------------------------------
    #[test]
    fn test_mw_water() {
        let m = mol("O");
        // H2O: 15.999 + 2*1.008 = 18.015
        assert!(pct2(molecular_weight(&m), 18.015), "water MW = {}", molecular_weight(&m));
    }

    // -- Test 3: ethanol molecular weight -----------------------------------
    #[test]
    fn test_mw_ethanol() {
        let m = mol("CCO");
        // C2H6O: 2*12.011 + 6*1.008 + 15.999 = 46.068
        assert!(pct2(molecular_weight(&m), 46.068), "ethanol MW = {}", molecular_weight(&m));
    }

    // -- Test 4: benzene molecular weight -----------------------------------
    #[test]
    fn test_mw_benzene() {
        let m = mol("c1ccccc1");
        // C6H6: 6*12.011 + 6*1.008 = 78.114
        assert!(pct2(molecular_weight(&m), 78.114), "benzene MW = {}", molecular_weight(&m));
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
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert_eq!(hba_count(&m), 4); // four O
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
}
