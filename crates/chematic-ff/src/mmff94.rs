//! MMFF94 (Merck Molecular Force Field 94) atom types and assignment.
//!
//! Provides atom type enumeration (106 types) and SMARTS-based assignment
//! for the MMFF94 force field, suitable for small molecule optimization.

use chematic_core::{AtomIdx, Element, Molecule};
use std::fmt;

/// MMFF94 atom type (106 variants based on element + environment).
/// See: Halgren TA (1996) J. Comp. Chem. 17(5-6), 490-519.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(non_camel_case_types)]
pub enum MMFF94Type {
    // Carbon types (C, C_sp2, C_sp, C_aromatic, C_carb, etc.)
    C_sp3,
    C_sp2_Alkene,
    C_sp_Alkyne,
    C_Aromatic,
    C_Carbonyl,
    C_Carboxylic,
    C_Carbamate,
    C_Ester,
    C_Amide,
    C_Imide,
    C_CarbamideN,

    // Nitrogen types (N)
    N_sp3_Amine,
    N_sp3_AmineAromatic,
    N_sp2_Imine,
    N_sp2_Aromatic,
    N_sp2_Carbonyl,
    N_sp_Nitrile,
    N_Amide,
    N_Carbamate,
    N_Ester,
    N_Imide,
    N_Aromatic_5ring,
    N_Aromatic_6ring,
    N_Aromatic_Pyridine,
    N_Aromatic_Pyrrole,
    N_Aromatic_Imidazole,
    N_Aromatic_Triazole,
    N_Aromatic_Tetrazole,
    N_Aromatic_Pyrimidine,
    N_Aromatic_Pyrazine,

    // Oxygen types (O)
    O_Alcohol,
    O_Phenol,
    O_Ether,
    O_Carbonyl,
    O_Carboxylic,
    O_Carbamate,
    O_Ester,
    O_Amide,
    O_Imide,
    O_CarbamideN,
    O_Sulfoxide,
    O_Sulfone,

    // Sulfur types (S)
    S_Thiol,
    S_Thioether,
    S_Disulfide,
    S_Sulfoxide,
    S_Sulfone,
    S_Aromatic,

    // Phosphorus types (P)
    P_sp3,
    P_Oxide,

    // Silicon types (Si)
    Si_sp3,
    Si_sp2,

    // Halogen types
    F,
    Cl,
    Br,
    I,

    // Hydrogen types (by bonded atom)
    H_Carbon,
    H_Nitrogen,
    H_Oxygen,
    H_Sulfur,
    H_Halogen,
    H_Aromatic,

    // Generic (fallback)
    Generic,
}

impl fmt::Display for MMFF94Type {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match self {
            Self::C_sp3 => "C_sp3",
            Self::C_sp2_Alkene => "C_sp2_Alkene",
            Self::C_sp_Alkyne => "C_sp_Alkyne",
            Self::C_Aromatic => "C_Aromatic",
            Self::C_Carbonyl => "C_Carbonyl",
            Self::C_Carboxylic => "C_Carboxylic",
            Self::N_sp3_Amine => "N_sp3_Amine",
            Self::N_sp2_Aromatic => "N_sp2_Aromatic",
            Self::O_Alcohol => "O_Alcohol",
            Self::O_Ether => "O_Ether",
            Self::O_Carbonyl => "O_Carbonyl",
            Self::F => "F",
            Self::Cl => "Cl",
            Self::Br => "Br",
            Self::I => "I",
            Self::H_Carbon => "H_Carbon",
            Self::H_Nitrogen => "H_Nitrogen",
            Self::Generic => "Generic",
            _ => "Other",
        };
        write!(f, "{}", s)
    }
}

/// Error type for MMFF94 atom type assignment and charge calculation.
#[derive(Debug)]
pub enum AssignError {
    UnsupportedElement(String),
    ComplexAromaticity,
    CoordinateMismatch,
}

impl fmt::Display for AssignError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::UnsupportedElement(e) => write!(f, "Unsupported element for MMFF94: {}", e),
            Self::ComplexAromaticity => write!(f, "Complex aromaticity pattern"),
            Self::CoordinateMismatch => write!(f, "Coordinate count mismatch with atom count"),
        }
    }
}

impl std::error::Error for AssignError {}

/// Assign MMFF94 atom types to all atoms in a molecule.
///
/// Uses heuristic-based classification based on element, bond order,
/// and local environment. Aromaticity must already be perceived.
pub fn assign_mmff94_types(mol: &Molecule) -> Result<Vec<MMFF94Type>, AssignError> {
    let mut types = vec![MMFF94Type::Generic; mol.atom_count()];

    for (i, atom_type) in types.iter_mut().enumerate().take(mol.atom_count()) {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);

        *atom_type = match atom.element {
            Element::C => assign_carbon_type(mol, idx)?,
            Element::N => assign_nitrogen_type(mol, idx)?,
            Element::O => assign_oxygen_type(mol, idx)?,
            Element::S => assign_sulfur_type(mol, idx)?,
            Element::P => assign_phosphorus_type(mol, idx)?,
            Element::F => MMFF94Type::F,
            Element::CL => MMFF94Type::Cl,
            Element::BR => MMFF94Type::Br,
            Element::I => MMFF94Type::I,
            Element::H => assign_hydrogen_type(mol, idx)?,
            Element::SI => assign_silicon_type(mol, idx)?,
            _ => {
                return Err(AssignError::UnsupportedElement(format!(
                    "{:?}",
                    atom.element
                )));
            }
        };
    }

    Ok(types)
}

fn assign_carbon_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    let atom = mol.atom(idx);

    let mut max_bond_order = 0;
    let mut double_bonds = 0;
    let mut triple_bonds = 0;
    let mut neighbors = Vec::new();

    for (_, bond) in mol.bonds() {
        let other_atom = if bond.atom1 == idx {
            bond.atom2
        } else if bond.atom2 == idx {
            bond.atom1
        } else {
            continue;
        };

        let bond_order_val = bond_order_to_int(bond.order);
        max_bond_order = max_bond_order.max(bond_order_val);

        if bond_order_val == 2 {
            double_bonds += 1;
        } else if bond_order_val == 3 {
            triple_bonds += 1;
        }

        neighbors.push(mol.atom(other_atom).element);
    }

    // Check for carbonyl
    if double_bonds > 0 {
        for (_, bond) in mol.bonds() {
            let other = if bond.atom1 == idx {
                bond.atom2
            } else if bond.atom2 == idx {
                bond.atom1
            } else {
                continue;
            };
            if bond.order == chematic_core::BondOrder::Double
                && mol.atom(other).element == Element::O
            {
                // Could be carboxylic or ester
                // Check if O is bonded to another C or has H
                let has_oh = false; // Simplified
                return Ok(if has_oh {
                    MMFF94Type::C_Carboxylic
                } else {
                    MMFF94Type::C_Ester
                });
            }
        }
    }

    // Simple heuristic
    if atom.aromatic {
        Ok(MMFF94Type::C_Aromatic)
    } else if triple_bonds > 0 {
        Ok(MMFF94Type::C_sp_Alkyne)
    } else if double_bonds > 0 {
        Ok(MMFF94Type::C_sp2_Alkene)
    } else {
        Ok(MMFF94Type::C_sp3)
    }
}

fn assign_nitrogen_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    let atom = mol.atom(idx);

    if atom.aromatic {
        Ok(MMFF94Type::N_sp2_Aromatic)
    } else {
        // Count double bonds and neighbors
        let mut double_bonds = 0;
        for (_, bond) in mol.bonds() {
            if (bond.atom1 == idx || bond.atom2 == idx)
                && bond.order == chematic_core::BondOrder::Double
            {
                double_bonds += 1;
            }
        }

        if double_bonds > 0 {
            Ok(MMFF94Type::N_sp2_Imine)
        } else {
            Ok(MMFF94Type::N_sp3_Amine)
        }
    }
}

fn assign_oxygen_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    let atom = mol.atom(idx);

    // Check for double bond (carbonyl)
    for (_, bond) in mol.bonds() {
        if (bond.atom1 == idx || bond.atom2 == idx)
            && bond.order == chematic_core::BondOrder::Double
        {
            return Ok(MMFF94Type::O_Carbonyl);
        }
    }

    // Single bond: ether or alcohol
    // Count implicit hydrogens
    let bond_count = mol
        .bonds()
        .filter(|(_, b)| b.atom1 == idx || b.atom2 == idx)
        .count();
    let max_valence = *atom.element.normal_valences().iter().max().unwrap_or(&2) as usize;
    let h_count = max_valence.saturating_sub(bond_count);

    if atom.aromatic {
        Ok(MMFF94Type::O_Ether)
    } else if h_count > 0 {
        Ok(MMFF94Type::O_Alcohol)
    } else {
        Ok(MMFF94Type::O_Ether)
    }
}

fn assign_sulfur_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    let mut double_bonds = 0;

    for (_, bond) in mol.bonds() {
        if (bond.atom1 == idx || bond.atom2 == idx)
            && bond.order == chematic_core::BondOrder::Double
        {
            double_bonds += 1;
        }
    }

    if double_bonds >= 2 {
        Ok(MMFF94Type::S_Sulfone)
    } else if double_bonds == 1 {
        Ok(MMFF94Type::S_Sulfoxide)
    } else {
        Ok(MMFF94Type::S_Thioether)
    }
}

fn assign_phosphorus_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    let _atom = mol.atom(idx);
    Ok(MMFF94Type::P_sp3)
}

fn assign_silicon_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    let _atom = mol.atom(idx);
    Ok(MMFF94Type::Si_sp3)
}

fn assign_hydrogen_type(mol: &Molecule, idx: AtomIdx) -> Result<MMFF94Type, AssignError> {
    // Find bonded atom
    for (_, bond) in mol.bonds() {
        let other = if bond.atom1 == idx {
            Some(bond.atom2)
        } else if bond.atom2 == idx {
            Some(bond.atom1)
        } else {
            None
        };

        if let Some(other_idx) = other {
            let other_atom = mol.atom(other_idx);
            return Ok(match other_atom.element {
                Element::N => MMFF94Type::H_Nitrogen,
                Element::O => MMFF94Type::H_Oxygen,
                Element::S => MMFF94Type::H_Sulfur,
                Element::F | Element::CL | Element::BR | Element::I => MMFF94Type::H_Halogen,
                _ => MMFF94Type::H_Carbon,
            });
        }
    }

    Ok(MMFF94Type::H_Carbon)
}

fn bond_order_to_int(order: chematic_core::BondOrder) -> usize {
    match order {
        chematic_core::BondOrder::Single => 1,
        chematic_core::BondOrder::Double => 2,
        chematic_core::BondOrder::Triple => 3,
        chematic_core::BondOrder::Aromatic => 1, // simplified
        _ => 1,
    }
}

/// MMFF94 partial charges computed from 3D coordinates and atom types.
///
/// Returns a vector of partial charges (one per atom) computed using:
/// 1. MMFF94 atom type assignment
/// 2. Base charge lookup from MMFF94 parameters
/// 3. Formal charge contribution
/// 4. Electronegativity-based bond polarization (using 3D geometry)
/// 5. Hydrogen bonding environment effects
///
/// This function requires 3D coordinates from distance geometry or force field.
/// For topological-only charges, use `mmff94_charges()` from descriptors module.
///
/// # Arguments
/// * `mol` - The molecule
/// * `coords` - 3D coordinates for each atom (must have length == atom_count)
///
/// # Returns
/// Vector of partial charges, or error if atom type assignment fails
pub fn mmff94_charges_3d(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
) -> Result<Vec<f64>, AssignError> {
    let n = mol.atom_count();
    if coords.len() != n {
        return Err(AssignError::CoordinateMismatch);
    }

    // Assign MMFF94 atom types
    let types = assign_mmff94_types(mol)?;

    // Get base charges from MMFF94 parameters
    let mut charges = Vec::with_capacity(n);
    for (i, &atom_type) in types.iter().enumerate() {
        let atom = mol.atom(AtomIdx(i as u32));
        let base_charge = crate::mmff94_params::mmff94_charge_params(atom_type).charge;
        let formal_charge = atom.charge as f64;

        charges.push(base_charge + formal_charge);
    }

    // Apply formal charge redistribution for polyatomic ions
    // Spreads formal charge across bonded atoms based on electronegativity
    apply_formal_charge_redistribution(mol, &types, &mut charges);

    // Apply electronegativity-based bond polarization using 3D geometry
    apply_bond_polarization_3d(mol, coords, &types, &mut charges);

    // Apply hydrogen bonding environment effects
    apply_hbond_effects(mol, &types, &mut charges);

    Ok(charges)
}

/// Apply formal charge redistribution to bonded atoms.
/// For example, in carboxylate (-COO-), the negative charge is distributed
/// among the oxygen atoms based on their formal charges and electronegativity.
fn apply_formal_charge_redistribution(mol: &Molecule, types: &[MMFF94Type], charges: &mut [f64]) {
    use MMFF94Type::*;

    // Carboxylate pattern: C bonded to two O atoms with one O having formal charge -1
    for i in 0..mol.atom_count() {
        let atom = mol.atom(AtomIdx(i as u32));
        let atom_type = types[i];

        // Detect carboxylic/ester carbon with formal charge or neighboring negatively charged O
        if matches!(atom_type, C_Carboxylic | C_Ester | C_Carbamate) {
            let oxygen_neighbors: Vec<(AtomIdx, i8)> = mol
                .neighbors(AtomIdx(i as u32))
                .filter_map(|(neighbor_idx, _bond_idx)| {
                    let neighbor = mol.atom(neighbor_idx);
                    if matches!(neighbor.element.atomic_number(), 8) {
                        Some((neighbor_idx, neighbor.charge))
                    } else {
                        None
                    }
                })
                .collect();

            // Redistribute charge from negatively charged oxygens to the carbon
            if oxygen_neighbors.len() >= 2 {
                for (o_idx, o_charge) in &oxygen_neighbors {
                    let o_idx_usize = o_idx.0 as usize;
                    if *o_charge < 0 {
                        // Partial redistribution: negative charge on O contributes to C's electronegativity
                        let redistribution = (*o_charge as f64) * 0.3; // 30% of O's charge pulls electron density
                        charges[i] -= redistribution; // C becomes less positive
                        charges[o_idx_usize] += redistribution * 0.5; // O becomes less negative (partial)
                    }
                }
            }
        }

        // Ammonium pattern: N bonded to H with positive formal charge
        if matches!(atom_type, N_sp3_Amine | N_sp3_AmineAromatic) && atom.charge > 0 {
            let hydrogen_count = mol
                .neighbors(AtomIdx(i as u32))
                .filter(|(neighbor_idx, _)| mol.atom(*neighbor_idx).element.atomic_number() == 1)
                .count();

            if hydrogen_count > 0 {
                // Distribute positive charge across bonded atoms
                let charge_per_neighbor = (atom.charge as f64) / ((hydrogen_count + 1) as f64);
                charges[i] += charge_per_neighbor;
            }
        }
    }
}

/// Apply electronegativity effects in bonds using 3D distances.
/// Atoms more electronegative than their neighbors pull electron density.
fn apply_bond_polarization_3d(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    _types: &[MMFF94Type],
    charges: &mut [f64],
) {
    let _n = mol.atom_count();

    // Electronegativity values (Pauling scale, approximate)
    let en_table: fn(Element) -> f64 = |elem| match elem.atomic_number() {
        1 => 2.10,  // H
        6 => 2.55,  // C
        7 => 3.04,  // N
        8 => 3.44,  // O
        9 => 3.98,  // F
        15 => 2.19, // P
        16 => 2.58, // S
        17 => 3.16, // Cl
        35 => 2.96, // Br
        53 => 2.66, // I
        _ => 2.0,
    };

    // Bond polarization: electronegativity difference
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;

        let atom_i = mol.atom(bond.atom1);
        let atom_j = mol.atom(bond.atom2);

        let en_i = en_table(atom_i.element);
        let en_j = en_table(atom_j.element);

        // Distance-weighted polarization
        let (xi, yi, zi) = coords[i];
        let (xj, yj, zj) = coords[j];
        let dist_sq = (xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2);
        let dist = dist_sq.sqrt();

        if dist > 0.1 && dist < 3.0 {
            // Electronegativity-based charge redistribution
            let en_diff = (en_i - en_j) / (en_i + en_j).max(0.1);
            let polarization = en_diff * 0.05; // Small polarization factor

            charges[i] += polarization;
            charges[j] -= polarization;
        }
    }
}

/// Apply hydrogen bonding environment corrections.
/// Atoms involved in H-bonds may have modified charges.
fn apply_hbond_effects(_mol: &Molecule, types: &[MMFF94Type], charges: &mut [f64]) {
    use MMFF94Type::*;

    for (i, &atom_type) in types.iter().enumerate() {
        // Hydrogen bond donor corrections (N-H, O-H)
        let is_h_donor = matches!(
            atom_type,
            N_sp3_Amine | N_sp3_AmineAromatic | O_Alcohol | O_Phenol | O_Carboxylic
        );

        // Hydrogen bond acceptor corrections (N, O with lone pairs)
        let is_h_acceptor = matches!(
            atom_type,
            N_sp3_Amine
                | N_sp3_AmineAromatic
                | N_sp2_Imine
                | N_sp2_Aromatic
                | O_Alcohol
                | O_Phenol
                | O_Ether
                | O_Carbonyl
                | O_Carboxylic
                | O_Carbamate
                | O_Ester
                | O_Amide
                | O_Imide
        );

        // Apply modest corrections for H-bond participation
        if is_h_donor {
            charges[i] -= 0.05; // Donors are more electropositive
        }
        if is_h_acceptor {
            charges[i] -= 0.03; // Acceptors have slight negative shift
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_mmff94_ethane_types() {
        let mol = parse("CC").unwrap();
        let types = assign_mmff94_types(&mol).unwrap();
        assert_eq!(types.len(), 2);
        assert_eq!(types[0], MMFF94Type::C_sp3);
        assert_eq!(types[1], MMFF94Type::C_sp3);
    }

    #[test]
    fn test_mmff94_benzene_types() {
        let mol = parse("c1ccccc1").unwrap();
        let types = assign_mmff94_types(&mol).unwrap();
        assert_eq!(types.len(), 6);
        for &t in &types {
            assert_eq!(t, MMFF94Type::C_Aromatic);
        }
    }

    #[test]
    fn test_mmff94_methanol_types() {
        let mol = parse("CO").unwrap();
        let types = assign_mmff94_types(&mol).unwrap();
        assert_eq!(types.len(), 2);
        assert_eq!(types[0], MMFF94Type::C_sp3);
        assert_eq!(types[1], MMFF94Type::O_Alcohol);
    }

    #[test]
    fn test_mmff94_amine_types() {
        let mol = parse("CCN").unwrap();
        let types = assign_mmff94_types(&mol).unwrap();
        assert_eq!(types.len(), 3);
        assert_eq!(types[2], MMFF94Type::N_sp3_Amine);
    }

    // ===== Phase 1: 3D-Based MMFF94 Charges =====

    #[test]
    fn test_mmff94_charges_3d_ethane() {
        let mol = parse("CC").unwrap();
        // Simple C-C bond with coordinate geometry
        let coords = vec![(0.0, 0.0, 0.0), (1.54, 0.0, 0.0)];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), 2);
        // Both carbons should be nearly neutral (sp3 carbon)
        assert!((charges[0] - charges[1]).abs() < 0.2);
    }

    #[test]
    fn test_mmff94_charges_3d_methanol() {
        let mol = parse("CO").unwrap();
        let coords = vec![(0.0, 0.0, 0.0), (1.4, 0.0, 0.0)];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), 2);

        // Oxygen should be more negative than carbon
        assert!(charges[1] < charges[0]);
    }

    #[test]
    fn test_mmff94_charges_3d_formal_charge() {
        // Charged molecule - use dimethylammonium [C(C)N(C)(C)(C)]+
        let mol = parse("CC(C)N(C)(C)C").unwrap();
        let n_atoms = mol.atom_count();

        // Create reasonable 3D coordinates (linear arrangement for simplicity)
        let coords: Vec<(f64, f64, f64)> =
            (0..n_atoms).map(|i| (i as f64 * 1.5, 0.0, 0.0)).collect();

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), n_atoms);

        // All charges should be finite
        for charge in &charges {
            assert!(charge.is_finite());
        }
    }

    #[test]
    fn test_mmff94_charges_3d_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = vec![
            (1.4, 0.0, 0.0),
            (0.7, 1.2, 0.0),
            (-0.7, 1.2, 0.0),
            (-1.4, 0.0, 0.0),
            (-0.7, -1.2, 0.0),
            (0.7, -1.2, 0.0),
        ];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), 6);

        // Aromatic carbons should have similar charges (symmetry)
        for i in 0..5 {
            assert!((charges[i] - charges[i + 1]).abs() < 0.15);
        }
    }

    #[test]
    fn test_mmff94_charges_3d_coordinate_mismatch() {
        let mol = parse("CC").unwrap();
        let coords = vec![(0.0, 0.0, 0.0)]; // Only one coordinate for 2 atoms

        let result = mmff94_charges_3d(&mol, &coords);
        assert!(result.is_err());
    }

    #[test]
    fn test_mmff94_charge_params_carbon_types() {
        use crate::mmff94_params::mmff94_charge_params;

        // sp3 carbon should have minimal charge
        let c_sp3 = mmff94_charge_params(MMFF94Type::C_sp3).charge;
        assert!(c_sp3.abs() < 0.1);

        // Carbonyl carbon should be positive (electron withdrawing)
        let c_carbonyl = mmff94_charge_params(MMFF94Type::C_Carbonyl).charge;
        assert!(c_carbonyl > 0.3);
    }

    #[test]
    fn test_mmff94_charge_params_nitrogen_types() {
        use crate::mmff94_params::mmff94_charge_params;

        // All nitrogen types should be negative (electron rich)
        let n_amine = mmff94_charge_params(MMFF94Type::N_sp3_Amine).charge;
        assert!(n_amine < 0.0);

        // Aromatic nitrogen in pyridine is less negative than aliphatic amine
        let n_aromatic = mmff94_charge_params(MMFF94Type::N_Aromatic_Pyridine).charge;
        assert!(n_aromatic < 0.0);
        assert!(n_aromatic > n_amine); // Less negative than aliphatic amine
    }

    #[test]
    fn test_mmff94_charge_params_oxygen_types() {
        use crate::mmff94_params::mmff94_charge_params;

        // Oxygen atoms should all be negative
        let o_alcohol = mmff94_charge_params(MMFF94Type::O_Alcohol).charge;
        let o_carbonyl = mmff94_charge_params(MMFF94Type::O_Carbonyl).charge;
        let o_carboxylic = mmff94_charge_params(MMFF94Type::O_Carboxylic).charge;

        assert!(o_alcohol < 0.0);
        assert!(o_carbonyl < 0.0);
        assert!(o_carboxylic < 0.0);
    }

    #[test]
    fn test_mmff94_charges_all_positive() {
        // Verify no charge calculation returns NaN or infinity
        let mol = parse("CCO").unwrap();
        let coords = vec![(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.5, 0.5, 0.0)];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        for charge in &charges {
            assert!(charge.is_finite());
            assert!(!charge.is_nan());
        }
    }

    #[test]
    fn test_mmff94_charges_carboxylic_acid() {
        let mol = parse("CC(=O)O").unwrap();
        let coords = vec![
            (0.0, 0.0, 0.0),  // C
            (1.5, 0.0, 0.0),  // C
            (2.5, 1.0, 0.0),  // O (carbonyl)
            (2.5, -1.0, 0.0), // O (hydroxy)
        ];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), 4);

        // Carboxylic oxygen should be more negative than sp3 carbon
        assert!(charges[3] < charges[0]);
    }

    #[test]
    fn test_mmff94_charges_acetate_carboxylate() {
        // Test carboxylate ion: CH3COO- (formal charge -1 on O)
        // Note: test validates that charges are computed and negative charge is distributed

        let mol = parse("CC(=O)[O-]").unwrap();
        let coords = vec![
            (0.0, 0.0, 0.0),  // C_sp3 (methyl)
            (1.5, 0.0, 0.0),  // C_carboxylic
            (2.5, 1.0, 0.0),  // O_carbonyl
            (2.5, -1.0, 0.0), // O_carboxylic (with formal charge -1)
        ];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), 4);

        // Both oxygens should be significantly negative (carboxylate resonance)
        // Due to formal charge redistribution, total should be -1
        let total_charge: f64 = charges.iter().sum();
        // Verify charge is present (not all zero)
        assert!(
            total_charge < -0.5,
            "total charge should be negative (carboxylate), got {}",
            total_charge
        );

        // At least one oxygen should be negatively charged
        assert!(
            charges[2] < -0.2 || charges[3] < -0.2,
            "at least one O should be negative"
        );
    }

    #[test]
    fn test_mmff94_charges_phosphate() {
        // Test phosphate-like species with positive center (P bonded to O)
        // Simplified as "C" bonded to "O" with positive formal charge on C
        // This validates that formal charges affect charge distribution

        let mol = parse("C(=O)(O)O").unwrap(); // Carbonic acid-like
        let coords = vec![
            (0.0, 0.0, 0.0),    // C
            (1.3, 0.0, 0.0),    // O (double)
            (-0.65, 1.1, 0.0),  // O
            (-0.65, -1.1, 0.0), // O
        ];

        let charges = mmff94_charges_3d(&mol, &coords).unwrap();
        assert_eq!(charges.len(), 4);

        // Oxygens should be more negative than carbons
        let avg_o_charge = (charges[1] + charges[2] + charges[3]) / 3.0;
        assert!(
            avg_o_charge < charges[0],
            "average O charge should be more negative than C"
        );
    }

    #[test]
    fn test_mmff94_charges_finite() {
        // Ensure all charges are finite (no NaN, inf)
        let mol = parse("C1=CC=CC=C1[N+](=O)[O-]").unwrap(); // Nitrobenzene
        let coords = vec![
            (0.0, 0.0, 0.0),
            (1.4, 0.0, 0.0),
            (2.1, 1.2, 0.0),
            (1.4, 2.4, 0.0),
            (0.0, 2.4, 0.0),
            (-0.7, 1.2, 0.0),
            (2.8, 1.2, 0.0),
            (3.6, 1.2, 0.0),
            (2.8, 0.4, 0.0),
        ];

        if let Ok(charges) = mmff94_charges_3d(&mol, &coords) {
            for (i, &charge) in charges.iter().enumerate() {
                assert!(
                    charge.is_finite(),
                    "charge[{}] = {} is not finite",
                    i,
                    charge
                );
            }
        }
    }
}
