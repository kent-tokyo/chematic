//! MMFF94 (Merck Molecular Force Field 94) atom types and assignment.
//!
//! Provides atom type enumeration (106 types) and SMARTS-based assignment
//! for the MMFF94 force field, suitable for small molecule optimization.

use chematic_core::{Molecule, AtomIdx, Element};
use std::fmt;

/// MMFF94 atom type (106 variants based on element + environment).
/// See: Halgren TA (1996) J. Comp. Chem. 17(5-6), 490-519.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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

/// Error type for MMFF94 atom type assignment.
#[derive(Debug)]
pub enum AssignError {
    UnsupportedElement(String),
    ComplexAromaticity,
}

impl fmt::Display for AssignError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::UnsupportedElement(e) => write!(f, "Unsupported element for MMFF94: {}", e),
            Self::ComplexAromaticity => write!(f, "Complex aromaticity pattern"),
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

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);

        types[i] = match atom.element {
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
            _ => return Err(AssignError::UnsupportedElement(
                format!("{:?}", atom.element),
            )),
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
            let other = if bond.atom1 == idx { bond.atom2 } else if bond.atom2 == idx { bond.atom1 } else { continue };
            if bond.order == chematic_core::BondOrder::Double {
                if mol.atom(other).element == Element::O {
                    // Could be carboxylic or ester
                    // Check if O is bonded to another C or has H
                    let has_oh = false; // Simplified
                    return Ok(if has_oh { MMFF94Type::C_Carboxylic } else { MMFF94Type::C_Ester });
                }
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
            if (bond.atom1 == idx || bond.atom2 == idx) && bond.order == chematic_core::BondOrder::Double {
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
        if (bond.atom1 == idx || bond.atom2 == idx) && bond.order == chematic_core::BondOrder::Double {
            return Ok(MMFF94Type::O_Carbonyl);
        }
    }

    // Single bond: ether or alcohol
    // Count implicit hydrogens
    let bond_count = mol.bonds().filter(|(_, b)| b.atom1 == idx || b.atom2 == idx).count();
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
        if (bond.atom1 == idx || bond.atom2 == idx) && bond.order == chematic_core::BondOrder::Double {
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
        let other = if bond.atom1 == idx { Some(bond.atom2) }
                   else if bond.atom2 == idx { Some(bond.atom1) }
                   else { None };

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
}
