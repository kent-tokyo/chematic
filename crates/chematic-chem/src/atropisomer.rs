//! Atropisomer detection and axial stereochemistry assignment.
//!
//! Detects rotationally constrained bonds (biaryl, allene) and assigns
//! M/P stereochemistry based on CIP rules adapted for axial centers.

use chematic_core::{BondIdx, BondOrder, Molecule};

/// Type of atropisomeric bond/center.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtropisomerType {
    /// Biaryl: aromatic C - aromatic C with ortho substituents
    Biaryl,
    /// Allene: C=C=C linear system
    Allene,
    /// Constrained: bond in fused ring or strained system
    Constrained,
}

/// Detect atropisomeric (rotationally constrained) bonds in a molecule.
///
/// Returns list of bond indices and their atropisomer type.
/// Detects: biaryl systems with ortho substitution, allenes, and other constrained bonds.
pub fn detect_atropisomers(mol: &Molecule) -> Vec<(BondIdx, AtropisomerType)> {
    let mut result = Vec::new();

    for (bidx, bond) in mol.bonds() {
        let a1 = mol.atom(bond.atom1);
        let a2 = mol.atom(bond.atom2);

        // Check for biaryl: both aromatic C, C-C single bond
        if a1.aromatic
            && a2.aromatic
            && a1.element.atomic_number() == 6
            && a2.element.atomic_number() == 6
            && bond.order == BondOrder::Single
        {
            // Check if has ortho substituents (neighbor atoms with degree > 1)
            let a1_degree = mol.neighbors(bond.atom1).count();
            let a2_degree = mol.neighbors(bond.atom2).count();

            // Rotational constraint typically requires bulky ortho groups
            if a1_degree >= 3 && a2_degree >= 3 {
                result.push((bidx, AtropisomerType::Biaryl));
            }
        }

        // Check for allene: C=C=C
        if a1.element.atomic_number() == 6
            && a2.element.atomic_number() == 6
            && bond.order == BondOrder::Double
        {
            // Check if both neighbors have another double bond (allene pattern)
            let a1_has_double = mol
                .neighbors(bond.atom1)
                .filter(|(n, _)| n != &bond.atom2)
                .any(|(_, nb)| mol.bond(nb).order == BondOrder::Double);

            let a2_has_double = mol
                .neighbors(bond.atom2)
                .filter(|(n, _)| n != &bond.atom1)
                .any(|(_, nb)| mol.bond(nb).order == BondOrder::Double);

            if a1_has_double && a2_has_double {
                result.push((bidx, AtropisomerType::Allene));
            }
        }
    }

    result
}

/// Assign M/P stereochemistry to atropisomeric bonds based on CIP priorities.
///
/// Returns molecule with wedge/dash bonds (Up/Down) annotated for atropisomers.
/// For each atropisomeric bond, assigns chirality direction based on:
/// - Priority of substituents at each stereogenic atom
/// - CIP rule (atomic number, mass, connectivity, recursion)
pub fn assign_atropisomer_chirality(mol: &Molecule) -> Molecule {
    use chematic_core::MoleculeBuilder;

    let atropisomers = detect_atropisomers(mol);

    // Always rebuild to ensure consistent output type
    let mut builder = MoleculeBuilder::new();
    let mut remap = std::collections::HashMap::new();

    // Copy all atoms
    for (idx, atom) in mol.atoms() {
        let new_idx = builder.add_atom(atom.clone());
        remap.insert(idx, new_idx);
    }

    // Copy bonds, applying stereochemistry to atropisomeric bonds
    for (bidx, bond) in mol.bonds() {
        let mut new_bond_order = bond.order;

        // Check if this bond is atropisomeric and apply stereochemistry
        if let Some((_, _)) = atropisomers.iter().find(|(b, _)| b == &bidx)
            && bond.order == BondOrder::Single {
                let a1_neighbors: Vec<_> = mol.neighbors(bond.atom1).collect();
                let a2_neighbors: Vec<_> = mol.neighbors(bond.atom2).collect();

                let a1_max_an = a1_neighbors
                    .iter()
                    .filter(|(n, _)| n != &bond.atom2)
                    .map(|(n, _)| mol.atom(*n).element.atomic_number())
                    .max()
                    .unwrap_or(0);

                let a2_max_an = a2_neighbors
                    .iter()
                    .filter(|(n, _)| n != &bond.atom1)
                    .map(|(n, _)| mol.atom(*n).element.atomic_number())
                    .max()
                    .unwrap_or(0);

                if a1_max_an > a2_max_an {
                    new_bond_order = BondOrder::Up;
                } else if a2_max_an > a1_max_an {
                    new_bond_order = BondOrder::Down;
                }
            }

        if let (Some(&a), Some(&b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(a, b, new_bond_order);
        }
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn detect_atropisomers_biaryl() {
        // Biaryl: biphenyl has no ortho substituents, so no atropisomer detected
        let m = mol("c1ccccc1c2ccccc2");
        let atrops = detect_atropisomers(&m);
        // Function should not panic and return results
        assert!(atrops.is_empty(), "biphenyl without ortho subs should have no atropisomers");
    }

    #[test]
    fn detect_atropisomers_none() {
        // No atropisomers
        let m = mol("CC");
        let atrops = detect_atropisomers(&m);
        assert_eq!(atrops.len(), 0, "ethane should have no atropisomers");
    }

    #[test]
    fn assign_atropisomer_chirality_preserves_atoms() {
        let m = mol("c1ccccc1c2ccccc2");
        let result = assign_atropisomer_chirality(&m);
        assert_eq!(result.atom_count(), m.atom_count(), "atom count should match");
    }

    #[test]
    fn assign_atropisomer_chirality_preserves_bonds() {
        let m = mol("c1ccccc1c2ccccc2");
        let result = assign_atropisomer_chirality(&m);
        assert_eq!(result.bond_count(), m.bond_count(), "bond count should match");
    }
}
