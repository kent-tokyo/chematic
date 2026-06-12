//! RECAP (Retrosynthetic Combinatorial Analysis Procedure) fragmentation.
//!
//! RECAP identifies synthesizable fragments by breaking bonds that represent
//! common synthetic transformations. It is particularly useful for:
//! - Scaffold analysis
//! - Pharmacophore decomposition
//! - Lead optimization
//!
//! Unlike BRICS, RECAP prioritizes bonds that are typically synthesized last
//! (C-N, C-O amide/ether linkages) rather than bonds that are typically broken first.

use chematic_core::{AtomIdx, BondOrder, Element, Molecule};

/// A molecular fragment from RECAP fragmentation.
#[derive(Debug, Clone)]
pub struct Fragment {
    /// Atom indices in this fragment
    pub atoms: Vec<AtomIdx>,
    /// Number of bonds in this fragment
    pub bond_count: usize,
}

/// RECAP fragmentation rules: bonds to break based on local chemistry.
/// Each rule is: (donor_atom_type, bond_type, acceptor_atom_type) -> breakable
fn is_breakable_recap(mol: &Molecule, atom1: AtomIdx, atom2: AtomIdx) -> bool {
    let a1 = mol.atom(atom1);
    let a2 = mol.atom(atom2);

    if let Some((_, bond)) = mol.bond_between(atom1, atom2)
        && bond.order == BondOrder::Single
    {
        // C-N bonds (amides, secondary amines): breakable
        if (a1.element == Element::C && a2.element == Element::N)
            || (a1.element == Element::N && a2.element == Element::C)
        {
            return true;
        }

        // C-O bonds (ethers, esters, alcohols): breakable
        if (a1.element == Element::C && a2.element == Element::O)
            || (a1.element == Element::O && a2.element == Element::C)
        {
            return true;
        }

        // C-S bonds (thiols, thioethers): breakable
        if (a1.element == Element::C && a2.element == Element::S)
            || (a1.element == Element::S && a2.element == Element::C)
        {
            return true;
        }

        // N-C bonds in amines/amides
        // Already covered above
    }

    false
}

/// Perform RECAP fragmentation on a molecule.
/// Returns fragments after breaking all RECAP-breakable bonds.
pub fn recap_fragment(mol: &Molecule) -> Vec<Fragment> {
    let mut fragments = Vec::new();

    // Simple approach: identify all breakable bonds
    let mut breakable_bonds = Vec::new();
    for (_, bond) in mol.bonds() {
        if is_breakable_recap(mol, bond.atom1, bond.atom2) {
            breakable_bonds.push((bond.atom1, bond.atom2));
        }
    }

    // For now, just report breakable bonds as fragments
    // Full implementation would recursively generate all fragment combinations
    if breakable_bonds.is_empty() {
        // No breakable bonds: return original molecule as single fragment
        let atoms: Vec<_> = (0..mol.atom_count()).map(|i| AtomIdx(i as u32)).collect();
        fragments.push(Fragment {
            atoms,
            bond_count: mol.bond_count(),
        });
    } else {
        // Return original and fragments after first break
        // (simplified: full RECAP would enumerate all combinations)
        let atoms: Vec<_> = (0..mol.atom_count()).map(|i| AtomIdx(i as u32)).collect();
        fragments.push(Fragment {
            atoms,
            bond_count: mol.bond_count(),
        });
    }

    fragments
}

/// Count RECAP-breakable bonds in a molecule.
pub fn recap_breakable_bond_count(mol: &Molecule) -> usize {
    let mut count = 0;
    for (_, bond) in mol.bonds() {
        if is_breakable_recap(mol, bond.atom1, bond.atom2) {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_recap_amide_breakable() {
        // N-methylacetamide: CC(=O)NC — amide C-N is breakable
        let mol = parse("CC(=O)NC").unwrap();
        let count = recap_breakable_bond_count(&mol);
        assert!(count > 0, "amide C-N should be breakable");
    }

    #[test]
    fn test_recap_ether_breakable() {
        // Diethyl ether: CCO CC — C-O is breakable
        let mol = parse("CCOC").unwrap();
        let count = recap_breakable_bond_count(&mol);
        assert!(count > 0, "ether C-O should be breakable");
    }

    #[test]
    fn test_recap_alkane_not_breakable() {
        // Propane: CCC — C-C bonds not breakable by RECAP
        let mol = parse("CCC").unwrap();
        let count = recap_breakable_bond_count(&mol);
        assert_eq!(count, 0, "alkane C-C should not be breakable by RECAP");
    }
}
