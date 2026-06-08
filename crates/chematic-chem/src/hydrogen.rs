//! Explicit hydrogen management.
//!
//! Converts between the compact hydrogen-implicit representation used
//! throughout chematic and a fully-explicit hydrogen graph where each
//! hydrogen is an atom node.

use std::collections::HashMap;

use chematic_core::{
    Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule, MoleculeBuilder, implicit_hcount,
};

/// Return a new molecule in which every implicit hydrogen is converted to an
/// explicit H atom node bonded to its parent heavy atom.
///
/// The heavy atoms in the returned molecule have `hydrogen_count = Some(0)`,
/// preventing further implicit-H generation.  All original bonds and atom
/// properties are preserved.
pub fn add_hydrogens(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Copy heavy atoms with hydrogen_count set to Some(0).
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();
        atom.hydrogen_count = Some(0);
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    // Copy all original bonds.
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(na, nb, bond.order);
        }
    }

    // Add explicit H atoms for each implicit hydrogen.
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let h_count = implicit_hcount(mol, old_idx);
        if h_count == 0 {
            continue;
        }
        let heavy_new = remap[&old_idx];
        for _ in 0..h_count {
            let h_atom = Atom::new(Element::H);
            let h_new = builder.add_atom(h_atom);
            let _ = builder.add_bond(heavy_new, h_new, BondOrder::Single);
        }
    }

    builder.build()
}

/// Return a new molecule in which explicit H atom nodes are removed and their
/// bonds are converted back to implicit hydrogens.
///
/// Only explicit hydrogen atoms (nodes with `element == H`) are removed.
/// Chirality annotations and other atom properties are preserved.
///
/// Heavy atoms that had explicit H neighbors will have `hydrogen_count` set
/// to `None` so that implicit H is recomputed from valence.
pub fn remove_hydrogens(mol: &Molecule) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        if mol.atom(old_idx).element == Element::H {
            continue;
        }
        let mut atom = mol.atom(old_idx).clone();
        // Restore implicit H computation for atoms that had explicit H set.
        if atom.hydrogen_count == Some(0) {
            atom.hydrogen_count = None;
        }
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    // Copy heavy–heavy bonds only.
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        let a1_is_h = mol.atom(bond.atom1).element == Element::H;
        let a2_is_h = mol.atom(bond.atom2).element == Element::H;
        if a1_is_h || a2_is_h {
            continue;
        }
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(na, nb, bond.order);
        }
    }

    builder.build()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn add_h_methane_atom_count() {
        // C → 1 C + 4 H = 5 atoms
        let m = add_hydrogens(&mol("C"));
        assert_eq!(m.atom_count(), 5, "methane + H should have 5 atoms");
    }

    #[test]
    fn add_h_methane_bond_count() {
        // 4 C-H bonds
        let m = add_hydrogens(&mol("C"));
        assert_eq!(m.bond_count(), 4, "methane + H should have 4 bonds");
    }

    #[test]
    fn add_h_ethane() {
        // CC → 2 C + 6 H = 8 atoms, 1 C-C + 6 C-H = 7 bonds
        let m = add_hydrogens(&mol("CC"));
        assert_eq!(m.atom_count(), 8, "ethane + H atoms");
        assert_eq!(m.bond_count(), 7, "ethane + H bonds");
    }

    #[test]
    fn add_h_benzene() {
        // c1ccccc1 → 6 C + 6 H = 12 atoms, 6 ring + 6 C-H = 12 bonds
        let m = add_hydrogens(&mol("c1ccccc1"));
        assert_eq!(m.atom_count(), 12, "benzene + H atoms");
        assert_eq!(m.bond_count(), 12, "benzene + H bonds");
    }

    #[test]
    fn add_remove_roundtrip_ethanol() {
        let orig = mol("CCO");
        let with_h = add_hydrogens(&orig);
        let restored = remove_hydrogens(&with_h);
        // Heavy-atom count and bond count should match original.
        assert_eq!(
            restored.atom_count(),
            orig.atom_count(),
            "roundtrip atom count"
        );
        assert_eq!(
            restored.bond_count(),
            orig.bond_count(),
            "roundtrip bond count"
        );
    }

    #[test]
    fn add_remove_roundtrip_aspirin() {
        let orig = mol("CC(=O)Oc1ccccc1C(=O)O");
        let with_h = add_hydrogens(&orig);
        let restored = remove_hydrogens(&with_h);
        assert_eq!(restored.atom_count(), orig.atom_count());
        assert_eq!(restored.bond_count(), orig.bond_count());
    }

    #[test]
    fn remove_h_no_h_atoms_unchanged() {
        // A molecule with no explicit H nodes: remove_hydrogens should be a no-op.
        let orig = mol("CC");
        let result = remove_hydrogens(&orig);
        assert_eq!(result.atom_count(), 2);
        assert_eq!(result.bond_count(), 1);
    }

    #[test]
    fn add_h_water() {
        // O → 1 O + 2 H = 3 atoms, 2 bonds
        let m = add_hydrogens(&mol("O"));
        assert_eq!(m.atom_count(), 3);
        assert_eq!(m.bond_count(), 2);
    }

    #[test]
    fn add_h_preserves_element_distribution() {
        // Aspirin: 9 C + 4 O = 13 heavy; 8 H added → 21 total
        let orig = mol("CC(=O)Oc1ccccc1C(=O)O");
        let with_h = add_hydrogens(&orig);
        let h_count = with_h
            .atoms()
            .filter(|(_, a)| a.element == Element::H)
            .count();
        assert_eq!(h_count, 8, "aspirin should gain 8 H atoms (C9H8O4)");
    }
}
