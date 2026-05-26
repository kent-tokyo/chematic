//! Murcko scaffold decomposition.
//!
//! Provides functions to extract the Murcko scaffold from a molecule:
//! - `murcko_scaffold`: ring atoms plus atoms on paths connecting ring systems.
//! - `generic_murcko_scaffold`: scaffold with all atoms replaced by C and all bonds by Single.

#![forbid(unsafe_code)]

use std::collections::{HashMap, HashSet};

use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};
use chematic_perception::find_sssr;

// ---------------------------------------------------------------------------
// Murcko scaffold
// ---------------------------------------------------------------------------

/// Extract the Murcko scaffold from `mol`.
///
/// The scaffold consists of:
/// - All ring atoms (atoms that participate in at least one ring).
/// - Linker atoms: non-ring atoms on paths that connect two ring systems.
///   Identified by an iterative expansion: a non-ring atom is a linker if it
///   has at least two heavy-atom neighbors that are ring atoms or already-selected
///   linker atoms. The expansion continues until no new linkers are found.
///
/// Returns an empty `Molecule` if the input contains no rings.
pub fn murcko_scaffold(mol: &Molecule) -> Molecule {
    let rings = find_sssr(mol);

    if rings.ring_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    // Collect ring atoms into a HashSet for O(1) lookup.
    let mut scaffold_atoms: HashSet<AtomIdx> = HashSet::new();
    for ring in rings.rings() {
        for &atom_idx in ring {
            scaffold_atoms.insert(atom_idx);
        }
    }

    // Iterative linker-atom expansion.
    // A non-ring atom is a linker if it has >= 2 heavy-atom neighbors that are
    // in scaffold_atoms (ring or already-found linkers). Repeat until stable.
    loop {
        let mut changed = false;

        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if scaffold_atoms.contains(&idx) {
                // Already in scaffold; skip.
                continue;
            }

            // Count how many heavy-atom neighbors are already in the scaffold.
            let scaffold_neighbor_count = mol
                .neighbors(idx)
                .filter(|(nb, _)| scaffold_atoms.contains(nb))
                .count();

            if scaffold_neighbor_count >= 2 {
                scaffold_atoms.insert(idx);
                changed = true;
            }
        }

        if !changed {
            break;
        }
    }

    // Rebuild the molecule with only scaffold atoms.
    build_subgraph(mol, &scaffold_atoms)
}

/// Generic Murcko scaffold: replace every atom with C and every bond with Single.
///
/// First extracts the Murcko scaffold, then maps all atoms to carbon and all
/// bond orders to single. Returns an empty `Molecule` if no rings exist.
pub fn generic_murcko_scaffold(mol: &Molecule) -> Molecule {
    let scaffold = murcko_scaffold(mol);

    if scaffold.atom_count() == 0 {
        return scaffold;
    }

    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Replace every atom with a plain carbon.
    for i in 0..scaffold.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let new_idx = builder.add_atom(Atom::organic(Element::C));
        remap.insert(old_idx, new_idx);
    }

    // Replace every bond with a single bond.
    for i in 0..scaffold.bond_count() {
        let bond = scaffold.bond(chematic_core::BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, BondOrder::Single);
        }
    }

    builder.build()
}

// ---------------------------------------------------------------------------
// Internal helper
// ---------------------------------------------------------------------------

/// Build a subgraph of `mol` containing only the atoms in `atom_set`,
/// preserving bonds whose both endpoints are in the set.
fn build_subgraph(mol: &Molecule, atom_set: &HashSet<AtomIdx>) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Add atoms in original index order to keep a stable mapping.
    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        if atom_set.contains(&old_idx) {
            let atom = mol.atom(old_idx).clone();
            let new_idx = builder.add_atom(atom);
            remap.insert(old_idx, new_idx);
        }
    }

    // Add bonds whose both endpoints are in the subgraph.
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }

    builder.build()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn murcko_benzene_preserves_all_atoms() {
        // Benzene: all 6 atoms are ring atoms; scaffold == entire molecule.
        let mol = parse("c1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 6);
    }

    #[test]
    fn murcko_toluene_removes_methyl() {
        // Toluene: benzene ring (6 atoms) + one methyl group (1 C); methyl is a side chain.
        let mol = parse("Cc1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 6, "methyl group should be removed");
    }

    #[test]
    fn murcko_ethylbenzene_removes_chain() {
        // Ethylbenzene: benzene ring + two-atom chain; chain is a side chain.
        let mol = parse("CCc1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 6, "ethyl chain should be removed");
    }

    #[test]
    fn murcko_acyclic_returns_empty() {
        // Ethane has no rings; scaffold should be empty.
        let mol = parse("CC").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert_eq!(scaffold.atom_count(), 0);
    }

    #[test]
    fn generic_murcko_benzene_all_carbon_single() {
        // Generic scaffold of benzene: 6 C atoms, all bonds Single.
        let mol = parse("c1ccccc1").unwrap();
        let generic = generic_murcko_scaffold(&mol);
        assert_eq!(generic.atom_count(), 6);
        for i in 0..generic.atom_count() {
            let atom = generic.atom(AtomIdx(i as u32));
            assert_eq!(atom.element, Element::C, "all atoms should be carbon");
        }
        for i in 0..generic.bond_count() {
            let bond = generic.bond(chematic_core::BondIdx(i as u32));
            assert_eq!(
                bond.order,
                BondOrder::Single,
                "all bonds should be Single"
            );
        }
    }

    #[test]
    fn murcko_biphenyl_keeps_all_ring_atoms() {
        // Biphenyl: two fused/connected phenyl rings, all atoms are ring atoms.
        let mol = parse("c1ccccc1c1ccccc1").unwrap();
        let scaffold = murcko_scaffold(&mol);
        assert!(
            scaffold.atom_count() >= 12,
            "biphenyl scaffold should have at least 12 atoms, got {}",
            scaffold.atom_count()
        );
    }
}
