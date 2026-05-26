//! Molecular standardization routines.
//!
//! Provides utilities for cleaning up molecular representations:
//! - Selecting the largest connected fragment.
//! - Neutralizing simple formal charges.

#![forbid(unsafe_code)]

use std::collections::{HashMap, VecDeque};

use chematic_core::{AtomIdx, Element, Molecule, MoleculeBuilder};

// ---------------------------------------------------------------------------
// Connected-component detection
// ---------------------------------------------------------------------------

/// Find all connected components of `mol` via BFS.
///
/// Returns a `Vec<Vec<AtomIdx>>` sorted in descending order by component size
/// (largest component first).
fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<AtomIdx>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }

        let mut component: Vec<AtomIdx> = Vec::new();
        let start_idx = AtomIdx(start as u32);
        visited[start] = true;
        let mut queue: VecDeque<AtomIdx> = VecDeque::new();
        queue.push_back(start_idx);

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for (neighbor, _bond_idx) in mol.neighbors(current) {
                let ni = neighbor.0 as usize;
                if !visited[ni] {
                    visited[ni] = true;
                    queue.push_back(neighbor);
                }
            }
        }

        components.push(component);
    }

    // Sort descending by component size so the largest fragment is first.
    components.sort_by(|a, b| b.len().cmp(&a.len()));
    components
}

// ---------------------------------------------------------------------------
// Largest fragment
// ---------------------------------------------------------------------------

/// Return a new `Molecule` containing only the largest connected fragment.
///
/// If the molecule has only one fragment (or is empty), a clone is returned.
pub fn largest_fragment(mol: &Molecule) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);

    // The first component is the largest after sorting.
    let largest = &components[0];

    // Build a remapping from old AtomIdx to new AtomIdx.
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut builder = MoleculeBuilder::new();

    for &old_idx in largest {
        let atom = mol.atom(old_idx).clone();
        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

    // Add bonds whose both endpoints are in the selected component.
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }

    builder.build()
}

// ---------------------------------------------------------------------------
// Charge neutralization
// ---------------------------------------------------------------------------

/// Neutralize simple formal charges in a molecule.
///
/// Rules applied:
/// - `[O-]` on a carbon neighbor: set charge to 0 and increment hydrogen_count by 1
///   (converts carboxylate to carboxylic acid).
/// - `[N+]` with at least one explicit H: remove one H and set charge to 0
///   (converts ammonium to amine).
/// - `[O+]` with at least one explicit H: remove one H and set charge to 0
///   (converts protonated ether to neutral ether).
///
/// Returns a new `Molecule` with modifications applied.
pub fn neutralize_charges(mol: &Molecule) -> Molecule {
    // Collect any modifications: (AtomIdx, new_charge, new_hydrogen_count)
    let mut modifications: HashMap<AtomIdx, (i8, Option<u8>)> = HashMap::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);

        match (atom.element, atom.charge) {
            // [O-]: neutralize if it has a carbon neighbor.
            (Element::O, -1) => {
                let has_carbon_neighbor = mol
                    .neighbors(idx)
                    .any(|(nb, _)| mol.atom(nb).element == Element::C);
                if has_carbon_neighbor {
                    let current_h = atom.hydrogen_count.unwrap_or(0);
                    modifications.insert(idx, (0, Some(current_h + 1)));
                }
            }
            // [N+]: neutralize if it carries at least one explicit H.
            (Element::N, 1) => {
                let current_h = atom.hydrogen_count.unwrap_or(0);
                if current_h > 0 {
                    modifications.insert(idx, (0, Some(current_h - 1)));
                }
            }
            // [O+]: neutralize if it carries at least one explicit H.
            (Element::O, 1) => {
                let current_h = atom.hydrogen_count.unwrap_or(0);
                if current_h > 0 {
                    modifications.insert(idx, (0, Some(current_h - 1)));
                }
            }
            _ => {}
        }
    }

    // Rebuild molecule, applying modifications where applicable.
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    for i in 0..mol.atom_count() {
        let old_idx = AtomIdx(i as u32);
        let mut atom = mol.atom(old_idx).clone();

        if let Some(&(new_charge, new_h)) = modifications.get(&old_idx) {
            atom.charge = new_charge;
            atom.hydrogen_count = new_h;
        }

        let new_idx = builder.add_atom(atom);
        remap.insert(old_idx, new_idx);
    }

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
    fn largest_fragment_two_fragments_picks_larger() {
        // "CC.CCC" — ethane (2 C) and propane (3 C)
        let mol = parse("CC.CCC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 3, "should keep propane (3 C)");
    }

    #[test]
    fn largest_fragment_single_fragment_unchanged() {
        // "CC" — ethane, only one fragment
        let mol = parse("CC").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 2);
    }

    #[test]
    fn largest_fragment_keeps_benzene_over_ethane() {
        // "CC.c1ccccc1" — ethane (2 C) vs benzene (6 C)
        let mol = parse("CC.c1ccccc1").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 6, "should keep benzene (6 atoms)");
    }

    #[test]
    fn largest_fragment_ionic_pair_keeps_one_atom() {
        // "[Na+].[Cl-]" — both fragments are single atoms; either is fine
        let mol = parse("[Na+].[Cl-]").unwrap();
        let result = largest_fragment(&mol);
        assert_eq!(result.atom_count(), 1);
    }

    #[test]
    fn neutralize_neutral_molecule_unchanged() {
        // "CC" is already neutral; no atom should gain/lose charge
        let mol = parse("CC").unwrap();
        let result = neutralize_charges(&mol);
        for i in 0..result.atom_count() {
            let atom = result.atom(AtomIdx(i as u32));
            assert_eq!(atom.charge, 0, "all atoms should remain neutral");
        }
    }

    #[test]
    fn neutralize_acetate_oxygen() {
        // "CC(=O)[O-]" — acetate; the [O-] should become neutral with H added
        let mol = parse("CC(=O)[O-]").unwrap();
        let result = neutralize_charges(&mol);

        // Find the oxygen that was originally [O-]: it should now have charge 0
        // and hydrogen_count == Some(1).
        let neutralized_o = (0..result.atom_count())
            .map(|i| result.atom(AtomIdx(i as u32)))
            .find(|a| a.element == Element::O && a.hydrogen_count == Some(1));

        assert!(
            neutralized_o.is_some(),
            "neutralized [O-] should have hydrogen_count == Some(1)"
        );
        assert_eq!(
            neutralized_o.unwrap().charge,
            0,
            "neutralized [O-] should have charge == 0"
        );
    }
}
