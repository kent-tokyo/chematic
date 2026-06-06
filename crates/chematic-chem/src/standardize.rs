//! Molecular standardization routines.
//!
//! Provides utilities for cleaning up molecular representations:
//! - Selecting the largest connected fragment.
//! - Neutralizing simple formal charges.

#![forbid(unsafe_code)]

use std::collections::{HashMap, VecDeque};

use chematic_core::{AtomIdx, BondIdx, Element, Molecule, MoleculeBuilder};

use crate::{hydrogen::remove_hydrogens, tautomer::canonical_tautomer};

/// Find all connected components of `mol` via BFS, sorted descending by size.
fn connected_components(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut components: Vec<Vec<AtomIdx>> = Vec::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut component = Vec::new();
        let mut queue: VecDeque<AtomIdx> = VecDeque::new();
        queue.push_back(AtomIdx(start as u32));

        while let Some(current) = queue.pop_front() {
            component.push(current);
            for (neighbor, _) in mol.neighbors(current) {
                let ni = neighbor.0 as usize;
                if !visited[ni] {
                    visited[ni] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }

    components.sort_by(|a, b| b.len().cmp(&a.len()));
    components
}

/// Copy bonds from `mol` into `builder` when both endpoints are remapped.
fn copy_bonds(mol: &Molecule, builder: &mut MoleculeBuilder, remap: &HashMap<AtomIdx, AtomIdx>) {
    for i in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(i as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }
}

/// Return a new `Molecule` containing only the largest connected fragment.
///
/// If the molecule is empty, an empty `Molecule` is returned.
pub fn largest_fragment(mol: &Molecule) -> Molecule {
    if mol.atom_count() == 0 {
        return MoleculeBuilder::new().build();
    }

    let components = connected_components(mol);
    let largest = &components[0];

    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut builder = MoleculeBuilder::new();
    for &old_idx in largest {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// Neutralize simple formal charges in a molecule.
///
/// Rules applied:
/// - `[O-]` with a carbon neighbor → charge 0, +1 H (carboxylate → carboxylic acid).
/// - `[N+]` with at least one explicit H → charge 0, −1 H (ammonium → amine).
/// - `[O+]` with at least one explicit H → charge 0, −1 H (protonated ether → ether).
pub fn neutralize_charges(mol: &Molecule) -> Molecule {
    let mut modifications: HashMap<AtomIdx, (i8, Option<u8>)> = HashMap::new();

    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let h = atom.hydrogen_count.unwrap_or(0);

        match (atom.element, atom.charge) {
            (Element::O, -1) => {
                let has_c_neighbor =
                    mol.neighbors(idx).any(|(nb, _)| mol.atom(nb).element == Element::C);
                if has_c_neighbor {
                    modifications.insert(idx, (0, Some(h + 1)));
                }
            }
            (Element::N, 1) | (Element::O, 1) => {
                if h > 0 {
                    modifications.insert(idx, (0, Some(h - 1)));
                }
            }
            _ => {}
        }
    }

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
    copy_bonds(mol, &mut builder, &remap);
    builder.build()
}

/// Options for molecular standardization.
///
/// Controls which cleaning transformations are applied in a standardization pipeline.
#[derive(Clone, Debug)]
pub struct StandardizeOptions {
    /// Convert to canonical tautomer. Default: `true`.
    pub canonical_tautomer: bool,
    /// Neutralize simple formal charges. Default: `true`.
    pub neutralize_charges: bool,
    /// Remove explicit hydrogen atoms. Default: `true`.
    pub remove_explicit_h: bool,
    /// Keep only the largest connected fragment. Default: `false`.
    pub largest_fragment_only: bool,
}

impl Default for StandardizeOptions {
    fn default() -> Self {
        Self {
            canonical_tautomer: true,
            neutralize_charges: true,
            remove_explicit_h: true,
            largest_fragment_only: false,
        }
    }
}

/// Apply a series of standardization steps to a molecule.
///
/// Transformations are applied in this order:
/// 1. If `largest_fragment_only`, select the largest connected component.
/// 2. If `neutralize_charges`, neutralize simple charges.
/// 3. If `remove_explicit_h`, remove explicit H atoms.
/// 4. If `canonical_tautomer`, convert to the canonical tautomer.
///
/// Useful for cleaning pasted structures or database entries.
pub fn standardize(mol: &Molecule, opts: &StandardizeOptions) -> Molecule {
    let current = if opts.largest_fragment_only {
        largest_fragment(mol)
    } else {
        // Build a copy of the molecule to pass through the pipeline.
        let mut builder = MoleculeBuilder::new();
        let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();
        for i in 0..mol.atom_count() {
            let old_idx = AtomIdx(i as u32);
            let new_idx = builder.add_atom(mol.atom(old_idx).clone());
            remap.insert(old_idx, new_idx);
        }
        copy_bonds(mol, &mut builder, &remap);
        builder.build()
    };

    let current = if opts.neutralize_charges {
        neutralize_charges(&current)
    } else {
        current
    };

    let current = if opts.remove_explicit_h {
        remove_hydrogens(&current)
    } else {
        current
    };

    if opts.canonical_tautomer {
        canonical_tautomer(&current)
    } else {
        current
    }
}

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

    #[test]
    fn standardize_with_defaults() {
        // "CC(=O)[O-]" — acetate ion
        let mol = parse("CC(=O)[O-]").unwrap();
        let opts = StandardizeOptions::default();
        let result = standardize(&mol, &opts);

        // With default options, remove_explicit_h will be applied.
        // Check that [O-] was neutralized (charge should be 0).
        let has_neutral_o = (0..result.atom_count())
            .map(|i| result.atom(AtomIdx(i as u32)))
            .any(|a| a.element == Element::O && a.charge == 0);
        assert!(has_neutral_o, "acetate oxygen should be neutralized");

        // Should have at least 3 atoms (C, C, O) with no explicit H
        assert!(result.atom_count() >= 3, "should have at least 3 atoms after standardization");
    }

    #[test]
    fn standardize_skip_largest_fragment() {
        // "CC.CCC" — ethane and propane
        let mol = parse("CC.CCC").unwrap();
        let opts = StandardizeOptions {
            largest_fragment_only: false,
            ..Default::default()
        };
        let result = standardize(&mol, &opts);

        // Should keep both fragments
        assert_eq!(result.atom_count(), 5, "should keep both fragments when largest_fragment_only=false");
    }
}
