//! Bounded RDKit Cleanup-compatible transformations.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder, implicit_hcount};

use super::{clone_molecule, is_metal, normalize_groups};
use crate::hydrogen::remove_hydrogens;

/// Apply the supported RDKit `MolStandardize.Cleanup` stages without changing
/// the broader [`super::standardize`] profile.
///
/// The input is stripped of removable explicit hydrogens, metal--ligand bonds
/// are disconnected with their bond-order charge transfer, and functional
/// groups are normalized. RDKit's acid/base-priority reionization is not yet
/// implemented here; molecules requiring that transfer can still differ.
pub fn cleanup(mol: &Molecule) -> Molecule {
    let without_h = remove_hydrogens(mol);
    let disconnected = disconnect_metals_for_cleanup(&without_h);
    let normalized = normalize_groups(&disconnected);
    let normalized = if disconnected.bond_count() < without_h.bond_count()
        && chematic_smiles::rdkit_smiles(
            &normalized,
            &chematic_smiles::RdkitSmilesParams::default(),
        )
        .is_err()
    {
        // Unsupported coordination and valence patterns must not turn a
        // readable input into an invalid result. Retain the original metal
        // bonds while keeping normalization of the other functional groups.
        normalize_groups(&without_h)
    } else {
        normalized
    };
    charge_originally_isolated_metals(&without_h, normalized)
}

/// RDKit Cleanup charges a bare Li/Na/K atom to +1 and a bare Mg/Ca atom
/// to +2. Restrict the rule to atoms already isolated before metal-bond
/// disconnection, so a newly detached coordinate donor is not reinterpreted.
fn charge_originally_isolated_metals(before: &Molecule, normalized: Molecule) -> Molecule {
    debug_assert_eq!(before.atom_count(), normalized.atom_count());
    let changes: Vec<_> = before
        .atoms()
        .filter_map(|(idx, original)| {
            let charge = match original.element.atomic_number() {
                3 | 11 | 19 => 1,
                12 | 20 => 2,
                _ => return None,
            };
            let atom = normalized.atom(idx);
            (before.degree(idx) == 0
                && original.charge == 0
                && atom.charge == 0
                && atom.hydrogen_count.unwrap_or(0) == 0)
                .then_some((idx, charge))
        })
        .collect();
    if changes.is_empty() {
        return normalized;
    }
    // RDKit pairs one newly charged isolated alkali atom with the sole
    // alcohol/phenol/carboxylic OH in a simple neutral mixture. Keep this
    // bounded: multiple metal atoms, charged groups, or competing OH sites
    // need an acid/base assignment that this Cleanup profile does not have.
    let balancing_oxygen = if changes.len() == 1
        && changes[0].1 == 1
        && before.atoms().all(|(_, atom)| atom.charge == 0)
        && matches!(
            before.atom(changes[0].0).element.atomic_number(),
            3 | 11 | 19
        ) {
        let mut candidates = normalized.atoms().filter_map(|(idx, atom)| {
            let (_, bond_idx) = normalized.neighbors(idx).next()?;
            (atom.element.atomic_number() == 8
                && atom.charge == 0
                && atom.hydrogen_count.is_none()
                && normalized.degree(idx) == 1
                && before.degree(idx) == 1
                && implicit_hcount(&normalized, idx) == 1
                && normalized.bond(bond_idx).order == BondOrder::Single
                && normalized
                    .neighbors(idx)
                    .next()
                    .is_some_and(|(neighbor, _)| {
                        normalized.atom(neighbor).element.atomic_number() == 6
                    }))
            .then_some(idx)
        });
        candidates.next().filter(|_| candidates.next().is_none())
    } else {
        None
    };
    let mut builder = MoleculeBuilder::from_molecule(&normalized);
    for (idx, charge) in changes {
        builder.set_atom_charge(idx, charge);
    }
    if let Some(oxygen) = balancing_oxygen {
        builder.set_atom_charge(oxygen, -1);
    }
    builder.build()
}

/// Disconnect a covalent metal--ligand bond by assigning its bond order to
/// opposite formal charges. Coordinate bonds contribute no charge transfer.
/// This is deliberately separate from the historical standardization profile,
/// which has different post-disconnection charge and hydrogen handling.
fn disconnect_metals_for_cleanup(mol: &Molecule) -> Molecule {
    let mut atoms: Vec<_> = mol.atoms().map(|(_, atom)| atom.clone()).collect();
    let mut removed = HashSet::new();
    // The source may be written in Kekule form. RDKit sanitizes it before
    // selecting bonds for the default metal disconnector.
    let aromaticity = mol
        .bonds()
        .any(|(_, bond)| {
            let first = mol.atom(bond.atom1).element;
            let second = mol.atom(bond.atom2).element;
            (is_metal(first) && second.atomic_number() == 6)
                || (is_metal(second) && first.atomic_number() == 6)
        })
        .then(|| chematic_perception::assign_aromaticity(mol));
    for (bond_idx, bond) in mol.bonds() {
        let first_is_metal = is_metal(atoms[bond.atom1.0 as usize].element);
        let second_is_metal = is_metal(atoms[bond.atom2.0 as usize].element);
        if first_is_metal == second_is_metal {
            continue;
        }
        let (metal, ligand) = if first_is_metal {
            (bond.atom1.0 as usize, bond.atom2.0 as usize)
        } else {
            (bond.atom2.0 as usize, bond.atom1.0 as usize)
        };
        let ligand_atom = &atoms[ligand];
        let metal_z = atoms[metal].element.atomic_number();
        let ligand_z = ligand_atom.element.atomic_number();
        if !default_metal_disconnector_match(metal_z, ligand_z)
            || (ligand_z == 6
                && (ligand_atom.aromatic
                    || aromaticity
                        .as_ref()
                        .is_some_and(|model| model.is_atom_aromatic(AtomIdx(ligand as u32)))))
        {
            continue;
        }
        removed.insert(bond_idx);
        transfer_metal_ligand_charge(mol, &mut atoms, metal, ligand, bond.order);
    }
    if removed.is_empty() {
        return clone_molecule(mol);
    }

    rebuild_after_metal_disconnection(mol, atoms, &removed)
}

/// RDKit 2026.09.1's default MetalNof and MetalNon SMARTS element sets.
/// Alkali/alkaline-earth metals only disconnect N/O/F, while the smaller
/// transition-metal/Al set also disconnects C, halogens, and other ligands.
fn default_metal_disconnector_match(metal_z: u8, ligand_z: u8) -> bool {
    let metal_nof =
        matches!(metal_z, 3 | 4 | 11 | 12 | 13 | 19..=31 | 37..=50 | 55 | 56 | 72..=83 | 87 | 88);
    let metal_non = matches!(metal_z, 13 | 21..=30 | 39..=48 | 72..=79);
    let ligand_nof = matches!(ligand_z, 7..=9);
    let ligand_non = matches!(ligand_z, 5 | 6 | 14..=17 | 33..=35 | 51..=53 | 85);
    (metal_nof && ligand_nof) || (metal_non && ligand_non)
}

fn transfer_metal_ligand_charge(
    mol: &Molecule,
    atoms: &mut [chematic_core::Atom],
    metal: usize,
    ligand: usize,
    bond_order: chematic_core::BondOrder,
) {
    use chematic_core::BondOrder;

    let ligand_atom = &atoms[ligand];
    // Sanitized RDKit regards saturated ammonia-to-metal bonds as dative
    // even when the source SMILES encodes a plain single bond.
    let ammonia_donor = ligand_atom.element.atomic_number() == 7
        && ligand_atom.charge == 0
        && ligand_atom.hydrogen_count == Some(3)
        && mol.degree(AtomIdx(ligand as u32)) == 1;
    let phosphorus_donor = ligand_atom.element.atomic_number() == 15;
    let lithium_carbonyl_donor = atoms[metal].element.atomic_number() == 3
        && ligand_atom.element.atomic_number() == 8
        && ligand_atom.charge == 0
        && mol
            .neighbors(AtomIdx(ligand as u32))
            .any(|(neighbor, adjacent)| {
                neighbor != AtomIdx(metal as u32) && mol.bond(adjacent).order == BondOrder::Double
            });
    if lithium_carbonyl_donor {
        // RDKit treats the overvalent carbonyl O -> Li bond as dative.
        // Once it is removed, neutral isolated Li has one implicit H.
        atoms[metal].hydrogen_count = Some(1);
    }
    if bond_order == BondOrder::Dative
        || ammonia_donor
        || lithium_carbonyl_donor
        || phosphorus_donor
    {
        return;
    }
    let order = i8::try_from(bond_order.order_int()).unwrap_or(0);
    atoms[metal].charge = atoms[metal].charge.saturating_add(order);
    atoms[ligand].charge = atoms[ligand].charge.saturating_sub(order);
}

fn rebuild_after_metal_disconnection(
    mol: &Molecule,
    atoms: Vec<chematic_core::Atom>,
    removed: &HashSet<BondIdx>,
) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for atom in atoms {
        builder.add_atom(atom);
    }
    for (bond_idx, bond) in mol.bonds() {
        if removed.contains(&bond_idx) {
            continue;
        }
        if let Ok(new_bond_idx) = builder.add_bond(bond.atom1, bond.atom2, bond.order)
            && let Some(direction) = mol.bond_direction(bond_idx)
        {
            builder.set_bond_direction(new_bond_idx, direction);
            if let Some(anchor) = mol.bond_direction_anchor(bond_idx) {
                builder.set_bond_direction_anchor(new_bond_idx, anchor);
            }
        }
    }
    builder.copy_stereo_from(mol);
    builder.copy_stereo_groups_from(mol);
    // The atom indices are unchanged, but a removed metal bond is no longer a
    // valid stereo neighbour. Keep the remaining order and implicit slots.
    for (idx, _) in mol.atoms() {
        if let Some(order) = mol.stereo_neighbor_order(idx) {
            let retained: Vec<_> = order
                .iter()
                .copied()
                .filter(|&neighbor| {
                    neighbor >= mol.atom_count() as u32
                        || mol
                            .bond_between(idx, AtomIdx(neighbor))
                            .is_none_or(|(bond_idx, _)| !removed.contains(&bond_idx))
                })
                .collect();
            if retained.len() != order.len() {
                builder.set_stereo_neighbor_order(idx, retained);
            }
        }
    }
    builder.build()
}

#[cfg(test)]
mod tests {
    use super::cleanup;
    use chematic_core::Molecule;
    use chematic_smiles::parse;

    #[test]
    fn preserves_charge_for_permanganate_and_copper_cyanide() {
        for (input, expected) in [
            (
                "O=[Mn](=O)(=O)[O-].[K+]",
                "[K+].[Mn+7].[O-2].[O-2].[O-2].[O-2]",
            ),
            ("C(#N)[Cu]", "[C-]#N.[Cu+]"),
        ] {
            let mol = parse(input).unwrap();
            let cleaned = cleanup(&mol);
            let total_charge = |mol: &Molecule| {
                mol.atoms()
                    .map(|(_, atom)| i32::from(atom.charge))
                    .sum::<i32>()
            };
            assert_eq!(total_charge(&cleaned), total_charge(&mol), "{input}");
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
                Ok(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn balances_one_isolated_alkali_with_one_oxygen_acid() {
        // RDKit 2026.09.1 Cleanup pairs one originally bare alkali metal
        // with the only OH site in these simple disconnected mixtures.
        for (input, expected) in [
            ("[Li].CO", "C[O-].[Li+]"),
            ("[Na].CCO", "CC[O-].[Na+]"),
            ("[K].Oc1ccccc1", "[K+].[O-]c1ccccc1"),
            ("[Na].CC(=O)O", "CC(=O)[O-].[Na+]"),
            ("[Na].O", "O.[Na+]"),
            ("[Na+].CCO", "CCO.[Na+]"),
            ("[Na].OCCO", "OCCO.[Na+]"),
        ] {
            let cleaned = cleanup(&parse(input).unwrap());
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
                Ok(expected),
                "{input}"
            );
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleanup(&cleaned)).as_deref(),
                Ok(expected),
                "second Cleanup: {input}"
            );
        }
    }
}
