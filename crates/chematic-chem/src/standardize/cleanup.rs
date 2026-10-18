//! Bounded RDKit Cleanup-compatible transformations.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, Molecule, MoleculeBuilder};

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
    if disconnected.bond_count() < without_h.bond_count()
        && chematic_smiles::rdkit_smiles(
            &normalized,
            &chematic_smiles::RdkitSmilesParams::default(),
        )
        .is_err()
    {
        // Unsupported coordination and valence patterns must not turn a
        // readable input into an invalid result. Retain the original metal
        // bonds while keeping normalization of the other functional groups.
        return normalize_groups(&without_h);
    }
    normalized
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
        // RDKit's default disconnector leaves aromatic metal-carbon bonds and
        // ordinary Sn-C bonds attached. Its other metal-carbon cases include
        // aliphatic organometallics and metal carbonyls.
        if ligand_atom.element.atomic_number() == 1
            || (ligand_atom.element.atomic_number() == 6
                && (ligand_atom.aromatic
                    || aromaticity
                        .as_ref()
                        .is_some_and(|model| model.is_atom_aromatic(AtomIdx(ligand as u32)))
                    || atoms[metal].element.atomic_number() == 50))
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
