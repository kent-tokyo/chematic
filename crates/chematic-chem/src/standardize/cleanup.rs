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
/// groups are normalized. A bounded existing-anion acid/base transfer follows;
/// RDKit's complete acid/base-priority table is not yet implemented, so
/// molecules requiring broader transfer can still differ.
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
    // Metal disconnection has its own charge-transfer ordering. Do not run
    // organic acid/base migration over the newly detached fragments.
    let reionized = if without_h.atoms().any(|(_, atom)| is_metal(atom.element)) {
        normalized
    } else {
        reionize_existing_anion(&normalized)
    };
    charge_originally_isolated_metals(&without_h, reionized)
}

/// Move an already-present negative charge to a stronger, neutral oxygen
/// acid. Unlike the public heuristic [`super::reionize`], this bounded
/// Cleanup stage never creates a new ion pair: it preserves total charge and
/// only consumes a C/N anion with an available proton. This covers the
/// RDKit Cleanup route for carboxylic, phosphorus, sulfur and phenolic acids
/// without changing ordinary neutral acids or bases.
fn reionize_existing_anion(mol: &Molecule) -> Molecule {
    let mut acids: Vec<(u8, AtomIdx)> = mol
        .atoms()
        .filter_map(|(idx, atom)| {
            if atom.element.atomic_number() != 8
                || atom.charge != 0
                || mol.degree(idx) != 1
                || implicit_hcount(mol, idx) != 1
            {
                return None;
            }
            let (neighbor, bond_idx) = mol.neighbors(idx).next()?;
            if mol.bond(bond_idx).order != BondOrder::Single {
                return None;
            }
            let center = mol.atom(neighbor);
            let has_terminal_oxo = mol.neighbors(neighbor).any(|(other, adjacent)| {
                other != idx
                    && mol.bond(adjacent).order == BondOrder::Double
                    && mol.atom(other).element.atomic_number() == 8
            });
            let priority = match center.element.atomic_number() {
                6 if has_terminal_oxo => 0,  // carboxylic acid
                15 if has_terminal_oxo => 1, // phosphorus acid
                16 if has_terminal_oxo => 2, // sulfur acid
                6 if center.aromatic => 3,   // phenol
                6 => 4,                      // alcohol
                _ => return None,
            };
            Some((priority, idx))
        })
        .collect();
    if acids.is_empty() {
        return clone_molecule(mol);
    }
    acids.sort_by_key(|&(priority, idx)| (priority, idx.0));
    let Some(donor) = mol.atoms().find_map(|(idx, atom)| {
        (atom.charge == -1
            && !atom.aromatic
            && implicit_hcount(mol, idx) > 0
            && matches!(atom.element.atomic_number(), 6 | 7))
        .then_some(idx)
    }) else {
        return clone_molecule(mol);
    };
    let target = acids[0].1;
    if donor == target {
        return clone_molecule(mol);
    }
    let mut atoms: Vec<_> = mol.atoms().map(|(_, atom)| atom.clone()).collect();
    let donor_atom = &mut atoms[donor.0 as usize];
    donor_atom.charge = donor_atom.charge.saturating_add(1);
    if let Some(hydrogens) = donor_atom.hydrogen_count {
        donor_atom.hydrogen_count = Some(hydrogens.saturating_add(1));
    }
    let target_atom = &mut atoms[target.0 as usize];
    target_atom.charge = target_atom.charge.saturating_sub(1);
    if let Some(hydrogens) = target_atom.hydrogen_count {
        target_atom.hydrogen_count = Some(hydrogens.saturating_sub(1));
    }
    rebuild_cleanup_graph(mol, atoms, &HashSet::new())
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

    // RDKit 2026.09.1's default metal disconnector represents permanganate
    // as Mn(+8) with four O(-2) atoms after removing the three double and one
    // single Mn--O bonds. Keep this narrow; it is a representation contract,
    // not an oxidation-state calculation.
    for (idx, original) in mol.atoms() {
        if original.element.atomic_number() != 25 || mol.degree(idx) != 4 {
            continue;
        }
        let mut doubles = 0usize;
        let mut singles = 0usize;
        let all_oxygen = mol.neighbors(idx).all(|(neighbor, bond_idx)| {
            match mol.bond(bond_idx).order {
                BondOrder::Double => doubles += 1,
                BondOrder::Single => singles += 1,
                _ => {}
            }
            mol.atom(neighbor).element.atomic_number() == 8
        });
        if all_oxygen && doubles == 3 && singles == 1 {
            atoms[idx.0 as usize].charge = 8;
        }
    }

    rebuild_cleanup_graph(mol, atoms, &removed)
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
    let copper_nitrile_donor = atoms[metal].element.atomic_number() == 29
        && ligand_atom.element.atomic_number() == 6
        && mol
            .neighbors(AtomIdx(ligand as u32))
            .any(|(neighbor, adjacent)| {
                neighbor != AtomIdx(metal as u32)
                    && mol.atom(neighbor).element.atomic_number() == 7
                    && mol.bond(adjacent).order == BondOrder::Triple
            });
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
    if bond_order == BondOrder::Dative && atoms[metal].element.atomic_number() == 12 {
        // RDKit's post-disconnection neutral magnesium is represented as
        // `[MgH2]`; record the implicit valence explicitly once its dative
        // oxygen donors have been removed.
        atoms[metal].hydrogen_count = Some(2);
    }
    if copper_nitrile_donor {
        // RDKit's default disconnector assigns cyanide to `[C-]#N` but keeps
        // the detached copper atom neutral in this representation.
        atoms[ligand].charge = atoms[ligand].charge.saturating_sub(1);
        return;
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

fn rebuild_cleanup_graph(
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
    fn matches_rdkit_charge_assignment_for_permanganate_and_copper_cyanide() {
        for (input, expected, expected_total_charge) in [
            (
                "O=[Mn](=O)(=O)[O-].[K+]",
                "[K+].[Mn+8].[O-2].[O-2].[O-2].[O-2]",
                1,
            ),
            ("C(#N)[Cu]", "[C-]#N.[Cu]", -1),
        ] {
            let mol = parse(input).unwrap();
            let cleaned = cleanup(&mol);
            let total_charge = |mol: &Molecule| {
                mol.atoms()
                    .map(|(_, atom)| i32::from(atom.charge))
                    .sum::<i32>()
            };
            assert_eq!(total_charge(&cleaned), expected_total_charge, "{input}");
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
                Ok(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn cleanup_reionizes_only_when_an_existing_anion_can_supply_a_proton() {
        for (input, expected) in [
            ("[CH2-]CC(=O)O", "CCC(=O)[O-]"),
            ("[NH-]C(C)C(=O)O", "CC(N)C(=O)[O-]"),
            ("[CH2-]CCO", "CCC[O-]"),
            ("CC(=O)O.CN", "CC(=O)O.CN"),
        ] {
            let cleaned = cleanup(&parse(input).unwrap());
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
                Ok(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn cleanup_keeps_dative_magnesium_valence_after_disconnection() {
        let cleaned = cleanup(&parse("O->[Mg](<-O)(<-O)<-O").unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok("O.O.O.O.[MgH2]")
        );
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
