//! RDKit Cleanup-compatible standardization transformations.

mod acid_base;

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder, implicit_hcount};

use super::{clone_molecule, is_metal, normalize_groups};
use crate::hydrogen::remove_hydrogens;

/// Apply the supported RDKit `MolStandardize.Cleanup` stages without changing
/// the broader [`super::standardize`] profile.
///
/// The input is stripped of removable explicit hydrogens, metal--ligand bonds
/// are disconnected with valence-aware charge adjustment, and functional
/// groups are normalized. The ordered acid/base catalog then moves protons
/// according to the pinned RDKit Cleanup profile.
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
    let normalized = charge_originally_isolated_metals(&without_h, normalized);
    acid_base::reionize_for_cleanup(&normalized)
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

    let cut_metals: HashSet<_> = mol
        .bonds()
        .filter(|(bond_idx, _)| removed.contains(bond_idx))
        .map(|(_, bond)| {
            if is_metal(mol.atom(bond.atom1).element) {
                bond.atom1
            } else {
                bond.atom2
            }
        })
        .collect();
    for metal in cut_metals {
        let idx = metal.0 as usize;
        let target: u8 = match (atoms[idx].element.atomic_number(), atoms[idx].charge) {
            (3 | 11 | 19, 0) => 1,
            (11, 3) => 2,
            (12 | 20, 0) => 2,
            (50, 0) => 2,
            (50, 1) => 3,
            _ => 0,
        };
        let remaining = mol
            .neighbors(metal)
            .filter(|(_, bond_idx)| !removed.contains(bond_idx))
            .map(|(_, bond_idx)| mol.bond(bond_idx).order.order_int())
            .sum::<u8>();
        atoms[idx].hydrogen_count = Some(target.saturating_sub(remaining));
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

    let ligand_atom = mol.atom(AtomIdx(ligand as u32));
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
    assign_valence_aware_metal_charge(mol, atoms, metal, ligand, bond_order);
}

fn assign_valence_aware_metal_charge(
    mol: &Molecule,
    atoms: &mut [chematic_core::Atom],
    metal: usize,
    ligand: usize,
    bond_order: BondOrder,
) {
    let ligand_atom = mol.atom(AtomIdx(ligand as u32));
    let order = i8::try_from(bond_order.order_int()).unwrap_or(0);
    let ligand_idx = AtomIdx(ligand as u32);
    let cut_metal_bonds = mol
        .neighbors(ligand_idx)
        .filter(|(neighbor, _)| is_metal(mol.atom(*neighbor).element))
        .count();
    let valences = ligand_atom.element.normal_valences();
    if cut_metal_bonds == 1 && !valences.is_empty() && order > 0 {
        // RDKit's MetalDisconnector determines the ligand charge from its
        // post-cut valence and radical count. A simple charge transfer by bond
        // order overcharges unsaturated C and initially charged O ligands.
        let original = mol.atom(ligand_idx);
        let h = i16::from(implicit_hcount(mol, ligand_idx));
        let valence_before = h + mol
            .neighbors(ligand_idx)
            .map(|(_, bond)| mol.bond(bond).order.order_int() as i16)
            .sum::<i16>();
        let radical_before = if original.hydrogen_count.is_some() {
            (i16::from(valences[0]) + i16::from(original.charge) - valence_before).max(0)
        } else {
            0
        };
        let valence_after = valence_before + radical_before - i16::from(order);
        let allowed = valences
            .iter()
            .map(|&value| i16::from(value))
            .find(|&value| value >= valence_after)
            .unwrap_or_else(|| i16::from(*valences.last().unwrap()));
        let mut ligand_charge = valence_after.min(allowed) - allowed;
        let mut lone_electrons = 0;
        if radical_before % 2 == 1 && ligand_charge < 0 {
            ligand_charge += 1;
            lone_electrons = 1;
        }
        let total_degree_after = i16::try_from(mol.degree(ligand_idx) - 1).unwrap_or(0) + h;
        if ligand_charge == -1 && valence_after == total_degree_after + 1 {
            ligand_charge = 0;
        }
        let metal_gain = if ligand_charge == -1
            && ligand_atom.element.atomic_number() == 6
            && valence_after == total_degree_after + 2
        {
            0
        } else {
            (-ligand_charge + lone_electrons).max(0)
        };
        atoms[ligand].charge = i8::try_from(ligand_charge).unwrap_or(original.charge);
        if matches!(ligand_atom.element.atomic_number(), 6 | 7)
            && ligand_charge == 0
            && radical_before == 0
            && ligand_atom.hydrogen_count == Some(0)
        {
            atoms[ligand].hydrogen_count = None;
        }
        atoms[metal].charge = atoms[metal]
            .charge
            .saturating_add(i8::try_from(metal_gain).unwrap_or(0));
        return;
    }
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
    use chematic_core::{Molecule, implicit_hcount};
    use chematic_smiles::parse;

    #[test]
    fn cleanup_removes_formally_double_bonded_hydrogen_like_rdkit() {
        for (input, expected) in [
            ("[H]=C", "C"),
            ("[H]=C(O)C", "CCO"),
            (
                "[H]=C(c1ccc(Cl)cc1)c1cc(N)ccc1N1CCC(C)CC1",
                "CC1CCN(c2ccc(N)cc2Cc2ccc(Cl)cc2)CC1",
            ),
            (
                "[H]=C(O)[C@H](CC#Cc1ccccc1F)NCP(=O)(O)O",
                "O=P(O)(O)CN[C@H](CO)CC#Cc1ccccc1F",
            ),
        ] {
            let result = cleanup(&parse(input).unwrap());
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&result).as_deref(),
                Ok(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn cleanup_preserves_bridgehead_nitrogen_stereo() {
        let input = "Cn1cc(C2=NC[C@@]3(C[N@@]4CC[C@@H]3C4)O2)c2ccccc21";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(input)
        );
    }

    #[test]
    fn cleanup_recomputes_hydrogens_on_disconnected_metals() {
        for (input, expected) in [
            ("O->[Mg](<-O)(<-O)<-O", "O.O.O.O.[MgH2]"),
            ("C[N]([SnH])C", "C[N-]C.[SnH3+]"),
        ] {
            let output = cleanup(&parse(input).unwrap());
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&output).as_deref(),
                Ok(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn matches_rdkit_metal_charge_repair_for_permanganate_and_copper_cyanide() {
        for (input, expected) in [
            (
                "O=[Mn](=O)(=O)[O-].[K+]",
                "[K+].[Mn+8].[O-2].[O-2].[O-2].[O-2]",
            ),
            ("C(#N)[Cu]", "[C-]#N.[Cu]"),
        ] {
            let mol = parse(input).unwrap();
            let cleaned = cleanup(&mol);
            // RDKit's metal charge repair can change the net formal charge
            // of an overvalent input; the pinned output is the contract here.
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

    #[test]
    fn transfers_one_acid_proton_to_simple_carbanions_and_alkoxide() {
        // Expected structures were pinned against RDKit 2026.09.1 Cleanup.
        for (input, expected) in [
            ("[CH3-].CC(=O)O", "C.CC(=O)[O-]"),
            ("[CH2-]C.CC(=O)O", "CC.CC(=O)[O-]"),
            ("[O-]C.CC(=O)O", "CC(=O)[O-].CO"),
            ("[O-]c1ccccc1.CC(=O)O", "CC(=O)[O-].Oc1ccccc1"),
            ("[CH3-].Oc1ccccc1", "C.[O-]c1ccccc1"),
            ("[CH3-].CCO", "C.CC[O-]"),
            ("[CH3-].Oc1ccccc1.CCO", "C.CCO.[O-]c1ccccc1"),
            ("[CH3-].CC(=O)O.Oc1ccccc1", "C.CC(=O)[O-].Oc1ccccc1"),
            ("[CH3-].O=S(=O)(O)C", "C.CS(=O)(=O)[O-]"),
            ("[CH3-].Oc1ccccc1.CS(=O)(=O)O", "C.CS(=O)(=O)[O-].Oc1ccccc1"),
        ] {
            let before = parse(input).unwrap();
            let after = cleanup(&before);
            let charge = |m: &Molecule| {
                m.atoms()
                    .map(|(_, atom)| i32::from(atom.charge))
                    .sum::<i32>()
            };
            let hydrogen = |m: &Molecule| {
                m.atoms()
                    .map(|(idx, _)| u32::from(implicit_hcount(m, idx)))
                    .sum::<u32>()
            };
            assert_eq!(charge(&after), charge(&before), "{input}");
            assert_eq!(hydrogen(&after), hydrogen(&before), "{input}");
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&after).as_deref(),
                Ok(expected),
                "{input}"
            );
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleanup(&after)).as_deref(),
                Ok(expected),
                "second Cleanup: {input}"
            );
        }
    }

    #[test]
    fn ordered_acid_base_pairs_match_rdkit() {
        for (input, expected) in [
            ("[CH3-].O", "O.[CH3-]"),
            ("[CH3-].OCCO", "C.[O-]CCO"),
            ("[CH3-].CC(=O)O.CC(=O)O", "C.CC(=O)O.CC(=O)[O-]"),
            ("[NH2-].CC(=O)O", "CC(=O)O.[NH2-]"),
            ("CC(=O)[O-].CC(=O)O", "CC(=O)O.CC(=O)[O-]"),
            ("[CH3-].CCP(C)(=O)O", "CCP(C)(=O)O.[CH3-]"),
            ("[CH3-].Oc1ccccc1.B(O)(O)C", "C.CB([O-])O.Oc1ccccc1"),
        ] {
            let after = cleanup(&parse(input).unwrap());
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&after).as_deref(),
                Ok(expected),
                "{input}"
            );
        }
    }
}
