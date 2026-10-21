//! RDKit Cleanup-compatible standardization transformations.

mod acid_base;

use std::collections::{HashSet, VecDeque};

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, MoleculeBuilder, SquarePlanarPermutation,
    implicit_hcount,
};

use super::{clone_molecule, is_metal, normalize_groups};
use crate::hydrogen::remove_hydrogens;

pub(super) fn reionize_with_catalog(mol: &Molecule) -> Molecule {
    acid_base::reionize_for_cleanup(mol)
}

/// Apply the supported RDKit `MolStandardize.Cleanup` stages without changing
/// the broader [`super::standardize`] profile.
///
/// The input is stripped of removable explicit hydrogens, metal--ligand bonds
/// are disconnected with valence-aware charge adjustment, and functional
/// groups are normalized. The ordered acid/base catalog then moves protons
/// according to the pinned RDKit Cleanup profile.
pub fn cleanup(mol: &Molecule) -> Molecule {
    let without_h = prepare_alkali_bridge_aromaticity(remove_hydrogens(mol));
    let disconnected = disconnect_metals_for_cleanup(&without_h);
    let normalized = normalize_for_cleanup(&disconnected);
    let normalized = if disconnected.bond_count() < without_h.bond_count()
        && chematic_smiles::rdkit_cleanup_smiles(
            &normalized,
            &chematic_smiles::RdkitSmilesParams::default(),
        )
        .is_err()
    {
        // Unsupported coordination and valence patterns must not turn a
        // readable input into an invalid result. Retain the original metal
        // bonds while keeping normalization of the other functional groups.
        normalize_for_cleanup(&without_h)
    } else {
        normalized
    };
    let normalized = charge_originally_isolated_metals(&without_h, normalized);
    acid_base::reionize_for_cleanup(&normalized)
}

fn prepare_alkali_bridge_aromaticity(mol: Molecule) -> Molecule {
    let shared_alkali_nitrogen = mol.atoms().any(|(idx, atom)| {
        atom.element.atomic_number() == 7
            && atom.aromatic
            && [3_u8, 11_u8].into_iter().all(|element| {
                mol.neighbors(idx)
                    .any(|(neighbor, _)| mol.atom(neighbor).element.atomic_number() == element)
            })
    });
    if !shared_alkali_nitrogen
        || mol.atom_count() > u16::MAX as usize
        || mol.atoms().any(|(_, atom)| atom.atom_map.is_some())
    {
        return mol;
    }

    // RDKit sanitizes before MetalDisconnector. A metal-bridged ring written
    // with lowercase source atoms can become nonaromatic at that point.
    // Preserve the sanitized bond orders by round-tripping mapped atoms,
    // then restoring the original indices before cutting metal--N bonds.
    let mut tagged_builder = MoleculeBuilder::new();
    for (idx, atom) in mol.atoms() {
        let mut tagged = atom.clone();
        tagged.atom_map = Some(idx.0 as u16 + 1);
        tagged_builder.add_atom(tagged);
    }
    for (_, bond) in mol.bonds() {
        let _ = tagged_builder.add_bond(bond.atom1, bond.atom2, bond.order);
    }
    tagged_builder.copy_stereo_from(&mol);
    tagged_builder.copy_bond_directions_from(&mol);
    tagged_builder.copy_stereo_groups_from(&mol);
    let tagged = tagged_builder.build();
    let params = chematic_smiles::RdkitSmilesParams {
        canonical: false,
        ..Default::default()
    };
    let Ok(text) = chematic_smiles::rdkit_smiles(&tagged, &params) else {
        return mol;
    };
    let Ok(sanitized) = chematic_smiles::parse(&text) else {
        return mol;
    };
    if sanitized.atom_count() != mol.atom_count() {
        return mol;
    }
    let mut by_original = vec![None; mol.atom_count()];
    for (idx, atom) in sanitized.atoms() {
        let Some(map) = atom.atom_map else { return mol };
        let Some(original_idx) = usize::from(map).checked_sub(1) else {
            return mol;
        };
        let Some(slot) = by_original.get_mut(original_idx) else {
            return mol;
        };
        *slot = Some(idx);
    }
    if by_original.iter().any(Option::is_none) {
        return mol;
    }
    let mut prepared = mol.clone();
    for (idx, _) in mol.atoms() {
        prepared.set_atom_aromatic(
            idx,
            sanitized
                .atom(by_original[idx.0 as usize].unwrap())
                .aromatic,
        );
    }
    for (bond_idx, bond) in mol.bonds() {
        let a = by_original[bond.atom1.0 as usize].unwrap();
        let b = by_original[bond.atom2.0 as usize].unwrap();
        let Some((_, sanitized_bond)) = sanitized.bond_between(a, b) else {
            return mol.clone();
        };
        prepared.set_bond_order(bond_idx, sanitized_bond.order);
    }
    prepared
}

/// The pinned RDKit `Cleanup` uses `NormalizeInPlace`. On N-rich fused
/// aromatic oxido-imidazolium systems that API recombines the O-/N+ pair,
/// whereas the copy-returning `Normalize` API leaves it unchanged. Keep this
/// behavior confined to Cleanup so the independent Normalize contract stays
/// aligned with RDKit's copy API.
fn normalize_for_cleanup(mol: &Molecule) -> Molecule {
    let mut normalized = normalize_groups(mol);
    let mut visited = HashSet::new();
    let assigned = chematic_perception::assign_aromaticity(mol);
    let aromatic_atom = |idx| mol.atom(idx).aromatic || assigned.is_atom_aromatic(idx);
    let aromatic_bond =
        |idx| mol.bond(idx).order == BondOrder::Aromatic || assigned.is_bond_aromatic(idx);
    for (start, _) in mol.atoms() {
        if !aromatic_atom(start) || !visited.insert(start) {
            continue;
        }
        let mut queue = VecDeque::from([start]);
        let mut component = Vec::new();
        while let Some(current) = queue.pop_front() {
            component.push(current);
            for (neighbor, bond) in mol.neighbors(current) {
                if aromatic_bond(bond) && aromatic_atom(neighbor) && visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
        if component.len() < 10
            || component
                .iter()
                .filter(|&&idx| mol.atom(idx).element.atomic_number() == 7)
                .count()
                < 4
        {
            continue;
        }
        let cations: Vec<_> = component
            .iter()
            .copied()
            .filter(|&idx| {
                let atom = mol.atom(idx);
                atom.element.atomic_number() == 7 && atom.charge == 1 && mol.degree(idx) == 3
            })
            .collect();
        let oxides: Vec<_> = component
            .iter()
            .flat_map(|&carbon| {
                mol.neighbors(carbon).filter_map(move |(oxygen, bond)| {
                    let atom = mol.atom(oxygen);
                    (mol.atom(carbon).element.atomic_number() == 6
                        && atom.element.atomic_number() == 8
                        && atom.charge == -1
                        && mol.degree(oxygen) == 1
                        && mol.bond(bond).order == BondOrder::Single)
                        .then_some((carbon, oxygen, bond))
                })
            })
            .collect();
        if let ([cation], [(carbon, oxygen, bond)]) = (cations.as_slice(), oxides.as_slice()) {
            for &idx in &component {
                normalized.set_atom_aromatic(idx, true);
            }
            for (idx, edge) in mol.bonds() {
                if component.contains(&edge.atom1)
                    && component.contains(&edge.atom2)
                    && aromatic_bond(idx)
                {
                    normalized.set_bond_order(idx, BondOrder::Aromatic);
                }
            }
            normalized.set_charge(*cation, 0);
            normalized.set_charge(*oxygen, 0);
            normalized.set_bond_order(*bond, BondOrder::Double);
            debug_assert_eq!(mol.atom(*carbon).element.atomic_number(), 6);
        }
    }
    normalized
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
        let aromatic_carbon = ligand_z == 6
            && aromaticity
                .as_ref()
                .is_some_and(|model| model.is_atom_aromatic(AtomIdx(ligand as u32)));
        let zinc_fused_aromatic_edge = metal_z == 30
            && ligand_z == 6
            && bond.order == BondOrder::Double
            && mol.neighbors(AtomIdx(metal as u32)).any(|(neighbor, _)| {
                neighbor.0 as usize != ligand
                    && mol.atom(neighbor).element.atomic_number() == 6
                    && (mol.atom(neighbor).aromatic
                        || aromaticity
                            .as_ref()
                            .is_some_and(|model| model.is_atom_aromatic(neighbor)))
            });
        if !default_metal_disconnector_match(metal_z, ligand_z)
            || aromatic_carbon
            || zinc_fused_aromatic_edge
        {
            continue;
        }
        removed.insert(bond_idx);
        if metal_z == 78
            && ligand_z == 6
            && ligand_atom.aromatic
            && !aromatic_carbon
            && bond.order == BondOrder::Single
        {
            // A hypervalent aromatic source token can sanitize to an
            // ordinary nonaromatic carbon. RDKit cuts its Pt bond without
            // charge transfer and lets the resulting carbon regain H.
            atoms[ligand].hydrogen_count = None;
            continue;
        }
        transfer_metal_ligand_charge(mol, &mut atoms, metal, ligand, bond.order);
        if metal_z == 30
            && ligand_z == 6
            && mol.degree(AtomIdx(ligand as u32)) == 4
            && mol.atom(AtomIdx(ligand as u32)).hydrogen_count == Some(1)
        {
            atoms[ligand].chirality = match atoms[ligand].chirality {
                Chirality::Clockwise => Chirality::CounterClockwise,
                Chirality::CounterClockwise => Chirality::Clockwise,
                other => other,
            };
        }
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
        if mol.degree(metal) == 3
            && mol
                .neighbors(metal)
                .all(|(_, bond_idx)| removed.contains(&bond_idx))
            && let Chirality::SquarePlanar(permutation) = atoms[idx].chirality
        {
            let has_incoming_ligand = mol.neighbors(metal).any(|(neighbor, _)| neighbor < metal);
            let remapped = match (has_incoming_ligand, permutation) {
                (true, SquarePlanarPermutation::SP1) => SquarePlanarPermutation::SP2,
                (true, SquarePlanarPermutation::SP2) => SquarePlanarPermutation::SP3,
                (true, SquarePlanarPermutation::SP3) => SquarePlanarPermutation::SP1,
                (false, SquarePlanarPermutation::SP1) => SquarePlanarPermutation::SP1,
                (false, SquarePlanarPermutation::SP2) => SquarePlanarPermutation::SP3,
                (false, SquarePlanarPermutation::SP3) => SquarePlanarPermutation::SP2,
            };
            atoms[idx].chirality = Chirality::SquarePlanar(remapped);
        }
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
        let metal_gain = if atoms[metal].element.atomic_number() == 30
            && ligand_atom.element.atomic_number() == 6
            && radical_before > 0
            && ligand_charge == 0
        {
            // Breaking a Zn--carbon bond to a carbon already carrying a
            // radical increases the carbon's radical count, not Zn charge.
            0
        } else if ligand_charge == -1
            && ligand_atom.element.atomic_number() == 6
            && valence_after == total_degree_after + 2
        {
            0
        } else {
            (-ligand_charge + lone_electrons).max(0)
        };
        atoms[ligand].charge = i8::try_from(ligand_charge).unwrap_or(original.charge);
        if ligand_atom.element.atomic_number() == 7
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
                let isolated_three_ligand_square_planar = mol.degree(idx) == 3
                    && mol.neighbors(idx).all(|(_, bond)| removed.contains(&bond))
                    && matches!(mol.atom(idx).chirality, Chirality::SquarePlanar(_));
                builder.set_stereo_neighbor_order(
                    idx,
                    if isolated_three_ligand_square_planar {
                        Vec::new()
                    } else {
                        retained
                    },
                );
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
    fn cleanup_matches_rdkit_in_place_normalization_on_fused_oxido_imidazolium() {
        for (input, expected) in [
            (
                "CC1=CN2C(=N1)N=C1C(=C2[O-])N(C)C=[N+]1C[200O]CCO",
                "Cc1cn2c(=O)c3c(nc2n1)n(C[200O]CCO)cn3C",
            ),
            (
                "OC[CH-]OC[n+]1c2nc3nc(cn3c([O-])c2n(C)c1)[CH]",
                "[CH]c1cn2c(=O)c3c(nc2n1)n(COCC[O-])cn3C",
            ),
            (
                "Cc1:C:n2c([O-])c3c(nc2n1)[n+](COCCO)cn3C",
                "Cc1cn2c(=O)c3c(nc2n1)n(COCCO)cn3C",
            ),
        ] {
            let mol = parse(input).unwrap();
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleanup(&mol)).as_deref(),
                Ok(expected),
                "{input}"
            );
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&super::normalize_groups(&mol))
                    .unwrap(),
                chematic_smiles::rdkit_cleanup_canonical_smiles(&mol).unwrap(),
                "copy-returning Normalize must retain the RDKit distinction: {input}"
            );
        }
    }

    #[test]
    fn cleanup_unbrackets_pentavalent_phosphorus_after_metal_disconnection() {
        let input = "O=C(O)[C@H]1C[C@@H](C/C=C/[P](=O)(O)[Fe])CCN1";
        let expected = "O=C(O)[C@H]1C[C@@H](C/C=C/P(=O)O)CCN1.[Fe]";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_uses_sanitized_aromaticity_for_metal_carbon_bonds() {
        for (input, expected) in [
            (
                "Cc1c[c]([Pt])c(Oc2ccc(O)cc2CC(=O)O)s(Cl)c1",
                "CC1=CC=C(Oc2ccc(O)cc2CC(=O)O)S(Cl)=C1.[Pt]",
            ),
            (
                "N(CCCN(C)CCCNC(C1C2C(C=C(Cl)C=1)=NC1C(=CC=CC=1)N=2)=O)C(C1C2[C](=[Zn][C]3=CC=CC=C3N=2)C=C(Cl)C=1)=O",
                "CN(CCCNC(=O)c1cc(Cl)c[c]2c1=Nc1cccc[c]1[Zn]=2)CCCNC(=O)c1cc(Cl)cc2nc3ccccc3nc12",
            ),
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
    fn cleanup_does_not_charge_zinc_for_detaching_an_existing_carbon_radical() {
        let input = "COc1cccc2c1C[C]([Zn]([CH3])[CH2]C)CC2.c1ccccc1";
        let expected = "C.CC.COc1cccc2c1[CH-][C]C[CH-]2.[Zn+2].c1ccccc1";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_remaps_three_ligand_square_planar_stereo_after_disconnection() {
        for (input, expected) in [
            ("F[Pt@SP1](Cl)Br", "[Br-].[Cl-].[F-].[Pt@SP2+3]"),
            ("F[Pt@SP2](Cl)Br", "[Br-].[Cl-].[F-].[Pt@SP3+3]"),
            ("F[Pt@SP3](Cl)Br", "[Br-].[Cl-].[F-].[Pt@SP1+3]"),
            ("[Pt@SP1](F)(Cl)Br", "[Br-].[Cl-].[F-].[Pt@SP1+3]"),
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
    fn cleanup_keeps_bridgehead_parity_after_dative_zinc_disconnection() {
        let input = "C1CC2CN3CCc4cc5c(OCO5)cc4[C@@H]3C[C@H]2([Zn])C[C@@H]1O";
        let expected = "O[C@@H]1CCC2CN3CCc4cc5c(cc4[C@@H]3C[C@@H]2C1)OCO5.[Zn]";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_reionizes_nonaromatic_imide_carbanion_after_zinc_cut() {
        let input = "CC1=C(CCOC2=CC=C(C[C]3([Zn])OC(=O)NC3=O)C=C2)N=C(C2=C[C+]=CO2)O1";
        let expected = "Cc1oc(C2=C[C+]=CO2)nc1[CH-]COc1ccc(CC2OC(=O)NC2=O)cc1.[Zn+]";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_detaches_alkali_metals_from_shared_amine_without_losing_carbon_bonds() {
        let input =
            "O=[C](O)[Li]1[c](-c2ccc3c(c2)OCO3)[c]2cc(OCc3ccccc3)c[cH][Na]2[n]1Cc1ccc2c(c1)OCO2";
        let expected = "O=[C](O)[Li+][C](=[C]1C=C(OCc2ccccc2)C=[CH][Na+]1)c1ccc2c(c1)OCO2.[N-2]Cc1ccc2c(c1)OCO2";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_preserves_nonaromatic_lithium_nitrogen_ring_representation() {
        let input = "ClC1=CC=C(C2=CC3=CC=CC=C3C3=NC[CH2][Li][N]23)C=C1";
        let expected = "[Li+][CH2]CN=C1[N-]C(c2ccc(Cl)cc2)=Cc2ccccc21";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_reionizes_cyclopentadienyl_anions_after_iron_cut() {
        let input = "C1=CC=CC1[Fe]C1C=CC=C1";
        let expected = "C1=CCC=C1.[C-]1=C[CH-]C=C1.[Fe+2]";
        let cleaned = cleanup(&parse(input).unwrap());
        assert_eq!(
            chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
            Ok(expected)
        );
    }

    #[test]
    fn cleanup_keeps_isolated_cyclopentadienyl_anion() {
        for input in ["C1=CC=C[CH-]1", "[CH-]1C=CC=C1"] {
            let cleaned = cleanup(&parse(input).unwrap());
            assert_eq!(
                chematic_smiles::rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
                Ok("c1cc[cH-]c1"),
                "{input}"
            );
        }
    }

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
