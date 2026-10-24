use super::*;

pub(super) fn write_cleanup_smiles(
    mol: &Molecule,
    params: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    let mut model = parse::from_chematic(mol)?;
    let in_place_oxide = has_in_place_oxide(mol);
    let lithium_nitrogen_ring = has_lithium_nitrogen_ring(mol);
    let iron_cyclopentadienyl_dianion = has_iron_cyclopentadienyl_dianion(mol);
    let before_sanitize =
        (in_place_oxide || lithium_nitrogen_ring || iron_cyclopentadienyl_dianion)
            .then(|| model.clone());
    let initial_atom_count = model.atoms.len();
    let isolated_tags = isolated_nontetrahedral_tags(mol, &model);
    let sanitized = if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut model)
    } else {
        sanitize::remove_hs_and_sanitize(&mut model)
    };
    if let Err(error) = sanitized {
        return write_unsanitized_fallback(before_sanitize, in_place_oxide, params, error);
    }
    if lithium_nitrogen_ring
        && let Some(raw) = before_sanitize.as_ref()
        && model.atoms.len() == raw.atoms.len()
    {
        restore_lithium_nitrogen_ring(&mut model, raw);
    }
    if iron_cyclopentadienyl_dianion
        && let Some(raw) = before_sanitize.as_ref()
        && model.atoms.len() == raw.atoms.len()
    {
        restore_aromatic_state(&mut model, raw);
    }
    if model.atoms.len() == initial_atom_count {
        clear_cleanup_radicals(mol, &mut model);
    }
    stereo::legacy_stereo_perception_unflagged(&mut model);
    if model.atoms.len() == initial_atom_count {
        for (idx, tag, permutation) in isolated_tags {
            model.atoms[idx].chiral = tag;
            model.atoms[idx].chiral_perm = permutation;
        }
    }
    write::mol_to_smiles_owned(model, params)
}

fn has_in_place_oxide(mol: &Molecule) -> bool {
    mol.atoms()
        .filter(|(_, atom)| atom.aromatic && atom.element.atomic_number() == 7)
        .count()
        >= 4
        && mol.bonds().any(|(_, bond)| {
            let (a, b) = (mol.atom(bond.atom1), mol.atom(bond.atom2));
            bond.order == chematic_core::BondOrder::Double
                && ((a.aromatic
                    && a.element.atomic_number() == 6
                    && b.element.atomic_number() == 8)
                    || (b.aromatic
                        && b.element.atomic_number() == 6
                        && a.element.atomic_number() == 8))
        })
}

fn has_lithium_nitrogen_ring(mol: &Molecule) -> bool {
    mol.atoms().any(|(idx, atom)| {
        atom.element.atomic_number() == 3 && atom.charge == 1 && mol.degree(idx) == 1
    }) && mol.atoms().any(|(idx, atom)| {
        atom.element.atomic_number() == 7
            && atom.charge == -1
            && !atom.aromatic
            && mol.degree(idx) == 2
    })
}

fn has_iron_cyclopentadienyl_dianion(mol: &Molecule) -> bool {
    mol.atom_count() == 11
        && mol
            .atoms()
            .any(|(_, atom)| atom.element.atomic_number() == 26 && atom.charge == 2)
        && mol
            .atoms()
            .filter(|(_, atom)| atom.element.atomic_number() == 6 && atom.charge == -1)
            .count()
            >= 2
        && mol
            .atoms()
            .all(|(_, atom)| matches!(atom.element.atomic_number(), 6 | 26))
}

fn isolated_nontetrahedral_tags(
    mol: &Molecule,
    model: &mol::Mol,
) -> Vec<(usize, mol::ChiralTag, u32)> {
    model
        .atoms
        .iter()
        .enumerate()
        .filter(|(idx, atom)| {
            model.degree(*idx) == 0
                && atom.chiral.nontet().is_some()
                && mol
                    .stereo_neighbor_order(chematic_core::AtomIdx(*idx as u32))
                    .is_some_and(|order| order.is_empty())
        })
        .map(|(idx, atom)| (idx, atom.chiral, atom.chiral_perm))
        .collect()
}

fn write_unsanitized_fallback(
    before_sanitize: Option<mol::Mol>,
    in_place_oxide: bool,
    params: &RdkitSmilesParams,
    error: RdkitSmilesError,
) -> Result<String, RdkitSmilesError> {
    let Some(mut raw) = before_sanitize else {
        return Err(error);
    };
    raw.find_rings()?;
    if in_place_oxide {
        for bond in &mut raw.bonds {
            if bond.aromatic {
                bond.bt = mol::BondType::Single;
            }
        }
    }
    write::mol_to_smiles_owned(raw, params)
}

fn restore_lithium_nitrogen_ring(model: &mut mol::Mol, raw: &mol::Mol) {
    let ring_info = model.ring_info().clone();
    for (nitrogen, atom) in raw.atoms.iter().enumerate() {
        if atom.anum != 7 || atom.charge != -1 || atom.aromatic {
            continue;
        }
        let Some(&ring_id) = ring_info
            .atom_members(nitrogen)
            .iter()
            .min_by_key(|&&ring| ring_info.atom_rings[ring].len())
        else {
            continue;
        };
        for &idx in &ring_info.atom_rings[ring_id] {
            let shared = ring_info.atom_members(idx).iter().any(|&other| {
                other != ring_id
                    && ring_info.atom_rings[other]
                        .iter()
                        .all(|&member| model.atoms[member].aromatic)
            });
            if !shared {
                model.atoms[idx].aromatic = raw.atoms[idx].aromatic;
            }
        }
        for &idx in &ring_info.bond_rings[ring_id] {
            let shared = ring_info.bond_members(idx).iter().any(|&other| {
                other != ring_id
                    && ring_info.atom_rings[other]
                        .iter()
                        .all(|&member| model.atoms[member].aromatic)
            });
            if !shared {
                model.bonds[idx].bt = raw.bonds[idx].bt;
                model.bonds[idx].aromatic = raw.bonds[idx].aromatic;
            }
        }
    }
}

fn restore_aromatic_state(model: &mut mol::Mol, raw: &mol::Mol) {
    for (current, original) in model.atoms.iter_mut().zip(&raw.atoms) {
        current.aromatic = original.aromatic;
    }
    for (current, original) in model.bonds.iter_mut().zip(&raw.bonds) {
        current.bt = original.bt;
        current.aromatic = original.aromatic;
    }
}

fn clear_cleanup_radicals(source: &Molecule, model: &mut mol::Mol) {
    let has_metal = source
        .atoms()
        .any(|(_, atom)| periodic::is_metal(u32::from(atom.element.atomic_number())));
    let aromaticity = has_metal.then(|| chematic_perception::assign_aromaticity(source));
    for (idx, atom) in source.atoms() {
        let carbon_after_cut = has_metal
            && atom.element.atomic_number() == 6
            && !atom.aromatic
            && !aromaticity
                .as_ref()
                .is_some_and(|view| view.is_atom_aromatic(idx))
            && atom.charge == 0
            && atom.hydrogen_count == Some(0)
            && source.degree(idx) == 2
            && source
                .neighbors(idx)
                .map(|(_, bond)| source.bond(bond).order.order_int())
                .sum::<u8>()
                == 3;
        let phosphorus_after_cut = atom.element.atomic_number() == 15
            && atom.charge == 0
            && atom.hydrogen_count == Some(0)
            && source.degree(idx) == 3
            && source.neighbors(idx).any(|(neighbor, bond)| {
                source.atom(neighbor).element.atomic_number() == 8
                    && source.bond(bond).order == chematic_core::BondOrder::Double
            });
        if (carbon_after_cut || phosphorus_after_cut) && model.atoms[idx.0 as usize].radicals == 1 {
            model.atoms[idx.0 as usize].radicals = 0;
        }
    }
}
