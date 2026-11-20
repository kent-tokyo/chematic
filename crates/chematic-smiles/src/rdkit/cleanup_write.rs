use super::*;

pub(super) fn write_cleanup_smiles(
    mol: &Molecule,
    params: &RdkitSmilesParams,
) -> Result<String, RdkitSmilesError> {
    let mut model = parse::from_chematic(mol)?;
    let lithium_nitrogen_ring = has_lithium_nitrogen_ring(mol);
    let iron_cyclopentadienyl_dianion = has_iron_cyclopentadienyl_dianion(mol);
    let before_sanitize =
        (lithium_nitrogen_ring || iron_cyclopentadienyl_dianion).then(|| model.clone());
    let initial_atom_count = model.atoms.len();
    let isolated_tags = isolated_nontetrahedral_tags(mol, &model);
    if let Err(error) = sanitize::sanitize_cleanup_output(&mut model, has_added_hydrogens(mol)) {
        return write_unsanitized_fallback(before_sanitize, params, error);
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
    stereo::legacy_stereo_perception_unflagged(&mut model);
    if model.atoms.len() == initial_atom_count {
        for (idx, tag, permutation) in isolated_tags {
            model.atoms[idx].chiral = tag;
            model.atoms[idx].chiral_perm = permutation;
        }
    }
    write::mol_to_smiles_owned(model, params)
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
    params: &RdkitSmilesParams,
    error: RdkitSmilesError,
) -> Result<String, RdkitSmilesError> {
    let Some(mut raw) = before_sanitize else {
        return Err(error);
    };
    raw.find_rings()?;
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
