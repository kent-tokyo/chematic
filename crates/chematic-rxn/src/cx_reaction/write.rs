use super::*;
use std::collections::{BTreeMap, HashSet};

fn participant(mol: &Molecule, canonical: bool) -> Result<(String, Vec<AtomIdx>), ReactionCxError> {
    if mol.atom_count() == 0 {
        return Err(ReactionCxError::Invalid(
            "empty participant molecule".into(),
        ));
    }
    let result = if canonical {
        chematic_smiles::canonical_smiles_with_atom_order(mol)
    } else {
        chematic_smiles::write_with_atom_order(mol)
    };
    if result.1.len() != mol.atom_count() {
        return Err(ReactionCxError::Invalid(
            "writer did not return every atom".into(),
        ));
    }
    Ok(result)
}

fn append_groups(
    mol: &Molecule,
    text: &str,
    order: &[AtomIdx],
    offset: usize,
    groups: &mut BTreeMap<GroupKey, Vec<usize>>,
) -> Result<(), ReactionCxError> {
    if mol.stereo_groups().is_empty() {
        return Ok(());
    }
    let parsed =
        chematic_smiles::parse(text).map_err(|e| ReactionCxError::Invalid(e.to_string()))?;
    if parsed.atom_count() != order.len() {
        return Err(ReactionCxError::Invalid(
            "written atom count differs from atom order".into(),
        ));
    }
    let mut inverse = vec![usize::MAX; mol.atom_count()];
    for (written, atom) in order.iter().enumerate() {
        if atom.0 as usize >= inverse.len() || inverse[atom.0 as usize] != usize::MAX {
            return Err(ReactionCxError::Invalid(
                "writer atom order is not a permutation".into(),
            ));
        }
        inverse[atom.0 as usize] = written;
    }
    let mut seen = HashSet::new();
    for (group_index, group) in mol.stereo_groups().iter().enumerate() {
        if mol
            .stereo_group_bond_indices(group_index)
            .is_some_and(|bonds| !bonds.is_empty())
        {
            return Err(ReactionCxError::Unsupported("bond stereo groups".into()));
        }
        if group.atom_indices.is_empty() {
            return Err(ReactionCxError::Invalid("empty stereo group".into()));
        }
        let target = groups.entry(group_key(&group.kind)?).or_default();
        for atom in &group.atom_indices {
            if atom.0 as usize >= inverse.len() || !seen.insert(*atom) {
                return Err(ReactionCxError::Invalid(
                    "invalid or overlapping group member".into(),
                ));
            }
            let written = inverse[atom.0 as usize];
            if !mol.atom(*atom).chirality.is_tetrahedral()
                || !parsed
                    .atom(AtomIdx(written as u32))
                    .chirality
                    .is_tetrahedral()
            {
                return Err(ReactionCxError::Invalid(
                    "group member lacks emitted tetrahedral stereo".into(),
                ));
            }
            target.push(offset + written);
        }
    }
    Ok(())
}

pub(super) fn write(reaction: &Reaction, canonical: bool) -> Result<String, ReactionCxError> {
    let mut groups = BTreeMap::new();
    let mut offset = 0usize;
    let mut roles = Vec::new();
    for molecules in [&reaction.reactants, &reaction.agents, &reaction.products] {
        let mut pieces = molecules
            .iter()
            .map(|mol| participant(mol, canonical).map(|(text, order)| (text, order, mol)))
            .collect::<Result<Vec<_>, _>>()?;
        if canonical {
            pieces.sort_by(|a, b| a.0.cmp(&b.0));
        }
        let mut texts = Vec::new();
        for (text, order, mol) in pieces {
            append_groups(mol, &text, &order, offset, &mut groups)?;
            offset += order.len();
            texts.push(text);
        }
        roles.push(texts.join("."));
    }
    let mut output = roles.join(">");
    if !groups.is_empty() {
        let mut fields = Vec::new();
        for (key, mut atoms) in groups {
            atoms.sort_unstable();
            fields.push(format!(
                "{}:{}",
                group_head(key),
                atoms
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ));
        }
        output.push_str(" |");
        output.push_str(&fields.join(","));
        output.push('|');
    }
    Ok(output)
}
