//! Checked concrete-core removal with an explicit query-to-source atom vector.
use chematic_core::{AtomIdx, BondIdx, Molecule};

use crate::bond_fragment::{fragment_single_bonds_checked, validate_graph};

/// One surviving side-chain attachment, labelled by the core atom index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoreAttachment {
    pub core_atom: AtomIdx,
    pub source_core_atom: AtomIdx,
    pub source_anchor: AtomIdx,
    pub source_bond: BondIdx,
    pub anchor: AtomIdx,
    pub dummy: AtomIdx,
}

/// New graph and order-preserving mapping of original atoms; removed atoms are None.
#[derive(Clone)]
pub struct CoreReplacement {
    pub molecule: Molecule,
    pub original_to_output: Vec<Option<AtomIdx>>,
    pub attachments: Vec<CoreAttachment>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoreReplacementError(pub String);

impl std::fmt::Display for CoreReplacementError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for CoreReplacementError {}

fn validate_match(
    mol: &Molecule,
    core: &Molecule,
    matching: &[AtomIdx],
) -> Result<Vec<Option<AtomIdx>>, CoreReplacementError> {
    if mol.atom_count() > 4096
        || mol.bond_count() > 8192
        || core.atom_count() == 0
        || core.atom_count() > 1024
        || core.bond_count() > 8192
    {
        return Err(CoreReplacementError("core removal requires 1..1024 core atoms, at most 4096 source atoms and 8192 bonds per graph".into()));
    }
    if matching.len() != core.atom_count() {
        return Err(CoreReplacementError(
            "match must contain one source atom per core atom".into(),
        ));
    }
    for graph in [mol, core] {
        validate_graph(graph).map_err(|e| CoreReplacementError(format!("core removal: {e}")))?;
    }
    let mut inverse = vec![None; mol.atom_count()];
    for (i, &source) in matching.iter().enumerate() {
        let Some(slot) = inverse.get_mut(source.0 as usize) else {
            return Err(CoreReplacementError(
                "matched source atom is out of range".into(),
            ));
        };
        if slot.is_some() {
            return Err(CoreReplacementError(
                "matched source atoms must be distinct".into(),
            ));
        }
        *slot = Some(AtomIdx(i as u32));
        let a = mol.atom(source);
        let b = core.atom(AtomIdx(i as u32));
        if a.element != b.element
            || a.aromatic != b.aromatic
            || a.charge != b.charge
            || a.isotope != b.isotope
        {
            return Err(CoreReplacementError(
                "matched atoms must agree in element, aromaticity, charge and isotope".into(),
            ));
        }
    }
    for (_, bond) in core.bonds() {
        if !mol
            .bond_between(
                matching[bond.atom1.0 as usize],
                matching[bond.atom2.0 as usize],
            )
            .is_some_and(|(_, source)| source.order == bond.order)
        {
            return Err(CoreReplacementError(
                "match does not preserve a core bond/order".into(),
            ));
        }
    }
    Ok(inverse)
}

/// Remove an explicitly mapped concrete core, labelling dummies by core index.
///
/// `matching[i]` is the source atom corresponding to core atom `i`. Only bridge
/// single bonds can cross the core/side-chain boundary, with at most 1024 cuts.
/// Both graphs must be concrete, standard-bond and free of stereo annotations.
/// Matching requires equal element/aromaticity/charge/isotope and core bond
/// orders, ignoring H counts and atom maps. No matching search is performed.
/// Source atoms not in the core retain order and attributes; returned mappings
/// identify their shifted indices. Dummy isotope zero is unlabelled. Unrelated
/// components are retained, and removing the whole source yields an empty Mol.
pub fn replace_core_checked(
    mol: &Molecule,
    core: &Molecule,
    matching: &[AtomIdx],
) -> Result<CoreReplacement, CoreReplacementError> {
    let inverse = validate_match(mol, core, matching)?;
    let mut cuts = Vec::new();
    let mut labels = Vec::new();
    for (index, bond) in mol.bonds() {
        let pair = (
            inverse[bond.atom1.0 as usize],
            inverse[bond.atom2.0 as usize],
        );
        if let (Some(core), None) | (None, Some(core)) = pair {
            cuts.push(index);
            labels.push((core.0 as u16, core.0 as u16));
        }
    }
    let fragmented = fragment_single_bonds_checked(mol, &cuts, Some(&labels))
        .map_err(|e| CoreReplacementError(format!("core removal: {e}")))?;
    let mut remove = vec![false; fragmented.molecule.atom_count()];
    for &index in matching {
        remove[index.0 as usize] = true;
    }
    for a in &fragmented.attachments {
        if inverse[a.anchor.0 as usize].is_some() {
            remove[a.dummy.0 as usize] = true;
        }
    }
    let mut next = 0;
    let remap: Vec<Option<AtomIdx>> = remove
        .iter()
        .map(|&removed| {
            if removed {
                None
            } else {
                let index = AtomIdx(next);
                next += 1;
                Some(index)
            }
        })
        .collect();
    let attachments = fragmented
        .attachments
        .iter()
        .filter(|a| inverse[a.anchor.0 as usize].is_none())
        .map(|a| {
            let source_core = mol.bond(a.source_bond).other(a.anchor).unwrap();
            CoreAttachment {
                core_atom: inverse[source_core.0 as usize].unwrap(),
                source_core_atom: source_core,
                source_anchor: a.anchor,
                source_bond: a.source_bond,
                anchor: remap[a.anchor.0 as usize].unwrap(),
                dummy: remap[a.dummy.0 as usize].unwrap(),
            }
        })
        .collect();
    let mut output = fragmented.molecule;
    for index in (0..remove.len()).rev() {
        if remove[index] {
            output.remove_atom(AtomIdx(index as u32));
        }
    }
    Ok(CoreReplacement {
        molecule: output,
        original_to_output: remap[..mol.atom_count()].to_vec(),
        attachments,
    })
}

#[cfg(test)]
#[path = "core_replacement_tests.rs"]
mod tests;
