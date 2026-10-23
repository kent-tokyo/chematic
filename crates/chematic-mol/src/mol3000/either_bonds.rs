//! Explicit V3000 either single bonds also make adjacent alkene stereo unknown.
use super::{AtomIdx, BondIdx, BondOrder, MolMetadata, Molecule};
use std::collections::HashSet;

pub(super) fn adjacent_unknown_double_bonds(
    mol: &Molecule,
    metadata: &MolMetadata,
) -> HashSet<BondIdx> {
    let mut endpoints = HashSet::<AtomIdx>::new();
    for (id, properties) in &metadata.v3000_bond_properties {
        if !properties.split_whitespace().any(|token| token == "CFG=2") {
            continue;
        }
        let Some(index) = id
            .checked_sub(1)
            .filter(|&i| (i as usize) < mol.bond_count())
        else {
            continue;
        };
        let bond = mol.bond(BondIdx(index));
        if bond.order == BondOrder::Single {
            endpoints.extend([bond.atom1, bond.atom2]);
        }
    }
    if endpoints.is_empty() {
        return HashSet::new();
    }
    mol.bonds()
        .filter_map(|(index, bond)| {
            (bond.order == BondOrder::Double
                && (endpoints.contains(&bond.atom1) || endpoints.contains(&bond.atom2)))
            .then_some(index)
        })
        .collect()
}
