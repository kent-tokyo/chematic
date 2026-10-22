//! Frontier expansion for the connected-growth MCS search.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondOrder, Molecule};

use super::{
    CartesianProduct, McsConfig, McsState, PartialMapping, atoms_compatible, bonds_compatible, grow,
};

/// Explore every compatible tuple for one frontier atom.
pub(super) fn extend_frontier_atom(
    state: &mut McsState<'_>,
    mapping: &mut PartialMapping,
    n0: AtomIdx,
) {
    let mol0 = state.mols[0];
    let mapped_parents: Vec<(usize, BondOrder)> = mol0
        .neighbors(n0)
        .filter_map(|(neighbor, bond_idx)| {
            mapping
                .query_idx_of(0, neighbor)
                .map(|query_idx| (query_idx, mol0.bond(bond_idx).order))
        })
        .collect();
    let mut seen_tuples = HashSet::new();

    for (parent_q, parent_bond_order) in mapped_parents {
        let candidates_per_mol = build_candidates(state, mapping, n0, parent_q, parent_bond_order);
        explore_candidate_tuples(state, mapping, n0, &candidates_per_mol, &mut seen_tuples);
        if state.timed_out || state.resource_limited {
            return;
        }
    }
}

fn explore_candidate_tuples(
    state: &mut McsState<'_>,
    mapping: &mut PartialMapping,
    n0: AtomIdx,
    candidates_per_mol: &[Vec<AtomIdx>],
    seen_tuples: &mut HashSet<Vec<AtomIdx>>,
) {
    for tuple in CartesianProduct::new(candidates_per_mol) {
        if state.timed_out || state.resource_limited {
            return;
        }
        if !seen_tuples.insert(tuple.clone()) {
            continue;
        }

        let mut atoms = Vec::with_capacity(state.mols.len());
        atoms.push(n0);
        atoms.extend_from_slice(&tuple);

        let extra_bonds = count_common_new_bonds(state.mols, mapping, n0, &tuple, state.config);
        debug_assert!(extra_bonds > 0);

        mapping.extend(&atoms, extra_bonds);
        grow(state, mapping);
        mapping.retract(&atoms, extra_bonds);
    }
}

/// For each molecule after the first, collect atoms that can extend the mapping
/// through the selected mapped parent.
fn build_candidates(
    state: &McsState<'_>,
    mapping: &PartialMapping,
    n0: AtomIdx,
    parent_q: usize,
    parent_bond_order: BondOrder,
) -> Vec<Vec<AtomIdx>> {
    let mols = state.mols;
    let config = state.config;
    let ring_sets = &state.ring_sets;
    let atom0 = mols[0].atom(n0);
    let n0_in_ring = !ring_sets.is_empty() && ring_sets[0].contains_atom(n0);
    let mut result = Vec::with_capacity(mols.len() - 1);

    for mi in 1..mols.len() {
        let parent_i = mapping.query_to_mol[parent_q][mi];
        let candidates = mols[mi]
            .atoms()
            .filter(|(ai, atom_i)| {
                !mapping.is_mapped(mi, *ai)
                    && atoms_compatible(
                        atom0,
                        atom_i,
                        &config.atom_compare,
                        state.rdkit_heavy,
                        config.match_chiral_tag,
                        config.match_charge,
                        config.match_isotope,
                    )
                    && (!config.ring_matches_ring_only
                        || n0_in_ring == ring_sets[mi].contains_atom(*ai))
                    && mols[mi]
                        .bond_between(*ai, parent_i)
                        .is_some_and(|(_, bond)| {
                            !config.match_bonds
                                || bonds_compatible(
                                    parent_bond_order,
                                    bond.order,
                                    &config.bond_compare,
                                )
                        })
            })
            .map(|(atom_idx, _)| atom_idx)
            .collect();
        result.push(candidates);
    }

    result
}

/// Count edges from a proposed new atom that are common to every molecule.
fn count_common_new_bonds(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    n0: AtomIdx,
    tuple: &[AtomIdx],
    config: &McsConfig,
) -> usize {
    mapping
        .query_to_mol
        .iter()
        .filter(|row| {
            let Some((_, bond0)) = mols[0].bond_between(n0, row[0]) else {
                return false;
            };
            tuple.iter().enumerate().all(|(offset, &candidate)| {
                let mi = offset + 1;
                mols[mi]
                    .bond_between(candidate, row[mi])
                    .is_some_and(|(_, bond)| {
                        !config.match_bonds
                            || bonds_compatible(bond0.order, bond.order, &config.bond_compare)
                    })
            })
        })
        .count()
}

pub(super) fn mapped_rows_share_common_bond(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    q1: usize,
    q2: usize,
    config: &McsConfig,
) -> bool {
    let row1 = &mapping.query_to_mol[q1];
    let row2 = &mapping.query_to_mol[q2];
    let Some((_, bond0)) = mols[0].bond_between(row1[0], row2[0]) else {
        return false;
    };
    mols.iter().enumerate().skip(1).all(|(mi, mol)| {
        mol.bond_between(row1[mi], row2[mi])
            .is_some_and(|(_, bond)| {
                !config.match_bonds
                    || bonds_compatible(bond0.order, bond.order, &config.bond_compare)
            })
    })
}
