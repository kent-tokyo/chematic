//! Project tetrahedral stereo onto the selected common bonds.

use chematic_core::{AtomIdx, BondIdx, Chirality, Molecule, STEREO_H_SENTINEL, StereoGroupKind};

use super::{
    Instant, McsConfig, McsState, PartialMapping, better_mapping, mapped_rows_share_common_bond,
    record_mapping, reserve_search_node,
};
use crate::match_chirality::{permutation_parity, tetrahedral_frame_parity};
use crate::query::{AtomPrimitive, AtomQuery, QueryMolecule};

pub(super) enum StereoConflict {
    IncompatibleTag,
    Orientation([BondIdx; 3]),
}

pub(super) fn mapped_stereo_conflict(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    config: &McsConfig,
) -> Option<StereoConflict> {
    mapping
        .query_to_mol
        .iter()
        .enumerate()
        .find_map(|(q, row)| {
            let source_atom = row[0];
            let source_tag = mols[0].atom(source_atom).chirality;
            if !source_tag.is_tetrahedral() {
                return None;
            }
            let source_allene = mols[0].is_allene_stereo(source_atom);
            let tag_mismatch = mols.iter().enumerate().skip(1).any(|(m, target)| {
                let target_tag = target.atom(row[m]).chirality;
                !target_tag.is_tetrahedral()
                    || source_allene != target.is_allene_stereo(row[m])
                    || (source_allene && source_tag != target_tag)
            });
            if tag_mismatch {
                return Some(StereoConflict::IncompatibleTag);
            }
            if source_allene {
                return None;
            }

            let ligands: Vec<_> = mols[0]
                .neighbors(source_atom)
                .filter_map(|(neighbor, bond)| Some((mapping.query_idx_of(0, neighbor)?, bond)))
                .filter(|&(other, _)| {
                    mapped_rows_share_common_bond(mols, mapping, q, other, config)
                })
                .collect();
            let [a, b, c, ..] = ligands.as_slice() else {
                return None;
            };
            let side = |m: usize| {
                let order = mols[m].stereo_neighbor_order(row[m])?;
                let projected: Vec<_> = ligands
                    .iter()
                    .map(|&(ligand, _)| mapping.query_to_mol[ligand][m].0)
                    .collect();
                Some(
                    (mols[m].atom(row[m]).chirality == Chirality::Clockwise)
                        == tetrahedral_frame_parity(&projected, order)?,
                )
            };
            let source_side = side(0);
            let inverted =
                mols.iter()
                    .enumerate()
                    .skip(1)
                    .any(|(m, target)| match (source_side, side(m)) {
                        (Some(source), Some(target)) => source != target,
                        _ => source_tag != target.atom(row[m]).chirality,
                    });
            // Any three common bonds determine the same parity agreement. For
            // an inverted tetrahedron, at least one of these must be excluded.
            inverted.then_some(StereoConflict::Orientation([a.1, b.1, c.1]))
        })
}

pub(super) fn consider_stereo_subgraphs(
    state: &mut McsState<'_>,
    mapping: &PartialMapping,
    initial: [BondIdx; 3],
) {
    if mapping.bond_count < mapping.size || !better_mapping(state, mapping) {
        // A connected tree cannot omit a bond and keep its vertex set.
        return;
    }
    let mut candidate = mapping.clone();
    let mut visited = rustc_hash::FxHashSet::default();
    let mut pending: Vec<Vec<BondIdx>> = initial.into_iter().map(|bond| vec![bond]).collect();
    while let Some(omitted) = pending.pop() {
        if !visited.insert(omitted.clone()) {
            continue;
        }
        if !reserve_search_node(state) {
            break;
        }
        if let Some(deadline) = state.deadline
            && Instant::now() >= deadline
        {
            state.timed_out = true;
            break;
        }
        candidate.bond_count = mapping.bond_count - omitted.len();
        candidate.omitted_bonds = omitted;
        if !better_mapping(state, &candidate) || !has_connected_common_bonds(state, &candidate) {
            continue;
        }
        match mapped_stereo_conflict(state.mols, &candidate, state.config) {
            None => {
                record_mapping(state, &candidate);
                if state.proven_optimal {
                    break;
                }
            }
            Some(StereoConflict::Orientation(bonds)) if candidate.bond_count >= candidate.size => {
                for bond in bonds {
                    let mut next = candidate.omitted_bonds.clone();
                    next.insert(
                        next.binary_search(&bond)
                            .expect_err("bond is still present"),
                        bond,
                    );
                    pending.push(next);
                }
            }
            _ => {}
        }
    }
}

fn has_connected_common_bonds(state: &McsState<'_>, mapping: &PartialMapping) -> bool {
    if mapping.bond_count + 1 < mapping.size {
        return false;
    }
    let mut seen = vec![false; mapping.size];
    let mut pending = vec![0];
    seen[0] = true;
    let mut count = 1;
    while let Some(q) = pending.pop() {
        let source = mapping.query_to_mol[q][0];
        for (neighbor, _) in state.mols[0].neighbors(source) {
            if let Some(other) = mapping.query_idx_of(0, neighbor)
                && !seen[other]
                && mapped_rows_share_common_bond(state.mols, mapping, q, other, state.config)
            {
                seen[other] = true;
                count += 1;
                pending.push(other);
            }
        }
    }
    count == mapping.size
}

pub(super) fn constrain_query_chirality(
    query: &mut QueryMolecule,
    mols: &[&Molecule],
    mapped: impl Fn(usize, usize) -> AtomIdx,
) {
    for (q, neighbors) in query.adj.iter().enumerate() {
        if mols[0].atom(mapped(q, 0)).chirality == Chirality::None {
            continue;
        }

        let has_bracket_h = neighbors.len() == 3
            && mols.iter().enumerate().all(|(m, mol)| {
                let center = mapped(q, m);
                mol.atom(center).hydrogen_count == Some(1)
                    && mol
                        .stereo_neighbor_order(center)
                        .is_some_and(|order| order.contains(&STEREO_H_SENTINEL))
            });
        let hydrogen_position = usize::from(neighbors.first().is_some_and(|&(_, n)| n < q));
        let mut common_kind = None;
        for (m, mol) in mols.iter().enumerate() {
            let center = mapped(q, m);
            let projections = neighbors.iter().map(|&(_, n)| mapped(n, m).0);
            let kind = projected_tetrahedral_kind(
                mol,
                center,
                projections,
                has_bracket_h,
                hydrogen_position,
            );
            if kind.is_none() || common_kind.is_some_and(|first| Some(first) != kind) {
                common_kind = None;
                break;
            }
            common_kind = kind;
        }

        let atom = &mut query.atoms[q].query;
        let constraint = |predicate| AtomQuery::Primitive(predicate);
        let and = |left, right| AtomQuery::And(Box::new(left), Box::new(right));
        // Keep unsupported stereo explicit so the SMARTS writer rejects it.
        *atom = and(
            atom.clone(),
            constraint(AtomPrimitive::Chirality(common_kind.unwrap_or(3))),
        );
        if common_kind.is_some() && has_bracket_h {
            *atom = and(atom.clone(), constraint(AtomPrimitive::HCount(1)));
        }
    }
}

fn projected_tetrahedral_kind(
    mol: &Molecule,
    center: AtomIdx,
    neighbors: impl Iterator<Item = u32>,
    has_bracket_h: bool,
    hydrogen_position: usize,
) -> Option<u8> {
    let atom = mol.atom(center);
    if !atom.chirality.is_tetrahedral()
        || mol.is_allene_stereo(center)
        || mol.stereo_groups().iter().any(|group| {
            group.kind != StereoGroupKind::Absolute && group.atom_indices.contains(&center)
        })
    {
        return None;
    }
    let order = mol.stereo_neighbor_order(center)?;
    if order.len() != 4 {
        return None;
    }
    let mut projected: Vec<u32> = neighbors.collect();
    if projected.len() < 3 {
        // SMARTS @/@@ with fewer than three query neighbors encodes the
        // presence of tetrahedral stereo and imposes no orientation.
        return Some(1);
    }
    let even = if has_bracket_h {
        projected.insert(hydrogen_position, STEREO_H_SENTINEL);
        permutation_parity(order, &projected)?
    } else {
        tetrahedral_frame_parity(&projected, order)?
    };
    let clockwise = (atom.chirality == Chirality::Clockwise) == even;
    Some(if clockwise { 2 } else { 1 })
}
