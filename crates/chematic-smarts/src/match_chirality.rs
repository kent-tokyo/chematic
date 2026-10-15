//! Tetrahedral SMARTS parity checks performed after a complete atom mapping.

use chematic_core::{AtomIdx, Chirality, Molecule, STEREO_H_SENTINEL};

use crate::query::{AtomPrimitive, AtomQuery, QueryMolecule};

const UNMAPPED: u32 = u32::MAX;

/// Evaluate every non-stereo term while leaving stereo as unknown until the
/// complete atom mapping makes neighbour-order parity available.
pub(crate) fn candidate_query_result(
    query: &AtomQuery,
    mut eval_primitive: impl FnMut(&AtomPrimitive) -> bool,
) -> Option<bool> {
    fn eval(
        query: &AtomQuery,
        eval_primitive: &mut impl FnMut(&AtomPrimitive) -> bool,
    ) -> Option<bool> {
        match query {
            AtomQuery::Primitive(AtomPrimitive::Chirality(_)) => None,
            AtomQuery::Primitive(primitive) => Some(eval_primitive(primitive)),
            AtomQuery::And(a, b) => match (eval(a, eval_primitive), eval(b, eval_primitive)) {
                (Some(false), _) | (_, Some(false)) => Some(false),
                (Some(true), Some(true)) => Some(true),
                _ => None,
            },
            AtomQuery::Or(a, b) => match (eval(a, eval_primitive), eval(b, eval_primitive)) {
                (Some(true), _) | (_, Some(true)) => Some(true),
                (Some(false), Some(false)) => Some(false),
                _ => None,
            },
            AtomQuery::Not(inner) => eval(inner, eval_primitive).map(|value| !value),
        }
    }
    eval(query, &mut eval_primitive)
}

pub(crate) fn mapped_chirality_compatible(
    query: &QueryMolecule,
    mapping: &[u32],
    mol: &Molecule,
    use_chirality: bool,
    mut eval_primitive: impl FnMut(&AtomPrimitive, AtomIdx) -> bool,
) -> bool {
    if !use_chirality {
        return true;
    }
    query.atoms.iter().enumerate().all(|(query_atom, atom)| {
        mapping
            .get(query_atom)
            .copied()
            .filter(|&target| target != UNMAPPED)
            .is_some_and(|target| {
                eval_mapped_query(
                    &atom.query,
                    query_atom,
                    AtomIdx(target),
                    query,
                    mapping,
                    mol,
                    use_chirality,
                    &mut eval_primitive,
                )
            })
    })
}

#[allow(clippy::too_many_arguments)]
fn eval_mapped_query(
    atom_query: &AtomQuery,
    query_atom: usize,
    target_atom: AtomIdx,
    query: &QueryMolecule,
    mapping: &[u32],
    mol: &Molecule,
    use_chirality: bool,
    eval_primitive: &mut impl FnMut(&AtomPrimitive, AtomIdx) -> bool,
) -> bool {
    match atom_query {
        AtomQuery::Primitive(AtomPrimitive::Chirality(kind)) => {
            mapped_chirality_primitive_with_map(
                *kind,
                query_atom,
                target_atom,
                query,
                mapping,
                mol,
                use_chirality,
            )
        }
        AtomQuery::Primitive(other) => eval_primitive(other, target_atom),
        AtomQuery::And(a, b) => {
            eval_mapped_query(
                a,
                query_atom,
                target_atom,
                query,
                mapping,
                mol,
                use_chirality,
                eval_primitive,
            ) && eval_mapped_query(
                b,
                query_atom,
                target_atom,
                query,
                mapping,
                mol,
                use_chirality,
                eval_primitive,
            )
        }
        AtomQuery::Or(a, b) => {
            eval_mapped_query(
                a,
                query_atom,
                target_atom,
                query,
                mapping,
                mol,
                use_chirality,
                eval_primitive,
            ) || eval_mapped_query(
                b,
                query_atom,
                target_atom,
                query,
                mapping,
                mol,
                use_chirality,
                eval_primitive,
            )
        }
        AtomQuery::Not(inner) => !eval_mapped_query(
            inner,
            query_atom,
            target_atom,
            query,
            mapping,
            mol,
            use_chirality,
            eval_primitive,
        ),
    }
}

pub(crate) fn mapped_chirality_primitive_with_map(
    kind: u8,
    query_atom: usize,
    target_atom: AtomIdx,
    query: &QueryMolecule,
    mapping: &[u32],
    mol: &Molecule,
    use_chirality: bool,
) -> bool {
    if !use_chirality {
        return true;
    }
    let target_chirality = mol.atom(target_atom).chirality;
    if target_chirality == Chirality::None {
        return false;
    }
    if !matches!(kind, 1 | 2) {
        return true;
    }

    let raw_matches = matches!(
        (kind, target_chirality),
        (1, Chirality::CounterClockwise) | (2, Chirality::Clockwise)
    );
    let Some(target_order) = mol.stereo_neighbor_order(target_atom) else {
        return raw_matches;
    };

    let mut mapped_order: Vec<u32> = query.adj[query_atom]
        .iter()
        .filter_map(|&(_, neighbor)| {
            mapping
                .get(neighbor)
                .copied()
                .filter(|&target| target != UNMAPPED)
        })
        .collect();
    if mapped_order.len() != query.adj[query_atom].len() {
        return raw_matches;
    }

    if target_order.len() == mapped_order.len() + 1
        && query_has_positive_h_count(&query.atoms[query_atom].query)
    {
        let missing: Vec<_> = target_order
            .iter()
            .copied()
            .filter(|value| !mapped_order.contains(value))
            .collect();
        if let [hydrogen] = missing.as_slice()
            && (*hydrogen == STEREO_H_SENTINEL
                || mol.atom(AtomIdx(*hydrogen)).element == chematic_core::Element::H)
        {
            let from_position = usize::from(
                query.adj[query_atom]
                    .first()
                    .is_some_and(|&(_, neighbor)| neighbor < query_atom),
            );
            mapped_order.insert(from_position, *hydrogen);
        }
    }

    let Some(even) = permutation_parity(&mapped_order, target_order) else {
        return raw_matches;
    };
    let query_chirality = if kind == 1 {
        Chirality::CounterClockwise
    } else {
        Chirality::Clockwise
    };
    (query_chirality == target_chirality) == even
}

fn query_has_positive_h_count(query: &AtomQuery) -> bool {
    match query {
        AtomQuery::Primitive(AtomPrimitive::HCount(count)) => *count > 0,
        AtomQuery::And(a, b) => query_has_positive_h_count(a) || query_has_positive_h_count(b),
        AtomQuery::Not(inner) => matches!(
            inner.as_ref(),
            AtomQuery::Primitive(AtomPrimitive::HCount(0))
        ),
        AtomQuery::Primitive(_) | AtomQuery::Or(_, _) => false,
    }
}

fn permutation_parity(from: &[u32], to: &[u32]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut positions = Vec::with_capacity(to.len());
    for value in to {
        let position = from.iter().position(|candidate| candidate == value)?;
        if positions.contains(&position) {
            return None;
        }
        positions.push(position);
    }
    let inversions = positions
        .iter()
        .enumerate()
        .map(|(i, value)| {
            positions[(i + 1)..]
                .iter()
                .filter(|next| value > next)
                .count()
        })
        .sum::<usize>();
    Some(inversions.is_multiple_of(2))
}
