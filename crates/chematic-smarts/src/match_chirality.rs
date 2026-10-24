//! Tetrahedral SMARTS parity checks performed after a complete atom mapping.

use chematic_core::ez_markers::EzFact;
use chematic_core::ez_markers::{direction_is_up, read_ez_facts};
use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL};
use std::sync::OnceLock;

use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

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
    target_ez_facts: &OnceLock<Option<(Vec<EzFact>, bool)>>,
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
    }) && mapped_alkene_compatible(query, mol, target_ez_facts, |query_atom| {
        mapping
            .get(query_atom)
            .copied()
            .filter(|&target| target != UNMAPPED)
            .map(AtomIdx)
    })
}

/// Check the relative sides encoded by two directional query bonds flanking
/// a double bond. An isolated `/` or `\` has no E/Z meaning, just as in RDKit.
/// Compare sides of the mapped substituents, so a reversed atom order or a
/// marker on the other substituent cannot invert the answer.
pub(crate) fn mapped_alkene_compatible(
    query: &QueryMolecule,
    mol: &Molecule,
    target_ez_facts: &OnceLock<Option<(Vec<EzFact>, bool)>>,
    mapped: impl Fn(usize) -> Option<AtomIdx>,
) -> bool {
    for bond in &query.bonds {
        if bond.query != BondQuery::Primitive(BondPrimitive::Double) {
            continue;
        }
        let (Some((left_ref, left_up)), Some((right_ref, right_up))) = (
            query_directional_end(query, bond.atom1, bond.atom2),
            query_directional_end(query, bond.atom2, bond.atom1),
        ) else {
            continue;
        };
        let (Some(left), Some(right)) = (mapped(bond.atom1), mapped(bond.atom2)) else {
            return false;
        };
        let Some((target_double, _)) = mol.bond_between(left, right) else {
            return false;
        };
        let facts = target_ez_facts.get_or_init(|| read_ez_facts(mol));
        let Some((target_facts, _)) = facts else {
            return false;
        };
        let Some(fact) = target_facts
            .iter()
            .find(|fact| fact.double == target_double)
        else {
            return false;
        };
        let (Some(left_ref), Some(right_ref)) = (mapped(left_ref), mapped(right_ref)) else {
            return false;
        };
        let Some((left_bond, _)) = mol.bond_between(left, left_ref) else {
            return false;
        };
        let Some((right_bond, _)) = mol.bond_between(right, right_ref) else {
            return false;
        };
        let side = |end: AtomIdx, reference: BondIdx| {
            fact.ends.iter().find(|part| part.atom == end).map(|part| {
                if part.reference == reference {
                    part.up
                } else {
                    !part.up
                }
            })
        };
        if side(left, left_bond)
            .zip(side(right, right_bond))
            .is_none_or(|(target_left, target_right)| {
                (left_up == right_up) != (target_left == target_right)
            })
        {
            return false;
        }
    }
    true
}

fn query_directional_end(query: &QueryMolecule, end: usize, other: usize) -> Option<(usize, bool)> {
    query.adj[end].iter().find_map(|&(bond_idx, reference)| {
        if reference == other {
            return None;
        }
        let bond = &query.bonds[bond_idx];
        let direction = match bond.query {
            BondQuery::Primitive(BondPrimitive::Up) => BondOrder::Up,
            BondQuery::Primitive(BondPrimitive::Down) => BondOrder::Down,
            _ => return None,
        };
        Some((
            reference,
            direction_is_up(direction, AtomIdx(bond.atom1 as u32), AtomIdx(end as u32)),
        ))
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

pub(crate) fn permutation_parity(from: &[u32], to: &[u32]) -> Option<bool> {
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
