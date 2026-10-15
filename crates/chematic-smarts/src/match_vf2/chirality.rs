//! Tetrahedral query parity after all neighbours have been mapped.
use super::{
    AtomIdx, AtomPrimitive, AtomQuery, EvalCtx, Molecule, QueryMolecule, UNMAPPED,
    eval_atom_primitive,
};
use chematic_core::STEREO_H_SENTINEL;

pub(super) fn eval_atom_query_candidate(
    q: &AtomQuery,
    idx: AtomIdx,
    ctx: &EvalCtx<'_>,
) -> Option<bool> {
    match q {
        AtomQuery::Primitive(AtomPrimitive::Chirality(_)) => None,
        AtomQuery::Primitive(p) => Some(eval_atom_primitive(p, idx, ctx)),
        AtomQuery::And(a, b) => match (
            eval_atom_query_candidate(a, idx, ctx),
            eval_atom_query_candidate(b, idx, ctx),
        ) {
            (Some(false), _) | (_, Some(false)) => Some(false),
            (Some(true), Some(true)) => Some(true),
            _ => None,
        },
        AtomQuery::Or(a, b) => match (
            eval_atom_query_candidate(a, idx, ctx),
            eval_atom_query_candidate(b, idx, ctx),
        ) {
            (Some(true), _) | (_, Some(true)) => Some(true),
            (Some(false), Some(false)) => Some(false),
            _ => None,
        },
        AtomQuery::Not(a) => eval_atom_query_candidate(a, idx, ctx).map(|value| !value),
    }
}

/// Re-evaluate an atom query after the complete embedding is known, replacing
/// the coarse chirality pre-filter with neighbour-order-aware tetrahedral
/// parity. All non-chirality primitives retain their ordinary meaning.
fn eval_atom_query_mapped(
    query: &AtomQuery,
    query_atom: usize,
    target_atom: AtomIdx,
    whole_query: &QueryMolecule,
    mapping: &[u32],
    ctx: &EvalCtx<'_>,
) -> bool {
    match query {
        AtomQuery::Primitive(AtomPrimitive::Chirality(kind)) => {
            mapped_chirality_primitive_with_map(
                *kind,
                query_atom,
                target_atom,
                whole_query,
                mapping,
                ctx.mol,
                ctx.config.use_chirality,
            )
        }
        AtomQuery::Primitive(other) => eval_atom_primitive(other, target_atom, ctx),
        AtomQuery::And(a, b) => {
            eval_atom_query_mapped(a, query_atom, target_atom, whole_query, mapping, ctx)
                && eval_atom_query_mapped(b, query_atom, target_atom, whole_query, mapping, ctx)
        }
        AtomQuery::Or(a, b) => {
            eval_atom_query_mapped(a, query_atom, target_atom, whole_query, mapping, ctx)
                || eval_atom_query_mapped(b, query_atom, target_atom, whole_query, mapping, ctx)
        }
        AtomQuery::Not(inner) => {
            !eval_atom_query_mapped(inner, query_atom, target_atom, whole_query, mapping, ctx)
        }
    }
}

pub(super) fn mapped_chirality_compatible(
    query: &QueryMolecule,
    mapping: &[u32],
    ctx: &EvalCtx<'_>,
) -> bool {
    if !ctx.config.use_chirality {
        return true;
    }
    query.atoms.iter().enumerate().all(|(q, atom)| {
        mapping
            .get(q)
            .copied()
            .filter(|&target| target != UNMAPPED)
            .is_some_and(|target| {
                eval_atom_query_mapped(&atom.query, q, AtomIdx(target), query, mapping, ctx)
            })
    })
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
    use chematic_core::Chirality;

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

    // Bracket `H` occupies the slot immediately after the atom from which the
    // centre was reached (or first for a root atom), matching the SMILES
    // parser's `stereo_neighbor_order` contract. Only add the virtual slot
    // when the target has one corresponding implicit or explicit H entry.
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
        // A partial chiral query such as `[C@H]` does not provide enough
        // neighbour identity to reconstruct absolute parity. Preserve the
        // historical raw-tag behavior for that explicitly bounded case.
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
        // An H constraint hidden behind OR does not establish a unique
        // virtual neighbour for a tetrahedral parity check.
        AtomQuery::Not(inner) => matches!(
            inner.as_ref(),
            AtomQuery::Primitive(AtomPrimitive::HCount(0))
        ),
        AtomQuery::Primitive(_) | AtomQuery::Or(_, _) => false,
    }
}

/// `Some(true)` for an even permutation, `Some(false)` for an odd one.
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

/// Cheap chirality pre-filter used while the mapping is incomplete.
///
/// `@`/`@@` is relative to neighbour encounter order, so the exact sense can
/// only be checked once all query neighbours have been mapped. During search
/// this therefore requires merely that the target carries a stereo tag;
/// [`mapped_chirality_compatible`] performs the parity-aware final check.
pub(super) fn eval_chirality(idx: AtomIdx, ctx: &EvalCtx<'_>, kind: u8) -> bool {
    if !ctx.config.use_chirality {
        return true;
    }
    use chematic_core::Chirality;
    let c = ctx.mol.atom(idx).chirality;
    let _ = kind;
    c != Chirality::None
}
