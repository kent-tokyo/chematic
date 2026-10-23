//! Necessary atom-count and element screens before VF2 state allocation.

use super::*;

/// Atomic number every target atom matching `q` must have, when `q` pins one
/// down (an element primitive at the top of an `And` chain, or both arms of an
/// `Or` pinning the same element). `None` means "no single element implied".
pub(super) fn implied_atomic_number(q: &AtomQuery) -> Option<u8> {
    match q {
        AtomQuery::Primitive(AtomPrimitive::AtomicNum(n)) => Some(*n),
        AtomQuery::Primitive(AtomPrimitive::Symbol(sym)) => {
            chematic_core::Element::from_symbol(sym)
                .filter(|e| e.symbol() == sym.as_str())
                .map(|e| e.atomic_number())
        }
        AtomQuery::And(a, b) => implied_atomic_number(a).or_else(|| implied_atomic_number(b)),
        AtomQuery::Or(a, b) => match (implied_atomic_number(a), implied_atomic_number(b)) {
            (Some(x), Some(y)) if x == y => Some(x),
            _ => None,
        },
        _ => None,
    }
}

/// Necessary condition for any embedding: the target has at least as many
/// atoms of each element as the query has atoms pinned to that element
/// (embeddings are injective). `false` proves there is no match.
pub(super) fn element_counts_allow_match(query: &QueryMolecule, mol: &Molecule) -> bool {
    let mut need = [0u16; 128];
    let mut any = false;
    for atom in &query.atoms {
        if let Some(z) = implied_atomic_number(&atom.query)
            && (z as usize) < need.len()
        {
            need[z as usize] += 1;
            any = true;
        }
    }
    if !any {
        return true;
    }
    for (_, atom) in mol.atoms() {
        let z = atom.element.atomic_number() as usize;
        if z < need.len() && need[z] > 0 {
            need[z] -= 1;
        }
    }
    need.iter().all(|&n| n == 0)
}

/// The screens [`has_match_with_config`] applies before searching: `false`
/// means no match. They read only atom counts and elements.
pub(crate) fn passes_size_and_element_screen(
    query: &QueryMolecule,
    mol: &Molecule,
    config: &MatchConfig,
) -> bool {
    !(query.atoms.is_empty()
        || query.atoms.len() > mol.atom_count()
        || (config.max_visit_budget.is_none() && !element_counts_allow_match(query, mol)))
}
