//! Lossless SMARTS serialization for the query subset produced by MCS search.

use std::fmt;

use crate::match_chirality::{permutation_parity, query_has_positive_h_count};
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

/// An MCS query cannot be represented faithfully as SMARTS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum McsSmartsError {
    InvalidGraph,
    UnsupportedAtomQuery,
    UnsupportedChirality,
    UnsupportedBondQuery,
    TooManyRingClosures,
}

impl fmt::Display for McsSmartsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if *self == Self::UnsupportedChirality {
            return f.write_str("MCS SMARTS can retain only complete absolute tetrahedral frames");
        }
        write!(f, "{self:?}")
    }
}

impl std::error::Error for McsSmartsError {}

fn atom_alternatives(query: &AtomQuery, out: &mut String) -> Result<(), McsSmartsError> {
    match query {
        AtomQuery::Primitive(AtomPrimitive::AtomicNum(number)) => {
            out.push('#');
            out.push_str(&number.to_string());
        }
        AtomQuery::Or(left, right) => {
            atom_alternatives(left, out)?;
            out.push(',');
            atom_alternatives(right, out)?;
        }
        _ => return Err(McsSmartsError::UnsupportedAtomQuery),
    }
    Ok(())
}

fn atom_expression(
    query: &AtomQuery,
    out: &mut String,
    chirality_swapped: &impl Fn() -> Result<bool, McsSmartsError>,
) -> Result<(), McsSmartsError> {
    let AtomQuery::And(left, right) = query else {
        return atom_alternatives(query, out);
    };
    atom_expression(left, out, chirality_swapped)?;
    // ';' binds more loosely than ',', so the predicate applies to every
    // observed element alternative rather than just to its last branch.
    match right.as_ref() {
        AtomQuery::Primitive(AtomPrimitive::Charge(charge)) => {
            out.push_str(&format!(";{charge:+}"));
        }
        AtomQuery::Primitive(AtomPrimitive::Isotope(mass)) => {
            out.push_str(&format!(";{mass}*"));
        }
        AtomQuery::Primitive(AtomPrimitive::Chirality(kind @ (1 | 2))) => {
            let clockwise = (*kind == 2) ^ chirality_swapped()?;
            out.push_str(if clockwise { ";@@" } else { ";@" });
        }
        AtomQuery::Primitive(AtomPrimitive::Chirality(_)) => {
            return Err(McsSmartsError::UnsupportedChirality);
        }
        AtomQuery::Primitive(AtomPrimitive::HCount(1)) => out.push_str(";H1"),
        _ => return Err(McsSmartsError::UnsupportedAtomQuery),
    }
    Ok(())
}

fn bond_expression(query: &BondQuery, out: &mut String) -> Result<(), McsSmartsError> {
    match query {
        BondQuery::Primitive(primitive) => out.push(match primitive {
            BondPrimitive::Single => '-',
            BondPrimitive::Double => '=',
            BondPrimitive::Triple => '#',
            BondPrimitive::Quadruple => '$',
            BondPrimitive::Aromatic => ':',
            BondPrimitive::Any => '~',
            _ => return Err(McsSmartsError::UnsupportedBondQuery),
        }),
        BondQuery::Or(left, right) => {
            bond_expression(left, out)?;
            out.push(',');
            bond_expression(right, out)?;
        }
        BondQuery::Any => {} // An unspecified SMARTS bond means single or aromatic.
        _ => return Err(McsSmartsError::UnsupportedBondQuery),
    }
    Ok(())
}

fn ring_label(number: usize, out: &mut String) {
    if number >= 10 {
        out.push('%');
    }
    out.push_str(&number.to_string());
}

/// Serialize an MCS `QueryMolecule` without collapsing atom or bond alternatives.
///
/// This accepts the `AtomicNum`/`Or` and bond-order/`Or` query subset emitted by
/// MCS, including trailing charge, isotope and complete absolute tetrahedral
/// constraints. Other predicates and malformed graphs return an error. An
/// empty query serializes to `""`.
pub fn mcs_query_to_smarts(query: &QueryMolecule) -> Result<String, McsSmartsError> {
    validate_graph(query)?;
    let layout = forest_layout(query)?;
    render_layout(query, &layout)
}

fn validate_graph(query: &QueryMolecule) -> Result<(), McsSmartsError> {
    let atom_count = query.atoms.len();
    if query.adj.len() != atom_count {
        return Err(McsSmartsError::InvalidGraph);
    }
    let mut adjacency_count = vec![0u8; query.bonds.len()];
    for (atom, neighbors) in query.adj.iter().enumerate() {
        for &(bond_idx, neighbor) in neighbors {
            let bond = query
                .bonds
                .get(bond_idx)
                .ok_or(McsSmartsError::InvalidGraph)?;
            if neighbor >= atom_count
                || atom == neighbor
                || !((bond.atom1 == atom && bond.atom2 == neighbor)
                    || (bond.atom2 == atom && bond.atom1 == neighbor))
            {
                return Err(McsSmartsError::InvalidGraph);
            }
            adjacency_count[bond_idx] = adjacency_count[bond_idx].saturating_add(1);
        }
    }
    if adjacency_count.iter().any(|&count| count != 2)
        || query.bonds.iter().any(|bond| {
            bond.atom1 >= atom_count || bond.atom2 >= atom_count || bond.atom1 == bond.atom2
        })
    {
        return Err(McsSmartsError::InvalidGraph);
    }

    Ok(())
}

struct Layout {
    parent: Vec<Option<usize>>,
    children: Vec<Vec<(usize, usize)>>,
    ring_marks: Vec<Vec<(usize, usize, bool)>>,
    roots: Vec<usize>,
}

fn forest_layout(query: &QueryMolecule) -> Result<Layout, McsSmartsError> {
    let atom_count = query.atoms.len();
    let mut parent = vec![None; atom_count];
    let mut children: Vec<Vec<(usize, usize)>> = vec![Vec::new(); atom_count];
    let mut ring_marks: Vec<Vec<(usize, usize, bool)>> = vec![Vec::new(); atom_count];
    let mut visited = vec![false; atom_count];
    let mut tree_edge = vec![false; query.bonds.len()];
    let mut ring_edge = vec![false; query.bonds.len()];
    let mut roots = Vec::new();
    let mut next_ring = 1;

    // With a preceding atom, a bracket H has the same parity convention
    // in the native and RDKit SMARTS readers.
    let has_bracket_h = |i: usize| query_has_positive_h_count(&query.atoms[i].query);
    let root_candidates = (0..atom_count)
        .filter(|&i| !has_bracket_h(i))
        .chain((0..atom_count).filter(|&i| has_bracket_h(i)));
    for root in root_candidates {
        if visited[root] {
            continue;
        }
        roots.push(root);
        visited[root] = true;
        let mut stack = vec![(root, 0usize)];
        while let Some((atom, cursor)) = stack.last_mut() {
            let current = *atom;
            if *cursor == query.adj[current].len() {
                stack.pop();
                continue;
            }
            let (bond_idx, neighbor) = query.adj[current][*cursor];
            *cursor += 1;
            if tree_edge[bond_idx] || ring_edge[bond_idx] {
                continue;
            }
            if !visited[neighbor] {
                visited[neighbor] = true;
                tree_edge[bond_idx] = true;
                parent[neighbor] = Some(current);
                children[current].push((bond_idx, neighbor));
                stack.push((neighbor, 0));
            } else {
                if next_ring > 99 {
                    return Err(McsSmartsError::TooManyRingClosures);
                }
                ring_edge[bond_idx] = true;
                ring_marks[neighbor].push((bond_idx, next_ring, true));
                ring_marks[current].push((bond_idx, next_ring, false));
                next_ring += 1;
            }
        }
    }

    Ok(Layout {
        parent,
        children,
        ring_marks,
        roots,
    })
}

fn chirality_swapped(
    query: &QueryMolecule,
    layout: &Layout,
    atom: usize,
) -> Result<bool, McsSmartsError> {
    let stored = &query.adj[atom];
    let bracket_h = query_has_positive_h_count(&query.atoms[atom].query);
    if stored.len() + usize::from(bracket_h) != 4 || (bracket_h && layout.parent[atom].is_none()) {
        return Err(McsSmartsError::UnsupportedChirality);
    }
    let stored_order: Vec<_> = stored.iter().map(|&(_, n)| n as u32).collect();
    let written_order: Vec<_> = layout.parent[atom]
        .into_iter()
        .chain(layout.ring_marks[atom].iter().map(|&(bond, _, _)| {
            let bond = &query.bonds[bond];
            if bond.atom1 == atom {
                bond.atom2
            } else {
                bond.atom1
            }
        }))
        .chain(layout.children[atom].iter().map(|&(_, n)| n))
        .map(|n| n as u32)
        .collect();
    let even =
        permutation_parity(&stored_order, &written_order).ok_or(McsSmartsError::InvalidGraph)?;
    let stored_has_predecessor = stored.first().is_some_and(|&(_, n)| n < atom);
    // A bracket H is first at a root and second after the preceding atom.
    let hydrogen_swapped = bracket_h && stored_has_predecessor != layout.parent[atom].is_some();
    Ok(!even ^ hydrogen_swapped)
}

fn render_layout(query: &QueryMolecule, layout: &Layout) -> Result<String, McsSmartsError> {
    enum Action {
        Atom(usize),
        Bond(usize),
        Open,
        Close,
    }
    let mut out = String::new();
    for &root in &layout.roots {
        if !out.is_empty() {
            out.push('.');
        }
        let mut actions = vec![Action::Atom(root)];
        while let Some(action) = actions.pop() {
            match action {
                Action::Open => out.push('('),
                Action::Close => out.push(')'),
                Action::Bond(bond_idx) => {
                    bond_expression(&query.bonds[bond_idx].query, &mut out)?;
                }
                Action::Atom(atom_idx) => {
                    let atom = &query.atoms[atom_idx];
                    if atom.atom_map.is_some() {
                        return Err(McsSmartsError::UnsupportedAtomQuery);
                    }
                    out.push('[');
                    atom_expression(&atom.query, &mut out, &|| {
                        chirality_swapped(query, layout, atom_idx)
                    })?;
                    out.push(']');
                    for &(bond_idx, ring, opening) in &layout.ring_marks[atom_idx] {
                        if opening {
                            bond_expression(&query.bonds[bond_idx].query, &mut out)?;
                        }
                        ring_label(ring, &mut out);
                    }
                    let last = layout.children[atom_idx].len().saturating_sub(1);
                    for (index, &(bond_idx, child)) in
                        layout.children[atom_idx].iter().enumerate().rev()
                    {
                        if index != last {
                            actions.push(Action::Close);
                        }
                        actions.push(Action::Atom(child));
                        actions.push(Action::Bond(bond_idx));
                        if index != last {
                            actions.push(Action::Open);
                        }
                    }
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{MatchConfig, McsConfig, find_matches, find_mcs_with_config, parse_smarts};

    #[test]
    fn complete_tetrahedral_query_round_trips_a_stereocenter_on_a_ring() {
        let config = McsConfig {
            match_chiral_tag: true,
            ..McsConfig::default()
        };
        let match_config = MatchConfig {
            use_chirality: true,
            ..MatchConfig::default()
        };
        let parse = |smiles| chematic_smiles::parse(smiles).unwrap();
        let left = parse("[C@H]1(N)CCOC1");
        let right = parse("N[C@H]1CCOC1");
        let inverted = parse("N[C@@H]1CCOC1");
        let mcs = find_mcs_with_config(&[&left, &right], &config);
        assert_eq!(mcs.atom_count(), left.atom_count());
        let smarts = mcs_query_to_smarts(&mcs).unwrap();
        let roundtrip = parse_smarts(&smarts).unwrap();
        for query in [&mcs, &roundtrip] {
            assert!(
                crate::has_match_with_config(query, &left, &match_config),
                "{smarts}"
            );
            assert!(
                crate::has_match_with_config(query, &right, &match_config),
                "{smarts}"
            );
            assert!(
                !crate::has_match_with_config(query, &inverted, &match_config),
                "{smarts}"
            );
        }
    }

    #[test]
    fn incomplete_stereocenters_do_not_serialize_as_complete_frames() {
        let config = McsConfig {
            match_chiral_tag: true,
            ..McsConfig::default()
        };
        let left = chematic_smiles::parse("C[C@H](O)Cl").unwrap();
        let right = chematic_smiles::parse("C[C@H](I)Cl").unwrap();
        let mcs = find_mcs_with_config(&[&left, &right], &config);
        assert_eq!(mcs.atom_count(), 3);
        assert_eq!(
            mcs_query_to_smarts(&mcs),
            Err(McsSmartsError::UnsupportedChirality)
        );
    }

    #[test]
    fn dummy_atom_in_mcs_serializes_with_its_effective_atomic_number() {
        let dummy = chematic_smiles::parse("[*]O").unwrap();
        let carbon = chematic_smiles::parse("CO").unwrap();
        let config = McsConfig::default();
        let common = find_mcs_with_config(&[&dummy, &carbon], &config);
        assert_eq!(common.atom_count(), 1);
        let common = find_mcs_with_config(&[&dummy, &dummy], &config);
        let smarts = mcs_query_to_smarts(&common).unwrap();
        assert!(smarts.contains("[#0]"), "{smarts}");
        let round_trip = parse_smarts(&smarts).unwrap();
        assert_eq!(find_matches(&round_trip, &dummy).len(), 1);
        assert!(find_matches(&round_trip, &carbon).is_empty());
    }

    #[test]
    fn round_trip_relaxed_atom_and_bond_alternatives() {
        for (left, right, config) in [
            (
                "CO",
                "NO",
                McsConfig {
                    atom_compare: crate::AtomCompare::Any,
                    ..McsConfig::default()
                },
            ),
            (
                "C=N",
                "C#N",
                McsConfig {
                    bond_compare: crate::BondCompare::Any,
                    ..McsConfig::default()
                },
            ),
            (
                "C1CC1",
                "N1CC1",
                McsConfig {
                    atom_compare: crate::AtomCompare::Any,
                    ..McsConfig::default()
                },
            ),
        ] {
            let left = chematic_smiles::parse(left).unwrap();
            let right = chematic_smiles::parse(right).unwrap();
            let original = find_mcs_with_config(&[&left, &right], &config);
            let smarts = mcs_query_to_smarts(&original).unwrap();
            assert!(smarts.contains(','), "alternatives were lost: {smarts}");
            let reparsed = parse_smarts(&smarts).unwrap();
            assert_eq!(reparsed.atom_count(), original.atom_count(), "{smarts}");
            assert_eq!(reparsed.bonds.len(), original.bonds.len(), "{smarts}");
            assert!(!find_matches(&reparsed, &left).is_empty(), "{smarts}");
            assert!(!find_matches(&reparsed, &right).is_empty(), "{smarts}");
        }
    }

    #[test]
    fn rejects_unsupported_predicates() {
        let mut query = QueryMolecule::new();
        query.add_atom(AtomQuery::Primitive(AtomPrimitive::Charge(1)));
        assert_eq!(
            mcs_query_to_smarts(&query),
            Err(McsSmartsError::UnsupportedAtomQuery)
        );
    }

    #[test]
    fn serialized_alternatives_exclude_unobserved_choices() {
        for (inputs, unseen, config) in [
            (
                ["CO", "NO"],
                "PO",
                McsConfig {
                    atom_compare: crate::AtomCompare::Any,
                    ..McsConfig::default()
                },
            ),
            (
                ["C=N", "CN"],
                "C#N",
                McsConfig {
                    bond_compare: crate::BondCompare::Any,
                    ..McsConfig::default()
                },
            ),
        ] {
            let inputs = inputs.map(|smiles| chematic_smiles::parse(smiles).unwrap());
            let unseen = chematic_smiles::parse(unseen).unwrap();
            let query = find_mcs_with_config(&[&inputs[0], &inputs[1]], &config);
            let smarts = mcs_query_to_smarts(&query).unwrap();
            let reparsed = parse_smarts(&smarts).unwrap();
            assert!(find_matches(&reparsed, &unseen).is_empty(), "{smarts}");
        }
    }
}
