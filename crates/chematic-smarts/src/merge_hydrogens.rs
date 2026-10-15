//! Conservative top-level query-H folding with tetrahedral parity preservation.
use crate::{AtomPrimitive as A, AtomQuery as Q, BondPrimitive as B, BondQuery, QueryMolecule};

/// A query cannot be folded without exceeding the supported representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryHydrogenError {
    InvalidGraph,
    ResourceLimit,
    RecursiveQuery,
    UnsupportedStereo { atom: usize },
}

impl core::fmt::Display for QueryHydrogenError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidGraph => write!(f, "query has inconsistent bonds or adjacency"),
            Self::ResourceLimit => {
                write!(f, "query-H folding exceeds its graph or expression limits")
            }
            Self::RecursiveQuery => write!(f, "query-H folding does not support recursive SMARTS"),
            Self::UnsupportedStereo { atom } => write!(
                f,
                "query atom {atom} needs one removable H and a complete tetrahedral frame"
            ),
        }
    }
}
impl std::error::Error for QueryHydrogenError {}

fn plain_h(q: &Q) -> bool {
    match q {
        Q::Primitive(A::AtomicNum(1)) => true,
        Q::And(a, b) => {
            matches!(a.as_ref(), Q::Primitive(A::Symbol(s)) if s == "H")
                && matches!(b.as_ref(), Q::Primitive(A::Aromatic(false)))
        }
        _ => false,
    }
}

fn definitely_h(q: &Q) -> bool {
    match q {
        Q::Primitive(A::Symbol(s)) => s == "H",
        Q::Primitive(A::AtomicNum(1)) => true,
        Q::And(a, b) => definitely_h(a) || definitely_h(b),
        Q::Or(a, b) => definitely_h(a) && definitely_h(b),
        _ => false,
    }
}

fn validate_graph(q: &QueryMolecule) -> Result<(), QueryHydrogenError> {
    if q.adj.len() != q.atoms.len() {
        return Err(QueryHydrogenError::InvalidGraph);
    }
    let mut occurrences = vec![0usize; q.bonds.len()];
    for (i, neighbors) in q.adj.iter().enumerate() {
        for &(bi, other) in neighbors {
            let b = q.bonds.get(bi).ok_or(QueryHydrogenError::InvalidGraph)?;
            if other >= q.atoms.len()
                || i == other
                || !((b.atom1 == i && b.atom2 == other) || (b.atom2 == i && b.atom1 == other))
                || neighbors.iter().filter(|&&(idx, _)| idx == bi).count() != 1
            {
                return Err(QueryHydrogenError::InvalidGraph);
            }
            occurrences[bi] += 1;
        }
    }
    if occurrences.iter().any(|&count| count != 2) {
        return Err(QueryHydrogenError::InvalidGraph);
    }
    Ok(())
}

fn validate(q: &QueryMolecule) -> Result<(), QueryHydrogenError> {
    if q.atoms.len() > 1024 || q.bonds.len() > 4096 {
        return Err(QueryHydrogenError::ResourceLimit);
    }
    validate_graph(q)?;
    let mut count = 0;
    for (atom_index, atom) in q.atoms.iter().enumerate() {
        let mut todo = vec![(&atom.query, 0usize)];
        while let Some((expr, depth)) = todo.pop() {
            count += 1;
            if depth > 32 || count > 16384 {
                return Err(QueryHydrogenError::ResourceLimit);
            }
            match expr {
                Q::Primitive(A::Recursive(_)) => return Err(QueryHydrogenError::RecursiveQuery),
                Q::Primitive(A::Chirality(kind)) if !matches!(kind, 1 | 2) => {
                    return Err(QueryHydrogenError::UnsupportedStereo { atom: atom_index });
                }
                Q::And(a, b) | Q::Or(a, b) => {
                    todo.push((a, depth + 1));
                    todo.push((b, depth + 1));
                }
                Q::Not(a) => todo.push((a, depth + 1)),
                Q::Primitive(_) => {}
            }
        }
    }
    for bond in &q.bonds {
        let mut todo = vec![(&bond.query, 0usize)];
        while let Some((expr, depth)) = todo.pop() {
            count += 1;
            if depth > 32 || count > 16384 {
                return Err(QueryHydrogenError::ResourceLimit);
            }
            match expr {
                BondQuery::And(a, b) | BondQuery::Or(a, b) => {
                    todo.push((a, depth + 1));
                    todo.push((b, depth + 1));
                }
                BondQuery::Not(a) => todo.push((a, depth + 1)),
                _ => {}
            }
        }
    }
    Ok(())
}

fn has_stereo(q: &Q) -> bool {
    match q {
        Q::Primitive(A::Chirality(_)) => true,
        Q::And(a, b) | Q::Or(a, b) => has_stereo(a) || has_stereo(b),
        Q::Not(a) => has_stereo(a),
        _ => false,
    }
}

fn flip(q: &mut Q) {
    match q {
        Q::Primitive(A::Chirality(kind)) => *kind = if *kind == 1 { 2 } else { 1 },
        Q::And(a, b) | Q::Or(a, b) => {
            flip(a);
            flip(b);
        }
        Q::Not(a) => flip(a),
        _ => {}
    }
}

/// Fold plain, unmapped terminal `[H]`/`[#1]` nodes on single/implicit bonds.
///
/// Adds a minimum *total* H count (`!H0`, also `!H1` for two H nodes, etc.),
/// preserving existing predicates and the original graph's adjacency order.
/// Isotopic, mapped, charged or otherwise constrained H nodes, H-H bonds, and
/// conditional/non-single bonds are retained. Targets are ordinary chemical
/// graphs with monovalent H; this is not a bond-query equivalence for malformed
/// multivalent H targets. Recursive SMARTS are rejected rather than projected.
///
/// Affected chiral atoms require one folded H in a four-neighbour frame. Their
/// `@`/`@@` sense is adjusted for the virtual-H position used by the matcher.
/// Partial/ambiguous frames are rejected. Input is immutable. Caps: 1024 atoms,
/// 4096 bonds, 16384 atom-expression nodes, expression depth 32.
/// Mappings returned for the result refer to its compacted query indices.
pub fn merge_query_hydrogens(input: &QueryMolecule) -> Result<QueryMolecule, QueryHydrogenError> {
    validate(input)?;
    let n = input.atoms.len();
    let mut removed = vec![false; n];
    let mut counts = vec![0usize; n];
    for (h, atom) in input.atoms.iter().enumerate() {
        if atom.atom_map.is_some() || !plain_h(&atom.query) || input.adj[h].len() != 1 {
            continue;
        }
        let (bi, neighbor) = input.adj[h][0];
        if definitely_h(&input.atoms[neighbor].query)
            || !matches!(
                input.bonds[bi].query,
                BondQuery::Any | BondQuery::Primitive(B::Single)
            )
        {
            continue;
        }
        removed[h] = true;
        counts[neighbor] += 1;
    }
    let mut output = QueryMolecule::new();
    let mut atom_map = vec![None; n];
    for (i, atom) in input.atoms.iter().enumerate() {
        if !removed[i] {
            atom_map[i] = Some(output.add_atom_with_map(atom.query.clone(), atom.atom_map));
        }
    }
    let mut bond_map = vec![None; input.bonds.len()];
    for (i, b) in input.bonds.iter().enumerate() {
        if let (Some(a), Some(c)) = (atom_map[b.atom1], atom_map[b.atom2]) {
            bond_map[i] = Some(output.bonds.len());
            output.add_bond(a, c, b.query.clone());
        }
    }
    for i in 0..n {
        let Some(new_i) = atom_map[i] else {
            continue;
        };
        output.adj[new_i] = input.adj[i]
            .iter()
            .filter_map(|&(bi, j)| Some((bond_map[bi]?, atom_map[j]?)))
            .collect();
        if counts[i] == 0 {
            continue;
        }
        if counts[i] > 255 {
            return Err(QueryHydrogenError::ResourceLimit);
        }
        if has_stereo(&input.atoms[i].query) {
            if counts[i] != 1 || input.adj[i].len() != 4 {
                return Err(QueryHydrogenError::UnsupportedStereo { atom: i });
            }
            let old_position = input.adj[i]
                .iter()
                .position(|&(_, j)| removed[j])
                .ok_or(QueryHydrogenError::InvalidGraph)?;
            let new_position =
                usize::from(output.adj[new_i].first().is_some_and(|&(_, j)| j < new_i));
            if old_position.abs_diff(new_position) % 2 == 1 {
                flip(&mut output.atoms[new_i].query);
            }
        }
        for h in 0..counts[i] {
            let old = core::mem::replace(&mut output.atoms[new_i].query, Q::Primitive(A::Wildcard));
            output.atoms[new_i].query = Q::And(
                Box::new(old),
                Box::new(Q::Not(Box::new(Q::Primitive(A::HCount(h as u8))))),
            );
        }
    }
    validate(&output)?;
    Ok(output)
}
