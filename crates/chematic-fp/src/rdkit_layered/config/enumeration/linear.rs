//! Rooted bond paths with RDKit's range-wide final-step closure convention.
use super::*;
use std::collections::BTreeSet;
type Adjacency = Vec<Vec<(usize, u32)>>;
pub(super) fn enumerate(
    mol: &Molecule,
    c: &RdkitLayeredConfig,
    budget: &mut Budget,
    visit: &mut Visitor<'_>,
) -> Result<(), RdkitLayeredError> {
    let mut adjacency: Adjacency = vec![Vec::new(); mol.atom_count()];
    for (i, b) in mol.bonds() {
        if mol.atom(b.atom1).element.atomic_number() == 1
            || mol.atom(b.atom2).element.atomic_number() == 1
        {
            continue;
        }
        adjacency[b.atom1.0 as usize].push((b.atom2.0 as usize, i.0));
        adjacency[b.atom2.0 as usize].push((b.atom1.0 as usize, i.0));
    }
    for a in &mut adjacency {
        a.sort_unstable();
    }
    for root in &c.from_atoms {
        let mut frontier = vec![vec![root.0 as usize]];
        for length in 1..=c.max_path + 1 {
            if length > c.min_path {
                emit(&frontier, &adjacency, budget, visit)?;
            }
            if length == c.max_path + 1 {
                break;
            }
            frontier = extend(frontier, &adjacency, c.max_path + 1, budget)?;
            if frontier.is_empty() {
                break;
            }
        }
    }
    Ok(())
}
fn emit(
    frontier: &[Vec<usize>],
    adj: &Adjacency,
    budget: &mut Budget,
    visit: &mut Visitor<'_>,
) -> Result<(), RdkitLayeredError> {
    let mut seen = BTreeSet::new();
    for atoms in frontier {
        budget.step()?;
        let edges = atoms
            .windows(2)
            .map(|p| {
                adj[p[0]]
                    .iter()
                    .find_map(|&(a, b)| (a == p[1]).then_some(b))
                    .ok_or_else(|| error("internal", "path adjacency is missing"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut key = edges.clone();
        key.sort_unstable();
        if seen.insert(key) {
            budget.record(1)?;
            visit(&edges, 1);
        }
    }
    Ok(())
}
fn extend(
    frontier: Vec<Vec<usize>>,
    adj: &Adjacency,
    max_atoms: usize,
    budget: &mut Budget,
) -> Result<Vec<Vec<usize>>, RdkitLayeredError> {
    let mut next = Vec::new();
    for path in frontier {
        // safe: every frontier is initialized with a root and only extended.
        let last = *path.last().expect("nonempty rooted walk");
        for &(atom, _) in &adj[last] {
            budget.step()?;
            let final_closure = max_atoms > 2
                && path.len() + 1 == max_atoms
                && path.get(path.len().saturating_sub(2)) != Some(&atom);
            if path.contains(&atom) && !final_closure {
                continue;
            }
            budget.frontier(next.len() + 1)?;
            let mut extended = path.clone();
            extended.push(atom);
            next.push(extended);
        }
    }
    Ok(next)
}
