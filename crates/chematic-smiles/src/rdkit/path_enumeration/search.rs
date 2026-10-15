use super::{Neighbors, RdkitPathLimits, RdkitPathParams, RdkitSmilesError, reserve, unsupported};
use std::collections::HashSet;

fn limit(detail: &str) -> RdkitSmilesError {
    unsupported(format!("path enumeration resource limit: {detail}"))
}

fn append(
    frontier: &mut Vec<Vec<usize>>,
    indices: &mut usize,
    prefix: &[usize],
    atom: usize,
    limits: &RdkitPathLimits,
) -> Result<(), RdkitSmilesError> {
    let len = prefix.len() + 1;
    if frontier.len() >= limits.max_frontier {
        return Err(limit("frontier count exceeded"));
    }
    let total = indices
        .checked_add(len)
        .ok_or_else(|| limit("frontier indices overflow"))?;
    if total > limits.max_path_indices {
        return Err(limit("frontier indices exceeded"));
    }
    let mut path = Vec::new();
    reserve(&mut path, len)?;
    path.extend_from_slice(prefix);
    path.push(atom);
    reserve(frontier, 1)?;
    frontier.push(path);
    *indices = total;
    Ok(())
}

pub(super) fn enumerate(
    neighbors: &Neighbors,
    target: usize,
    params: &RdkitPathParams,
    limits: &RdkitPathLimits,
) -> Result<Vec<Vec<usize>>, RdkitSmilesError> {
    let mut frontier = Vec::new();
    let mut indices = 0;
    let n = neighbors.len();
    for atom in 0..n {
        if params.rooted_at_atom.is_none_or(|root| root == atom) {
            append(&mut frontier, &mut indices, &[], atom, limits)?;
        }
    }
    let distances = if params.only_shortest_paths {
        Some(distances(neighbors)?)
    } else {
        None
    };
    let mut visits = 0;
    for size in 1..target {
        if frontier.is_empty() {
            break;
        }
        let mut next = Vec::new();
        indices = 0;
        for path in &frontier {
            for &(atom, _) in &neighbors[path[size - 1]] {
                if visits >= limits.max_visits {
                    return Err(limit("visit count exceeded"));
                }
                visits += 1;
                if distances
                    .as_ref()
                    .is_some_and(|d| usize::from(d[path[0] * n + atom]) < size)
                {
                    continue;
                }
                let repeat = path.contains(&atom);
                let closure = target > 2 && size == target - 1 && atom != path[size - 2];
                if !repeat || closure {
                    append(&mut next, &mut indices, path, atom, limits)?;
                }
            }
        }
        frontier = next;
    }
    deduplicate(frontier, neighbors, params.use_bonds, limits)
}

// Unweighted all-pairs BFS is equivalent to RDKit's default distance matrix.
// Storage is bounded by the already checked 1024-atom model and allocated fallibly.
fn distances(neighbors: &Neighbors) -> Result<Vec<u16>, RdkitSmilesError> {
    let n = neighbors.len();
    let mut out = Vec::new();
    reserve(&mut out, n * n)?;
    out.resize(n * n, u16::MAX);
    let mut queue = Vec::new();
    reserve(&mut queue, n)?;
    for root in 0..n {
        queue.clear();
        queue.push(root);
        out[root * n + root] = 0;
        let mut cursor = 0;
        while cursor < queue.len() {
            let atom = queue[cursor];
            cursor += 1;
            for &(other, _) in &neighbors[atom] {
                if out[root * n + other] == u16::MAX {
                    out[root * n + other] = out[root * n + atom] + 1;
                    queue.push(other);
                }
            }
        }
    }
    Ok(out)
}

fn deduplicate(
    frontier: Vec<Vec<usize>>,
    neighbors: &Neighbors,
    use_bonds: bool,
    limits: &RdkitPathLimits,
) -> Result<Vec<Vec<usize>>, RdkitSmilesError> {
    let mut result = Vec::new();
    let mut seen = HashSet::new();
    let mut indices = 0usize;
    for path in frontier {
        let mut bonds = Vec::new();
        reserve(&mut bonds, path.len().saturating_sub(1))?;
        for pair in path.windows(2) {
            let edge = neighbors[pair[0]]
                .binary_search_by_key(&pair[1], |(atom, _)| *atom)
                .map_err(|_| unsupported("path edge missing from sanitized model"))?;
            bonds.push(neighbors[pair[0]][edge].1);
        }
        let mut key = Vec::new();
        reserve(&mut key, bonds.len())?;
        key.extend_from_slice(&bonds);
        key.sort_unstable();
        // Single-atom paths are individually retained, including isolated Hs.
        if !key.is_empty() {
            if seen.contains(&key) {
                continue;
            }
            seen.try_reserve(1)
                .map_err(|e| unsupported(format!("path dedup allocation failed: {e}")))?;
            seen.insert(key);
        }
        let output = if use_bonds { bonds } else { path };
        let total = indices
            .checked_add(output.len())
            .ok_or_else(|| limit("result indices overflow"))?;
        if result.len() >= limits.max_results {
            return Err(limit("result count exceeded"));
        }
        if total > limits.max_path_indices {
            return Err(limit("result indices exceeded"));
        }
        reserve(&mut result, 1)?;
        result.push(output);
        indices = total;
    }
    Ok(result)
}
