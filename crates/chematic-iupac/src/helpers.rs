//! Shared graph/naming helpers used by every compound-class handler module.

use chematic_core::{AtomIdx, BondOrder, Molecule};
use std::collections::{HashSet, VecDeque};

pub(crate) fn atoms_of(mol: &Molecule, atomic_num: u8) -> Vec<AtomIdx> {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() == atomic_num)
        .map(|(i, _)| i)
        .collect()
}

/// BFS count of C atoms reachable from `start` without crossing `blocked`.
pub(crate) fn count_c_chain(mol: &Molecule, start: AtomIdx, blocked: AtomIdx) -> usize {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    visited.insert(start);
    queue.push_back(start);
    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(cur) {
            if nb == blocked {
                continue;
            }
            if mol.atom(nb).element.atomic_number() == 6 && visited.insert(nb) {
                queue.push_back(nb);
            }
        }
    }
    visited.len()
}

/// Find the longest carbon chain in a C-subgraph using two-pass BFS.
///
/// Returns the sequence of AtomIdx forming the longest simple path.
/// For branched alkanes this gives the principal chain (IUPAC rule: longest chain).
pub(crate) fn find_longest_c_chain(mol: &Molecule, carbons: &[AtomIdx]) -> Vec<AtomIdx> {
    if carbons.is_empty() {
        return Vec::new();
    }

    let c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();

    // BFS to find the farthest atom from a given start, returning (farthest, parents).
    let bfs_far = |start: AtomIdx| -> (AtomIdx, std::collections::HashMap<AtomIdx, AtomIdx>) {
        let mut parent: std::collections::HashMap<AtomIdx, AtomIdx> =
            std::collections::HashMap::new();
        let mut visited: std::collections::HashSet<AtomIdx> = std::collections::HashSet::new();
        let mut queue = VecDeque::new();
        let mut farthest = start;
        visited.insert(start);
        queue.push_back(start);
        while let Some(cur) = queue.pop_front() {
            farthest = cur;
            for (nb, _) in mol.neighbors(cur) {
                if c_set.contains(&nb) && visited.insert(nb) {
                    parent.insert(nb, cur);
                    queue.push_back(nb);
                }
            }
        }
        (farthest, parent)
    };

    let reconstruct = |end: AtomIdx,
                       start: AtomIdx,
                       parents: &std::collections::HashMap<AtomIdx, AtomIdx>|
     -> Vec<AtomIdx> {
        let mut path = vec![end];
        let mut cur = end;
        while cur != start {
            cur = parents[&cur];
            path.push(cur);
        }
        path.reverse();
        path
    };

    // Pass 1: BFS from first carbon to find one endpoint of the longest chain.
    let (end1, _) = bfs_far(carbons[0]);
    // Pass 2: BFS from end1 to find the other endpoint.
    let (end2, parents) = bfs_far(end1);

    reconstruct(end2, end1, &parents)
}

/// Maximum number of tied-longest-chain candidates
/// [`find_longest_c_chain_candidates`] will return. Guards against
/// pathological blowup on highly symmetric/branched trees (same style as
/// `find_bridge_sizes`'s cap); real molecules with more than a handful of
/// chemically distinct equal-length arms off one branch point are not a
/// case this project has seen in practice.
const MAX_CHAIN_CANDIDATES: usize = 8;

/// Find every longest carbon chain tied for maximum length, not just one.
///
/// [`find_longest_c_chain`]'s two-pass BFS only ever returns a single
/// diameter path -- when a branch point has multiple arms of equal length
/// that are NOT graph-automorphic (chemically distinguishable, e.g. two
/// ethyl arms + one isopropyl arm, all reaching the same BFS depth), IUPAC
/// rule P-44.3 requires picking the chain with the most substituents among
/// those tied for length, which needs the full candidate set, not just one
/// arbitrary pick.
///
/// Finds one endpoint `A` via a single BFS pass from `carbons[0]` (valid
/// for any starting vertex in a tree: the farthest node from any vertex is
/// guaranteed to be an endpoint of *some* longest path), then does a full
/// BFS from `A` and returns the reconstructed path to every node tied for
/// maximum depth. This is deliberately NOT an exhaustive enumeration of
/// every longest path in the tree (a longest path not touching `A` at all
/// is possible in principle) -- it covers the realistic branch-point-tie
/// shape this rule exists for, bounded by [`MAX_CHAIN_CANDIDATES`].
///
/// Only used by [`crate::acyclic::Namer::name_branched_alkane`] --
/// [`find_longest_c_chain`] itself is untouched and still used by every
/// other caller (functional-group naming, where chain choice is
/// constrained by the functional group's position, a different problem).
pub(crate) fn find_longest_c_chain_candidates(
    mol: &Molecule,
    carbons: &[AtomIdx],
) -> Vec<Vec<AtomIdx>> {
    if carbons.is_empty() {
        return Vec::new();
    }
    if carbons.len() == 1 {
        return vec![vec![carbons[0]]];
    }

    let c_set: std::collections::HashSet<AtomIdx> = carbons.iter().copied().collect();

    let bfs_far = |start: AtomIdx| -> AtomIdx {
        let mut visited: std::collections::HashSet<AtomIdx> = std::collections::HashSet::new();
        let mut queue = VecDeque::new();
        let mut farthest = start;
        visited.insert(start);
        queue.push_back(start);
        while let Some(cur) = queue.pop_front() {
            farthest = cur;
            for (nb, _) in mol.neighbors(cur) {
                if c_set.contains(&nb) && visited.insert(nb) {
                    queue.push_back(nb);
                }
            }
        }
        farthest
    };

    // Pass 1: any endpoint of some longest chain.
    let anchor = bfs_far(carbons[0]);

    // Pass 2: full BFS from anchor, tracking depth + parent for every node.
    let mut parent: std::collections::HashMap<AtomIdx, AtomIdx> = std::collections::HashMap::new();
    let mut depth: std::collections::HashMap<AtomIdx, usize> = std::collections::HashMap::new();
    let mut visited: std::collections::HashSet<AtomIdx> = std::collections::HashSet::new();
    let mut queue = VecDeque::new();
    visited.insert(anchor);
    depth.insert(anchor, 0);
    queue.push_back(anchor);
    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(cur) {
            if c_set.contains(&nb) && visited.insert(nb) {
                parent.insert(nb, cur);
                depth.insert(nb, depth[&cur] + 1);
                queue.push_back(nb);
            }
        }
    }

    let max_depth = *depth.values().max().unwrap_or(&0);
    let mut endpoints: Vec<AtomIdx> = depth
        .iter()
        .filter(|&(_, &d)| d == max_depth)
        .map(|(&a, _)| a)
        .collect();
    // Deterministic order (not spelling-invariant -- see MAX_CHAIN_CANDIDATES
    // doc comment and the caller's own further tie-break for the residual gap).
    endpoints.sort_unstable();
    endpoints.truncate(MAX_CHAIN_CANDIDATES);

    endpoints
        .into_iter()
        .map(|end| {
            let mut path = vec![end];
            let mut cur = end;
            while cur != anchor {
                cur = parent[&cur];
                path.push(cur);
            }
            path.reverse();
            path
        })
        .collect()
}

/// Format substituents as an IUPAC prefix string ("2-methyl", "2,2-dimethyl", etc.).
pub(crate) fn format_substituents(subs: &[(usize, usize)]) -> String {
    // Group by alkyl name; sort alphabetically.
    let mut groups: std::collections::BTreeMap<&str, Vec<usize>> =
        std::collections::BTreeMap::new();
    for &(pos, len) in subs {
        let alkyl = match len {
            1 => "methyl",
            2 => "ethyl",
            3 => "propyl",
            4 => "butyl",
            _ => continue,
        };
        groups.entry(alkyl).or_default().push(pos);
    }

    let mut parts: Vec<String> = Vec::new();
    for (alkyl, mut positions) in groups {
        positions.sort_unstable();
        let locants = positions
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let mult = match positions.len() {
            1 => String::new(),
            2 => "di".to_string(),
            3 => "tri".to_string(),
            _ => "?".to_string(),
        };
        parts.push(format!("{}-{}{}", locants, mult, alkyl));
    }
    parts.join("-")
}

/// BFS chain anchored at `anchor` (always at index 0 = IUPAC position 1 in result).
pub(crate) fn chain_from_anchor(
    mol: &Molecule,
    c_set: &HashSet<AtomIdx>,
    anchor: AtomIdx,
) -> Vec<AtomIdx> {
    let mut parent: std::collections::HashMap<AtomIdx, AtomIdx> = std::collections::HashMap::new();
    let mut visited: HashSet<AtomIdx> = HashSet::new();
    let mut queue = VecDeque::new();
    let mut farthest = anchor;
    visited.insert(anchor);
    queue.push_back(anchor);
    while let Some(cur) = queue.pop_front() {
        farthest = cur;
        for (nb, _) in mol.neighbors(cur) {
            if c_set.contains(&nb) && visited.insert(nb) {
                parent.insert(nb, cur);
                queue.push_back(nb);
            }
        }
    }
    let mut path = vec![farthest];
    let mut cur = farthest;
    while cur != anchor {
        cur = parent[&cur];
        path.push(cur);
    }
    path.reverse();
    path
}

/// Return the IUPAC locant (1-based, lowest) of a double or triple bond on the chain.
pub(crate) fn unsaturation_locant(mol: &Molecule, carbons: &[AtomIdx], order: BondOrder) -> usize {
    let chain = find_longest_c_chain(mol, carbons);
    let n = chain.len();
    for (_, b) in mol.bonds() {
        if b.order == order
            && let (Some(p1), Some(p2)) = (
                chain.iter().position(|&c| c == b.atom1),
                chain.iter().position(|&c| c == b.atom2),
            )
        {
            let fwd = p1.min(p2) + 1; // 1-based lower position in forward direction
            let rev = n - p1.max(p2); // 1-based lower position in reversed direction
            return fwd.min(rev);
        }
    }
    1
}

/// Return ring atoms in cyclic traversal order.
pub(crate) fn ring_order_traversal(mol: &Molecule, ring_atoms: &HashSet<AtomIdx>) -> Vec<AtomIdx> {
    if ring_atoms.is_empty() {
        return Vec::new();
    }
    let start = *ring_atoms.iter().next().unwrap();
    let mut order = vec![start];
    let first_nb = mol
        .neighbors(start)
        .find(|(nb, _)| ring_atoms.contains(nb))
        .map(|(nb, _)| nb);
    let mut cur = match first_nb {
        Some(nb) => nb,
        None => return order,
    };
    let mut prev = start;
    while cur != start {
        order.push(cur);
        let next = mol
            .neighbors(cur)
            .find(|(nb, _)| ring_atoms.contains(nb) && *nb != prev)
            .map(|(nb, _)| nb);
        prev = cur;
        match next {
            Some(nb) => cur = nb,
            None => break,
        }
    }
    order
}

/// Find the minimum IUPAC locant assignment for `attach_points` on a ring.
/// Returns sorted `(locant, attach_atom)` pairs.
pub(crate) fn best_benzene_locants(
    mol: &Molecule,
    ring_atoms: &HashSet<AtomIdx>,
    attach_points: &[AtomIdx],
) -> Vec<(usize, AtomIdx)> {
    let ring_order = ring_order_traversal(mol, ring_atoms);
    let ring_n = ring_order.len();
    if ring_n == 0 {
        return Vec::new();
    }
    let n = attach_points.len();
    let pos_of: Vec<usize> = attach_points
        .iter()
        .map(|a| ring_order.iter().position(|r| r == a).unwrap_or(0))
        .collect();
    let mut best_locs: Option<Vec<usize>> = None;
    let mut best_assignment: Vec<(usize, AtomIdx)> = Vec::new();
    for start in 0..n {
        for &reverse in &[false, true] {
            let mut assignment: Vec<(usize, AtomIdx)> = Vec::new();
            for k in 0..n {
                let idx = (start + k) % n;
                let pos = if !reverse {
                    (pos_of[idx] + ring_n - pos_of[start]) % ring_n
                } else {
                    (pos_of[start] + ring_n - pos_of[idx]) % ring_n
                };
                assignment.push((pos + 1, attach_points[idx]));
            }
            assignment.sort_by_key(|&(l, _)| l);
            let locs: Vec<usize> = assignment.iter().map(|&(l, _)| l).collect();
            let is_better = best_locs.as_ref().is_none_or(|b| locs < *b);
            if is_better {
                best_locs = Some(locs);
                best_assignment = assignment;
            }
        }
    }
    best_assignment
}

pub(crate) fn count_components(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    if n == 0 {
        return 0;
    }
    let mut visited = vec![false; n];
    let mut count = 0;
    for start in 0..n {
        if visited[start] {
            continue;
        }
        count += 1;
        let mut queue = VecDeque::new();
        queue.push_back(AtomIdx(start as u32));
        visited[start] = true;
        while let Some(cur) = queue.pop_front() {
            for (nb, _) in mol.neighbors(cur) {
                if !visited[nb.0 as usize] {
                    visited[nb.0 as usize] = true;
                    queue.push_back(nb);
                }
            }
        }
    }
    count
}

/// Return the alkyl group name prefix (e.g. 1 → "methyl", 2 → "ethyl").
pub(crate) fn alkyl_prefix(n: usize) -> String {
    format!("{}yl", alkane_stem(n))
}

pub(crate) fn alkane_stem(n: usize) -> &'static str {
    match n {
        1 => "meth",
        2 => "eth",
        3 => "prop",
        4 => "but",
        5 => "pent",
        6 => "hex",
        7 => "hept",
        8 => "oct",
        9 => "non",
        10 => "dec",
        _ => "long",
    }
}

/// Find the lengths of all simple paths between two bridgehead atoms.
///
/// Returns a Vec of bridge lengths (number of atoms in bridge, excluding bridgeheads).
/// Uses DFS to enumerate paths, stopping early if more than 3 are found.
pub(crate) fn find_bridge_sizes(
    mol: &Molecule,
    bh0: AtomIdx,
    bh1: AtomIdx,
    ring_atoms: &HashSet<AtomIdx>,
) -> Vec<usize> {
    let mut bridges: Vec<usize> = Vec::new();
    let mut stack: Vec<(AtomIdx, Vec<AtomIdx>)> = vec![(bh0, vec![bh0])];

    while let Some((curr, path)) = stack.pop() {
        if bridges.len() > 4 {
            break; // Guard against complex polycyclics
        }
        for (nb, _) in mol.neighbors(curr) {
            if !ring_atoms.contains(&nb) {
                continue;
            }
            if nb == bh1 {
                // Found a bridge; length = path len - 1 (excludes bh0, counts intermediate atoms)
                bridges.push(path.len() - 1);
                continue;
            }
            if nb == bh0 || path.contains(&nb) {
                continue; // Avoid revisiting start or cycling
            }
            let mut new_path = path.clone();
            new_path.push(nb);
            stack.push((nb, new_path));
        }
    }
    bridges
}

/// Stem with "an" appended — base for most suffix compounds.
pub(crate) fn alkane_base(n: usize) -> String {
    format!("{}an", alkane_stem(n))
}

pub(crate) fn alkane_suffix(n: usize) -> String {
    match n {
        1 => "methane".into(),
        2 => "ethane".into(),
        3 => "propane".into(),
        4 => "butane".into(),
        5 => "pentane".into(),
        6 => "hexane".into(),
        7 => "heptane".into(),
        8 => "octane".into(),
        9 => "nonane".into(),
        10 => "decane".into(),
        11 => "undecane".into(),
        12 => "dodecane".into(),
        13 => "tridecane".into(),
        14 => "tetradecane".into(),
        15 => "pentadecane".into(),
        16 => "hexadecane".into(),
        17 => "heptadecane".into(),
        18 => "octadecane".into(),
        19 => "nonadecane".into(),
        20 => "icosane".into(),
        _ => format!("{n}alkane"),
    }
}

pub(crate) fn alkene_suffix(n: usize) -> String {
    alkane_suffix(n).replace("ane", "ene")
}
pub(crate) fn alkyne_suffix(n: usize) -> String {
    alkane_suffix(n).replace("ane", "yne")
}
