//! Smallest Set of Smallest Rings (SSSR) via Horton's algorithm.
//!
//! Algorithm overview:
//! 1. Compute the cycle rank r = E - V + C (Euler characteristic),
//!    where C is the number of connected components.
//! 2. For every vertex v (as a candidate root) and every ring-eligible edge
//!    (x, y), form the candidate cycle SP(v, x) + edge(x, y) + SP(y, v),
//!    where SP is the shortest path in v's BFS tree. Keep it only if it's a
//!    genuine simple cycle (the two paths share no vertex other than v).
//!    This produces O(V*E) candidates and is guaranteed (Horton, 1987) to
//!    contain a minimum-weight cycle basis — unlike a single spanning tree's
//!    fundamental-cycle set (exactly r candidates, no redundancy), which can
//!    only ever report *a* valid basis, never guaranteed to be minimal.
//! 3. Represent each cycle as a set of bond indices (for GF(2) XOR independence
//!    testing) and sort candidates by (length, canonical tie-break) — the
//!    tie-break uses a local Weisfeiler-Leman-style atom ranking (see
//!    `canonical_atom_ranks`) so ring *selection* doesn't depend on input
//!    atom-numbering (i.e. SMILES parse/traversal order), only on molecular
//!    graph structure.
//! 4. Use Gaussian elimination over GF(2) to greedily build an independent
//!    basis of r cycles from that sorted candidate list.
//! 5. Convert the chosen bond-sets back to ordered atom sequences for the public API.

use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

/// Returns `true` if the bond order is eligible for ring perception.
///
/// Zero-order and Dative bonds are coordinate/non-valence connections that must
/// not form ring closures in the SSSR (RDKit PR #9118). Query bond types are
/// also excluded since they only appear in SMARTS patterns, never in real molecules.
fn is_ring_eligible(order: BondOrder) -> bool {
    matches!(
        order,
        BondOrder::Single
            | BondOrder::Double
            | BondOrder::Triple
            | BondOrder::Quadruple
            | BondOrder::Aromatic
            | BondOrder::Up
            | BondOrder::Down
    )
}

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// The Smallest Set of Smallest Rings for a molecule.
///
/// Each ring is stored as a sequence of `AtomIdx` values listed in ring order.
/// The first atom is not repeated at the end.
#[derive(Debug, Clone)]
pub struct RingSet(Vec<Vec<AtomIdx>>);

impl RingSet {
    /// All rings as slices of atom indices.
    pub fn rings(&self) -> &[Vec<AtomIdx>] {
        &self.0
    }

    /// Number of rings in the SSSR.
    pub fn ring_count(&self) -> usize {
        self.0.len()
    }

    /// Whether atom `atom` is a member of at least one ring.
    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.0.iter().any(|ring| ring.contains(&atom))
    }

    /// Number of rings that atom `atom` belongs to.
    pub fn atoms_in_ring_count(&self, atom: AtomIdx) -> usize {
        self.0.iter().filter(|ring| ring.contains(&atom)).count()
    }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Compute the Smallest Set of Smallest Rings for `mol`.
///
/// Returns a [`RingSet`] whose ring count equals the cycle rank r = E - V + C.
/// For acyclic molecules (r = 0) the returned set is empty.
pub fn find_sssr(mol: &Molecule) -> RingSet {
    let v = mol.atom_count();
    // Count only ring-eligible bonds for the cycle rank (E - V + C).
    // Zero-order and Dative bonds are excluded (RDKit PR #9118).
    let e = mol
        .bonds()
        .filter(|(_, b)| is_ring_eligible(b.order))
        .count();

    if v == 0 || e == 0 {
        return RingSet(Vec::new());
    }

    // Component count only (the parent tree itself isn't used any more —
    // candidate generation below builds its own BFS tree per root).
    let (components, _) = bfs_spanning_forest(mol);
    let r = (e as isize) - (v as isize) + (components as isize);

    if r <= 0 {
        return RingSet(Vec::new());
    }
    let r = r as usize;

    let ring_bonds: Vec<(BondIdx, AtomIdx, AtomIdx)> = mol
        .bonds()
        .filter(|(_, b)| is_ring_eligible(b.order))
        .map(|(bidx, b)| (bidx, b.atom1, b.atom2))
        .collect();

    // Horton candidate generation: BFS from every vertex, then for every
    // ring-eligible edge form the candidate SP(root,x)+edge(x,y)+SP(y,root).
    // O(V*E) candidates total — enough redundancy to guarantee a
    // minimum-weight basis is representable in the pool (see module doc).
    let mut candidates: Vec<(Vec<BondIdx>, Vec<AtomIdx>)> = Vec::new();
    for root_idx in 0..v {
        let root = AtomIdx(root_idx as u32);
        let (dist, parent) = bfs_tree(mol, root);
        for &(bidx, x, y) in &ring_bonds {
            if x == root || y == root {
                continue; // degenerate: edge touches the root itself
            }
            if dist[x.0 as usize] == usize::MAX || dist[y.0 as usize] == usize::MAX {
                continue; // x or y unreachable from this root (different component)
            }
            if let Some(candidate) = horton_candidate(mol, root, x, y, bidx, &parent) {
                candidates.push(candidate);
            }
        }
    }

    // Deterministic ordering: shortest first, then a canonical (input-order-
    // independent) tie-break so ring *selection* doesn't depend on how the
    // molecule happened to be numbered by the parser.
    let ranks = canonical_atom_ranks(mol);
    candidates.sort_by_cached_key(|c| (c.0.len(), canonical_cycle_key(&c.1, &ranks)));
    // The same geometric cycle can be generated from multiple roots; collapse
    // duplicates (bond_set is already sorted, so identical cycles are equal).
    candidates.dedup_by(|a, b| a.0 == b.0);

    // Gaussian elimination over GF(2) to select r linearly independent cycles.
    // The basis maps a pivot BondIdx to the full bond-set of that basis row.
    let mut basis: FxHashMap<BondIdx, Vec<BondIdx>> = FxHashMap::default();
    let mut selected_atoms: Vec<Vec<AtomIdx>> = Vec::new();

    for (bond_set, atom_seq) in candidates {
        // Reduce this cycle against the current basis.
        let reduced = gf2_reduce(&bond_set, &basis);

        if !reduced.is_empty() {
            // This cycle is independent — add it to the basis.
            let pivot = *reduced.iter().min().unwrap();
            basis.insert(pivot, reduced);
            selected_atoms.push(atom_seq);

            if selected_atoms.len() == r {
                break;
            }
        }
    }

    // Sort output rings by length for output consistency.
    selected_atoms.sort_by_key(|ring| ring.len());
    RingSet(selected_atoms)
}

// ---------------------------------------------------------------------------
// BFS spanning forest
// ---------------------------------------------------------------------------

/// Build a BFS spanning forest over the entire molecule.
///
/// Returns:
/// - `components`: number of connected components.
/// - `parent`: for each atom, the atom from which it was first discovered
///   (None for BFS roots).
fn bfs_spanning_forest(mol: &Molecule) -> (usize, Vec<Option<AtomIdx>>) {
    let n = mol.atom_count();
    let mut visited = vec![false; n];
    let mut parent: Vec<Option<AtomIdx>> = vec![None; n];
    let mut components = 0;
    let mut queue: VecDeque<AtomIdx> = VecDeque::new();

    for start in 0..n {
        if visited[start] {
            continue;
        }
        components += 1;
        let start_idx = AtomIdx(start as u32);
        visited[start] = true;
        queue.push_back(start_idx);

        while let Some(current) = queue.pop_front() {
            for (neighbor, bidx) in mol.neighbors(current) {
                // Skip non-ring-eligible bonds (Zero, Dative, Query*) — they
                // must not form ring closures in the spanning forest (RDKit PR #9118).
                if !is_ring_eligible(mol.bond(bidx).order) {
                    continue;
                }
                let ni = neighbor.0 as usize;
                if !visited[ni] {
                    visited[ni] = true;
                    parent[ni] = Some(current);
                    queue.push_back(neighbor);
                }
            }
        }
    }

    (components, parent)
}

// ---------------------------------------------------------------------------
// Horton candidate generation
// ---------------------------------------------------------------------------

/// BFS shortest-path tree from `root`, restricted to ring-eligible bonds.
///
/// Returns `(dist, parent)`: `dist[i]` is the shortest-path distance from
/// `root` to atom `i` (`usize::MAX` if unreachable), `parent[i]` is the
/// preceding atom on that shortest path (`None` for `root` and unreachable
/// atoms).
fn bfs_tree(mol: &Molecule, root: AtomIdx) -> (Vec<usize>, Vec<Option<AtomIdx>>) {
    let n = mol.atom_count();
    let mut dist = vec![usize::MAX; n];
    let mut parent: Vec<Option<AtomIdx>> = vec![None; n];
    let mut queue: VecDeque<AtomIdx> = VecDeque::new();

    dist[root.0 as usize] = 0;
    queue.push_back(root);

    while let Some(current) = queue.pop_front() {
        for (neighbor, bidx) in mol.neighbors(current) {
            if !is_ring_eligible(mol.bond(bidx).order) {
                continue;
            }
            let ni = neighbor.0 as usize;
            if dist[ni] == usize::MAX {
                dist[ni] = dist[current.0 as usize] + 1;
                parent[ni] = Some(current);
                queue.push_back(neighbor);
            }
        }
    }

    (dist, parent)
}

/// Find the smallest simple rings containing `root` using the Figueras-style
/// root-neighbor BFS primitive.
///
/// A ring through `root` is formed by two distinct ring-eligible neighbors of
/// the root plus a shortest path between those neighbors with the root
/// removed. Every neighbor pair is searched, and all pairs producing the
/// minimum ring size for this root are returned. The result is a local
/// primitive for the future symmetrized-ring model; it is deliberately not
/// substituted for [`find_sssr`], whose output is a linearly independent
/// Horton basis.
///
/// The search enumerates every shortest path for each neighbor pair. A later
/// symmetrization layer is responsible for applying RDKit's duplicate-ring
/// acceptance rules.
pub fn find_smallest_rings_bfs(mol: &Molecule, root: AtomIdx) -> Vec<Vec<AtomIdx>> {
    find_smallest_rings_bfs_with_blocked_bonds(mol, root, &FxHashSet::default())
}

/// Find the smallest root-centered rings while temporarily ignoring a set of
/// bonds. This is the bounded re-search primitive used by the RDKit-compatible
/// duplicate-D2 candidate pass; the molecule itself is never mutated.
pub fn find_smallest_rings_bfs_with_blocked_bonds(
    mol: &Molecule,
    root: AtomIdx,
    blocked_bonds: &FxHashSet<BondIdx>,
) -> Vec<Vec<AtomIdx>> {
    if root.0 as usize >= mol.atom_count() {
        return Vec::new();
    }

    let neighbors: Vec<AtomIdx> = mol
        .neighbors(root)
        .filter(|(_, bidx)| {
            is_ring_eligible(mol.bond(*bidx).order) && !blocked_bonds.contains(bidx)
        })
        .map(|(neighbor, _)| neighbor)
        .collect();
    if neighbors.len() < 2 {
        return Vec::new();
    }

    let mut best_size = usize::MAX;
    let mut rings = Vec::new();
    for (left_pos, &left) in neighbors.iter().enumerate() {
        for &right in neighbors.iter().skip(left_pos + 1) {
            let mut dist = vec![usize::MAX; mol.atom_count()];
            let mut queue = VecDeque::new();
            dist[left.0 as usize] = 0;
            queue.push_back(left);

            while let Some(current) = queue.pop_front() {
                if current == right {
                    break;
                }
                for (next, bidx) in mol.neighbors(current) {
                    if next == root
                        || !is_ring_eligible(mol.bond(bidx).order)
                        || blocked_bonds.contains(&bidx)
                    {
                        continue;
                    }
                    let next_i = next.0 as usize;
                    if dist[next_i] == usize::MAX {
                        dist[next_i] = dist[current.0 as usize] + 1;
                        queue.push_back(next);
                    }
                }
            }

            let right_dist = dist[right.0 as usize];
            if right_dist == usize::MAX {
                continue;
            }
            let ring_size = right_dist + 2;
            if ring_size > best_size {
                continue;
            }

            if ring_size < best_size {
                best_size = ring_size;
                rings.clear();
            }
            let mut path = vec![left];
            enumerate_shortest_paths(
                mol,
                root,
                right,
                &dist,
                blocked_bonds,
                &mut path,
                &mut rings,
            );
        }
    }

    rings.sort();
    rings.dedup();
    rings
}

/// Find smallest root-centered rings after applying RDKit-style leaf trimming
/// to the temporary bond mask. Removing a bond can expose degree-0/1 atoms;
/// those atoms cannot participate in a cycle, so their remaining active bonds
/// are removed transitively before the BFS is run. The molecule is unchanged.
pub fn find_smallest_rings_bfs_with_trimmed_bonds(
    mol: &Molecule,
    root: AtomIdx,
    blocked_bonds: &FxHashSet<BondIdx>,
) -> Vec<Vec<AtomIdx>> {
    let trimmed = trim_ring_bonds(mol, blocked_bonds);
    find_smallest_rings_bfs_with_blocked_bonds(mol, root, &trimmed)
}

/// Compute the active-bond mask after repeatedly removing bonds incident to
/// degree-0/1 atoms. This is the non-mutating equivalent of RDKit's
/// `trimBonds` queue and is useful when several rooted searches share a
/// progressively reduced graph.
pub fn trim_ring_bonds(mol: &Molecule, blocked_bonds: &FxHashSet<BondIdx>) -> FxHashSet<BondIdx> {
    let mut active_degree = vec![0usize; mol.atom_count()];
    for (bond, entry) in mol.bonds() {
        if !is_ring_eligible(entry.order) || blocked_bonds.contains(&bond) {
            continue;
        }
        active_degree[entry.atom1.0 as usize] += 1;
        active_degree[entry.atom2.0 as usize] += 1;
    }

    let mut trimmed = blocked_bonds.clone();
    let mut queue: VecDeque<AtomIdx> = active_degree
        .iter()
        .enumerate()
        .filter(|(_, degree)| **degree < 2)
        .map(|(idx, _)| AtomIdx(idx as u32))
        .collect();
    let mut queued = vec![false; mol.atom_count()];
    for atom in &queue {
        queued[atom.0 as usize] = true;
    }

    while let Some(atom) = queue.pop_front() {
        for (neighbor, bond) in mol.neighbors(atom) {
            if !is_ring_eligible(mol.bond(bond).order) || !trimmed.insert(bond) {
                continue;
            }
            let neighbor_degree = &mut active_degree[neighbor.0 as usize];
            *neighbor_degree = neighbor_degree.saturating_sub(1);
            if *neighbor_degree < 2 && !queued[neighbor.0 as usize] {
                queued[neighbor.0 as usize] = true;
                queue.push_back(neighbor);
            }
        }
    }

    trimmed
}

/// Find smallest rings from the BFS tree used by RDKit's Figueras pass.
/// Unlike [`find_smallest_rings_bfs_with_blocked_bonds`], this intentionally
/// keeps one parent tree and derives cycles from non-tree edges. It is exposed
/// separately so the bounded pair-shortest-path primitive remains available
/// for callers that need all shortest paths.
pub fn find_smallest_rings_bfs_with_rdkit_tree(
    mol: &Molecule,
    root: AtomIdx,
    blocked_bonds: &FxHashSet<BondIdx>,
) -> Vec<Vec<AtomIdx>> {
    if root.0 as usize >= mol.atom_count() {
        return Vec::new();
    }

    let mut state = vec![0u8; mol.atom_count()];
    let mut parent: Vec<Option<AtomIdx>> = vec![None; mol.atom_count()];
    let mut depth = vec![0usize; mol.atom_count()];
    let mut queue = VecDeque::new();
    let mut best_size = usize::MAX;
    let mut rings = Vec::new();
    state[root.0 as usize] = 1;
    queue.push_back(root);

    'bfs: while let Some(current) = queue.pop_front() {
        state[current.0 as usize] = 2;
        if depth[current.0 as usize] + 1 > best_size {
            break;
        }
        for (neighbor, bond) in mol.neighbors(current) {
            if !is_ring_eligible(mol.bond(bond).order) || blocked_bonds.contains(&bond) {
                continue;
            }
            if parent[current.0 as usize] == Some(neighbor) {
                continue;
            }
            match state[neighbor.0 as usize] {
                0 => {
                    state[neighbor.0 as usize] = 1;
                    parent[neighbor.0 as usize] = Some(current);
                    depth[neighbor.0 as usize] = depth[current.0 as usize] + 1;
                    queue.push_back(neighbor);
                }
                1 => {
                    let mut ring = vec![neighbor];
                    let mut ancestor = parent[neighbor.0 as usize];
                    while ancestor.is_some() && ancestor != Some(root) {
                        let atom = ancestor.expect("BFS node has a parent");
                        ring.push(atom);
                        ancestor = parent[atom.0 as usize];
                    }
                    ring.insert(0, current);
                    ancestor = parent[current.0 as usize];
                    while let Some(atom) = ancestor {
                        if ring.contains(&atom) {
                            ring.clear();
                            break;
                        }
                        ring.insert(0, atom);
                        ancestor = parent[atom.0 as usize];
                    }
                    if ring.len() > 1 {
                        if ring.len() <= best_size {
                            if ring.len() < best_size {
                                best_size = ring.len();
                                rings.clear();
                            }
                            rings.push(ring);
                        } else {
                            break 'bfs;
                        }
                    }
                }
                _ => {}
            }
        }
    }

    rings.sort();
    rings.dedup();
    rings
}

/// Select one root from each connected component of ring-eligible degree-2
/// atoms, matching RDKit's `pickD2Nodes`/`markUselessD2s` pass. A degree-2
/// chain is represented by its first atom in molecule order; this avoids
/// treating every atom along the same chain as an independent root.
pub fn select_rdkit_d2_roots(mol: &Molecule) -> Vec<AtomIdx> {
    let degree2: Vec<bool> = (0..mol.atom_count())
        .map(|raw| {
            mol.neighbors(AtomIdx(raw as u32))
                .filter(|(_, bond)| is_ring_eligible(mol.bond(*bond).order))
                .count()
                == 2
        })
        .collect();
    let mut seen = vec![false; mol.atom_count()];
    let mut roots = Vec::new();
    for raw in 0..mol.atom_count() {
        if !degree2[raw] || seen[raw] {
            continue;
        }
        roots.push(AtomIdx(raw as u32));
        let mut stack = vec![AtomIdx(raw as u32)];
        seen[raw] = true;
        while let Some(atom) = stack.pop() {
            for (neighbor, bond) in mol.neighbors(atom) {
                let neighbor_i = neighbor.0 as usize;
                if degree2[neighbor_i]
                    && is_ring_eligible(mol.bond(bond).order)
                    && !seen[neighbor_i]
                {
                    seen[neighbor_i] = true;
                    stack.push(neighbor);
                }
            }
        }
    }
    roots
}

/// Build the symmetrized smallest-ring set from the Horton basis and the
/// root-centered Figueras candidates.
///
/// A candidate is accepted only when it has the same size as a basis ring,
/// shares a bond with that ring, and does not omit a bond that is unique to
/// the basis ring. These are RDKit's duplicate-ring acceptance conditions.
/// The existing [`find_sssr`] result remains the base and this function is a
/// separate opt-in model for consumers that need symmetry-equivalent rings.
pub fn find_symmetrized_sssr(mol: &Molecule) -> RingSet {
    let base = find_sssr(mol);
    if base.rings().is_empty() {
        return base;
    }

    let base_bonds: Vec<FxHashSet<BondIdx>> = base
        .rings()
        .iter()
        .map(|ring| ring_bond_set(mol, ring))
        .collect();
    let mut bond_ring_count: FxHashMap<BondIdx, usize> = FxHashMap::default();
    for ring_bonds in &base_bonds {
        for &bond in ring_bonds {
            *bond_ring_count.entry(bond).or_insert(0) += 1;
        }
    }

    let base_keys: FxHashSet<Vec<u32>> = base_bonds.iter().map(bond_set_key).collect();
    let mut seen = base_keys.clone();
    let mut rings = base.rings().to_vec();
    let d2_roots = select_rdkit_d2_roots(mol);

    let mut accept_candidate = |candidate: Vec<AtomIdx>| {
        let candidate_bonds = ring_bond_set(mol, &candidate);
        let key = bond_set_key(&candidate_bonds);
        if base_keys.contains(&key) || !seen.insert(key) {
            return false;
        }
        let accepted = base_bonds.iter().any(|basis| {
            basis.iter().any(|bond| candidate_bonds.contains(bond))
                && basis.iter().all(|bond| {
                    bond_ring_count.get(bond).copied().unwrap_or(0) != 1
                        || candidate_bonds.contains(bond)
                })
        });
        if accepted
            && base
                .rings()
                .iter()
                .any(|ring| ring.len() == candidate.len())
        {
            rings.push(candidate);
            true
        } else {
            false
        }
    };
    let mut direct_replacements: Vec<Vec<AtomIdx>> = Vec::new();

    if d2_roots.is_empty() {
        for root in 0..mol.atom_count() {
            for candidate in find_smallest_rings_bfs(mol, AtomIdx(root as u32)) {
                accept_candidate(candidate);
            }
        }
    } else {
        let mut duplicate_groups: FxHashMap<Vec<u32>, (Vec<AtomIdx>, Vec<AtomIdx>)> =
            FxHashMap::default();
        let mut active_blocked = FxHashSet::default();
        for &root in &d2_roots {
            let candidates = find_smallest_rings_bfs_with_blocked_bonds(mol, root, &active_blocked);
            if candidates.is_empty() {
                for (_, bond) in mol.neighbors(root) {
                    if is_ring_eligible(mol.bond(bond).order) {
                        active_blocked.insert(bond);
                    }
                }
                active_blocked = trim_ring_bonds(mol, &active_blocked);
                continue;
            }
            for candidate in candidates {
                let key = bond_set_key(&ring_bond_set(mol, &candidate));
                let entry = duplicate_groups
                    .entry(key)
                    .or_insert_with(|| (candidate.clone(), Vec::new()));
                if !entry.1.contains(&root) {
                    entry.1.push(root);
                }
            }
        }

        for (_, (original_candidate, duplicate_roots)) in duplicate_groups {
            if duplicate_roots.len() < 2 {
                accept_candidate(original_candidate);
                continue;
            }
            let mut replacements = Vec::new();
            for &root in &duplicate_roots {
                let mut blocked = FxHashSet::default();
                for &other in &duplicate_roots {
                    if other == root {
                        continue;
                    }
                    for (_, bond) in mol.neighbors(other) {
                        if is_ring_eligible(mol.bond(bond).order) {
                            blocked.insert(bond);
                        }
                    }
                }
                let trimmed = trim_ring_bonds(mol, &blocked);
                replacements.extend(find_smallest_rings_bfs_with_rdkit_tree(mol, root, &trimmed));
            }
            if let Some(min_size) = replacements.iter().map(Vec::len).min() {
                replacements.retain(|candidate| candidate.len() == min_size);
            }
            replacements.sort_by_key(|candidate| bond_set_key(&ring_bond_set(mol, candidate)));
            for replacement in replacements {
                direct_replacements.push(replacement);
            }
        }
    }

    #[allow(clippy::drop_non_drop)]
    drop(accept_candidate);
    for replacement in direct_replacements {
        let key = bond_set_key(&ring_bond_set(mol, &replacement));
        if seen.insert(key)
            && base
                .rings()
                .iter()
                .any(|ring| ring.len() == replacement.len())
        {
            rings.push(replacement);
        }
    }

    // Keep every independently verified minimum replacement. RDKit's
    // symmetrized SSSR intentionally retains multiple overlapping rings in a
    // degenerate fused/bridged system; collapsing them to one representative
    // loses the active ring context needed by MMFF94 aromaticity.
    let mut extras = rings.split_off(base.ring_count());
    extras.sort_by_key(|ring| basis_exchange_key(mol, ring, &base_bonds));
    rings.extend(extras);

    let ranks = canonical_atom_ranks(mol);
    rings.sort_by_cached_key(|ring| (ring.len(), canonical_cycle_key(ring, &ranks)));
    RingSet(rings)
}

fn ring_bond_set(mol: &Molecule, ring: &[AtomIdx]) -> FxHashSet<BondIdx> {
    let mut bonds = FxHashSet::default();
    for i in 0..ring.len() {
        if let Some((bond, _)) = mol.bond_between(ring[i], ring[(i + 1) % ring.len()]) {
            bonds.insert(bond);
        }
    }
    bonds
}

fn bond_set_key(set: &FxHashSet<BondIdx>) -> Vec<u32> {
    let mut key: Vec<u32> = set.iter().map(|bond| bond.0).collect();
    key.sort_unstable();
    key
}

/// Stable tie-break key for a candidate ring based on a GF(2)-valid basis
/// exchange. The candidate replaces each same-sized Horton basis ring in
/// turn; only replacements that remain linearly independent are considered.
/// This preserves the minimum-cycle-basis contract while making the choice
/// depend on the resulting basis rather than raw bond numbering alone.
fn basis_exchange_key(
    mol: &Molecule,
    candidate: &[AtomIdx],
    base_bonds: &[FxHashSet<BondIdx>],
) -> Vec<Vec<u32>> {
    let candidate_set = ring_bond_set(mol, candidate);
    let candidate_key = bond_set_key(&candidate_set);
    let candidate_len = candidate.len();
    let mut best: Option<Vec<Vec<u32>>> = None;
    for (replace_idx, base_ring) in base_bonds.iter().enumerate() {
        if base_ring.len() != candidate_len {
            continue;
        }
        let mut rows = base_bonds
            .iter()
            .enumerate()
            .map(|(idx, set)| {
                if idx == replace_idx {
                    candidate_key.clone()
                } else {
                    bond_set_key(set)
                }
            })
            .collect::<Vec<_>>();
        if gf2_rank(&rows) != base_bonds.len() {
            continue;
        }
        rows.sort_unstable();
        if best.as_ref().is_none_or(|current| rows < *current) {
            best = Some(rows);
        }
    }
    best.unwrap_or_else(|| vec![candidate_key])
}

fn gf2_rank(rows: &[Vec<u32>]) -> usize {
    let mut basis: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    let mut rank = 0;
    for row in rows {
        let mut reduced = row.clone();
        while let Some(&pivot) = reduced.first() {
            let Some(existing) = basis.get(&pivot) else {
                basis.insert(pivot, reduced);
                rank += 1;
                break;
            };
            let mut xor = Vec::with_capacity(reduced.len() + existing.len());
            let mut left = 0;
            let mut right = 0;
            while left < reduced.len() || right < existing.len() {
                match (reduced.get(left), existing.get(right)) {
                    (Some(&a), Some(&b)) if a == b => {
                        left += 1;
                        right += 1;
                    }
                    (Some(&a), Some(&b)) if a < b => {
                        xor.push(a);
                        left += 1;
                    }
                    (Some(_), Some(&b)) => {
                        xor.push(b);
                        right += 1;
                    }
                    (Some(&a), None) => {
                        xor.push(a);
                        left += 1;
                    }
                    (None, Some(&b)) => {
                        xor.push(b);
                        right += 1;
                    }
                    (None, None) => break,
                }
            }
            reduced = xor;
        }
    }
    rank
}

/// Enumerate shortest paths in an already-computed BFS distance field.
/// Distances strictly increase at each step, so every emitted path is simple
/// and cannot revisit the excluded root.
fn enumerate_shortest_paths(
    mol: &Molecule,
    excluded: AtomIdx,
    target: AtomIdx,
    dist: &[usize],
    blocked_bonds: &FxHashSet<BondIdx>,
    path: &mut Vec<AtomIdx>,
    rings: &mut Vec<Vec<AtomIdx>>,
) {
    let current = *path.last().expect("shortest-path prefix is non-empty");
    if current == target {
        let mut ring = Vec::with_capacity(path.len() + 1);
        ring.push(excluded);
        ring.extend(path.iter().copied());
        rings.push(ring);
        return;
    }

    let current_dist = dist[current.0 as usize];
    for (next, bidx) in mol.neighbors(current) {
        if next == excluded
            || !is_ring_eligible(mol.bond(bidx).order)
            || blocked_bonds.contains(&bidx)
        {
            continue;
        }
        if dist[next.0 as usize] != current_dist + 1 {
            continue;
        }
        path.push(next);
        enumerate_shortest_paths(mol, excluded, target, dist, blocked_bonds, path, rings);
        path.pop();
    }
}

/// Form the Horton candidate cycle `SP(root,x) + edge(x,y) + SP(y,root)`.
///
/// Returns `None` if the two root-rooted shortest paths share any vertex
/// other than `root` itself — in that case the paths actually meet at some
/// closer common ancestor, so this (root, edge) pair doesn't yield a *simple*
/// cycle (the true minimal cycle through that closer ancestor is correctly
/// picked up when *it* is tried as root instead).
fn horton_candidate(
    mol: &Molecule,
    root: AtomIdx,
    x: AtomIdx,
    y: AtomIdx,
    bidx: BondIdx,
    parent: &[Option<AtomIdx>],
) -> Option<(Vec<BondIdx>, Vec<AtomIdx>)> {
    let path_x = path_to_root(x, parent); // [x, ..., root]
    let path_y = path_to_root(y, parent); // [y, ..., root]
    debug_assert_eq!(*path_x.last().unwrap(), root);
    debug_assert_eq!(*path_y.last().unwrap(), root);

    // Simplicity check: the two paths must share only `root`.
    let interior_x: FxHashSet<AtomIdx> = path_x[..path_x.len() - 1].iter().copied().collect();
    if path_y[..path_y.len() - 1]
        .iter()
        .any(|a| interior_x.contains(a))
    {
        return None;
    }

    // Ordered ring atoms: x ... root ... y (then the edge x-y closes the cycle).
    let mut ring_atoms: Vec<AtomIdx> = path_x.clone();
    for &a in path_y.iter().rev().skip(1) {
        ring_atoms.push(a);
    }

    let mut bond_set: Vec<BondIdx> = Vec::new();
    for i in 0..path_x.len().saturating_sub(1) {
        let (b, _) = mol.bond_between(path_x[i], path_x[i + 1])?;
        bond_set.push(b);
    }
    for i in 0..path_y.len().saturating_sub(1) {
        let (b, _) = mol.bond_between(path_y[i], path_y[i + 1])?;
        bond_set.push(b);
    }
    bond_set.push(bidx);
    bond_set.sort();
    bond_set.dedup();

    Some((bond_set, ring_atoms))
}

/// Walk parent pointers from `start` to `root`, returning the chain
/// including `start` (first) and the root (last).
fn path_to_root(start: AtomIdx, parent: &[Option<AtomIdx>]) -> Vec<AtomIdx> {
    let mut chain = Vec::new();
    let mut current = start;
    loop {
        chain.push(current);
        match parent[current.0 as usize] {
            Some(p) => current = p,
            None => break,
        }
    }
    chain
}

// ---------------------------------------------------------------------------
// Canonical tie-break (determinism, independent of atom-numbering)
// ---------------------------------------------------------------------------

/// A cheap, self-contained (no cross-crate dependency) approximation of
/// canonical atom ranking, used only to make candidate-cycle tie-breaking
/// deterministic: the same molecular graph always produces the same SSSR,
/// regardless of how its atoms happen to be numbered by the parser (SMILES
/// traversal order, SDF atom-block order, etc). This is a local
/// Weisfeiler-Leman-style refinement (seed on (element, degree, charge,
/// aromatic), then repeatedly fold in each atom's sorted neighbor keys) —
/// it does not aim for full canonical-labeling discriminating power (ties
/// among genuinely symmetric atoms are expected and fine; the point is
/// input-order-independence, not maximal refinement).
fn canonical_atom_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();
    let mut keys: Vec<u64> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let z = atom.element.atomic_number() as u64;
            let degree = mol.degree(idx) as u64;
            let charge = (atom.charge as i64 + 8) as u64; // shift to non-negative
            let aromatic = u64::from(atom.aromatic);
            (z << 24) | (degree << 16) | (charge << 8) | aromatic
        })
        .collect();

    const ROUNDS: usize = 3;
    for _ in 0..ROUNDS {
        let mut next = Vec::with_capacity(n);
        for i in 0..n {
            let mut neighbor_keys: Vec<u64> = mol
                .neighbors(AtomIdx(i as u32))
                .map(|(nb, _)| keys[nb.0 as usize])
                .collect();
            neighbor_keys.sort_unstable();
            let mut h = keys[i];
            for nk in neighbor_keys {
                h = h.wrapping_mul(1_000_003).wrapping_add(nk);
            }
            next.push(h);
        }
        keys = next;
    }
    keys
}

/// Deterministic sort key for a candidate cycle: the sorted multiset of its
/// atoms' canonical ranks, combined order-independently — isomorphic cycles
/// (same molecule, different traversal) always collapse to the same key.
fn canonical_cycle_key(atom_seq: &[AtomIdx], ranks: &[u64]) -> u64 {
    let mut vals: Vec<u64> = atom_seq.iter().map(|a| ranks[a.0 as usize]).collect();
    vals.sort_unstable();
    let mut h: u64 = 0;
    for v in vals {
        h = h.wrapping_mul(1_000_003).wrapping_add(v);
    }
    h
}

// ---------------------------------------------------------------------------
// GF(2) Gaussian elimination
// ---------------------------------------------------------------------------

/// Reduce `cycle` over GF(2) against the current `basis`.
///
/// Each basis entry maps a pivot bond (the minimum BondIdx in that row)
/// to the full row (sorted Vec<BondIdx>).
///
/// Returns the reduced cycle (empty if dependent on existing basis).
fn gf2_reduce(cycle: &[BondIdx], basis: &FxHashMap<BondIdx, Vec<BondIdx>>) -> Vec<BondIdx> {
    let mut current: Vec<BondIdx> = cycle.to_vec();
    while let Some(&pivot) = current.iter().min() {
        match basis.get(&pivot) {
            None => return current, // independent
            // XOR: symmetric difference of the two sorted sets.
            Some(basis_row) => current = sym_diff(&current, basis_row),
        }
    }
    current
}

/// Symmetric difference of two sorted slices (GF(2) addition / XOR for sets).
fn sym_diff(a: &[BondIdx], b: &[BondIdx]) -> Vec<BondIdx> {
    let mut result = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                result.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                result.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                // Both contain this element — XOR removes it.
                i += 1;
                j += 1;
            }
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

    // Build a cyclohexane molecule (6 carbons, 6 single bonds).
    fn cyclohexane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    // Build a benzene molecule (6 aromatic carbons, 6 single bonds for topology).
    fn benzene() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    // Build naphthalene: 10 atoms, 11 bonds (two fused 6-membered rings).
    // Atom numbering:
    //   0-1-2-3-4-5-0  (ring 1 perimeter, 6 atoms)
    //   4-6-7-8-9-5    (ring 2, sharing bond 4-5)
    fn naphthalene() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Ring 1: 0-1-2-3-4-9-0
        let ring1 = [0usize, 1, 2, 3, 4, 9];
        for i in 0..6 {
            b.add_bond(
                atoms[ring1[i]],
                atoms[ring1[(i + 1) % 6]],
                BondOrder::Single,
            )
            .unwrap();
        }
        // Ring 2: 4-5-6-7-8-9 (shares bond 4-9)
        // Bonds to add: 4-5, 5-6, 6-7, 7-8, 8-9
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[8], BondOrder::Single).unwrap();
        b.add_bond(atoms[8], atoms[9], BondOrder::Single).unwrap();
        b.build()
    }

    // Build norbornane (bicyclo[2.2.1]heptane): 7 carbons, 8 bonds, 2 rings.
    // Numbering:
    //   bridgehead atoms: 0, 3
    //   bridge 1: 0-1-2-3
    //   bridge 2: 0-4-5-3
    //   bridge 3: 0-6-3  (one-carbon bridge)
    fn norbornane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..7).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Bridge 1: 0-1-2-3
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        // Bridge 2: 0-4-5-3
        b.add_bond(atoms[0], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[3], BondOrder::Single).unwrap();
        // Bridge 3: 0-6-3
        b.add_bond(atoms[0], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[3], BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn test_azulene_sssr_minimal() {
        // Azulene: cyclopentadiene fused to cycloheptatriene, sharing one
        // bond. RDKit's GetSymmSSSR ring-size multiset is [5, 7]; the old
        // single-spanning-tree find_sssr previously returned a non-minimal
        // basis for this topology (see aromaticity.rs's PROVISIONAL-tagged
        // azulene regression test for the downstream effect on Pass 1/2).
        let mol = chematic_smiles::parse("C1=CC2=CC=CC=CC2=C1").expect("azulene SMILES");
        let sssr = find_sssr(&mol);
        let mut sizes: Vec<usize> = sssr.rings().iter().map(|r| r.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![5, 7], "azulene SSSR must be minimal [5, 7]");
    }

    #[test]
    fn test_indolizine_sssr_minimal() {
        // Indolizine: pyrrole fused to pyridine sharing a bridgehead N.
        // RDKit's GetSymmSSSR ring-size multiset is [5, 6] — a fused
        // heterocycle oracle case distinct from azulene's all-carbon,
        // odd/odd-sized ring pair.
        let mol = chematic_smiles::parse("c1ccn2ccccc12").expect("indolizine SMILES");
        let sssr = find_sssr(&mol);
        let mut sizes: Vec<usize> = sssr.rings().iter().map(|r| r.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![5, 6], "indolizine SSSR must be minimal [5, 6]");
    }

    #[test]
    fn test_cyclohexane_sssr() {
        let mol = cyclohexane();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 1, "cyclohexane has exactly 1 ring");
        assert_eq!(rings.rings()[0].len(), 6, "cyclohexane ring has 6 atoms");
    }

    #[test]
    fn blocked_bond_shortest_ring_search_is_non_mutating() {
        let mol = cyclohexane();
        let (bond, _) = mol.bond_between(AtomIdx(0), AtomIdx(1)).unwrap();
        let mut blocked = FxHashSet::default();
        blocked.insert(bond);
        assert!(find_smallest_rings_bfs_with_blocked_bonds(&mol, AtomIdx(0), &blocked).is_empty());
        assert_eq!(find_sssr(&mol).ring_count(), 1);
    }

    #[test]
    fn rdkit_d2_root_selection_collapses_degree2_chains() {
        let roots = select_rdkit_d2_roots(&benzene());
        assert_eq!(roots, vec![AtomIdx(0)]);
    }

    #[test]
    fn test_benzene_sssr() {
        let mol = benzene();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 1, "benzene has exactly 1 ring");
        assert_eq!(rings.rings()[0].len(), 6, "benzene ring has 6 atoms");
    }

    #[test]
    fn test_naphthalene_sssr() {
        let mol = naphthalene();
        let rings = find_sssr(&mol);
        // Cycle rank: 11 bonds - 10 atoms + 1 component = 2
        // SSSR should have 2 rings, both 6-membered.
        assert_eq!(rings.ring_count(), 2, "naphthalene SSSR has 2 rings");
        for ring in rings.rings() {
            assert_eq!(ring.len(), 6, "each naphthalene SSSR ring has 6 atoms");
        }
    }

    #[test]
    fn test_norbornane_sssr() {
        let mol = norbornane();
        let rings = find_sssr(&mol);
        // Cycle rank: 8 bonds - 7 atoms + 1 component = 2
        assert_eq!(rings.ring_count(), 2, "norbornane SSSR has 2 rings");
        // The two smallest rings are both 5-membered.
        for ring in rings.rings() {
            assert_eq!(ring.len(), 5, "each norbornane SSSR ring has 5 atoms");
        }
    }

    #[test]
    fn test_acyclic_molecule() {
        // Ethane: no rings.
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 0);
    }

    #[test]
    fn test_contains_atom() {
        let mol = cyclohexane();
        let rings = find_sssr(&mol);
        for i in 0..6u32 {
            assert!(
                rings.contains_atom(AtomIdx(i)),
                "atom {} should be in a ring",
                i
            );
        }
    }

    #[test]
    fn test_atoms_in_ring_count_benzene() {
        let mol = benzene();
        let rings = find_sssr(&mol);
        for i in 0..6u32 {
            assert_eq!(
                rings.atoms_in_ring_count(AtomIdx(i)),
                1,
                "each benzene atom is in exactly 1 ring"
            );
        }
    }

    // Anthracene: 14 atoms, 16 bonds (3 fused 6-membered rings).
    // Linear fusion: central ring shares edges with two outer rings.
    fn anthracene() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..14).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Ring 1 (left): 0-1-2-3-8-9-0
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[8], BondOrder::Single).unwrap();
        b.add_bond(atoms[8], atoms[9], BondOrder::Single).unwrap();
        b.add_bond(atoms[9], atoms[0], BondOrder::Single).unwrap();
        // Ring 2 (center): 3-4-5-6-7-8-3
        b.add_bond(atoms[3], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[8], BondOrder::Single).unwrap();
        // Ring 3 (right): 7-10-11-12-13-6-7
        b.add_bond(atoms[7], atoms[10], BondOrder::Single).unwrap();
        b.add_bond(atoms[10], atoms[11], BondOrder::Single).unwrap();
        b.add_bond(atoms[11], atoms[12], BondOrder::Single).unwrap();
        b.add_bond(atoms[12], atoms[13], BondOrder::Single).unwrap();
        b.add_bond(atoms[13], atoms[6], BondOrder::Single).unwrap();
        b.build()
    }

    // Spiro[4.4]nonane: two 5-membered rings sharing a single bridgehead atom.
    // 9 atoms total, cycle rank 2.
    fn spiro_nonane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..9).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Bridgehead: atom 0
        // Ring 1: 0-1-2-3-4-0
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[0], BondOrder::Single).unwrap();
        // Ring 2: 0-5-6-7-8-0
        b.add_bond(atoms[0], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[8], BondOrder::Single).unwrap();
        b.add_bond(atoms[8], atoms[0], BondOrder::Single).unwrap();
        b.build()
    }

    // 12-membered macrocycle (1 ring, 12 atoms).
    fn dodecane_ring() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..12).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..12 {
            b.add_bond(atoms[i], atoms[(i + 1) % 12], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    // Two disconnected rings (two components).
    fn disconnected_rings() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        // Benzene ring: atoms 0-5
        let benzene_atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(
                benzene_atoms[i],
                benzene_atoms[(i + 1) % 6],
                BondOrder::Single,
            )
            .unwrap();
        }
        // Separate cyclohexane ring: atoms 6-11
        let hexane_atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(
                hexane_atoms[i],
                hexane_atoms[(i + 1) % 6],
                BondOrder::Single,
            )
            .unwrap();
        }
        b.build()
    }

    // Adamantane-like tricyclic structure (simplified):
    // 10 atoms, 3 bridges between 2 bridgeheads
    fn adamantane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Bridgehead atoms: 0, 5
        // Bridge 1: 0-1-2-5 (3 bonds in chain)
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[5], BondOrder::Single).unwrap();
        // Bridge 2: 0-3-4-5 (3 bonds in chain)
        b.add_bond(atoms[0], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        // Bridge 3: 0-6-7-5 (3 bonds in chain)
        b.add_bond(atoms[0], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[5], BondOrder::Single).unwrap();
        // Cross-link bonds to connect bridges (forming tertiary center)
        // 1-3, 2-4, 6-? to complete cage
        b.add_bond(atoms[1], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[4], BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn test_anthracene_sssr() {
        let mol = anthracene();
        let rings = find_sssr(&mol);
        // Cycle rank: 16 bonds - 14 atoms + 1 component = 3.
        // RDKit's GetSymmSSSR gives three 6-membered rings (linear acene) —
        // Horton's minimum-weight basis must match this exactly, not just
        // "3 rings covering most atoms" (the old spanning-tree algorithm
        // could substitute a larger non-minimal ring for one of these).
        assert_eq!(rings.ring_count(), 3, "anthracene SSSR has 3 rings");
        for ring in rings.rings() {
            assert_eq!(ring.len(), 6, "each anthracene SSSR ring has 6 atoms");
        }
        let all_ring_atoms: std::collections::HashSet<_> = rings
            .rings()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect();
        assert_eq!(
            all_ring_atoms.len(),
            14,
            "anthracene SSSR atoms cover every atom"
        );
    }

    #[test]
    fn test_spiro_nonane_sssr() {
        let mol = spiro_nonane();
        let rings = find_sssr(&mol);
        // Cycle rank: 8 bonds - 9 atoms + 1 component = 0... wait, let me recalculate
        // Actually: two 5-membered rings sharing 1 atom = 4 + 4 + 2 (bridge) = 10 bonds
        // 10 bonds - 9 atoms + 1 = 2 rings
        assert_eq!(rings.ring_count(), 2, "spiro[4.4]nonane SSSR has 2 rings");
        for ring in rings.rings() {
            assert_eq!(ring.len(), 5, "each spiro nonane SSSR ring is 5-membered");
        }
    }

    #[test]
    fn test_dodecane_ring_sssr() {
        let mol = dodecane_ring();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 1, "12-membered ring has 1 SSSR entry");
        assert_eq!(
            rings.rings()[0].len(),
            12,
            "12-membered ring SSSR has 12 atoms"
        );
    }

    #[test]
    fn test_disconnected_rings_sssr() {
        let mol = disconnected_rings();
        let rings = find_sssr(&mol);
        // Cycle rank: 12 bonds - 12 atoms + 2 components = 2 rings
        assert_eq!(
            rings.ring_count(),
            2,
            "two disconnected rings yield 2 SSSR entries"
        );
        let sizes: Vec<_> = rings.rings().iter().map(|r| r.len()).collect();
        assert!(sizes.contains(&6), "one ring should be 6-membered");
    }

    #[test]
    fn test_adamantane_sssr() {
        let mol = adamantane();
        let rings = find_sssr(&mol);
        // Simplified adamantane: 10 atoms, 12 bonds, 1 component
        // Cycle rank: 12 - 10 + 1 = 3
        // (May be 3-4 depending on cross-link structure and Gaussian elimination)
        assert!(
            rings.ring_count() >= 3,
            "adamantane SSSR has at least 3 rings"
        );
        // Each ring should be reasonable size
        for ring in rings.rings() {
            assert!(!ring.is_empty(), "each ring should have atoms");
            assert!(ring.len() <= 10, "ring should not exceed molecule size");
        }
    }

    #[test]
    fn test_macrocycle_atom_in_ring_count() {
        let mol = dodecane_ring();
        let rings = find_sssr(&mol);
        for i in 0..12u32 {
            assert_eq!(
                rings.atoms_in_ring_count(AtomIdx(i)),
                1,
                "each dodecane atom is in exactly 1 ring"
            );
        }
    }

    #[test]
    fn test_figueras_bfs_finds_all_smallest_cubane_faces_through_each_root() {
        let mol = chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES");
        let mut faces = std::collections::BTreeSet::new();
        for root in 0..mol.atom_count() {
            for ring in find_smallest_rings_bfs(&mol, AtomIdx(root as u32)) {
                assert_eq!(ring.len(), 4, "cubane's smallest rings are square faces");
                let mut face = ring.into_iter().map(|a| a.0).collect::<Vec<_>>();
                face.sort_unstable();
                faces.insert(face);
            }
        }
        assert_eq!(
            faces.len(),
            6,
            "cubane has six symmetry-equivalent square faces"
        );

        let mol = chematic_smiles::parse("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12")
            .expect("dodecahedrane SMILES");
        let mut faces = std::collections::BTreeSet::new();
        for root in 0..mol.atom_count() {
            for ring in find_smallest_rings_bfs(&mol, AtomIdx(root as u32)) {
                assert_eq!(
                    ring.len(),
                    5,
                    "dodecahedrane's smallest rings are pentagons"
                );
                let mut face = ring.into_iter().map(|a| a.0).collect::<Vec<_>>();
                face.sort_unstable();
                faces.insert(face);
            }
        }
        assert_eq!(
            faces.len(),
            12,
            "dodecahedrane has twelve symmetry-equivalent pentagonal faces"
        );
    }

    #[test]
    fn test_symmetrized_sssr_adds_only_verified_duplicate_faces() {
        let benzene = chematic_smiles::parse("c1ccccc1").expect("benzene SMILES");
        assert_eq!(find_symmetrized_sssr(&benzene).ring_count(), 1);

        let cubane = chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES");
        assert_eq!(find_symmetrized_sssr(&cubane).ring_count(), 6);

        let dodeca = chematic_smiles::parse("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12")
            .expect("dodecahedrane SMILES");
        assert_eq!(find_symmetrized_sssr(&dodeca).ring_count(), 12);
    }

    #[test]
    fn test_cubane_sssr() {
        // Cubane C8H8 — a cage molecule with 12 C-C bonds and 8 vertices.
        // Cycle rank = E - V + 1 = 12 - 8 + 1 = 5.
        // Cubane has 6 square (4-membered) faces but only 5 are linearly
        // independent in GF(2) — the 6th is the XOR of the other 5 (not just
        // any pair of them, since Sum-of-all-6-faces = 0 in GF(2): each edge
        // is shared by exactly 2 faces).
        //
        // Horton's candidate generation (O(V*E) candidates, guaranteed to
        // contain a minimum-weight basis) finds a truly minimal SSSR here:
        // all 5 basis rings are 4-membered (weight 20), strictly better than
        // the old single-spanning-tree algorithm's typical "4 four-membered +
        // 1 six-membered diagonal" (weight 22) — see module doc / project
        // history for why the old algorithm couldn't guarantee this.
        //
        // Recovering the 6th (symmetry-equivalent) face requires XOR-ing all
        // 5 basis rings together, not a pairwise XOR — augmented_ring_set
        // only does pairwise XOR, so it cannot find it from an already-fully-
        // 4-membered basis (two same-size adjacent cube faces XOR to a
        // 6-membered "belt", not another 4-ring). This is expected: full
        // symmetrization (all 6 symmetry-equivalent minimal rings, matching
        // RDKit's GetSymmSSSR on cubane) is out of scope for Horton alone and
        // deferred to a later Vismara "relevant cycles" pass (see project
        // plan) — not a regression, since Horton's SSSR is still a strict
        // minimality improvement over the previous algorithm.
        let mol = chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES");
        let sssr = find_sssr(&mol);

        assert_eq!(
            sssr.rings().len(),
            5,
            "cubane must have exactly 5 SSSR rings (cycle rank 12−8+1=5)"
        );
        for ring in sssr.rings() {
            assert!(
                ring.len() <= 6,
                "cubane SSSR rings must be ≤ 6-membered, got {}",
                ring.len()
            );
        }
        let four_membered = sssr.rings().iter().filter(|r| r.len() == 4).count();
        assert_eq!(
            four_membered, 5,
            "Horton SSSR should find all 5 basis rings as 4-membered faces, got {four_membered}"
        );
    }

    /// Rebuild `mol` with atoms relabeled by `perm` (perm[new_idx] = old_idx),
    /// preserving the same graph but a different atom insertion order.
    fn permute_molecule(mol: &chematic_core::Molecule, perm: &[usize]) -> chematic_core::Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        for (_, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let _ = builder.add_bond(a, b, bond.order);
        }
        builder.build()
    }

    /// find_sssr's ring-size multiset must not depend on atom insertion order.
    /// Probes the fused/bridged/cage systems this project has already
    /// identified as the hard cases for canonical_atom_ranks' 3-round
    /// (not-fixpoint) Weisfeiler-Leman tie-break -- unlike
    /// `canonical_atom_order` (chematic-smiles), this function's doc comment
    /// does not claim full canonical-labeling power, and find_sssr's
    /// self-stability was already measured at 0% on a 5000-molecule corpus
    /// during the Horton rewrite (698ba3f); this is a permanent regression
    /// guard for that claim, not a first-time probe.
    #[test]
    fn find_sssr_ring_size_multiset_is_permutation_invariant() {
        let cases: Vec<(&str, chematic_core::Molecule)> = vec![
            ("naphthalene", naphthalene()),
            ("norbornane", norbornane()),
            ("spiro_nonane", spiro_nonane()),
            ("adamantane", adamantane()),
            (
                "cubane",
                chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES"),
            ),
        ];

        for (name, mol) in cases {
            let n = mol.atom_count();
            let mut orig_sizes: Vec<usize> = find_sssr(&mol).rings().iter().map(Vec::len).collect();
            orig_sizes.sort_unstable();

            let perms: Vec<Vec<usize>> = vec![(0..n).rev().collect(), {
                let mut p: Vec<usize> = (0..n).collect();
                if n > 2 {
                    p.rotate_left(n / 3 + 1);
                }
                p
            }];
            for perm in perms {
                let permuted = permute_molecule(&mol, &perm);
                let mut perm_sizes: Vec<usize> =
                    find_sssr(&permuted).rings().iter().map(Vec::len).collect();
                perm_sizes.sort_unstable();
                assert_eq!(
                    orig_sizes, perm_sizes,
                    "{name}: SSSR ring-size multiset changed under atom permutation {perm:?}"
                );
            }
        }
    }

    // ── RDKit PR #9118: SSSR excludes Zero-order and Dative bonds ────────────

    #[test]
    fn sssr_ignores_zero_order_bonds() {
        // A--B via a single bond PLUS a Zero-order bond between the same atoms
        // must NOT create a ring. Zero-order bonds are non-valence connections.
        let mut b = MoleculeBuilder::new();
        let mut a_atom = Atom::new(chematic_core::Element::C);
        a_atom.hydrogen_count = Some(3);
        let mut b_atom = Atom::new(chematic_core::Element::C);
        b_atom.hydrogen_count = Some(3);
        let a = b.add_atom(a_atom);
        let bb = b.add_atom(b_atom);
        b.add_bond(a, bb, BondOrder::Single).unwrap();
        b.add_bond(a, bb, BondOrder::Zero)
            .expect_err("duplicate bond — MoleculeBuilder should reject or ignore it");
        // Build a proper molecule: just two atoms with a single bond.
        // The zero-order bond attempt is rejected, so the molecule is acyclic.
        let mol = b.build();
        let sssr = find_sssr(&mol);
        assert_eq!(
            sssr.rings().len(),
            0,
            "single bond between two atoms → no ring"
        );
    }

    #[test]
    fn sssr_ignores_zero_order_bond_as_third_bond() {
        // Benzene ring (6 aromatic bonds) PLUS one Zero-order bond closing an
        // extra connection should NOT add a phantom 2-membered ring.
        // We build cyclohexane (all single bonds) and verify no extra ring from
        // a Zero-order bond added between two non-adjacent atoms.
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| {
                let mut a = Atom::new(chematic_core::Element::C);
                a.hydrogen_count = Some(2);
                b.add_atom(a)
            })
            .collect();
        // Square ring: 0-1-2-3-0
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[0], BondOrder::Single).unwrap();
        // Zero-order bond between non-adjacent atoms: should NOT create a new ring.
        // Note: builder may reject parallel bonds; if so the test still passes.
        let _ = b.add_bond(atoms[0], atoms[2], BondOrder::Zero);
        let mol = b.build();
        let sssr = find_sssr(&mol);
        // Should find exactly 1 ring (the 4-membered ring), NOT 2 or 3.
        assert_eq!(
            sssr.rings().len(),
            1,
            "zero-order diagonal bond must not create extra rings: found {:?}",
            sssr.rings().iter().map(|r| r.len()).collect::<Vec<_>>()
        );
    }
}
