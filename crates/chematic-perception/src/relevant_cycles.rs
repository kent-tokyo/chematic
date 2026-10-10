//! Relevant cycles: the union of all minimum cycle bases (Vismara 1997).
//!
//! A cycle is relevant when it is not the GF(2) sum of strictly shorter
//! cycles. Unlike an SSSR, the set does not depend on atom order or on a
//! choice between equally short rings, so every atom's count of relevant
//! cycles is an invariant. RDKit 2026.09.1 counts rings this way for
//! `[R<n>]` (RingDecomposerLib): on the 310,000-cell exposed-10k SMARTS grid
//! its 12 changed cells (six bis-quinolinium macrocycles) are exactly the
//! atoms whose relevant-cycle count differs from RDKit 2026.03.6's ring list.
//!
//! Vismara's algorithm: for every vertex `r`, a breadth-first search over
//! the vertices not after `r` (atom index order) gives shortest paths from
//! `r`; an edge `(y, z)` with `d(y) = d(z)` closes an odd candidate and two
//! predecessors `p`, `q` of a vertex `y` close an even one, when the two
//! shortest paths meet only at `r`. Candidates are taken by length; one is
//! relevant when it is independent of all shorter candidates. Each relevant
//! candidate stands for its family: every cycle formed by any pair of
//! shortest paths to the same two vertices that meet only at `r`.

use chematic_core::{AtomIdx, BondOrder, Molecule};

/// More relevant cycles than the caller's cap (families of equally short
/// rings in large fused or cage systems can grow combinatorially).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RelevantCyclesTooMany {
    /// The cap that was exceeded.
    pub cap: usize,
}

impl std::fmt::Display for RelevantCyclesTooMany {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "more than {} relevant cycles", self.cap)
    }
}

impl std::error::Error for RelevantCyclesTooMany {}

#[derive(Clone, Copy)]
enum Closure {
    /// Odd cycle: paths to `y` and `z`, edge `y`-`z`.
    Odd { y: usize, z: usize },
    /// Even cycle: paths to `p` and `q`, both bonded to `y`.
    Even { p: usize, y: usize, q: usize },
}

struct Candidate {
    r: usize,
    closure: Closure,
    len: usize,
    edges: Vec<u64>,
}

/// Per-root shortest-path data.
struct RootSearch {
    dist: Vec<usize>,
    preds: Vec<Vec<usize>>,
}

const UNREACHED: usize = usize::MAX;

struct CycleGraph {
    adj: Vec<Vec<(usize, usize)>>,
    alive: Vec<bool>,
    edge_count: usize,
}

impl CycleGraph {
    fn from_molecule(mol: &Molecule) -> Self {
        let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); mol.atom_count()];
        let mut edge_count = 0usize;
        for (_, bond) in mol.bonds() {
            if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
                continue;
            }
            let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            if a == b {
                continue;
            }
            adj[a].push((b, edge_count));
            adj[b].push((a, edge_count));
            edge_count += 1;
        }

        // Retain the 2-core: only atoms on a cycle (or on a path between
        // cycles) can contribute a relevant cycle.
        let mut alive = vec![true; mol.atom_count()];
        let mut degree: Vec<usize> = adj.iter().map(Vec::len).collect();
        let mut stack: Vec<usize> = (0..mol.atom_count())
            .filter(|&vertex| degree[vertex] < 2)
            .collect();
        while let Some(vertex) = stack.pop() {
            if !alive[vertex] {
                continue;
            }
            alive[vertex] = false;
            for &(other, _) in &adj[vertex] {
                if alive[other] {
                    degree[other] -= 1;
                    if degree[other] < 2 {
                        stack.push(other);
                    }
                }
            }
        }

        Self {
            adj,
            alive,
            edge_count,
        }
    }

    fn has_cycles(&self) -> bool {
        self.alive.iter().any(|&alive| alive)
    }

    fn words(&self) -> usize {
        self.edge_count.div_ceil(64)
    }

    fn edge_between(&self, a: usize, b: usize) -> usize {
        self.adj[a]
            .iter()
            .find(|&&(other, _)| other == b)
            .map(|&(_, edge)| edge)
            .expect("bonded")
    }

    fn search_from(&self, root: usize) -> RootSearch {
        let mut dist = vec![UNREACHED; self.adj.len()];
        let mut preds: Vec<Vec<usize>> = vec![Vec::new(); self.adj.len()];
        let mut queue = std::collections::VecDeque::new();
        dist[root] = 0;
        queue.push_back(root);
        while let Some(vertex) = queue.pop_front() {
            for &(other, _) in &self.adj[vertex] {
                if other > root || !self.alive[other] {
                    continue;
                }
                if dist[other] == UNREACHED {
                    dist[other] = dist[vertex] + 1;
                    queue.push_back(other);
                }
                if dist[other] == dist[vertex] + 1 {
                    preds[other].push(vertex);
                }
            }
        }
        for predecessors in &mut preds {
            predecessors.sort_unstable();
        }
        RootSearch { dist, preds }
    }

    fn fixed_path(search: &RootSearch, mut vertex: usize) -> Vec<usize> {
        let mut path = vec![vertex];
        while search.dist[vertex] > 0 {
            vertex = search.preds[vertex][0];
            path.push(vertex);
        }
        path
    }

    fn path_edges(&self, path: &[usize], bits: &mut [u64]) {
        for pair in path.windows(2) {
            let edge = self.edge_between(pair[0], pair[1]);
            bits[edge / 64] ^= 1 << (edge % 64);
        }
    }

    fn collect_candidates(&self) -> (Vec<Candidate>, Vec<Option<RootSearch>>) {
        let mut candidates = Vec::new();
        let mut searches: Vec<Option<RootSearch>> = (0..self.adj.len()).map(|_| None).collect();
        for (root, search_slot) in searches.iter_mut().enumerate() {
            if !self.alive[root] {
                continue;
            }
            let search = self.search_from(root);
            self.collect_root_candidates(root, &search, &mut candidates);
            *search_slot = Some(search);
        }
        candidates.sort_by_key(|candidate| candidate.len);
        (candidates, searches)
    }

    fn collect_root_candidates(
        &self,
        root: usize,
        search: &RootSearch,
        candidates: &mut Vec<Candidate>,
    ) {
        for y in 0..=root {
            if y == root || !self.alive[y] || search.dist[y] == UNREACHED {
                continue;
            }
            let path_y = Self::fixed_path(search, y);
            self.collect_odd_candidates(root, y, search, &path_y, candidates);
            self.collect_even_candidates(root, y, search, candidates);
        }
    }

    fn collect_odd_candidates(
        &self,
        root: usize,
        y: usize,
        search: &RootSearch,
        path_y: &[usize],
        candidates: &mut Vec<Candidate>,
    ) {
        for &(z, edge) in &self.adj[y] {
            if z >= y || z > root || !self.alive[z] || search.dist[z] != search.dist[y] {
                continue;
            }
            let path_z = Self::fixed_path(search, z);
            if !meet_only_at_root(path_y, &path_z) {
                continue;
            }
            let mut edges = vec![0u64; self.words()];
            self.path_edges(path_y, &mut edges);
            self.path_edges(&path_z, &mut edges);
            edges[edge / 64] ^= 1 << (edge % 64);
            candidates.push(Candidate {
                r: root,
                closure: Closure::Odd { y, z },
                len: 2 * search.dist[y] + 1,
                edges,
            });
        }
    }

    fn collect_even_candidates(
        &self,
        root: usize,
        y: usize,
        search: &RootSearch,
        candidates: &mut Vec<Candidate>,
    ) {
        let predecessors = &search.preds[y];
        for left in 0..predecessors.len() {
            for right in (left + 1)..predecessors.len() {
                let (p, q) = (predecessors[left], predecessors[right]);
                let path_p = Self::fixed_path(search, p);
                let path_q = Self::fixed_path(search, q);
                if !meet_only_at_root(&path_p, &path_q) {
                    continue;
                }
                let mut edges = vec![0u64; self.words()];
                self.path_edges(&path_p, &mut edges);
                self.path_edges(&path_q, &mut edges);
                for edge in [self.edge_between(p, y), self.edge_between(q, y)] {
                    edges[edge / 64] ^= 1 << (edge % 64);
                }
                candidates.push(Candidate {
                    r: root,
                    closure: Closure::Even { p, y, q },
                    len: 2 * search.dist[y],
                    edges,
                });
            }
        }
    }
}

fn meet_only_at_root(a: &[usize], b: &[usize]) -> bool {
    let (a, b) = (&a[..a.len() - 1], &b[..b.len() - 1]);
    a.iter().all(|vertex| !b.contains(vertex))
}

fn reduce_by_basis(basis: &[(usize, Vec<u64>)], vector: &mut [u64]) {
    for (pivot, row) in basis {
        if vector[pivot / 64] >> (pivot % 64) & 1 == 1 {
            for (value, basis_value) in vector.iter_mut().zip(row) {
                *value ^= basis_value;
            }
        }
    }
}

fn leading_bit(vector: &[u64]) -> Option<usize> {
    vector
        .iter()
        .enumerate()
        .find(|(_, word)| **word != 0)
        .map(|(index, word)| index * 64 + word.trailing_zeros() as usize)
}

fn relevant_candidate_indices(candidates: &[Candidate]) -> Vec<usize> {
    let mut basis: Vec<(usize, Vec<u64>)> = Vec::new();
    let mut relevant = Vec::new();
    let mut first = 0;
    while first < candidates.len() {
        let len = candidates[first].len;
        let mut end = first;
        while end < candidates.len() && candidates[end].len == len {
            end += 1;
        }

        let shorter_basis_len = basis.len();
        for (index, candidate) in candidates.iter().enumerate().take(end).skip(first) {
            let mut vector = candidate.edges.clone();
            reduce_by_basis(&basis[..shorter_basis_len], &mut vector);
            if leading_bit(&vector).is_some() {
                relevant.push(index);
            }
        }

        for candidate in &candidates[first..end] {
            let mut vector = candidate.edges.clone();
            reduce_by_basis(&basis, &mut vector);
            if let Some(pivot) = leading_bit(&vector) {
                for (_, row) in &mut basis {
                    if row[pivot / 64] >> (pivot % 64) & 1 == 1 {
                        for (value, candidate_value) in row.iter_mut().zip(&vector) {
                            *value ^= candidate_value;
                        }
                    }
                }
                basis.push((pivot, vector));
            }
        }
        first = end;
    }
    relevant
}

fn expand_relevant_families(
    graph: &CycleGraph,
    candidates: &[Candidate],
    searches: &[Option<RootSearch>],
    relevant: &[usize],
    max_cycles: usize,
) -> Result<Vec<Vec<AtomIdx>>, RelevantCyclesTooMany> {
    let mut seen = std::collections::HashSet::new();
    let mut rings = Vec::new();
    for &candidate_index in relevant {
        let candidate = &candidates[candidate_index];
        let search = searches[candidate.r].as_ref().expect("searched");
        let (a, b, middle) = match candidate.closure {
            Closure::Odd { y, z } => (y, z, None),
            Closure::Even { p, y, q } => (p, q, Some(y)),
        };
        let paths_a = all_shortest_paths(search, a, max_cycles)
            .ok_or(RelevantCyclesTooMany { cap: max_cycles })?;
        let paths_b = all_shortest_paths(search, b, max_cycles)
            .ok_or(RelevantCyclesTooMany { cap: max_cycles })?;
        for path_a in &paths_a {
            for path_b in &paths_b {
                if !meet_only_at_root(path_a, path_b) {
                    continue;
                }
                let mut ring = path_a.clone();
                ring.extend(path_b.iter().rev().skip(1));
                if let Some(middle) = middle {
                    ring.push(middle);
                }
                let mut bits = vec![0u64; graph.words()];
                for index in 0..ring.len() {
                    let edge = graph.edge_between(ring[index], ring[(index + 1) % ring.len()]);
                    bits[edge / 64] ^= 1 << (edge % 64);
                }
                if seen.insert(bits) {
                    if rings.len() >= max_cycles {
                        return Err(RelevantCyclesTooMany { cap: max_cycles });
                    }
                    rings.push(
                        ring.into_iter()
                            .map(|vertex| AtomIdx(vertex as u32))
                            .collect(),
                    );
                }
            }
        }
    }
    Ok(rings)
}

/// The relevant cycles of `mol` as atom rings (each in ring order), at most
/// `max_cycles` of them. Bonds of order zero and dative bonds are not part
/// of any ring; bonds to metals otherwise are (as in RDKit's ring finding).
pub fn relevant_cycles(
    mol: &Molecule,
    max_cycles: usize,
) -> Result<Vec<Vec<AtomIdx>>, RelevantCyclesTooMany> {
    let graph = CycleGraph::from_molecule(mol);
    if !graph.has_cycles() {
        return Ok(Vec::new());
    }
    let (candidates, searches) = graph.collect_candidates();
    let relevant = relevant_candidate_indices(&candidates);
    expand_relevant_families(&graph, &candidates, &searches, &relevant, max_cycles)
}

/// Every shortest path from `v` back to the search root (root last), or
/// `None` past `cap` paths.
fn all_shortest_paths(s: &RootSearch, v: usize, cap: usize) -> Option<Vec<Vec<usize>>> {
    let mut done: Vec<Vec<usize>> = Vec::new();
    let mut stack: Vec<Vec<usize>> = vec![vec![v]];
    while let Some(path) = stack.pop() {
        let last = *path.last().expect("non-empty");
        if s.dist[last] == 0 {
            done.push(path);
            if done.len() > cap {
                return None;
            }
            continue;
        }
        for &p in &s.preds[last] {
            let mut next = path.clone();
            next.push(p);
            stack.push(next);
        }
    }
    done.sort();
    Some(done)
}

/// Per atom, the number of relevant cycles through it (see
/// [`relevant_cycles`]).
pub fn relevant_cycle_counts(
    mol: &Molecule,
    max_cycles: usize,
) -> Result<Vec<usize>, RelevantCyclesTooMany> {
    let mut counts = vec![0usize; mol.atom_count()];
    for ring in relevant_cycles(mol, max_cycles)? {
        for a in ring {
            counts[a.0 as usize] += 1;
        }
    }
    Ok(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(smiles: &str) -> Vec<usize> {
        let mol = chematic_smiles_free_parse(smiles);
        relevant_cycle_counts(&mol, 10_000).unwrap()
    }

    // chematic-perception cannot depend on chematic-smiles; build graphs.
    fn chematic_smiles_free_parse(spec: &str) -> Molecule {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        // spec: "n;a-b,c-d,..."
        let (n, edges) = spec.split_once(';').unwrap();
        let n: usize = n.parse().unwrap();
        let mut b = MoleculeBuilder::new();
        for _ in 0..n {
            b.add_atom(Atom::new(Element::C));
        }
        for e in edges.split(',') {
            let (x, y) = e.split_once('-').unwrap();
            b.add_bond(
                AtomIdx(x.parse().unwrap()),
                AtomIdx(y.parse().unwrap()),
                BondOrder::Single,
            )
            .unwrap();
        }
        b.build()
    }

    #[test]
    fn cyclohexane_and_naphthalene_skeleton() {
        assert_eq!(counts("6;0-1,1-2,2-3,3-4,4-5,5-0"), vec![1; 6]);
        // decalin skeleton: two 6-rings sharing 0-5; the 10-ring is the sum.
        let c = counts("10;0-1,1-2,2-3,3-4,4-5,5-0,5-6,6-7,7-8,8-9,9-0");
        assert_eq!(c[0], 2);
        assert_eq!(c[5], 2);
        assert_eq!(c[1], 1);
        assert_eq!(c[7], 1);
    }

    #[test]
    fn cubane_has_six_relevant_four_rings() {
        // Cube: SSSR has 5 rings, all six faces are relevant.
        let c = counts("8;0-1,1-2,2-3,3-0,4-5,5-6,6-7,7-4,0-4,1-5,2-6,3-7");
        assert_eq!(c, vec![3; 8]);
    }

    #[test]
    fn bicyclo_222_has_three_relevant_six_rings() {
        // Bridgeheads 0 and 1, three two-atom bridges.
        let c = counts("8;0-2,2-3,3-1,0-4,4-5,5-1,0-6,6-7,7-1");
        assert_eq!(c[0], 3);
        assert_eq!(c[1], 3);
        assert_eq!(c[2], 2);
    }

    #[test]
    fn atom_order_does_not_change_counts() {
        // K4 relabelled: every triangle is relevant (4), each atom in 3.
        assert_eq!(counts("4;0-1,0-2,0-3,1-2,1-3,2-3"), vec![3; 4]);
        assert_eq!(counts("4;3-2,3-1,3-0,2-1,2-0,1-0"), vec![3; 4]);
    }

    #[test]
    fn cap_is_a_typed_error() {
        let mol = chematic_smiles_free_parse("8;0-1,1-2,2-3,3-0,4-5,5-6,6-7,7-4,0-4,1-5,2-6,3-7");
        assert_eq!(
            relevant_cycles(&mol, 3),
            Err(RelevantCyclesTooMany { cap: 3 })
        );
    }
}
