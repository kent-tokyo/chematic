use super::*;

// Frozen pre-optimization collector: compare every candidate, including its
// closure, root, length and edge vector, before relevance reduction/expansion.
impl CycleGraph {
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

    fn legacy_collect_root_candidates(
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
            self.legacy_collect_odd_candidates(root, y, search, &path_y, candidates);
            self.legacy_collect_even_candidates(root, y, search, candidates);
        }
    }

    fn legacy_collect_odd_candidates(
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

    fn legacy_collect_even_candidates(
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

fn candidate_key(c: &Candidate) -> (usize, usize, (usize, usize, usize, usize), Vec<u64>) {
    let closure = match c.closure {
        Closure::Odd { y, z } => (0, y, z, 0),
        Closure::Even { p, y, q } => (1, p, y, q),
    };
    (c.r, c.len, closure, c.edges.clone())
}

fn compare_collectors(graph: &CycleGraph) {
    let mut expected = Vec::new();
    for root in 0..graph.adj.len() {
        if !graph.alive[root] {
            continue;
        }
        let search = graph.search_from(root);
        let mut legacy = Vec::new();
        graph.legacy_collect_root_candidates(root, &search, &mut legacy);
        let mut optimized = Vec::new();
        graph.collect_root_candidates(root, &search, &mut optimized);
        assert_eq!(
            optimized.iter().map(candidate_key).collect::<Vec<_>>(),
            legacy.iter().map(candidate_key).collect::<Vec<_>>(),
            "root {root}"
        );
        expected.extend(legacy);
    }
    expected.sort_by_key(|candidate| candidate.len);
    let (actual, searches) = graph.collect_candidates();
    assert_eq!(
        actual.iter().map(candidate_key).collect::<Vec<_>>(),
        expected.iter().map(candidate_key).collect::<Vec<_>>()
    );
    for candidate in actual {
        assert!(searches[candidate.r].is_some());
    }
}

#[test]
fn all_five_vertex_graphs_retain_candidate_sequence() {
    let edges: Vec<_> = (0..5)
        .flat_map(|a| ((a + 1)..5).map(move |b| (a, b, true)))
        .collect();
    for mask in 0u32..(1 << edges.len()) {
        let bonds: Vec<_> = edges
            .iter()
            .enumerate()
            .filter_map(|(i, &e)| ((mask >> i) & 1 == 1).then_some(e))
            .collect();
        compare_collectors(&CycleGraph::from_bonds(5, &bonds));
    }
}

#[test]
fn multiword_vectors_disconnected_graphs_and_ineligible_edges() {
    // K12 has 66 edges: prefix copies and XOR must cross a u64 boundary.
    let bonds: Vec<_> = (0..12)
        .flat_map(|a| ((a + 1)..12).map(move |b| (a, b, true)))
        .collect();
    compare_collectors(&CycleGraph::from_bonds(12, &bonds));
    let mut joined = bonds;
    joined.extend([
        (12, 13, true),
        (13, 14, true),
        (14, 12, true),
        (15, 16, true),
        (16, 17, true),
        (17, 15, false),
        (1, 12, false),
    ]);
    compare_collectors(&CycleGraph::from_bonds(19, &joined));
}

#[test]
fn relabelled_sparse_and_dense_graphs_retain_candidate_sequence() {
    let mut state = 0x5eed_u64;
    for sample in 0..256 {
        let mut bonds = Vec::new();
        for a in 0..14 {
            for b in (a + 1)..14 {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                if state % 8 < (sample % 7 + 1) {
                    bonds.push(((a * 5 + 3) % 14, (b * 5 + 3) % 14, true));
                }
            }
        }
        compare_collectors(&CycleGraph::from_bonds(14, &bonds));
    }
}

#[test]
fn long_paths_cross_word_boundaries_without_changing_candidates() {
    let mut bonds = Vec::new();
    // Three 40-edge bridges between 0 and 1, with distinct internal atoms.
    for branch in 0..3 {
        let mut last = 0;
        for step in 0..39 {
            let vertex = 2 + branch * 39 + step;
            bonds.push((last, vertex, true));
            last = vertex;
        }
        bonds.push((last, 1, true));
    }
    compare_collectors(&CycleGraph::from_bonds(119, &bonds));
    bonds.reverse();
    for (a, b, _) in &mut bonds {
        *a = 118 - *a;
        *b = 118 - *b;
    }
    compare_collectors(&CycleGraph::from_bonds(119, &bonds));
}

#[test]
fn fixed_path_rows_exclude_peeled_tails_and_other_components() {
    let mut bonds = vec![
        (0, 1, true),
        (1, 2, true),
        (2, 0, true),
        (1001, 1002, true),
        (1002, 1003, true),
        (1003, 1001, true),
    ];
    for vertex in 3..1001 {
        bonds.push((vertex - 1, vertex, true));
    }
    let graph = CycleGraph::from_bonds(1004, &bonds);
    let search = graph.search_from(2);
    let paths = graph.fixed_tree_paths(2, &search);
    assert_eq!(paths.edges.len(), 3 * graph.words());
    assert!(paths.rows[3..].iter().all(|&row| row == UNREACHED));
    compare_collectors(&graph);
}
