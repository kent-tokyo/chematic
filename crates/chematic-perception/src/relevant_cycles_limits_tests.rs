use super::*;

// Frozen pre-change expansion: its single cap applies to both paths and rings.
fn frozen_expand_relevant_families(
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
                if !frozen_meet_only_at_root(path_a, path_b) {
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

fn frozen_meet_only_at_root(a: &[usize], b: &[usize]) -> bool {
    let (a, b) = (&a[..a.len() - 1], &b[..b.len() - 1]);
    a.iter().all(|vertex| !b.contains(vertex))
}

fn compare_limits(graph: &CycleGraph) {
    for root in 0..graph.adj.len() {
        let search = graph.search_from(root);
        let paths: Vec<_> = (0..=root)
            .filter(|&v| search.dist[v] != UNREACHED)
            .flat_map(|v| all_shortest_paths(&search, v, 100_000).unwrap())
            .collect();
        for a in &paths {
            for b in &paths {
                assert_eq!(meet_only_at_root(a, b), frozen_meet_only_at_root(a, b));
            }
        }
    }
    let (candidates, searches) = graph.collect_candidates();
    let relevant = relevant_candidate_indices(&candidates);
    let expected =
        frozen_expand_relevant_families(graph, &candidates, &searches, &relevant, 100_000).unwrap();
    for cap in 0..=expected.len() + 1 {
        let result =
            expand_relevant_families(graph, &candidates, &searches, &relevant, cap, 100_000);
        if expected.len() > cap {
            assert_eq!(result, Err(RelevantCyclesTooMany { cap }));
        } else {
            assert_eq!(result.unwrap(), expected);
        }
    }
    for cap in [0, 1, 2, 3, 8, 100_000] {
        assert_eq!(
            expand_relevant_families(graph, &candidates, &searches, &relevant, cap, cap),
            frozen_expand_relevant_families(graph, &candidates, &searches, &relevant, cap)
        );
    }
}

#[test]
fn separate_limits_preserve_ring_sequence_and_legacy_errors_on_all_five_vertex_graphs() {
    let edges: Vec<_> = (0..5)
        .flat_map(|a| ((a + 1)..5).map(move |b| (a, b, true)))
        .collect();
    for mask in 0u32..(1 << edges.len()) {
        let bonds: Vec<_> = edges
            .iter()
            .enumerate()
            .filter_map(|(i, &e)| ((mask >> i) & 1 == 1).then_some(e))
            .collect();
        compare_limits(&CycleGraph::from_bonds(5, &bonds));
    }
}

#[test]
fn output_limit_does_not_shorten_endpoint_path_search() {
    // Six paths to 8 overlap the sole path to 9 at vertex 1: no output ring.
    // A zero output limit must not preempt the separate path limit.
    let mut bonds = vec![(0, 1, true), (1, 9, true)];
    for v in 2..8 {
        bonds.extend([(1, v, true), (v, 8, true)]);
    }
    let graph = CycleGraph::from_bonds(10, &bonds);
    let mut preds = vec![Vec::new(); 10];
    preds[1] = vec![0];
    for item in &mut preds[2..8] {
        *item = vec![1];
    }
    preds[8] = (2..8).collect();
    preds[9] = vec![1];
    let searches = vec![Some(RootSearch {
        dist: vec![0, 1, 2, 2, 2, 2, 2, 2, 3, 2],
        preds,
    })];
    let candidates = vec![Candidate {
        r: 0,
        closure: Closure::Odd { y: 8, z: 9 },
        len: 6,
        edges: Vec::new(),
    }];
    assert_eq!(
        expand_relevant_families(&graph, &candidates, &searches, &[0], 0, 6),
        Ok(Vec::new())
    );
    assert_eq!(
        expand_relevant_families(&graph, &candidates, &searches, &[0], 10, 5),
        Err(RelevantCyclesTooMany { cap: 5 })
    );
}

#[test]
fn macrocycles_cages_and_disconnected_graphs_preserve_capped_outputs() {
    let mut bonds: Vec<_> = (0..8)
        .flat_map(|a| ((a + 1)..8).map(move |b| (a, b, true)))
        .collect();
    bonds.extend([(8, 9, true), (9, 10, true), (10, 8, true), (0, 8, false)]);
    compare_limits(&CycleGraph::from_bonds(12, &bonds));
    let mut theta = Vec::new();
    for branch in 0..3 {
        let mut last = 0;
        for step in 0..39 {
            let v = 2 + branch * 39 + step;
            theta.push((last, v, true));
            last = v;
        }
        theta.push((last, 1, true));
    }
    compare_limits(&CycleGraph::from_bonds(120, &theta));
}
