//! Ordered-ring oracle frozen from main 235df7a before this optimization.
use super::*;

#[allow(clippy::too_many_arguments)]
fn reference_smallest_rings_bfs_in(
    g: &Graph,
    root: usize,
    active: &[bool],
    forbidden: &[usize],
    done: &mut Vec<u8>,
    parents: &mut Vec<isize>,
    depths: &mut Vec<usize>,
    queue: &mut VecDeque<usize>,
) -> Result<Vec<Vec<usize>>, TooBig> {
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;
    let n = g.adj.len();
    done.clear();
    done.resize(n, WHITE);
    for &f in forbidden {
        done[f] = BLACK;
    }
    parents.clear();
    parents.resize(n, -1);
    depths.clear();
    depths.resize(n, 0);
    queue.clear();
    queue.push_back(root);
    let mut rings: Vec<Vec<usize>> = Vec::new();
    let mut cur_size = usize::MAX;
    while let Some(curr) = queue.pop_front() {
        if queue.len() + 1 >= MAX_BFSQ_SIZE {
            return Err(TooBig);
        }
        done[curr] = BLACK;
        let depth = depths[curr] + 1;
        if depth > cur_size {
            break;
        }
        for &(nbr, bi) in &g.adj[curr] {
            if !active[bi] {
                continue;
            }
            if done[nbr] == BLACK || parents[curr] == nbr as isize {
                continue;
            }
            if done[nbr] == WHITE {
                parents[nbr] = curr as isize;
                done[nbr] = GRAY;
                depths[nbr] = depth;
                queue.push_back(nbr);
            } else {
                let mut ring = Vec::with_capacity(depths[nbr] + depth + 1);
                ring.push(nbr);
                let mut parent = parents[nbr];
                while parent != -1 && parent != root as isize {
                    ring.push(parent as usize);
                    parent = parents[parent as usize];
                }
                ring.insert(0, curr);
                parent = parents[curr];
                while parent != -1 {
                    if ring.contains(&(parent as usize)) {
                        ring.clear();
                        break;
                    }
                    ring.insert(0, parent as usize);
                    parent = parents[parent as usize];
                }
                if ring.len() > 1 {
                    if ring.len() <= cur_size {
                        cur_size = ring.len();
                        rings.push(ring);
                    } else {
                        return Ok(rings);
                    }
                }
            }
        }
    }
    Ok(rings)
}

fn compare_search(g: &Graph, root: usize, active: &[bool], forbidden: &[usize]) -> usize {
    let mut done = Vec::new();
    let mut parents = Vec::new();
    let mut depths = Vec::new();
    let mut queue = VecDeque::new();
    let expected = reference_smallest_rings_bfs_in(
        g,
        root,
        active,
        forbidden,
        &mut done,
        &mut parents,
        &mut depths,
        &mut queue,
    );
    let actual = smallest_rings_bfs(g, root, active, forbidden);
    match (expected, actual) {
        (Ok(expected), Ok(actual)) => {
            assert_eq!(actual, expected, "root={root}, forbidden={forbidden:?}");
            actual.len()
        }
        (Err(_), Err(_)) => 0,
        _ => panic!("search limit differs: root={root}, forbidden={forbidden:?}"),
    }
}

fn graph(n: usize, ends: Vec<(usize, usize)>) -> Graph {
    let mut adj = vec![SmallVec::new(); n];
    for (index, &(left, right)) in ends.iter().enumerate() {
        adj[left].push((right, index));
        adj[right].push((left, index));
    }
    let orders = vec![BondOrder::Single; ends.len()];
    Graph { adj, ends, orders }
}

#[test]
fn every_small_graph_keeps_ordered_rings_for_masks_and_forbidden_atoms() {
    let n = 5;
    let edges: Vec<_> = (0..n)
        .flat_map(|a| ((a + 1)..n).map(move |b| (a, b)))
        .collect();
    let mut queries = 0usize;
    let mut rings = 0usize;
    for mask in 0usize..(1 << edges.len()) {
        let ends: Vec<_> = edges
            .iter()
            .enumerate()
            .filter_map(|(bit, &edge)| ((mask >> bit) & 1 != 0).then_some(edge))
            .collect();
        for reversed in [false, true] {
            let mut ends = ends.clone();
            if reversed {
                ends.reverse();
            }
            let g = graph(n, ends);
            for root in 0..n {
                let mut active = vec![true; g.ends.len()];
                rings += compare_search(&g, root, &active, &[]);
                rings += compare_search(&g, root, &active, &[(root + 1) % n]);
                queries += 2;
                for removed in 0..active.len() {
                    active[removed] = false;
                    rings += compare_search(&g, root, &active, &[]);
                    queries += 1;
                    active[removed] = true;
                }
            }
        }
    }
    assert!(rings > 0);
    eprintln!("small-graph ordered BFS oracle: {queries} queries, {rings} rings");
}

#[test]
fn pinned_corpora_keep_every_ordered_search_result() {
    let corpora = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"
        )),
    ];
    let mut rows = 0usize;
    let mut queries = 0usize;
    let mut rings = 0usize;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let smiles = line.split_whitespace().next().unwrap();
            let mol = chematic_smiles::parse_template(smiles).unwrap();
            let g = Graph::new(&mol);
            let mut active: Vec<_> = g.orders.iter().map(|&order| ring_eligible(order)).collect();
            for root in 0..g.adj.len() {
                rings += compare_search(&g, root, &active, &[]);
                rings += compare_search(&g, root, &active, &[(root + 1) % g.adj.len()]);
                queries += 2;
                if !active.is_empty() {
                    let removed = root % active.len();
                    let previous = active[removed];
                    active[removed] = false;
                    rings += compare_search(&g, root, &active, &[]);
                    queries += 1;
                    active[removed] = previous;
                }
            }
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert!(queries > 500_000 && rings > 0);
    eprintln!("pinned-corpus ordered BFS oracle: {rows} rows, {queries} queries, {rings} rings");
}

#[test]
fn long_paths_and_parallel_edges_keep_the_original_traversal_order() {
    for n in [3usize, 6, 32, 65, 128, 257, 1024] {
        let mut edges: Vec<_> = (0..n).map(|a| (a, (a + 1) % n)).collect();
        for extra in [false, true] {
            if extra {
                edges.push((0, n / 2));
            }
            let g = graph(n, edges.clone());
            for root in [0, 1, n / 2, n - 1] {
                compare_search(&g, root, &vec![true; g.ends.len()], &[]);
                compare_search(&g, root, &vec![true; g.ends.len()], &[n / 3]);
            }
        }
    }
    let g = graph(4, vec![(0, 1), (1, 2), (1, 2), (2, 3), (3, 0)]);
    for root in 0..4 {
        compare_search(&g, root, &[true; 5], &[]);
    }
}
