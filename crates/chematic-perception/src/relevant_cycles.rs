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

/// The relevant cycles of `mol` as atom rings (each in ring order), at most
/// `max_cycles` of them. Bonds of order zero and dative bonds are not part
/// of any ring; bonds to metals otherwise are (as in RDKit's ring finding).
pub fn relevant_cycles(
    mol: &Molecule,
    max_cycles: usize,
) -> Result<Vec<Vec<AtomIdx>>, RelevantCyclesTooMany> {
    let n = mol.atom_count();
    let mut adj: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    let mut n_edges = 0usize;
    for (_, bond) in mol.bonds() {
        if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
            continue;
        }
        let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if a == b {
            continue;
        }
        adj[a].push((b, n_edges));
        adj[b].push((a, n_edges));
        n_edges += 1;
    }
    // 2-core: only atoms on a cycle (or on a path between cycles) remain.
    let mut alive = vec![true; n];
    let mut degree: Vec<usize> = adj.iter().map(Vec::len).collect();
    let mut stack: Vec<usize> = (0..n).filter(|&v| degree[v] < 2).collect();
    while let Some(v) = stack.pop() {
        if !alive[v] {
            continue;
        }
        alive[v] = false;
        for &(w, _) in &adj[v] {
            if alive[w] {
                degree[w] -= 1;
                if degree[w] < 2 {
                    stack.push(w);
                }
            }
        }
    }
    if !alive.iter().any(|&a| a) {
        return Ok(Vec::new());
    }
    let words = n_edges.div_ceil(64);
    let edge_between = |a: usize, b: usize| -> usize {
        adj[a]
            .iter()
            .find(|&&(w, _)| w == b)
            .map(|&(_, e)| e)
            .expect("bonded")
    };

    let search = |r: usize| -> RootSearch {
        let mut dist = vec![UNREACHED; n];
        let mut preds: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut queue = std::collections::VecDeque::new();
        dist[r] = 0;
        queue.push_back(r);
        while let Some(v) = queue.pop_front() {
            for &(w, _) in &adj[v] {
                if w > r || !alive[w] {
                    continue;
                }
                if dist[w] == UNREACHED {
                    dist[w] = dist[v] + 1;
                    queue.push_back(w);
                }
                if dist[w] == dist[v] + 1 {
                    preds[w].push(v);
                }
            }
        }
        for p in &mut preds {
            p.sort_unstable();
        }
        RootSearch { dist, preds }
    };

    // The fixed shortest path from `v` back to the root (lowest-index
    // predecessor at each step), root last.
    let fixed_path = |s: &RootSearch, mut v: usize| -> Vec<usize> {
        let mut path = vec![v];
        while s.dist[v] > 0 {
            v = s.preds[v][0];
            path.push(v);
        }
        path
    };
    let meet_only_at_root = |a: &[usize], b: &[usize]| -> bool {
        // Both end at the root; compare the rest.
        let (a, b) = (&a[..a.len() - 1], &b[..b.len() - 1]);
        a.iter().all(|x| !b.contains(x))
    };
    let path_edges = |path: &[usize], bits: &mut [u64]| {
        for w in path.windows(2) {
            let e = edge_between(w[0], w[1]);
            bits[e / 64] ^= 1 << (e % 64);
        }
    };

    let mut candidates: Vec<Candidate> = Vec::new();
    let mut searches: Vec<Option<RootSearch>> = (0..n).map(|_| None).collect();
    for r in 0..n {
        if !alive[r] {
            continue;
        }
        let s = search(r);
        for y in 0..=r {
            if y == r || !alive[y] || s.dist[y] == UNREACHED {
                continue;
            }
            let py = fixed_path(&s, y);
            for &(z, e) in &adj[y] {
                if z < y && z <= r && alive[z] && s.dist[z] == s.dist[y] {
                    let pz = fixed_path(&s, z);
                    if meet_only_at_root(&py, &pz) {
                        let mut edges = vec![0u64; words];
                        path_edges(&py, &mut edges);
                        path_edges(&pz, &mut edges);
                        edges[e / 64] ^= 1 << (e % 64);
                        candidates.push(Candidate {
                            r,
                            closure: Closure::Odd { y, z },
                            len: 2 * s.dist[y] + 1,
                            edges,
                        });
                    }
                }
            }
            let preds = &s.preds[y];
            for i in 0..preds.len() {
                for j in (i + 1)..preds.len() {
                    let (p, q) = (preds[i], preds[j]);
                    let pp = fixed_path(&s, p);
                    let pq = fixed_path(&s, q);
                    if meet_only_at_root(&pp, &pq) {
                        let mut edges = vec![0u64; words];
                        path_edges(&pp, &mut edges);
                        path_edges(&pq, &mut edges);
                        let (e1, e2) = (edge_between(p, y), edge_between(q, y));
                        edges[e1 / 64] ^= 1 << (e1 % 64);
                        edges[e2 / 64] ^= 1 << (e2 % 64);
                        candidates.push(Candidate {
                            r,
                            closure: Closure::Even { p, y, q },
                            len: 2 * s.dist[y],
                            edges,
                        });
                    }
                }
            }
        }
        searches[r] = Some(s);
    }
    candidates.sort_by_key(|c| c.len);

    // GF(2) elimination by length class.
    let mut basis: Vec<(usize, Vec<u64>)> = Vec::new(); // (pivot bit, row)
    let reduce = |basis: &[(usize, Vec<u64>)], v: &mut Vec<u64>| {
        for (pivot, row) in basis {
            if v[pivot / 64] >> (pivot % 64) & 1 == 1 {
                for (a, b) in v.iter_mut().zip(row) {
                    *a ^= b;
                }
            }
        }
    };
    let leading = |v: &[u64]| -> Option<usize> {
        v.iter()
            .enumerate()
            .find(|(_, w)| **w != 0)
            .map(|(i, w)| i * 64 + w.trailing_zeros() as usize)
    };
    let mut relevant: Vec<usize> = Vec::new();
    let mut i = 0;
    while i < candidates.len() {
        let len = candidates[i].len;
        let mut j = i;
        while j < candidates.len() && candidates[j].len == len {
            j += 1;
        }
        let shorter = basis.len();
        for (k, cand) in candidates.iter().enumerate().take(j).skip(i) {
            let mut v = cand.edges.clone();
            reduce(&basis[..shorter], &mut v);
            if leading(&v).is_some() {
                relevant.push(k);
            }
        }
        for cand in &candidates[i..j] {
            let mut v = cand.edges.clone();
            reduce(&basis, &mut v);
            if let Some(pivot) = leading(&v) {
                // Keep rows reduced at their pivots so later reductions see
                // a triangular basis.
                for (_, row) in basis.iter_mut() {
                    if row[pivot / 64] >> (pivot % 64) & 1 == 1 {
                        for (a, b) in row.iter_mut().zip(&v) {
                            *a ^= b;
                        }
                    }
                }
                basis.push((pivot, v));
            }
        }
        i = j;
    }

    // Expand each relevant prototype into its family.
    let mut seen: std::collections::HashSet<Vec<u64>> = std::collections::HashSet::new();
    let mut out: Vec<Vec<AtomIdx>> = Vec::new();
    for k in relevant {
        let cand = &candidates[k];
        let s = searches[cand.r].as_ref().expect("searched");
        let (a, b, middle) = match cand.closure {
            Closure::Odd { y, z } => (y, z, None),
            Closure::Even { p, y, q } => (p, q, Some(y)),
        };
        let paths_a = all_shortest_paths(s, a, max_cycles)
            .ok_or(RelevantCyclesTooMany { cap: max_cycles })?;
        let paths_b = all_shortest_paths(s, b, max_cycles)
            .ok_or(RelevantCyclesTooMany { cap: max_cycles })?;
        for pa in &paths_a {
            for pb in &paths_b {
                if !meet_only_at_root(pa, pb) {
                    continue;
                }
                // Ring order: a ... r ... b (, y).
                let mut ring: Vec<usize> = pa.clone();
                ring.extend(pb.iter().rev().skip(1));
                if let Some(y) = middle {
                    ring.push(y);
                }
                let mut bits = vec![0u64; words];
                for w in 0..ring.len() {
                    let e = edge_between(ring[w], ring[(w + 1) % ring.len()]);
                    bits[e / 64] ^= 1 << (e % 64);
                }
                if seen.insert(bits) {
                    if out.len() >= max_cycles {
                        return Err(RelevantCyclesTooMany { cap: max_cycles });
                    }
                    out.push(ring.into_iter().map(|v| AtomIdx(v as u32)).collect());
                }
            }
        }
    }
    Ok(out)
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
