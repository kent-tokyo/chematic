//! The order in which RDKit lists a molecule's rings.
//!
//! A port of RDKit 2026.03.6 `FindRings.cpp` (`findSSSR` followed by
//! `symmetrizeSSSR`, the ring perception of `SanitizeMol`). chematic's own
//! [`crate::find_sssr`] returns the same ring *sets* in a canonical order;
//! RDKit's order follows its Figueras-style search (rings at degree-two atoms
//! first, then at the first degree-three atom, ...) and decides results of
//! RDKit algorithms that walk `RingInfo` in order and stop early, such as
//! MMFF94's aromaticity pass.

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use smallvec::SmallVec;

/// RDKit `RingUtils::MAX_BFSQ_SIZE`.
const MAX_BFSQ_SIZE: usize = 200_000;

struct Graph {
    /// Per atom: (neighbour, bond index) in RDKit's bond creation order.
    adj: Vec<SmallVec<[(usize, usize); 4]>>,
    /// Per bond (RDKit's index): its two atoms.
    ends: Vec<(usize, usize)>,
    /// Per bond (RDKit's index): its order.
    orders: Vec<BondOrder>,
}

/// RDKit's largest allowed valence per atomic number (`-1`: any).
const MAX_VALENCE: [i8; 119] = [
    -1, 1, 0, -1, 2, 3, 4, 3, 2, 1, 0, -1, -1, 3, 4, 5, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, 3, 4, 5, 6, 1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, 4, 5, 6,
    5, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, 4, 5, 6, 5, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
];

/// RDKit `QueryOps::isMetal`.
fn is_metal(z: u8) -> bool {
    !matches!(
        z,
        0 | 1
            | 2
            | 5
            | 6
            | 7
            | 8
            | 9
            | 10
            | 14
            | 15
            | 16
            | 17
            | 18
            | 33
            | 34
            | 35
            | 36
            | 52
            | 53
            | 54
            | 85
            | 86
    )
}

/// The single bonds `SanitizeMol`'s `cleanUpOrganometallics` turns into
/// dative bonds (which RDKit's ring perception skips): for each non-metal,
/// in canonical-rank order, whose explicit valence is over the largest
/// valence of its charge-shifted element (or equal to it on an aromatic
/// atom of total degree four), the single bond to the bonded metal with the
/// fewest dative bonds so far (higher canonical rank first on a tie).
fn organometallic_dative_bonds(mol: &Molecule) -> Vec<bool> {
    let n = mol.atom_count();
    let mut dative = vec![false; mol.bond_count()];
    let z = |a: usize| -> u8 {
        let atom = mol.atom(AtomIdx(a as u32));
        if atom.wildcard {
            0
        } else {
            atom.element.atomic_number()
        }
    };
    if !(0..n).any(|a| is_metal(z(a))) {
        return dative;
    }
    let explicit_valence = |a: usize, dative: &[bool]| -> i32 {
        let idx = AtomIdx(a as u32);
        let mut accum = 0.0f64;
        for (_, bi) in mol.neighbors(idx) {
            let b = mol.bond(bi);
            accum += if dative[bi.0 as usize] {
                if b.atom2 == idx { 1.0 } else { 0.0 }
            } else {
                match b.order {
                    BondOrder::Aromatic => 1.5,
                    BondOrder::Dative => {
                        if b.atom2 == idx {
                            1.0
                        } else {
                            0.0
                        }
                    }
                    o => f64::from(o.order_int()),
                }
            };
        }
        accum += f64::from(mol.atom(idx).hydrogen_count.unwrap_or(0));
        (accum + 0.1).round() as i32
    };
    let hypervalent = |a: usize, dative: &[bool]| -> bool {
        let za = z(a);
        if is_metal(za) || matches!(za, 1 | 2 | 9 | 10) {
            return false;
        }
        let atom = mol.atom(AtomIdx(a as u32));
        let eff = i32::from(za) - i32::from(atom.charge);
        if eff <= 0 || eff > 118 {
            return false;
        }
        let max_v = i32::from(MAX_VALENCE[eff as usize]);
        let ev = explicit_valence(a, dative);
        let total_degree = mol.degree(AtomIdx(a as u32))
            + usize::from(chematic_core::implicit_hcount(mol, AtomIdx(a as u32)));
        max_v > 0 && (ev > max_v || (ev == max_v && atom.aromatic && total_degree == 4))
    };
    let single_metal = |a: usize, bi: BondIdx| -> Option<usize> {
        let b = mol.bond(bi);
        let other = if b.atom1.0 as usize == a {
            b.atom2.0 as usize
        } else {
            b.atom1.0 as usize
        };
        (matches!(b.order, BondOrder::Single | BondOrder::Up | BondOrder::Down)
            && is_metal(z(other)))
        .then_some(other)
    };
    let needs = (0..n).any(|a| {
        hypervalent(a, &dative)
            && mol
                .neighbors(AtomIdx(a as u32))
                .any(|(_, bi)| single_metal(a, bi).is_some())
    });
    if !needs {
        return dative;
    }
    // RDKit ranks here before any ring perception (`fastFindRings`); ring
    // membership only breaks ties, so chematic's SSSR stands in for it.
    let rings: Vec<Vec<usize>> = crate::find_sssr(mol)
        .rings()
        .iter()
        .map(|r| r.iter().map(|a| a.0 as usize).collect())
        .collect();
    let ranks = crate::rdkit_canon::rank_with_rings(mol, &rings, &[]);
    let mut by_rank: Vec<usize> = (0..n).collect();
    by_rank.sort_by_key(|&a| ranks[a]);
    for a in by_rank {
        if !hypervalent(a, &dative) {
            continue;
        }
        let mut metals: Vec<(usize, BondIdx)> = mol
            .neighbors(AtomIdx(a as u32))
            .filter_map(|(_, bi)| {
                if dative[bi.0 as usize] {
                    None
                } else {
                    single_metal(a, bi).map(|m| (m, bi))
                }
            })
            .collect();
        if metals.is_empty() {
            continue;
        }
        let n_dative = |m: usize, dative: &[bool]| {
            mol.neighbors(AtomIdx(m as u32))
                .filter(|(_, bi)| dative[bi.0 as usize] || mol.bond(*bi).order == BondOrder::Dative)
                .count()
        };
        metals.sort_by(|x, y| {
            let (dx, dy) = (n_dative(x.0, &dative), n_dative(y.0, &dative));
            dx.cmp(&dy).then(ranks[y.0].cmp(&ranks[x.0]))
        });
        dative[metals[0].1.0 as usize] = true;
    }
    dative
}

impl Graph {
    /// Bonds are numbered as RDKit numbers them ([`Molecule::rdkit_bond_order`]:
    /// a SMILES's ring-closure bonds come last), so neighbour lists and the
    /// searches that walk them visit atoms in RDKit's order.
    fn new(mol: &Molecule) -> Self {
        let mut adj = vec![SmallVec::new(); mol.atom_count()];
        let mut ends = Vec::with_capacity(mol.bond_count());
        let mut orders = Vec::with_capacity(mol.bond_count());
        let dative = organometallic_dative_bonds(mol);
        for bidx in mol.rdkit_bond_order() {
            let bond = mol.bond(bidx);
            let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            let k = ends.len();
            ends.push((a, b));
            orders.push(if dative[bidx.0 as usize] {
                BondOrder::Dative
            } else {
                bond.order
            });
            adj[a].push((b, k));
            adj[b].push((a, k));
        }
        Self { adj, ends, orders }
    }

    fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.adj[a]
            .iter()
            .find(|&&(nb, _)| nb == b)
            .map(|&(_, bi)| bi)
    }
}

/// Ring invariant: the ring's atoms as a sorted list (RDKit uses a bitset).
fn invariant(ring: &[usize]) -> Vec<usize> {
    let mut v = ring.to_vec();
    v.sort_unstable();
    v
}

/// `boost::dynamic_bitset` order of two atom sets of one molecule: compared
/// from the highest atom index down.
fn bitset_cmp(a: &[usize], b: &[usize]) -> std::cmp::Ordering {
    let (mut i, mut j) = (a.len(), b.len());
    while i > 0 && j > 0 {
        let (x, y) = (a[i - 1], b[j - 1]);
        if x != y {
            return x.cmp(&y);
        }
        i -= 1;
        j -= 1;
    }
    i.cmp(&j)
}

#[derive(PartialEq, Eq)]
struct Invariant(Vec<usize>);

impl PartialOrd for Invariant {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Invariant {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        bitset_cmp(&self.0, &other.0)
    }
}

fn trim_bonds(
    g: &Graph,
    cand: usize,
    changed: &mut VecDeque<usize>,
    degrees: &mut [i32],
    active: &mut [bool],
) {
    for &(other, bi) in &g.adj[cand] {
        if !active[bi] {
            continue;
        }
        if degrees[other] <= 2 {
            changed.push_back(other);
        }
        active[bi] = false;
        degrees[other] -= 1;
        degrees[cand] -= 1;
    }
}

struct TooBig;

/// RDKit `smallestRingsBfs`: every smallest ring through `root`.
fn smallest_rings_bfs(
    g: &Graph,
    root: usize,
    active: &[bool],
    forbidden: &[usize],
) -> Result<Vec<Vec<usize>>, TooBig> {
    BFS_SCRATCH.with(|scratch| {
        let mut scratch = scratch.borrow_mut();
        let BfsScratch {
            done,
            parents,
            depths,
            queue,
        } = &mut *scratch;
        smallest_rings_bfs_in(g, root, active, forbidden, done, parents, depths, queue)
    })
}

/// Per-thread buffers for [`smallest_rings_bfs`], which runs once per ring
/// candidate (thousands of times on a large ring system).
#[derive(Default)]
struct BfsScratch {
    done: Vec<u8>,
    parents: Vec<isize>,
    depths: Vec<usize>,
    queue: VecDeque<usize>,
}

thread_local! {
    static BFS_SCRATCH: std::cell::RefCell<BfsScratch> = std::cell::RefCell::new(BfsScratch::default());
}

#[allow(clippy::too_many_arguments)]
fn smallest_rings_bfs_in(
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
                let mut ring = vec![nbr];
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

fn pick_d2_nodes(g: &Graph, frag: &[usize], degrees: &[i32], active: &[bool]) -> Vec<usize> {
    let mut d2 = Vec::new();
    let mut forb = vec![false; g.adj.len()];
    while let Some(&root) = frag.iter().find(|&&a| degrees[a] == 2 && !forb[a]) {
        d2.push(root);
        forb[root] = true;
        // markUselessD2s: degree-two atoms reachable through degree-two
        // atoms are represented by `root`.
        let mut stack = vec![root];
        while let Some(at) = stack.pop() {
            for &(other, bi) in &g.adj[at] {
                if active[bi] && !forb[other] && degrees[other] == 2 {
                    forb[other] = true;
                    stack.push(other);
                }
            }
        }
    }
    d2
}

struct Search<'g> {
    g: &'g Graph,
    invars: HashSet<Vec<usize>>,
    ring_atoms: Vec<bool>,
    ring_bonds: Vec<bool>,
}

impl Search<'_> {
    fn add_new(&mut self, res: &mut Vec<Vec<usize>>, ring: Vec<usize>) -> bool {
        if self.invars.insert(invariant(&ring)) {
            res.push(ring);
            true
        } else {
            false
        }
    }

    fn find_rings_d2_nodes(
        &mut self,
        res: &mut Vec<Vec<usize>>,
        d2nodes: &[usize],
        degrees: &mut [i32],
        active: &mut [bool],
    ) -> Result<(), TooBig> {
        let g = self.g;
        let mut dup_d2_cands: BTreeMap<Invariant, Vec<usize>> = BTreeMap::new();
        let mut dup_map: HashMap<usize, Vec<usize>> = HashMap::new();
        for &cand in d2nodes {
            let srings = smallest_rings_bfs(g, cand, active, &[])?;
            for nring in &srings {
                let invr = invariant(nring);
                let dups = dup_d2_cands.entry(Invariant(invr.clone())).or_default();
                if !self.invars.contains(&invr) {
                    res.push(nring.clone());
                    self.invars.insert(invr);
                    for w in nring.windows(2) {
                        if let Some(bi) = g.bond_between(w[0], w[1]) {
                            self.ring_bonds[bi] = true;
                        }
                        self.ring_atoms[w[0]] = true;
                    }
                    if let Some(bi) = g.bond_between(nring[0], nring[nring.len() - 1]) {
                        self.ring_bonds[bi] = true;
                    }
                    self.ring_atoms[nring[nring.len() - 1]] = true;
                } else {
                    for &other in dups.iter() {
                        dup_map.entry(cand).or_default().push(other);
                        dup_map.entry(other).or_default().push(cand);
                    }
                }
                dups.push(cand);
            }
            if srings.is_empty() {
                let mut changed = VecDeque::from([cand]);
                while let Some(local) = changed.pop_front() {
                    trim_bonds(g, local, &mut changed, degrees, active);
                }
            }
        }
        // findSSSRforDupCands
        for dup_cands in dup_d2_cands.values() {
            if dup_cands.len() <= 1 {
                continue;
            }
            let mut nrings: Vec<Vec<usize>> = Vec::new();
            let mut min_size = usize::MAX;
            for &dup in dup_cands {
                let mut degrees_copy = degrees.to_vec();
                let mut active_copy = active.to_vec();
                let mut changed = VecDeque::new();
                if let Some(others) = dup_map.get(&dup) {
                    for &dni in others {
                        trim_bonds(g, dni, &mut changed, &mut degrees_copy, &mut active_copy);
                    }
                }
                for sring in smallest_rings_bfs(g, dup, &active_copy, &[])? {
                    min_size = min_size.min(sring.len());
                    nrings.push(sring);
                }
            }
            for nring in nrings {
                if nring.len() == min_size {
                    self.add_new(res, nring);
                }
            }
        }
        Ok(())
    }

    fn find_rings_d3_node(
        &mut self,
        res: &mut Vec<Vec<usize>>,
        cand: usize,
        active: &[bool],
    ) -> Result<(), TooBig> {
        let g = self.g;
        let srings = smallest_rings_bfs(g, cand, active, &[])?;
        let nsmall = srings.len();
        for nring in &srings {
            self.add_new(res, nring.clone());
        }
        if nsmall >= 3 {
            return Ok(());
        }
        let nbrs: Vec<usize> = g.adj[cand]
            .iter()
            .filter(|&&(_, bi)| active[bi])
            .map(|&(nb, _)| nb)
            .take(3)
            .collect();
        if nbrs.len() < 3 {
            return Ok(());
        }
        let (n1, n2, n3) = (nbrs[0], nbrs[1], nbrs[2]);
        if nsmall == 2 {
            let in_both = |x: usize| srings[0].contains(&x) && srings[1].contains(&x);
            let f = if in_both(n1) {
                n1
            } else if in_both(n2) {
                n2
            } else if in_both(n3) {
                n3
            } else {
                return Ok(());
            };
            for nring in smallest_rings_bfs(g, cand, active, &[f])? {
                self.add_new(res, nring);
            }
        } else if nsmall == 1 {
            let (f1, f2) = if !srings[0].contains(&n1) {
                (n2, n3)
            } else if !srings[0].contains(&n2) {
                (n1, n3)
            } else if !srings[0].contains(&n3) {
                (n1, n2)
            } else {
                return Ok(());
            };
            for nring in smallest_rings_bfs(g, cand, active, &[f2])? {
                self.add_new(res, nring);
            }
            for nring in smallest_rings_bfs(g, cand, active, &[f1])? {
                self.add_new(res, nring);
            }
        }
        Ok(())
    }
}

/// libstdc++'s `std::sort` (introsort, median-of-three pivot, heap sort past
/// the depth limit, final insertion sort over 16-element runs), which RDKit's
/// Linux builds use: for more than 16 elements the order of equal elements
/// is this algorithm's, not the input's.
fn libstdcxx_sort<T>(v: &mut [T], less: impl Fn(&T, &T) -> bool + Copy) {
    const THRESHOLD: usize = 16;
    fn move_median_to_first<T>(
        v: &mut [T],
        result: usize,
        a: usize,
        b: usize,
        c: usize,
        less: impl Fn(&T, &T) -> bool,
    ) {
        if less(&v[a], &v[b]) {
            if less(&v[b], &v[c]) {
                v.swap(result, b);
            } else if less(&v[a], &v[c]) {
                v.swap(result, c);
            } else {
                v.swap(result, a);
            }
        } else if less(&v[a], &v[c]) {
            v.swap(result, a);
        } else if less(&v[b], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, b);
        }
    }
    fn unguarded_partition<T>(
        v: &mut [T],
        mut first: usize,
        mut last: usize,
        pivot: usize,
        less: impl Fn(&T, &T) -> bool,
    ) -> usize {
        loop {
            while less(&v[first], &v[pivot]) {
                first += 1;
            }
            last -= 1;
            while less(&v[pivot], &v[last]) {
                last -= 1;
            }
            if first >= last {
                return first;
            }
            v.swap(first, last);
            first += 1;
        }
    }
    fn adjust_heap<T>(
        v: &mut [T],
        first: usize,
        mut hole: usize,
        len: usize,
        less: impl Fn(&T, &T) -> bool,
    ) {
        // libstdc++ __adjust_heap followed by __push_heap, on indices.
        let top = hole;
        let mut child = hole;
        while child < (len - 1) / 2 {
            child = 2 * (child + 1);
            if less(&v[first + child], &v[first + child - 1]) {
                child -= 1;
            }
            v.swap(first + hole, first + child);
            hole = child;
        }
        if len.is_multiple_of(2) && child == (len - 2) / 2 {
            child = 2 * (child + 1);
            v.swap(first + hole, first + child - 1);
            hole = child - 1;
        }
        // The value now at `hole` is the one that started at `top`.
        let mut parent = if hole > 0 { (hole - 1) / 2 } else { 0 };
        while hole > top && less(&v[first + parent], &v[first + hole]) {
            v.swap(first + hole, first + parent);
            hole = parent;
            parent = if hole > 0 { (hole - 1) / 2 } else { 0 };
        }
    }
    fn heap_sort<T>(v: &mut [T], first: usize, last: usize, less: impl Fn(&T, &T) -> bool + Copy) {
        let len = last - first;
        if len < 2 {
            return;
        }
        // make_heap
        let mut parent = (len - 2) / 2;
        loop {
            adjust_heap(v, first, parent, len, less);
            if parent == 0 {
                break;
            }
            parent -= 1;
        }
        // sort_heap (pop_heap repeatedly)
        let mut end = last;
        while end - first > 1 {
            end -= 1;
            v.swap(first, end);
            adjust_heap(v, first, 0, end - first, less);
        }
    }
    fn introsort_loop<T>(
        v: &mut [T],
        first: usize,
        mut last: usize,
        mut depth: usize,
        less: impl Fn(&T, &T) -> bool + Copy,
    ) {
        while last - first > THRESHOLD {
            if depth == 0 {
                heap_sort(v, first, last, less);
                return;
            }
            depth -= 1;
            let mid = first + (last - first) / 2;
            move_median_to_first(v, first, first + 1, mid, last - 1, less);
            let cut = unguarded_partition(v, first + 1, last, first, less);
            introsort_loop(v, cut, last, depth, less);
            last = cut;
        }
    }
    fn insertion_sort<T>(v: &mut [T], first: usize, last: usize, less: impl Fn(&T, &T) -> bool) {
        for i in first + 1..last {
            if less(&v[i], &v[first]) {
                v[first..=i].rotate_right(1);
            } else {
                let mut j = i;
                while less(&v[j], &v[j - 1]) {
                    v.swap(j, j - 1);
                    j -= 1;
                }
            }
        }
    }
    let n = v.len();
    if n < 2 {
        return;
    }
    let lg = usize::BITS as usize - 1 - n.leading_zeros() as usize;
    introsort_loop(v, 0, n, 2 * lg, less);
    if n > THRESHOLD {
        insertion_sort(v, 0, THRESHOLD, less);
        for i in THRESHOLD..n {
            let mut j = i;
            while less(&v[j], &v[j - 1]) {
                v.swap(j, j - 1);
                j -= 1;
            }
        }
    } else {
        insertion_sort(v, 0, n, less);
    }
}

fn ring_bond_set(g: &Graph, ring: &[usize]) -> Vec<usize> {
    let mut bonds: Vec<usize> = (0..ring.len())
        .filter_map(|i| g.bond_between(ring[i], ring[(i + 1) % ring.len()]))
        .collect();
    bonds.sort_unstable();
    bonds
}

/// RDKit `removeExtraRings`: keeps, among rings sorted by size, those that
/// add bonds; returns (kept, extras).
fn remove_extra_rings(g: &Graph, mut res: Vec<Vec<usize>>) -> (Vec<Vec<usize>>, Vec<Vec<usize>>) {
    // `std::sort` on size: not stable, so equal sizes keep libstdc++'s order.
    libstdcxx_sort(&mut res, |a, b| a.len() < b.len());
    let brings: Vec<Vec<usize>> = res.iter().map(|r| ring_bond_set(g, r)).collect();
    let nb = g.ends.len();
    let bits: Vec<Vec<bool>> = brings
        .iter()
        .map(|b| {
            let mut v = vec![false; nb];
            for &i in b {
                v[i] = true;
            }
            v
        })
        .collect();
    let subset = |a: &[bool], u: &[bool]| a.iter().zip(u).all(|(&x, &y)| !x || y);
    let mut avail = vec![true; res.len()];
    let mut keep = vec![false; res.len()];
    let mut munion = vec![false; nb];
    for i in 0..res.len() {
        if subset(&bits[i], &munion) {
            avail[i] = false;
        }
        if !avail[i] {
            continue;
        }
        for (u, &x) in munion.iter_mut().zip(&bits[i]) {
            *u |= x;
        }
        keep[i] = true;
        let mut consider: Vec<bool> = (0..res.len())
            .map(|j| j > i && avail[j] && brings[j].len() == brings[i].len())
            .collect();
        while consider.iter().any(|&c| c) {
            let mut best_j = i + 1;
            let mut best_overlap: isize = -1;
            let mut j = i + 1;
            while j < res.len() && brings[j].len() == brings[i].len() {
                if consider[j] && avail[j] {
                    let overlap = bits[j]
                        .iter()
                        .zip(&munion)
                        .filter(|&(&x, &y)| x && y)
                        .count() as isize;
                    if overlap > best_overlap {
                        best_overlap = overlap;
                        best_j = j;
                    }
                }
                j += 1;
            }
            consider[best_j] = false;
            if !subset(&bits[best_j], &munion) {
                keep[best_j] = true;
                for (u, &x) in munion.iter_mut().zip(&bits[best_j]) {
                    *u |= x;
                }
            }
            avail[best_j] = false;
        }
    }
    let mut kept = Vec::new();
    let mut extras = Vec::new();
    for (i, ring) in res.into_iter().enumerate() {
        if keep[i] {
            kept.push(ring);
        } else {
            extras.push(ring);
        }
    }
    (kept, extras)
}

fn ring_eligible(order: BondOrder) -> bool {
    !matches!(order, BondOrder::Zero | BondOrder::Dative)
}

/// The rings RDKit's `SanitizeMol` stores for `mol` (symmetrized SSSR), in
/// RDKit's order, each starting where RDKit starts it. `None` where RDKit
/// itself would fall back to its approximate ring finder (or a search grows
/// past RDKit's limit), so callers can keep their own order.
pub fn rdkit_sssr_ring_order(mol: &Molecule) -> Option<Vec<Vec<AtomIdx>>> {
    let g = Graph::new(mol);
    let n = mol.atom_count();
    let mut active: Vec<bool> = g.orders.iter().map(|&o| ring_eligible(o)).collect();
    let mut degrees: Vec<i32> = (0..n)
        .map(|a| g.adj[a].iter().filter(|&&(_, bi)| active[bi]).count() as i32)
        .collect();
    let degrees_all: Vec<i32> = (0..n).map(|a| g.adj[a].len() as i32).collect();

    // Fragments (all bonds), atoms ascending, ordered by their lowest atom.
    let mut frag_of = vec![usize::MAX; n];
    let mut frags: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if frag_of[start] != usize::MAX {
            continue;
        }
        let id = frags.len();
        let mut members = vec![start];
        frag_of[start] = id;
        let mut i = 0;
        while i < members.len() {
            let a = members[i];
            i += 1;
            for &(nb, _) in &g.adj[a] {
                if frag_of[nb] == usize::MAX {
                    frag_of[nb] = id;
                    members.push(nb);
                }
            }
        }
        members.sort_unstable();
        frags.push(members);
    }

    let mut search = Search {
        g: &g,
        invars: HashSet::new(),
        ring_atoms: vec![false; n],
        ring_bonds: vec![false; g.ends.len()],
    };
    let mut res: Vec<Vec<usize>> = Vec::new();
    let mut extras_all: Vec<Vec<usize>> = Vec::new();
    for frag in &frags {
        if frag.len() < 3 {
            continue;
        }
        let mut changed = VecDeque::new();
        let mut bonds_with_zero: i64 = 0;
        let mut nbnds: i64 = 0;
        for &a in frag {
            bonds_with_zero += i64::from(degrees_all[a]);
            nbnds += i64::from(degrees[a]);
            if degrees[a] < 2 {
                changed.push_back(a);
            }
        }
        let possible = bonds_with_zero / 2 - frag.len() as i64 + 1;
        if possible < 1 {
            continue;
        }
        let nbnds = nbnds / 2;
        let mut done = vec![false; n];
        let mut n_done = 0usize;
        let mut frag_res: Vec<Vec<usize>> = Vec::new();
        while n_done + 3 <= frag.len() {
            while let Some(cand) = changed.pop_front() {
                if !done[cand] {
                    done[cand] = true;
                    n_done += 1;
                    trim_bonds(&g, cand, &mut changed, &mut degrees, &mut active);
                }
            }
            let d2nodes = pick_d2_nodes(&g, frag, &degrees, &active);
            if !d2nodes.is_empty() {
                search
                    .find_rings_d2_nodes(&mut frag_res, &d2nodes, &mut degrees, &mut active)
                    .ok()?;
                for &d2 in &d2nodes {
                    done[d2] = true;
                    n_done += 1;
                    trim_bonds(&g, d2, &mut changed, &mut degrees, &mut active);
                }
            } else if n_done + 3 <= frag.len() {
                let Some(&cand) = frag.iter().find(|&&a| degrees[a] == 3) else {
                    break;
                };
                search
                    .find_rings_d3_node(&mut frag_res, cand, &active)
                    .ok()?;
                done[cand] = true;
                n_done += 1;
                trim_bonds(&g, cand, &mut changed, &mut degrees, &mut active);
            }
        }
        let expected = nbnds - frag.len() as i64 + 1;
        let found = frag_res.len() as i64;
        if found < expected {
            // RDKit's Issue 3514824 search and approximate fallback.
            return None;
        }
        if found > expected {
            let (kept, extras) = remove_extra_rings(&g, frag_res);
            frag_res = kept;
            extras_all.extend(extras);
        }
        res.extend(frag_res);
    }

    // symmetrizeSSSR: an extra ring that can stand in for one SSSR ring of
    // the same size without dropping a bond only that ring provides.
    let bond_rings: Vec<Vec<usize>> = res.iter().map(|r| ring_bond_set(&g, r)).collect();
    let mut bond_counts = vec![0usize; g.ends.len()];
    for r in &bond_rings {
        for &b in r {
            bond_counts[b] += 1;
        }
    }
    let mut out = res.clone();
    for extra in &extras_all {
        let extra_bonds = ring_bond_set(&g, extra);
        for ring in &bond_rings {
            if ring.len() != extra_bonds.len() {
                continue;
            }
            let mut share = false;
            let mut replaces_unique = true;
            for &b in ring {
                if bond_counts[b] == 1 || !share {
                    if extra_bonds.binary_search(&b).is_ok() {
                        share = true;
                    } else if bond_counts[b] == 1 {
                        replaces_unique = false;
                    }
                }
            }
            if share && replaces_unique {
                out.push(extra.clone());
                break;
            }
        }
    }
    Some(
        out.into_iter()
            .map(|r| r.into_iter().map(|a| AtomIdx(a as u32)).collect())
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order(smiles: &str) -> Vec<Vec<u32>> {
        let mol = chematic_smiles::parse(smiles).unwrap();
        rdkit_sssr_ring_order(&mol)
            .unwrap()
            .into_iter()
            .map(|r| r.into_iter().map(|a| a.0).collect())
            .collect()
    }

    #[test]
    fn ring_closure_bonds_follow_rdkit_creation_order() {
        // RDKit closes SMILES rings after parsing, so ring-closure bonds come
        // last in its bond list and its SSSR search meets them last; the
        // order of the last two rings depends on it (RDKit 2026.03.6).
        assert_eq!(
            order("CC(=O)OC1CCC2(C)C3(C1)C=CC4(C5CCC(C(C)=O)C5(C)CC6OC246)C7C3C(=O)OC7=O"),
            vec![
                vec![25, 24, 26],
                vec![15, 14, 21, 17, 16],
                vec![29, 28, 27, 32, 31],
                vec![4, 10, 9, 7, 6, 5],
                vec![11, 12, 13, 27, 28, 9],
                vec![23, 24, 26, 13, 14, 21],
                vec![7, 26, 13, 27, 28, 9],
                vec![11, 12, 13, 26, 7, 9]
            ]
        );
    }

    #[test]
    fn hypervalent_carbanion_metal_bond_is_dative_like_rdkit() {
        // `cleanUpOrganometallics` makes the four-bonded C- to Fe bond dative,
        // so RDKit finds nine rings, one of them four-membered.
        assert_eq!(
            order("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10"),
            vec![
                vec![5, 6, 9],
                vec![8, 7, 9],
                vec![6, 7, 9],
                vec![10, 9, 11],
                vec![10, 9, 14],
                vec![11, 12, 9],
                vec![14, 13, 9],
                vec![9, 12, 13],
                vec![4, 8, 9, 5]
            ]
        );
    }

    #[test]
    fn matches_rdkit_ring_info_order() {
        // RDKit 2026.03.6 `GetRingInfo().AtomRings()` after `MolFromSmiles`.
        assert_eq!(
            order("O=C1NC(=O)c2c1c1c3ccccc3n3c1c1c2c2ccccc2n1[C@H]1CC[C@@H]3O1"),
            vec![
                vec![1, 2, 3, 5, 6],
                vec![9, 10, 11, 12, 13, 8],
                vec![19, 20, 21, 22, 23, 18],
                vec![26, 25, 29, 28, 27],
                vec![29, 28, 14, 15, 16, 24, 25],
                vec![5, 17, 16, 15, 7, 6],
                vec![8, 7, 15, 14, 13],
                vec![18, 17, 16, 24, 23],
            ]
        );
        assert_eq!(
            order("c1ccc2ccccc2c1"),
            vec![vec![0, 9, 8, 3, 2, 1], vec![4, 5, 6, 7, 8, 3]]
        );
    }
}
