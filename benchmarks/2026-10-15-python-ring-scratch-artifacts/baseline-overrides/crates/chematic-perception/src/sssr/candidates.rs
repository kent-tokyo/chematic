/// Horton candidate pool, restricted to the cyclic subgraph and collapsed to
/// the first generation-order occurrence of each geometric cycle.
///
/// This is an allocation-light, output-identical replacement for the
/// historical "BFS from every atom over every ring-eligible bond" loop kept in
/// [`find_sssr_horton_reference`]. The restriction is exact:
///
/// * An atom that lies on no cycle can never be the root of a *simple*
///   candidate cycle (a simple cycle through the root would make it a ring
///   atom), so non-ring roots contribute nothing.
/// * A bridge bond can never lie on a simple cycle, so only cyclic
///   (non-bridge) bonds can close a candidate.
/// * Every shortest path between two atoms of the same cyclic component stays
///   inside the cyclic subgraph, and restricting the adjacency lists (in their
///   original order) preserves BFS discovery order among those atoms, hence the
///   same parent tree. Atoms in a different cyclic component are only
///   reachable through a shared bridge path, which the simplicity check always
///   rejected; here they are simply unreachable.
fn horton_candidates_fast(mol: &Molecule, keep: Option<&[bool]>) -> Vec<CycleCandidate> {
    let graph = CyclicGraph::new(mol);
    let mut search = HortonCandidateSearch::new(mol);
    for root in 0..mol.atom_count() {
        if !graph.ring_atom[root]
            || graph.skip_root[root]
            || keep.is_some_and(|mask| !mask[root])
        {
            continue;
        }
        search.breadth_first(root, &graph);
        let component = graph.component[root] as usize;
        for &(bond, left, right) in graph.component_bonds(component) {
            search.inspect_closing_edge(root as u32, bond, left, right);
        }
    }
    search.candidates
}

struct CyclicGraph {
    adjacency_start: Vec<u32>,
    adjacency: Vec<(u32, BondIdx)>,
    ring_atom: Vec<bool>,
    skip_root: Vec<bool>,
    component: Vec<u32>,
    group_start: Vec<u32>,
    bonds: Vec<(BondIdx, u32, u32)>,
}

impl CyclicGraph {
    fn new(mol: &Molecule) -> Self {
        let cyclic = ring_bond_flags_shared(mol);
        let (adjacency_start, adjacency, ring_atom) = cyclic_adjacency(mol, &cyclic);
        let (component, skip_root, component_count) =
            cyclic_components(&adjacency_start, &adjacency, &ring_atom);
        let (group_start, bonds) = component_bonds(mol, &cyclic, &component, component_count);
        Self {
            adjacency_start,
            adjacency,
            ring_atom,
            skip_root,
            component,
            group_start,
            bonds,
        }
    }

    fn neighbors(&self, atom: usize) -> &[(u32, BondIdx)] {
        &self.adjacency
            [self.adjacency_start[atom] as usize..self.adjacency_start[atom + 1] as usize]
    }

    fn component_bonds(&self, component: usize) -> &[(BondIdx, u32, u32)] {
        &self.bonds[self.group_start[component] as usize..self.group_start[component + 1] as usize]
    }
}

fn cyclic_adjacency(
    mol: &Molecule,
    cyclic: &[bool],
) -> (Vec<u32>, Vec<(u32, BondIdx)>, Vec<bool>) {
    let mut start = vec![0; mol.atom_count() + 1];
    let mut adjacency = Vec::with_capacity(2 * mol.bond_count());
    let mut ring_atom = vec![false; mol.atom_count()];
    for atom in 0..mol.atom_count() {
        start[atom] = adjacency.len() as u32;
        for (neighbor, bond) in mol.neighbors(AtomIdx(atom as u32)) {
            if cyclic[bond.0 as usize] {
                adjacency.push((neighbor.0, bond));
                ring_atom[atom] = true;
            }
        }
    }
    start[mol.atom_count()] = adjacency.len() as u32;
    (start, adjacency, ring_atom)
}

fn cyclic_components(
    start: &[u32],
    adjacency: &[(u32, BondIdx)],
    ring_atom: &[bool],
) -> (Vec<u32>, Vec<bool>, u32) {
    let mut component = vec![u32::MAX; ring_atom.len()];
    let mut skip_root = vec![false; ring_atom.len()];
    let mut labelled = vec![false; ring_atom.len()];
    let mut component_count = 0;
    for first in 0..ring_atom.len() {
        if !ring_atom[first] || labelled[first] {
            continue;
        }
        let mut stack = vec![first as u32];
        let mut members = Vec::new();
        let mut simple = true;
        labelled[first] = true;
        while let Some(current) = stack.pop() {
            let index = current as usize;
            members.push(current);
            component[index] = component_count;
            simple &= start[index + 1] - start[index] == 2;
            for &(neighbor, _) in &adjacency[start[index] as usize..start[index + 1] as usize] {
                if !labelled[neighbor as usize] {
                    labelled[neighbor as usize] = true;
                    stack.push(neighbor);
                }
            }
        }
        if simple {
            for member in members.into_iter().filter(|&member| member as usize != first) {
                skip_root[member as usize] = true;
            }
        }
        component_count += 1;
    }
    (component, skip_root, component_count)
}

fn component_bonds(
    mol: &Molecule,
    cyclic: &[bool],
    component: &[u32],
    component_count: u32,
) -> (Vec<u32>, Vec<(BondIdx, u32, u32)>) {
    let mut start = vec![0; component_count as usize + 1];
    for (index, bond) in mol.bonds() {
        if cyclic[index.0 as usize] {
            start[component[bond.atom1.0 as usize] as usize + 1] += 1;
        }
    }
    for index in 0..component_count as usize {
        start[index + 1] += start[index];
    }
    let mut bonds = vec![(BondIdx(0), 0, 0); start[component_count as usize] as usize];
    let mut fill = start.clone();
    for (index, bond) in mol.bonds() {
        if cyclic[index.0 as usize] {
            let group = component[bond.atom1.0 as usize] as usize;
            bonds[fill[group] as usize] = (index, bond.atom1.0, bond.atom2.0);
            fill[group] += 1;
        }
    }
    (start, bonds)
}

struct HortonCandidateSearch {
    distance: Vec<u32>,
    parent: Vec<u32>,
    parent_bond: Vec<BondIdx>,
    queue: Vec<u32>,
    mark: Vec<u32>,
    stamp: u32,
    seen_bonds: FxHashSet<Vec<BondIdx>>,
    seen_masks: FxHashSet<u128>,
    path_mask: Vec<u128>,
    use_mask: bool,
    candidates: Vec<CycleCandidate>,
}

impl HortonCandidateSearch {
    const UNSEEN: u32 = u32::MAX;

    fn new(mol: &Molecule) -> Self {
        let use_mask = mol.bond_count() <= 128;
        Self {
            distance: vec![Self::UNSEEN; mol.atom_count()],
            parent: vec![Self::UNSEEN; mol.atom_count()],
            parent_bond: vec![BondIdx(0); mol.atom_count()],
            queue: Vec::with_capacity(mol.atom_count()),
            mark: vec![0; mol.atom_count()],
            stamp: 0,
            seen_bonds: FxHashSet::default(),
            seen_masks: FxHashSet::default(),
            path_mask: if use_mask { vec![0; mol.atom_count()] } else { Vec::new() },
            use_mask,
            candidates: Vec::new(),
        }
    }

    fn breadth_first(&mut self, root: usize, graph: &CyclicGraph) {
        for &atom in &self.queue {
            self.distance[atom as usize] = Self::UNSEEN;
            self.parent[atom as usize] = Self::UNSEEN;
        }
        self.queue.clear();
        self.distance[root] = 0;
        if self.use_mask {
            self.path_mask[root] = 0;
        }
        self.queue.push(root as u32);
        let mut head = 0;
        while head < self.queue.len() {
            let current = self.queue[head] as usize;
            head += 1;
            for &(neighbor, bond) in graph.neighbors(current) {
                let neighbor_index = neighbor as usize;
                if self.distance[neighbor_index] != Self::UNSEEN {
                    continue;
                }
                self.distance[neighbor_index] = self.distance[current] + 1;
                self.parent[neighbor_index] = current as u32;
                self.parent_bond[neighbor_index] = bond;
                if self.use_mask {
                    self.path_mask[neighbor_index] = self.path_mask[current] | (1u128 << bond.0);
                }
                self.queue.push(neighbor);
            }
        }
    }

    fn inspect_closing_edge(&mut self, root: u32, bond: BondIdx, left: u32, right: u32) {
        if left == root
            || right == root
            || self.distance[left as usize] == Self::UNSEEN
            || self.distance[right as usize] == Self::UNSEEN
        {
            return;
        }
        if self.use_mask {
            self.add_mask_candidate(bond, left, right);
        } else if self.paths_are_disjoint(root, left, right) {
            self.add_large_candidate(root, bond, left, right);
        }
    }

    fn add_mask_candidate(&mut self, bond: BondIdx, left: u32, right: u32) {
        let (left_path, right_path) = (self.path_mask[left as usize], self.path_mask[right as usize]);
        if left_path & right_path != 0 {
            return;
        }
        let mask = left_path | right_path | (1u128 << bond.0);
        if self.seen_masks.insert(mask) {
            self.candidates.push(CycleCandidate {
                mask,
                bonds: Vec::new(),
                atoms: Vec::new(),
                len: self.distance[left as usize] + self.distance[right as usize] + 1,
                x: left,
                closing: bond,
            });
        }
    }

    fn paths_are_disjoint(&mut self, root: u32, left: u32, mut right: u32) -> bool {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.mark.fill(0);
            self.stamp = 1;
        }
        let mut atom = left;
        while atom != root {
            self.mark[atom as usize] = self.stamp;
            atom = self.parent[atom as usize];
        }
        while right != root {
            if self.mark[right as usize] == self.stamp {
                return false;
            }
            right = self.parent[right as usize];
        }
        true
    }

    fn add_large_candidate(&mut self, root: u32, bond: BondIdx, left: u32, right: u32) {
        let mut bonds = self.path_bonds(root, left);
        bonds.extend(self.path_bonds(root, right));
        bonds.push(bond);
        bonds.sort_unstable();
        if self.seen_bonds.contains(&bonds) {
            return;
        }
        let atoms = self.ordered_ring_atoms(root, left, right);
        self.seen_bonds.insert(bonds.clone());
        self.candidates.push(CycleCandidate {
            mask: 0,
            bonds,
            len: atoms.len() as u32,
            atoms,
            x: left,
            closing: bond,
        });
    }

    fn path_bonds(&self, root: u32, mut atom: u32) -> Vec<BondIdx> {
        let mut bonds = Vec::new();
        while atom != root {
            bonds.push(self.parent_bond[atom as usize]);
            atom = self.parent[atom as usize];
        }
        bonds
    }

    fn ordered_ring_atoms(&self, root: u32, left: u32, mut right: u32) -> Vec<AtomIdx> {
        let mut atoms = Vec::new();
        let mut atom = left;
        loop {
            atoms.push(AtomIdx(atom));
            if atom == root {
                break;
            }
            atom = self.parent[atom as usize];
        }
        let tail = atoms.len();
        while right != root {
            atoms.push(AtomIdx(right));
            right = self.parent[right as usize];
        }
        atoms[tail..].reverse();
        atoms
    }
}

/// Unoptimized reference for [`canonical_atom_ranks`] (differential oracle).
#[doc(hidden)]
pub fn canonical_atom_ranks_reference(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();
    let mut keys: Vec<u64> = (0..n)
        .map(|i| {
            let idx = AtomIdx(i as u32);
            let atom = mol.atom(idx);
            let z = atom.element.atomic_number() as u64;
            // SSSR selection is a property of the cyclic heavy-atom graph.
            // Non-ring-eligible attachments (most importantly explicit H)
            // must not change ranks, otherwise adding/removing explicit H can
            // select a different fused-ring basis for the same heavy graph.
            let degree = mol
                .neighbors(idx)
                .filter(|(_, bond)| is_ring_eligible(mol.bond(*bond).order))
                .count() as u64;
            let charge = (atom.charge as i64 + 8) as u64; // shift to non-negative
            let aromatic = u64::from(atom.aromatic);
            (z << 24) | (degree << 16) | (charge << 8) | aromatic
        })
        .collect();

    const ROUNDS: usize = 8;
    for _ in 0..ROUNDS {
        let mut next = Vec::with_capacity(n);
        for i in 0..n {
            let mut neighbor_keys: Vec<u64> = mol
                .neighbors(AtomIdx(i as u32))
                .filter(|(_, bidx)| is_ring_eligible(mol.bond(*bidx).order))
                .map(|(nb, bidx)| {
                    keys[nb.0 as usize]
                        .wrapping_mul(257)
                        .wrapping_add(mol.bond(bidx).order.order_int() as u64)
                })
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

/// Unoptimized reference for [`canonical_cycle_order_key`] (differential oracle).
#[doc(hidden)]
pub fn canonical_cycle_order_key_reference(
    mol: &Molecule,
    atom_seq: &[AtomIdx],
    ranks: &[u64],
) -> Vec<(u64, u8)> {
    if atom_seq.is_empty() {
        return Vec::new();
    }
    let n = atom_seq.len();
    let mut sequence = Vec::with_capacity(n);
    for i in 0..n {
        let current = atom_seq[i];
        let next = atom_seq[(i + 1) % n];
        let order = mol
            .bond_between(current, next)
            .map(|(_, bond)| bond.order.order_int())
            .unwrap_or_default();
        sequence.push((ranks[current.0 as usize], order));
    }

    let mut best: Option<Vec<(u64, u8)>> = None;
    for reversed in [false, true] {
        for offset in 0..n {
            let candidate = (0..n)
                .map(|step| {
                    let index = if reversed {
                        (offset + n - (step % n)) % n
                    } else {
                        (offset + step) % n
                    };
                    if reversed {
                        // Reverse traversal uses the bond entering the atom,
                        // so preserve the same pair semantics as forward.
                        let previous = (index + n - 1) % n;
                        (ranks[atom_seq[index].0 as usize], sequence[previous].1)
                    } else {
                        sequence[index]
                    }
                })
                .collect::<Vec<_>>();
            if best.as_ref().is_none_or(|current| candidate < *current) {
                best = Some(candidate);
            }
        }
    }
    best.expect("non-empty ring has a canonical traversal")
}

/// The historical, unoptimized Horton SSSR implementation.
///
/// Retained verbatim as the differential-testing oracle for [`find_sssr`];
/// not intended for production use.
#[doc(hidden)]
pub fn find_sssr_horton_reference(mol: &Molecule) -> RingSet {
    let v = mol.atom_count();
    let e = mol
        .bonds()
        .filter(|(_, b)| is_ring_eligible(b.order))
        .count();
    if v == 0 || e == 0 {
        return RingSet::empty();
    }
    let (components, _) = bfs_spanning_forest(mol);
    let r = (e as isize) - (v as isize) + (components as isize);
    if r <= 0 {
        return RingSet::empty();
    }
    let r = r as usize;
    if r == 1 {
        return single_cycle_sssr(mol);
    }
    let ring_bonds: Vec<(BondIdx, AtomIdx, AtomIdx)> = mol
        .bonds()
        .filter(|(_, b)| is_ring_eligible(b.order))
        .map(|(bidx, b)| (bidx, b.atom1, b.atom2))
        .collect();
    let mut candidates: Vec<(Vec<BondIdx>, Vec<AtomIdx>)> = Vec::new();
    for root_idx in 0..v {
        let root = AtomIdx(root_idx as u32);
        let (dist, parent) = bfs_tree(mol, root);
        for &(bidx, x, y) in &ring_bonds {
            if x == root || y == root {
                continue;
            }
            if dist[x.0 as usize] == usize::MAX || dist[y.0 as usize] == usize::MAX {
                continue;
            }
            if let Some(candidate) = horton_candidate(mol, root, x, y, bidx, &parent) {
                candidates.push(candidate);
            }
        }
    }
    let ranks = canonical_atom_ranks_reference(mol);
    candidates.sort_by_cached_key(|c| {
        (
            c.0.len(),
            canonical_cycle_order_key_reference(mol, &c.1, &ranks),
            canonical_cycle_key(&c.1, &ranks),
        )
    });
    candidates.dedup_by(|a, b| a.0 == b.0);
    let mut basis: FxHashMap<BondIdx, Vec<BondIdx>> = FxHashMap::default();
    let mut selected_atoms: Vec<Vec<AtomIdx>> = Vec::new();
    for (bond_set, atom_seq) in candidates {
        let reduced = gf2_reduce(&bond_set, &basis);
        if !reduced.is_empty() {
            let pivot = *reduced.iter().min().unwrap();
            basis.insert(pivot, reduced);
            selected_atoms.push(atom_seq);
            if selected_atoms.len() == r {
                break;
            }
        }
    }
    selected_atoms.sort_by_key(|ring| ring.len());
    RingSet::from_rings(selected_atoms, v)
}

/// Return one flag per atom indicating whether it belongs to any eligible
/// cycle, without constructing an SSSR.
///
/// Morgan connectivity invariants need only RDKit's `isRingAtom` predicate,
/// not the selected rings themselves. An eligible bond is in a cycle exactly
/// when it is not a bridge; therefore both endpoints of every non-bridge bond
/// are ring atoms. This is equivalent to [`find_sssr`]'s atom membership:
/// a cycle basis must contain every non-bridge edge in at least one selected
/// cycle, while bridge edges cannot occur in any cycle.
///
/// The DFS is iterative so a malformed or deliberately very long acyclic
/// graph cannot exhaust the Rust call stack merely while determining that it
/// has no ring atoms.
/// Return an index-aligned flag for every bond that belongs to at least one
/// cycle in the ring-eligible molecular graph.
///
/// This is a linear-time bridge computation, not an SSSR construction.  A
/// ring-eligible edge is cyclic exactly when it is not a bridge, so callers
/// that only need ring membership can avoid the substantially heavier cycle
/// basis machinery.
pub fn ring_bond_flags(mol: &Molecule) -> Vec<bool> {
    (*ring_bond_flags_shared(mol)).clone()
}

/// Shared, memoized [`ring_bond_flags`] (no clone of the flag vector).
pub fn ring_bond_flags_shared(mol: &Molecule) -> std::sync::Arc<Vec<bool>> {
    // Memoized on `mol`; see `chematic_core::DerivedSlot`.
    if let Some(flags) = mol.derived_if_computed(chematic_core::DerivedSlot::RingBondFlags) {
        return flags;
    }
    compute_ring_flags_and_components(mol).0
}

/// Connected components of the cyclic subgraph (the 2-edge-connected
/// components of the ring-eligible graph).
#[derive(Debug)]
#[non_exhaustive]
pub struct RingComponents {
    /// Component label of every atom (atoms on no cycle are singletons).
    pub label: Vec<u32>,
    /// Number of labels.
    pub count: u32,
}

/// Shared, memoized [`RingComponents`] of `mol`, computed together with
/// [`ring_bond_flags_shared`]. Labels identify connected components after
/// removing bridges and ring-ineligible bonds; they are not stable atom IDs.
/// Graph mutations invalidate both cached values. A graph with externally
/// seeded bond flags may compute its component labels on this first call.
pub fn ring_components_shared(mol: &Molecule) -> std::sync::Arc<RingComponents> {
    if let Some(c) = mol.derived_if_computed(chematic_core::DerivedSlot::RingComponents) {
        return c;
    }
    compute_ring_flags_and_components(mol).1
}

/// Compute (once) and memoize both the cyclic-bond flags and the component
/// labels; a slot already filled keeps its value.
fn compute_ring_flags_and_components(
    mol: &Molecule,
) -> (std::sync::Arc<Vec<bool>>, std::sync::Arc<RingComponents>) {
    let (flags, components) = ring_bond_flags_uncached(mol);
    let flags = mol.derived(chematic_core::DerivedSlot::RingBondFlags, || flags);
    let components = mol.derived(chematic_core::DerivedSlot::RingComponents, || components);
    (flags, components)
}

/// Seed `copy`'s cyclic-bond flags and components from `mol`'s (same graph,
/// same ring-eligible bonds).
pub(crate) fn seed_ring_data_from(copy: &Molecule, mol: &Molecule) {
    copy.seed_derived(
        chematic_core::DerivedSlot::RingBondFlags,
        ring_bond_flags_shared(mol),
    );
    copy.seed_derived(
        chematic_core::DerivedSlot::RingComponents,
        ring_components_shared(mol),
    );
}

/// Bridge flags (cyclic = eligible non-bridge bonds) and 2-edge-connected
/// component labels from one iterative low-link DFS: a vertex whose low-link
/// equals its discovery time closes a component (the vertices above it on
/// the visit stack), exactly when the tree edge into it is a bridge.
fn ring_bond_flags_uncached(mol: &Molecule) -> (Vec<bool>, RingComponents) {
    let atom_count = mol.atom_count();
    let bond_count = mol.bond_count();
    let mut flags = vec![false; bond_count];
    if atom_count == 0 || bond_count == 0 {
        return (
            flags,
            RingComponents {
                label: (0..atom_count as u32).collect(),
                count: atom_count as u32,
            },
        );
    }
    // Bridges are a property of the graph, so the traversal order is free:
    // walk the molecule's own adjacency lists, skipping ineligible bonds.
    let all_eligible = mol.bonds().all(|(_, b)| is_ring_eligible(b.order));

    const UNSEEN: u32 = u32::MAX;
    const NO_BOND: u32 = u32::MAX;
    // discovery[i], low[i] interleaved.
    let mut dl = vec![(UNSEEN, UNSEEN); atom_count];
    let mut time = 0u32;
    // `(atom, parent edge, next neighbor offset)` frames.
    let mut stack: Vec<(u32, u32, u32)> = Vec::with_capacity(atom_count.min(64));
    let mut visit: Vec<u32> = Vec::with_capacity(atom_count);
    let mut label = vec![0u32; atom_count];
    let mut count = 0u32;
    // Non-bridge eligible bonds become `true`; bridges stay `false`.
    let mut any_cycle_edge = false;

    for root in 0..atom_count {
        if dl[root].0 != UNSEEN {
            continue;
        }
        dl[root] = (time, time);
        time += 1;
        visit.push(root as u32);
        stack.push((root as u32, NO_BOND, 0));

        while let Some(&mut (atom, parent_bond, ref mut next)) = stack.last_mut() {
            let a = atom as usize;
            let nbrs = mol.neighbor_slice(AtomIdx(atom));
            if (*next as usize) < nbrs.len() {
                let (neighbor, bond) = nbrs[*next as usize];
                *next += 1;
                if bond.0 == parent_bond
                    || (!all_eligible && !is_ring_eligible(mol.bond(bond).order))
                {
                    continue;
                }
                let nb = neighbor.0 as usize;
                if dl[nb].0 == UNSEEN {
                    dl[nb] = (time, time);
                    time += 1;
                    visit.push(neighbor.0);
                    stack.push((neighbor.0, bond.0, 0));
                } else {
                    // Back (or already-finished forward) edge: always on a cycle.
                    if dl[nb].0 < dl[a].1 {
                        dl[a].1 = dl[nb].0;
                    }
                    flags[bond.0 as usize] = true;
                    any_cycle_edge = true;
                }
                continue;
            }
            stack.pop();
            if dl[a].1 == dl[a].0 {
                // `a` roots a 2-edge-connected component.
                while let Some(v) = visit.pop() {
                    label[v as usize] = count;
                    if v == atom {
                        break;
                    }
                }
                count += 1;
            }
            if parent_bond != NO_BOND {
                // Its parent is the preceding stack frame: this DFS only
                // pushes a child immediately after its parent frame.
                let parent = stack
                    .last()
                    .expect("non-root DFS frame must retain its parent")
                    .0 as usize;
                let low_a = dl[a].1;
                if low_a < dl[parent].1 {
                    dl[parent].1 = low_a;
                }
                // Tree edge: a bridge exactly when the child's subtree has
                // no back edge reaching the parent or above.
                if low_a <= dl[parent].0 {
                    flags[parent_bond as usize] = true;
                }
            }
        }
    }
    let _ = any_cycle_edge;
    (flags, RingComponents { label, count })
}

/// Return an index-aligned flag for every atom that belongs to at least one
/// cycle. Uses [`ring_bond_flags`] so this remains linear in graph size.
pub fn ring_atom_flags(mol: &Molecule) -> Vec<bool> {
    if let Some(flags) = single_smiles_cycle_atom_flags(mol) {
        return flags;
    }
    let ring_bonds = ring_bond_flags_shared(mol);
    let mut flags = vec![false; mol.atom_count()];
    for (bond_idx, bond) in mol.bonds() {
        if ring_bonds[bond_idx.0 as usize] {
            flags[bond.atom1.0 as usize] = true;
            flags[bond.atom2.0 as usize] = true;
        }
    }
    flags
}

/// Fast path for an unmodified SMILES graph with exactly one closure bond.
/// Removing that bond leaves a forest, so its endpoints have one unique path;
/// that path plus the closure is the molecule's only cycle. This computes the
/// same atom predicate as the bridge algorithm with fewer allocations on the
/// common monocyclic case. Zero closures is deliberately not a shortcut: it
/// is also how non-SMILES molecules represent unknown provenance.
fn single_smiles_cycle_atom_flags(mol: &Molecule) -> Option<Vec<bool>> {
    if mol.smiles_ring_closure_bond_count() != 1 {
        return None;
    }
    let (closure_idx, closure) = mol
        .bonds()
        .find(|(idx, _)| mol.is_smiles_ring_closure(*idx))?;
    if !is_ring_eligible(closure.order) {
        return None;
    }

    const UNSEEN: u32 = u32::MAX;
    let n = mol.atom_count();
    let start = closure.atom1;
    let end = closure.atom2;
    let mut predecessor: SmallVec<[u32; 64]> = SmallVec::from_elem(UNSEEN, n);
    let mut stack: SmallVec<[AtomIdx; 64]> = SmallVec::new();
    predecessor[start.0 as usize] = start.0;
    stack.push(start);
    while let Some(atom) = stack.pop() {
        if atom == end {
            break;
        }
        for (neighbor, bond_idx) in mol.neighbors(atom) {
            if bond_idx == closure_idx || !is_ring_eligible(mol.bond(bond_idx).order) {
                continue;
            }
            let slot = &mut predecessor[neighbor.0 as usize];
            if *slot == UNSEEN {
                *slot = atom.0;
                stack.push(neighbor);
            }
        }
    }
    if predecessor[end.0 as usize] == UNSEEN {
        return None;
    }

    let mut flags = vec![false; n];
    let mut atom = end.0;
    loop {
        flags[atom as usize] = true;
        if atom == start.0 {
            break;
        }
        atom = predecessor[atom as usize];
    }
    Some(flags)
}

fn single_cycle_sssr(mol: &Molecule) -> RingSet {
    let n = mol.atom_count();
    let mut degree = vec![0usize; n];
    for (_, bond) in mol.bonds() {
        if is_ring_eligible(bond.order) {
            degree[bond.atom1.0 as usize] += 1;
            degree[bond.atom2.0 as usize] += 1;
        }
    }

    let mut queue = VecDeque::new();
    let mut removed = vec![false; n];
    for (idx, &d) in degree.iter().enumerate() {
        if d < 2 {
            queue.push_back(idx);
        }
    }
    while let Some(idx) = queue.pop_front() {
        if removed[idx] {
            continue;
        }
        removed[idx] = true;
        for (neighbor, bond_idx) in mol.neighbors(AtomIdx(idx as u32)) {
            if !is_ring_eligible(mol.bond(bond_idx).order) {
                continue;
            }
            let neighbor = neighbor.0 as usize;
            if !removed[neighbor] {
                degree[neighbor] = degree[neighbor].saturating_sub(1);
                if degree[neighbor] < 2 {
                    queue.push_back(neighbor);
                }
            }
        }
    }

    let Some(start) = (0..n).find(|&idx| !removed[idx]) else {
        return RingSet::empty();
    };
    let mut ring = Vec::new();
    let mut current = start;
    let mut previous = None;
    loop {
        ring.push(AtomIdx(current as u32));
        let next = mol
            .neighbors(AtomIdx(current as u32))
            .filter(|&(_, bond_idx)| is_ring_eligible(mol.bond(bond_idx).order))
            .map(|(neighbor, _)| neighbor.0 as usize)
            .filter(|&neighbor| !removed[neighbor] && Some(neighbor) != previous)
            .find(|&neighbor| neighbor == start || !ring.contains(&AtomIdx(neighbor as u32)));
        let Some(next) = next else {
            return RingSet::empty();
        };
        if next == start {
            break;
        }
        previous = Some(current);
        current = next;
    }
    RingSet::from_rings(vec![ring], n)
}
