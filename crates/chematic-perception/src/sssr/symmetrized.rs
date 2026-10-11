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

/// Select D2 roots for the symmetrized-ring pass without using the first raw
/// atom in a component as a chemical tie-break.  All minimum-rank roots are
/// retained when a component is genuinely symmetric; callers must therefore
/// treat this as a candidate set, not as a single representative.
fn select_permutation_invariant_d2_roots(mol: &Molecule) -> Vec<AtomIdx> {
    let ranks = canonical_atom_ranks(mol);
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
        let mut stack = vec![AtomIdx(raw as u32)];
        let mut component = Vec::new();
        seen[raw] = true;
        while let Some(atom) = stack.pop() {
            component.push(atom);
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
        let min_rank = component
            .iter()
            .map(|atom| ranks[atom.0 as usize])
            .min()
            .expect("degree-2 component is non-empty");
        roots.extend(
            component
                .into_iter()
                .filter(|atom| ranks[atom.0 as usize] == min_rank),
        );
    }
    roots.sort_by_key(|atom| (ranks[atom.0 as usize], atom.0));
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
pub fn find_symmetrized_sssr_with_diagnostics(mol: &Molecule) -> SymmetrizedSssrResult {
    find_symmetrized_sssr_with_diagnostics_bounded(mol, None)
}

/// Bounded variant of [`find_symmetrized_sssr_with_diagnostics`]. When
/// `max_candidates` is `Some`, candidate generation stops as soon as the
/// measured number of inspected candidate rings exceeds that cap and returns
/// the complete Horton basis with [`SymmetrizedSssrStatus::CapExhausted`].
/// The optional bound is intentionally separate from the historical
/// extra-ring cap, so callers can propagate their own resource budget without
/// changing the default perception behavior.
pub fn find_symmetrized_sssr_with_diagnostics_bounded(
    mol: &Molecule,
    max_candidates: Option<usize>,
) -> SymmetrizedSssrResult {
    // Explicit hydrogen nodes are not part of the cyclic heavy-atom graph,
    // but they can perturb the symmetrized candidate search through molecule
    // size, atom metadata, and tie-breaking. Normalize them away at this
    // boundary, then map the resulting rings back to the caller's indices.
    // This makes SymmSSSR invariant under compact-vs-explicit-H representation
    // while preserving the public index contract.
    if mol.atoms().any(|(_, atom)| atom.element == Element::H) {
        return symmetrized_heavy_atom_view(mol, max_candidates);
    }
    let base = find_sssr(mol);
    if base.rings().is_empty() {
        return SymmetrizedSssrResult {
            ring_set: base,
            status: SymmetrizedSssrStatus::Complete,
            candidates_examined: 0,
        };
    }

    let d2_roots = select_permutation_invariant_d2_roots(mol);
    let mut expansion = SymmetrizedExpansion::new(mol, base, max_candidates);
    if d2_roots.is_empty() {
        expansion.collect_all_roots();
    } else {
        expansion.collect_duplicate_roots(&d2_roots);
    }
    expansion.finish()
}

const MAX_SYMMETRIZED_EXTRA_RINGS: usize = 256;

fn symmetrized_heavy_atom_view(
    mol: &Molecule,
    max_candidates: Option<usize>,
) -> SymmetrizedSssrResult {
    let (heavy_mol, heavy_to_original) = heavy_atom_ring_view(mol);
    let result = find_symmetrized_sssr_with_diagnostics_bounded(&heavy_mol, max_candidates);
    let rings = result
        .ring_set
        .rings
        .into_iter()
        .map(|ring| {
            ring.into_iter()
                .map(|atom| heavy_to_original[atom.0 as usize])
                .collect()
        })
        .collect();
    SymmetrizedSssrResult {
        ring_set: RingSet::from_rings(rings, mol.atom_count()),
        status: result.status,
        candidates_examined: result.candidates_examined,
    }
}

struct SymmetrizedExpansion<'a> {
    mol: &'a Molecule,
    base: RingSet,
    base_bonds: Vec<FxHashSet<BondIdx>>,
    bond_ring_count: FxHashMap<BondIdx, usize>,
    ranks: Vec<u64>,
    base_keys: FxHashSet<Vec<u32>>,
    seen: FxHashSet<Vec<u32>>,
    rings: Vec<Vec<AtomIdx>>,
    direct_replacements: Vec<Vec<AtomIdx>>,
    max_candidates: Option<usize>,
    candidates_examined: usize,
    cap_exhausted: bool,
}

impl<'a> SymmetrizedExpansion<'a> {
    fn new(mol: &'a Molecule, base: RingSet, max_candidates: Option<usize>) -> Self {
        let base_bonds: Vec<_> = base
            .rings()
            .iter()
            .map(|ring| ring_bond_set(mol, ring))
            .collect();
        let mut bond_ring_count = FxHashMap::default();
        for ring in &base_bonds {
            for &bond in ring {
                *bond_ring_count.entry(bond).or_insert(0) += 1;
            }
        }
        let base_keys: FxHashSet<_> = base_bonds.iter().map(bond_set_key).collect();
        Self {
            mol,
            ranks: canonical_atom_ranks(mol),
            seen: base_keys.clone(),
            rings: base.rings().to_vec(),
            base,
            base_bonds,
            bond_ring_count,
            base_keys,
            direct_replacements: Vec::new(),
            max_candidates,
            candidates_examined: 0,
            cap_exhausted: false,
        }
    }

    fn collect_all_roots(&mut self) {
        for root in 0..self.mol.atom_count() {
            for candidate in find_smallest_rings_bfs(self.mol, AtomIdx(root as u32)) {
                if !self.inspect_candidate() {
                    return;
                }
                self.accept_candidate(candidate);
            }
        }
    }

    fn collect_duplicate_roots(&mut self, roots: &[AtomIdx]) {
        let mut groups: FxHashMap<Vec<u32>, (Vec<AtomIdx>, Vec<AtomIdx>)> = FxHashMap::default();
        for &root in roots {
            for candidate in find_smallest_rings_bfs(self.mol, root) {
                if !self.inspect_candidate() {
                    return;
                }
                let key = bond_set_key(&ring_bond_set(self.mol, &candidate));
                let entry = groups
                    .entry(key)
                    .or_insert_with(|| (candidate.clone(), Vec::new()));
                if !entry.1.contains(&root) {
                    entry.1.push(root);
                }
            }
        }
        let mut groups: Vec<_> = groups.into_iter().collect();
        groups.sort_by(|(left, _), (right, _)| self.compare_group_keys(left, right));
        for (_, (candidate, mut duplicate_roots)) in groups {
            duplicate_roots.sort_by_key(|root| (self.ranks[root.0 as usize], root.0));
            if duplicate_roots.len() < 2 {
                self.accept_candidate(candidate);
            } else {
                self.collect_replacements(&duplicate_roots);
            }
            if self.cap_exhausted {
                return;
            }
        }
    }

    fn collect_replacements(&mut self, roots: &[AtomIdx]) {
        let mut replacements = Vec::new();
        for &root in roots {
            let blocked = self.blocked_other_roots(root, roots);
            let trimmed = trim_ring_bonds(self.mol, &blocked);
            for replacement in find_smallest_rings_bfs_with_rdkit_tree(self.mol, root, &trimmed) {
                if !self.inspect_candidate() {
                    return;
                }
                replacements.push(replacement);
            }
        }
        if let Some(minimum) = replacements.iter().map(Vec::len).min() {
            replacements.retain(|candidate| candidate.len() == minimum);
        }
        replacements.sort_by_key(|candidate| {
            (
                canonical_cycle_order_key(self.mol, candidate, &self.ranks),
                canonical_bond_set_key(
                    self.mol,
                    &ring_bond_set(self.mol, candidate),
                    &self.ranks,
                ),
            )
        });
        self.direct_replacements.extend(replacements);
    }

    fn blocked_other_roots(&self, root: AtomIdx, roots: &[AtomIdx]) -> FxHashSet<BondIdx> {
        let mut blocked = FxHashSet::default();
        for &other in roots {
            if other == root {
                continue;
            }
            for (_, bond) in self.mol.neighbors(other) {
                if is_ring_eligible(self.mol.bond(bond).order) {
                    blocked.insert(bond);
                }
            }
        }
        blocked
    }

    fn inspect_candidate(&mut self) -> bool {
        self.candidates_examined = self.candidates_examined.saturating_add(1);
        if self
            .max_candidates
            .is_some_and(|cap| self.candidates_examined > cap)
        {
            self.cap_exhausted = true;
            false
        } else {
            true
        }
    }

    fn accept_candidate(&mut self, candidate: Vec<AtomIdx>) -> bool {
        if self.rings.len().saturating_sub(self.base.ring_count())
            >= MAX_SYMMETRIZED_EXTRA_RINGS
        {
            self.cap_exhausted = true;
            return false;
        }
        let bonds = ring_bond_set(self.mol, &candidate);
        let key = bond_set_key(&bonds);
        if self.base_keys.contains(&key) || !self.seen.insert(key) {
            return false;
        }
        let shares_required_bonds = self.base_bonds.iter().any(|basis| {
            basis.iter().any(|bond| bonds.contains(bond))
                && basis.iter().all(|bond| {
                    self.bond_ring_count.get(bond).copied().unwrap_or(0) != 1
                        || bonds.contains(bond)
                })
        });
        if shares_required_bonds
            && self.preserves_basis_independence(&bonds)
            && self.same_size_as_base(&candidate)
        {
            self.rings.push(candidate);
            true
        } else {
            false
        }
    }

    fn preserves_basis_independence(&self, candidate: &FxHashSet<BondIdx>) -> bool {
        self.base_bonds
            .iter()
            .enumerate()
            .any(|(replace_index, base_ring)| {
                base_ring.len() == candidate.len()
                    && gf2_rank(
                        &self
                            .base_bonds
                            .iter()
                            .enumerate()
                            .map(|(index, ring)| {
                                bond_set_key(if index == replace_index { candidate } else { ring })
                            })
                            .collect::<Vec<_>>(),
                    ) == self.base_bonds.len()
            })
    }

    fn same_size_as_base(&self, candidate: &[AtomIdx]) -> bool {
        self.base
            .rings()
            .iter()
            .any(|ring| ring.len() == candidate.len())
    }

    fn compare_group_keys(&self, left: &[u32], right: &[u32]) -> std::cmp::Ordering {
        let canonical = |key: &[u32]| {
            canonical_bond_set_key(
                self.mol,
                &key.iter().map(|&raw| BondIdx(raw)).collect(),
                &self.ranks,
            )
        };
        canonical(left)
            .cmp(&canonical(right))
            .then_with(|| left.cmp(right))
    }

    fn finish(mut self) -> SymmetrizedSssrResult {
        if self.cap_exhausted
            || self
                .rings
                .len()
                .saturating_sub(self.base.ring_count())
                .saturating_add(self.direct_replacements.len())
                > MAX_SYMMETRIZED_EXTRA_RINGS
        {
            return self.cap_result();
        }
        self.add_direct_replacements();
        self.remove_replaced_macrocycle();
        self.sort_rings();
        SymmetrizedSssrResult {
            ring_set: RingSet::from_rings(self.rings, self.mol.atom_count()),
            status: SymmetrizedSssrStatus::Complete,
            candidates_examined: self.candidates_examined,
        }
    }

    fn add_direct_replacements(&mut self) {
        for replacement in std::mem::take(&mut self.direct_replacements) {
            let bonds = ring_bond_set(self.mol, &replacement);
            let key = bond_set_key(&bonds);
            if self.preserves_basis_independence(&bonds)
                && self.seen.insert(key)
                && self.same_size_as_base(&replacement)
            {
                self.rings.push(replacement);
            }
        }
    }

    fn remove_replaced_macrocycle(&mut self) {
        let macrocycles: Vec<_> = self
            .base
            .rings()
            .iter()
            .enumerate()
            .filter(|(_, ring)| ring.len() >= 20)
            .collect();
        if macrocycles.len() != 1 {
            return;
        }
        let (base_index, base_ring) = macrocycles[0];
        let matching_extras = self
            .rings
            .iter()
            .skip(self.base.ring_count())
            .filter(|ring| ring.len() == base_ring.len())
            .count();
        let has_double = ring_bond_set(self.mol, base_ring)
            .iter()
            .any(|&bond| self.mol.bond(bond).order == BondOrder::Double);
        if matching_extras == 3 && has_double {
            self.rings.remove(base_index);
        }
    }

    fn sort_rings(&mut self) {
        let mut extras = self.rings.split_off(self.base.ring_count());
        extras.sort_by_key(|ring| {
            basis_exchange_key(self.mol, ring, &self.base_bonds, &self.ranks)
        });
        self.rings.extend(extras);
        self.rings.sort_by_cached_key(|ring| {
            (
                ring.len(),
                canonical_cycle_order_key(self.mol, ring, &self.ranks),
                canonical_cycle_key(ring, &self.ranks),
            )
        });
    }

    fn cap_result(self) -> SymmetrizedSssrResult {
        SymmetrizedSssrResult {
            ring_set: self.base,
            status: SymmetrizedSssrStatus::CapExhausted,
            candidates_examined: self.candidates_examined,
        }
    }
}

/// Build the graph used by ring perception while retaining a mapping to the
/// original atom indices. Hydrogen neighbors are represented in the retained
/// heavy atom's hydrogen count so aromatic kekulization still sees the same
/// donor information as a compact representation.
fn heavy_atom_ring_view(mol: &Molecule) -> (Molecule, Vec<AtomIdx>) {
    let mut builder = MoleculeBuilder::new();
    let mut original_to_heavy = vec![None; mol.atom_count()];
    let mut heavy_to_original = Vec::new();

    for (original_idx, atom) in mol.atoms() {
        if atom.element == Element::H {
            continue;
        }
        let mut retained = atom.clone();
        let explicit_h_count = mol
            .neighbors(original_idx)
            .filter(|(neighbor, _)| mol.atom(*neighbor).element == Element::H)
            .count();
        if explicit_h_count > 0 {
            retained.hydrogen_count = Some(explicit_h_count as u8);
        } else if retained.hydrogen_count == Some(0) {
            // `add_hydrogens` uses Some(0) as a sentinel after moving all
            // implicit H to explicit nodes. Restore the compact form for the
            // ring-only copy; no real bracket-H information was lost above.
            retained.hydrogen_count = None;
        }
        let heavy_idx = builder.add_atom(retained);
        original_to_heavy[original_idx.0 as usize] = Some(heavy_idx);
        heavy_to_original.push(original_idx);
    }

    for (_, bond) in mol.bonds() {
        let (Some(a), Some(b)) = (
            original_to_heavy[bond.atom1.0 as usize],
            original_to_heavy[bond.atom2.0 as usize],
        ) else {
            continue;
        };
        builder
            .add_bond(a, b, bond.order)
            .expect("heavy-atom ring view must preserve unique bonds");
    }
    (builder.build(), heavy_to_original)
}

/// Compute the bounded symmetrized SSSR, preserving the historical ring-only
/// API. Call [`find_symmetrized_sssr_with_diagnostics`] when callers need to
/// distinguish a complete family from a fail-closed Horton-basis fallback.
pub fn find_symmetrized_sssr(mol: &Molecule) -> RingSet {
    find_symmetrized_sssr_with_diagnostics(mol).into_ring_set()
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
    ranks: &[u64],
) -> Vec<Vec<(u64, u64, u8)>> {
    let candidate_set = ring_bond_set(mol, candidate);
    let candidate_key = canonical_bond_set_key(mol, &candidate_set, ranks);
    let candidate_len = candidate.len();
    let mut best: Option<Vec<Vec<(u64, u64, u8)>>> = None;
    for (replace_idx, base_ring) in base_bonds.iter().enumerate() {
        if base_ring.len() != candidate_len {
            continue;
        }
        let rank_rows = base_bonds
            .iter()
            .enumerate()
            .map(|(idx, set)| {
                let source = if idx == replace_idx {
                    &candidate_set
                } else {
                    set
                };
                bond_set_key(source)
            })
            .collect::<Vec<_>>();
        if gf2_rank(&rank_rows) != base_bonds.len() {
            continue;
        }
        let mut rows = base_bonds
            .iter()
            .enumerate()
            .map(|(idx, set)| {
                if idx == replace_idx {
                    candidate_key.clone()
                } else {
                    canonical_bond_set_key(mol, set, ranks)
                }
            })
            .collect::<Vec<_>>();
        rows.sort_unstable();
        if best.as_ref().is_none_or(|current| rows < *current) {
            best = Some(rows);
        }
    }
    best.unwrap_or_else(|| vec![candidate_key])
}

/// Canonical, input-order-independent representation of a bond set.
///
/// This is used only for ordering alternative basis exchanges. Raw bond
/// indices are intentionally excluded: they are storage artifacts and would
/// make symmetrized-ring selection depend on atom/bond insertion order.
fn canonical_bond_set_key(
    mol: &Molecule,
    bonds: &FxHashSet<BondIdx>,
    ranks: &[u64],
) -> Vec<(u64, u64, u8)> {
    let mut key = bonds
        .iter()
        .map(|&bidx| {
            let bond = mol.bond(bidx);
            let a = ranks[bond.atom1.0 as usize];
            let b = ranks[bond.atom2.0 as usize];
            (a.min(b), a.max(b), bond.order.order_int())
        })
        .collect::<Vec<_>>();
    key.sort_unstable();
    key
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
