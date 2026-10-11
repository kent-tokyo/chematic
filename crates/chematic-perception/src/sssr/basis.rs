/// A Horton candidate cycle: its bond set as a bitmask (molecules with at most
/// 128 bonds) or as a sorted bond list (larger molecules; the mask is then 0),
/// its length, and its ordered ring atoms `x ... root ... y` (the bond
/// `closing` = y-x closes the ring).
///
/// In mask mode the atoms are left empty at generation time and rebuilt on
/// demand by [`CycleCandidate::materialize`]: most candidates are longer than
/// every ring the basis needs and are never looked at again.
struct CycleCandidate {
    mask: u128,
    bonds: Vec<BondIdx>,
    atoms: Vec<AtomIdx>,
    len: u32,
    x: u32,
    closing: BondIdx,
}

impl CycleCandidate {
    /// Fill `atoms` (mask mode). The cycle is simple, so walking it from `x`
    /// away from the closing bond visits exactly `x ... root ... y`, the
    /// order the BFS parent paths give.
    fn materialize(&mut self, mol: &Molecule) {
        if !self.atoms.is_empty() {
            return;
        }
        let mut atoms = Vec::with_capacity(self.len as usize);
        let mut cur = AtomIdx(self.x);
        let mut came = self.closing;
        for _ in 0..self.len {
            atoms.push(cur);
            let Some((nb, b)) = mol
                .neighbors(cur)
                .find(|&(_, b)| b != came && b.0 < 128 && self.mask & (1u128 << b.0) != 0)
            else {
                break;
            };
            came = b;
            cur = nb;
        }
        debug_assert_eq!(atoms.len(), self.len as usize);
        self.atoms = atoms;
    }
}

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
pub struct RingSet {
    rings: Vec<Vec<AtomIdx>>,
    atom_membership: Vec<bool>,
}

impl RingSet {
    fn from_rings(rings: Vec<Vec<AtomIdx>>, atom_count: usize) -> Self {
        let mut atom_membership = vec![false; atom_count];
        for ring in &rings {
            for &atom in ring {
                if let Some(member) = atom_membership.get_mut(atom.0 as usize) {
                    *member = true;
                }
            }
        }
        Self {
            rings,
            atom_membership,
        }
    }

    fn empty() -> Self {
        Self {
            rings: Vec::new(),
            atom_membership: Vec::new(),
        }
    }

    /// All rings as slices of atom indices.
    pub fn rings(&self) -> &[Vec<AtomIdx>] {
        &self.rings
    }

    /// Number of rings in the SSSR.
    pub fn ring_count(&self) -> usize {
        self.rings.len()
    }

    /// Whether atom `atom` is a member of at least one ring.
    pub fn contains_atom(&self, atom: AtomIdx) -> bool {
        self.atom_membership
            .get(atom.0 as usize)
            .copied()
            .unwrap_or(false)
    }

    /// Number of rings that atom `atom` belongs to.
    pub fn atoms_in_ring_count(&self, atom: AtomIdx) -> usize {
        self.rings
            .iter()
            .filter(|ring| ring.contains(&atom))
            .count()
    }
}

/// Outcome of the bounded symmetrized-ring expansion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetrizedSssrStatus {
    /// The candidate search completed without exceeding its safety cap.
    Complete,
    /// Candidate generation exceeded the safety cap and the returned rings
    /// are the complete Horton basis, not a partial symmetrized family.
    CapExhausted,
}

/// Symmetrized SSSR together with the bounded-search outcome.
#[derive(Debug, Clone)]
pub struct SymmetrizedSssrResult {
    ring_set: RingSet,
    status: SymmetrizedSssrStatus,
    candidates_examined: usize,
}

impl SymmetrizedSssrResult {
    /// Return the selected ring set.
    pub fn rings(&self) -> &RingSet {
        &self.ring_set
    }

    /// Return the bounded-search outcome.
    pub fn status(&self) -> SymmetrizedSssrStatus {
        self.status
    }

    /// Number of candidate rings inspected before completion or cap
    /// exhaustion. This is diagnostic evidence for callers that enforce a
    /// caller-owned resource budget; it is never inferred from the cap.
    pub fn candidates_examined(&self) -> usize {
        self.candidates_examined
    }

    /// Consume the result and return its ring set.
    pub fn into_ring_set(self) -> RingSet {
        self.ring_set
    }
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

/// Compute the Smallest Set of Smallest Rings for `mol`.
///
/// Returns a [`RingSet`] whose ring count equals the cycle rank r = E - V + C.
/// For acyclic molecules (r = 0) the returned set is empty.
///
/// The result is memoized on `mol` (see [`find_sssr_shared`]); repeated calls
/// on an unmodified molecule are a cheap clone.
pub fn find_sssr(mol: &Molecule) -> RingSet {
    (*find_sssr_shared(mol)).clone()
}

/// Shared, memoized [`find_sssr`] result. Prefer this in hot paths to avoid
/// cloning the ring vectors.
pub fn find_sssr_shared(mol: &Molecule) -> std::sync::Arc<RingSet> {
    mol.derived(chematic_core::DerivedSlot::Sssr, || find_sssr_uncached(mol))
}

/// Number of rings [`find_sssr`] returns, without constructing them.
///
/// The SSSR is a cycle basis, so its size is the cycle rank `E - V + C` of the
/// ring-eligible graph (zero for acyclic molecules). Linear time.
pub fn sssr_ring_count(mol: &Molecule) -> usize {
    let v = mol.atom_count();
    if v == 0 {
        return 0;
    }
    // Every ring-closure bond written by a valid SMILES adds exactly one to
    // the graph's cycle rank. The parser already records those bonds for
    // RDKit-compatible bond ordering, so parsed cyclic molecules need no
    // union-find pass here. Edits clear the provenance together with the
    // derived caches, making this shortcut safe for mutable molecules.
    let parsed_cycle_rank = mol.smiles_ring_closure_bond_count();
    if parsed_cycle_rank != 0 && mol.bonds().all(|(_, bond)| is_ring_eligible(bond.order)) {
        return parsed_cycle_rank;
    }
    // Acyclic SMILES seed exact all-false flags. Non-SMILES molecules do not,
    // so an empty/missing value cannot accidentally classify them as trees.
    if let Some(flags) =
        mol.derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RingBondFlags)
        && !flags.iter().any(|&is_ring| is_ring)
    {
        return 0;
    }
    // E - V + C is zero for a forest, so no acyclicity pre-check is needed;
    // one union-find pass counts both E and C.
    let mut parent: Vec<u32> = (0..v as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let p = parent[x as usize];
            parent[x as usize] = parent[p as usize];
            x = p;
        }
        x
    }
    let mut cycles = 0usize; // edges joining atoms already connected
    for (_, b) in mol.bonds() {
        if !is_ring_eligible(b.order) {
            continue;
        }
        let a = find(&mut parent, b.atom1.0);
        let c = find(&mut parent, b.atom2.0);
        if a == c {
            cycles += 1;
        } else {
            parent[a as usize] = c;
        }
    }
    cycles
}

fn find_sssr_uncached(mol: &Molecule) -> RingSet {
    let v = mol.atom_count();
    if v == 0 {
        return RingSet::empty();
    }
    // Cycle rank r = E - V + C over the ring-eligible graph. Zero-order and
    // Dative bonds are excluded (RDKit PR #9118). Bridges and acyclic atoms
    // do not change r, so an acyclic graph is detected from the (memoized,
    // linear-time) bridge flags without counting components.
    let cyclic = ring_bond_flags_shared(mol);
    if !cyclic.iter().any(|&c| c) {
        return RingSet::empty();
    }
    let e = mol
        .bonds()
        .filter(|(_, b)| is_ring_eligible(b.order))
        .count();
    let components = ring_eligible_component_count(mol);
    let r = (e as isize) - (v as isize) + (components as isize);

    if r <= 0 {
        return RingSet::empty();
    }
    let r = r as usize;

    // With cycle rank one there is exactly one cycle in the graph. Strip
    // acyclic tails and walk the remaining degree-2 cycle directly. This is
    // equivalent to the Horton basis result, but avoids running a BFS from
    // every atom for the common single-ring case.
    if r == 1 {
        return single_cycle_sssr(mol);
    }

    let candidates = horton_candidates_fast(mol, None);
    let rings = select_horton_basis(mol, &cyclic, candidates, r, true, true)
        .expect("ranked selection always completes");
    RingSet::from_rings(rings, v)
}

/// The rings [`find_sssr`] returns that lie in the cyclic-subgraph components
/// selected by `keep` (one flag per atom, set for every atom of a kept
/// component), in the same relative order, without computing the rings of
/// the other components.
///
/// Horton candidates are cycles through their BFS root, and a root's BFS never
/// leaves its cyclic component, so the candidates of a kept component are
/// generated identically and in the same relative order with or without the
/// other roots. Cycles of different components have disjoint bond sets, so
/// GF(2) independence of a kept component's candidates does not depend on
/// other components' rings; the full selection saturates every component
/// (it stops only at the total cycle rank), and the candidates of a kept
/// component are visited in the same relative order here (same length
/// groups, same canonical keys from the same whole-molecule ranks, stable
/// sorts). Hence exactly the same rings are selected, and the final stable
/// sort by length keeps their relative order.
#[cfg(test)]
pub(crate) fn find_sssr_in_components(mol: &Molecule, keep: &[bool]) -> Vec<Vec<AtomIdx>> {
    sssr_rings_in_components(mol, keep, true, true).expect("ranked selection always completes")
}

/// The *set* of rings [`find_sssr_in_components`] returns, in an unspecified
/// (deterministic) order; cheaper when same-length rings need no canonical
/// ordering (see [`select_horton_basis`]). `keep` may select every component.
pub(crate) fn find_sssr_ring_set_in_components(mol: &Molecule, keep: &[bool]) -> Vec<Vec<AtomIdx>> {
    if let Some(full) = mol.derived_if_computed::<RingSet>(chematic_core::DerivedSlot::Sssr) {
        return full
            .rings()
            .iter()
            .filter(|ring| keep[ring[0].0 as usize])
            .cloned()
            .collect();
    }
    sssr_rings_in_components(mol, keep, false, true).expect("ranked selection always completes")
}

/// [`find_sssr_ring_set_in_components`] when it can be decided from the graph
/// alone (no same-length tie needs the canonical atom ranks, which read
/// aromatic flags and bond orders); `None` otherwise. The set is then the
/// same for every molecule with this graph and adjacency order, e.g. for a
/// Kekulé form of `mol`.
pub(crate) fn find_sssr_ring_set_in_components_unranked(
    mol: &Molecule,
    keep: &[bool],
) -> Option<Vec<Vec<AtomIdx>>> {
    sssr_rings_in_components(mol, keep, false, false)
}

fn sssr_rings_in_components(
    mol: &Molecule,
    keep: &[bool],
    order_needed: bool,
    allow_ranks: bool,
) -> Option<Vec<Vec<AtomIdx>>> {
    let cyclic = ring_bond_flags_shared(mol);
    // Cycle rank of the kept components: E - V + C over their cyclic bonds.
    let n = mol.atom_count();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let p = parent[x as usize];
            parent[x as usize] = parent[p as usize];
            x = p;
        }
        x
    }
    let mut r = 0usize;
    for (bidx, b) in mol.bonds() {
        if !cyclic[bidx.0 as usize] || !keep[b.atom1.0 as usize] {
            continue;
        }
        let a = find(&mut parent, b.atom1.0);
        let c = find(&mut parent, b.atom2.0);
        if a == c {
            r += 1;
        } else {
            parent[a as usize] = c;
        }
    }
    if r == 0 {
        return Some(Vec::new());
    }
    if sssr_ring_count(mol) == 1 {
        // The whole molecule takes the single-cycle path (graph-only); `r ==
        // 1` here means that cycle is in a kept component.
        return Some(if allow_ranks {
            find_sssr_shared(mol).rings().to_vec()
        } else {
            single_cycle_sssr(mol).rings().to_vec()
        });
    }
    if order_needed {
        let candidates = horton_candidates_fast(mol, Some(keep));
        return select_horton_basis(mol, &cyclic, candidates, r, order_needed, allow_ranks);
    }
    // Components whose minimum cycle basis is unique and read off the graph
    // directly (see `unique_basis_of_small_component`); the others go through
    // the Horton selection. Components never share cycles, and per component
    // the selected set does not depend on the other components' candidates.
    let mut edges = vec![0u32; n];
    let mut atoms = vec![0u32; n];
    let mut is_ring_atom = vec![false; n];
    for (bidx, b) in mol.bonds() {
        if cyclic[bidx.0 as usize] && keep[b.atom1.0 as usize] {
            edges[find(&mut parent, b.atom1.0) as usize] += 1;
            is_ring_atom[b.atom1.0 as usize] = true;
            is_ring_atom[b.atom2.0 as usize] = true;
        }
    }
    let mut root_of = vec![u32::MAX; n];
    for a in 0..n {
        if is_ring_atom[a] {
            let root = find(&mut parent, a as u32);
            root_of[a] = root;
            atoms[root as usize] += 1;
        }
    }
    let mut rings: Vec<Vec<AtomIdx>> = Vec::new();
    let mut general = vec![false; n];
    let mut r_rest = 0usize;
    for root in 0..n {
        if root_of[root] != root as u32 {
            continue;
        }
        let rank = (edges[root] + 1).saturating_sub(atoms[root]) as usize;
        let small = rank <= 2
            && unique_basis_of_small_component(mol, &cyclic, &root_of, root as u32, &mut rings);
        if !small {
            general[root] = true;
            r_rest += rank;
        }
    }
    if r_rest > 0 {
        let keep_rest: Vec<bool> = (0..n)
            .map(|a| is_ring_atom[a] && general[root_of[a] as usize])
            .collect();
        let candidates = horton_candidates_fast(mol, Some(&keep_rest));
        let rest = select_horton_basis(mol, &cyclic, candidates, r_rest, false, allow_ranks)?;
        rings.extend(rest);
    }
    Some(rings)
}

/// Appends the minimum cycle basis of the cyclic component rooted at `root`
/// (atoms with `root_of[a] == root`) when its cycle rank is at most two and
/// the basis is unique, returning `true`; returns `false` (appending
/// nothing) otherwise.
///
/// * rank 1: the component is one simple cycle, its only cycle;
/// * rank 2 with one atom of cyclic degree 4 (spiro): two edge-disjoint
///   cycles through it, the only simple cycles there are;
/// * rank 2 with two atoms of cyclic degree 3 (theta graph): three internally
///   disjoint paths of lengths `a <= b <= c` between them, whose pairs form
///   the only simple cycles (`a+b <= a+c <= b+c`). Any two of them are a
///   basis, so the minimum basis is `{a+b, a+c}`, unique exactly when
///   `a < b` (with `a == b` the cycles `a+c` and `b+c` tie).
///
/// A unique minimum basis is what every SSSR selection returns, whatever its
/// tie-break, so the rings do not depend on canonical ranks.
pub(crate) fn unique_basis_of_small_component(
    mol: &Molecule,
    cyclic: &[bool],
    root_of: &[u32],
    root: u32,
    rings: &mut Vec<Vec<AtomIdx>>,
) -> bool {
    let cyclic_nbrs = |a: u32| {
        mol.neighbor_slice(AtomIdx(a))
            .iter()
            .filter(|&&(_, b)| cyclic[b.0 as usize])
            .map(|&(nb, _)| nb.0)
    };
    let members = || (0..root_of.len() as u32).filter(|&a| root_of[a as usize] == root);
    let mut branch: [u32; 2] = [u32::MAX; 2];
    let mut n_branch = 0usize;
    let mut spiro = u32::MAX;
    for a in members() {
        match cyclic_nbrs(a).count() {
            2 => {}
            3 if n_branch < 2 => {
                branch[n_branch] = a;
                n_branch += 1;
            }
            4 if spiro == u32::MAX => spiro = a,
            _ => return false,
        }
    }
    // Walk from `start` through `first` along degree-2 atoms until an atom
    // with `stop(atom)` is reached; returns the interior atoms and the end.
    let walk = |start: u32, first: u32, stop: &dyn Fn(u32) -> bool| -> (Vec<AtomIdx>, u32) {
        let mut interior = Vec::new();
        let (mut prev, mut cur) = (start, first);
        while !stop(cur) {
            interior.push(AtomIdx(cur));
            let next = cyclic_nbrs(cur)
                .find(|&x| x != prev)
                .expect("degree-2 cycle atom");
            prev = cur;
            cur = next;
        }
        (interior, prev)
    };
    match (n_branch, spiro != u32::MAX) {
        (0, false) => {
            // One simple cycle.
            let Some(start) = members().next() else {
                return false;
            };
            let first = cyclic_nbrs(start)
                .next()
                .expect("cycle atom has neighbours");
            let (interior, _) = walk(start, first, &|x| x == start);
            let mut ring = Vec::with_capacity(interior.len() + 1);
            ring.push(AtomIdx(start));
            ring.extend(interior);
            rings.push(ring);
            true
        }
        (0, true) => {
            let mut used: Vec<u32> = Vec::with_capacity(4);
            let nbrs: Vec<u32> = cyclic_nbrs(spiro).collect();
            for first in nbrs {
                if used.contains(&first) {
                    continue;
                }
                let (interior, last) = walk(spiro, first, &|x| x == spiro);
                used.push(first);
                used.push(last);
                let mut ring = Vec::with_capacity(interior.len() + 1);
                ring.push(AtomIdx(spiro));
                ring.extend(interior);
                rings.push(ring);
            }
            true
        }
        (2, false) => {
            let (u, v) = (branch[0], branch[1]);
            let mut paths: Vec<Vec<AtomIdx>> = cyclic_nbrs(u)
                .map(|first| {
                    if first == v {
                        Vec::new()
                    } else {
                        walk(u, first, &|x| x == v).0
                    }
                })
                .collect();
            paths.sort_by_key(|p| p.len()); // stable
            if paths.len() != 3 || paths[0].len() == paths[1].len() {
                return false;
            }
            for other in [&paths[1], &paths[2]] {
                let mut ring = Vec::with_capacity(paths[0].len() + other.len() + 2);
                ring.push(AtomIdx(u));
                ring.extend(paths[0].iter().copied());
                ring.push(AtomIdx(v));
                ring.extend(other.iter().rev().copied());
                rings.push(ring);
            }
            true
        }
        _ => false,
    }
}

/// Greedy GF(2) selection of `r` independent candidates in canonical order,
/// returned sorted by length (stable).
///
/// With `order_needed == false` only the returned *set* of rings is
/// guaranteed (it is the same set). A length group whose candidates are all
/// independent of the basis and of each other is then taken whole without
/// the canonical sort (greedy selection takes every candidate of such a group
/// in any order), as is a group of which none is independent; any other group
/// is sorted canonically as usual.
///
/// With `allow_ranks == false` (only meaningful with `order_needed ==
/// false`), a group that would need the canonical sort makes the selection
/// return `None` instead: the result then never depends on the canonical
/// atom ranks, i.e. on anything but the graph.
fn select_horton_basis(
    mol: &Molecule,
    cyclic: &[bool],
    mut candidates: Vec<CycleCandidate>,
    r: usize,
    order_needed: bool,
    allow_ranks: bool,
) -> Option<Vec<Vec<AtomIdx>>> {
    // Deterministic ordering: shortest first, then a canonical (input-order-
    // independent) tie-break so ring *selection* doesn't depend on how the
    // molecule happened to be numbered by the parser.
    //
    // `horton_candidates_fast` already collapsed geometric duplicates to their
    // first generation-order occurrence. Duplicates share an identical sort
    // key (both keys are rotation/reflection invariant), and
    // `sort_by_cached_key` is stable, so the relative order of first
    // occurrences — the only thing the GF(2) selection below can observe —
    // is exactly the order the historical sort+adjacent-dedup produced.
    //
    // Candidates are processed one ring length at a time: the global sort key
    // is (length, order key, multiset key), so sorting each length group by
    // the two canonical keys (stably, in generation order) visits exactly the
    // same sequence, while the canonical keys of lengths beyond the point
    // where `r` independent cycles are found are never computed. A group with
    // a single candidate needs no keys at all, and the canonical atom ranks
    // are only computed once some group actually has to be sorted.
    candidates.sort_by_key(|c| c.len); // stable: keeps generation order
    let mut ranks: Option<Vec<u64>> = None;

    // Gaussian elimination over GF(2) to select r linearly independent cycles.
    // Rows are bitsets over bond indices when they fit in 128 bits (the pivot
    // is the lowest set bit == the minimum BondIdx, as in the Vec form).
    let use_bits = mol.bond_count() <= 128;
    let mut basis: FxHashMap<BondIdx, Vec<BondIdx>> = FxHashMap::default();
    let mut basis_bits: [u128; 128] = [0; 128];
    let mut selected_atoms: Vec<Vec<AtomIdx>> = Vec::with_capacity(r);

    let mut rest: &mut [CycleCandidate] = &mut candidates;
    'lengths: while !rest.is_empty() {
        let len = rest[0].len;
        let group_len = rest.iter().take_while(|c| c.len == len).count();
        let (group, tail) = std::mem::take(&mut rest).split_at_mut(group_len);
        rest = tail;
        if group.len() > 1 && !order_needed && use_bits {
            match take_independent_group(mol, group, &mut basis_bits, &mut selected_atoms, r) {
                GroupShortcut::Complete => break 'lengths,
                GroupShortcut::Taken => continue,
                GroupShortcut::Dependent => continue,
                GroupShortcut::NeedsOrdering => {}
            }
        }
        if group.len() > 1 && !allow_ranks {
            return None;
        }
        order_candidate_group(mol, cyclic, group, &mut ranks);
        for candidate in group.iter_mut() {
            if candidate_is_independent(candidate, use_bits, &mut basis_bits, &mut basis) {
                candidate.materialize(mol);
                selected_atoms.push(std::mem::take(&mut candidate.atoms));
                if selected_atoms.len() == r {
                    break 'lengths;
                }
            }
        }
    }

    // Sort output rings by length for output consistency.
    selected_atoms.sort_by_key(|ring| ring.len());
    Some(selected_atoms)
}

enum GroupShortcut {
    Complete,
    Taken,
    Dependent,
    NeedsOrdering,
}

fn take_independent_group(
    mol: &Molecule,
    group: &mut [CycleCandidate],
    basis: &mut [u128; 128],
    selected: &mut Vec<Vec<AtomIdx>>,
    target: usize,
) -> GroupShortcut {
    let mut filled = [0u8; 128];
    let mut filled_count = 0;
    for candidate in group.iter() {
        let mut row = candidate.mask;
        while row != 0 {
            let pivot = row.trailing_zeros() as usize;
            if basis[pivot] == 0 {
                basis[pivot] = row;
                filled[filled_count] = pivot as u8;
                filled_count += 1;
                break;
            }
            row ^= basis[pivot];
        }
    }
    if filled_count == group.len() {
        for candidate in group {
            candidate.materialize(mol);
            selected.push(std::mem::take(&mut candidate.atoms));
        }
        return if selected.len() == target {
            GroupShortcut::Complete
        } else {
            GroupShortcut::Taken
        };
    }
    if filled_count == 0 {
        return GroupShortcut::Dependent;
    }
    for &pivot in &filled[..filled_count] {
        basis[pivot as usize] = 0;
    }
    GroupShortcut::NeedsOrdering
}

fn order_candidate_group(
    mol: &Molecule,
    cyclic: &[bool],
    group: &mut [CycleCandidate],
    ranks: &mut Option<Vec<u64>>,
) {
    if group.len() <= 1 {
        return;
    }
    let ranks = ranks.get_or_insert_with(|| canonical_ring_atom_ranks(mol, cyclic));
    for candidate in group.iter_mut() {
        candidate.materialize(mol);
    }
    group.sort_by_cached_key(|candidate| {
        (
            canonical_cycle_order_key(mol, &candidate.atoms, ranks),
            canonical_cycle_key(&candidate.atoms, ranks),
        )
    });
}

fn candidate_is_independent(
    candidate: &CycleCandidate,
    use_bits: bool,
    bit_basis: &mut [u128; 128],
    bond_basis: &mut FxHashMap<BondIdx, Vec<BondIdx>>,
) -> bool {
    if !use_bits {
        let reduced = gf2_reduce(&candidate.bonds, bond_basis);
        if reduced.is_empty() {
            return false;
        }
        let pivot = *reduced.iter().min().expect("non-empty reduced cycle");
        bond_basis.insert(pivot, reduced);
        return true;
    }
    let mut row = candidate.mask;
    loop {
        if row == 0 {
            return false;
        }
        let pivot = row.trailing_zeros() as usize;
        if bit_basis[pivot] == 0 {
            bit_basis[pivot] = row;
            return true;
        }
        row ^= bit_basis[pivot];
    }
}

/// Number of connected components of the ring-eligible graph (every atom,
/// including isolated ones, is its own component). Union-find; equivalent to
/// the component count of [`bfs_spanning_forest`].
fn ring_eligible_component_count(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let p = parent[x as usize];
            parent[x as usize] = parent[p as usize];
            x = p;
        }
        x
    }
    let mut components = n;
    for (_, b) in mol.bonds() {
        if !is_ring_eligible(b.order) {
            continue;
        }
        let a = find(&mut parent, b.atom1.0);
        let c = find(&mut parent, b.atom2.0);
        if a != c {
            parent[a as usize] = c;
            components -= 1;
        }
    }
    components
}
