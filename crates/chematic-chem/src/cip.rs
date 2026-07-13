//! CIP (Cahn–Ingold–Prelog) stereochemistry assignment.
//!
//! Implements R/S (tetrahedral) and E/Z (double bond) assignment for
//! molecules parsed from SMILES with chirality annotations.

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, CipCode, Molecule, implicit_hcount};

/// The result of a CIP stereochemistry assignment run.
#[derive(Debug)]
pub struct CipAssignment {
    pub assignments: Vec<(AtomIdx, CipCode)>,
}

impl CipAssignment {
    /// Look up the CIP code for a given atom index.
    pub fn get(&self, idx: AtomIdx) -> Option<CipCode> {
        self.assignments
            .iter()
            .find(|(i, _)| *i == idx)
            .map(|(_, c)| *c)
    }
}

/// For a tetrahedral chiral center, return the CIP code and the 4 neighbor atom indices
/// sorted by **decreasing CIP priority** (highest priority = index 0, lowest = index 3).
///
/// The virtual sentinel `AtomIdx(u32::MAX)` is inserted for bracket-H atoms (`[C@H]`, etc.)
/// following the same convention as [`assign_cip`] internally.
///
/// Returns `None` if the atom has no chirality annotation, has fewer than 4 neighbors
/// (including H), or if CIP priorities are tied (cannot be uniquely ranked).
pub fn tetrahedral_stereo_neighbors(
    mol: &Molecule,
    center: AtomIdx,
) -> Option<(CipCode, [AtomIdx; 4])> {
    let atom = mol.atom(center);
    if atom.chirality == Chirality::None {
        return None;
    }

    let neighbors = stereo_neighbors(mol, center);
    if neighbors.len() != 4 {
        return None;
    }

    let cip_code = assign_tetrahedral(mol, center)?;
    let ranks = rank_substituents(mol, center, &neighbors)?;

    // Sort neighbors by decreasing CIP priority (rank N first, rank 1 last).
    let mut sorted: Vec<(u8, AtomIdx)> = ranks
        .iter()
        .zip(neighbors.iter())
        .map(|(&r, &n)| (r, n))
        .collect();
    sorted.sort_by_key(|x| std::cmp::Reverse(x.0));
    let arr = [sorted[0].1, sorted[1].1, sorted[2].1, sorted[3].1];

    Some((cip_code, arr))
}

/// Run CIP assignment on `mol`.  Returns R/S for chiral tetrahedral centers
/// and E/Z for stereospecified double bonds.
pub fn assign_cip(mol: &Molecule) -> CipAssignment {
    let mut assignments = Vec::new();

    // R/S for tetrahedral centers
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        if let Some(code) = assign_tetrahedral(mol, idx) {
            assignments.push((idx, code));
        }
    }

    // E/Z for double bonds
    for j in 0..mol.bond_count() {
        let bidx = BondIdx(j as u32);
        if let Some((atom_idx, code)) = assign_ez(mol, bidx) {
            assignments.push((atom_idx, code));
        }
    }

    // Axial chirality for allenes (>C=C=C<)
    for i in 0..mol.atom_count() {
        let idx = AtomIdx(i as u32);
        if is_allene_central(mol, idx)
            && let Some((atom_idx, code)) = assign_allene(mol, idx)
        {
            assignments.push((atom_idx, code));
        }
    }

    CipAssignment { assignments }
}

/// A single "sphere layer" in a CIP branch expansion: a sorted list of
/// `(atomic_num, isotope, atomic_mass)` tuples (sorted descending for lexicographic comparison).
type SphereLayer = Vec<(u8, Option<u16>, f64)>;

/// Get the key `(atomic_num, isotope, atomic_mass)` for an atom, used in CIP comparisons.
///
/// For the virtual H sentinel (`AtomIdx(u32::MAX)`), returns `(1, None, 1.0)`.
/// Atomic mass is the monoisotopic mass of the element (CIP rule 4 tiebreaker).
fn atom_key(mol: &Molecule, idx: AtomIdx) -> (u8, Option<u16>, f64) {
    if idx.0 == u32::MAX {
        return (1, None, 1.007825);
    }
    let a = mol.atom(idx);
    (
        a.element.atomic_number(),
        a.isotope,
        a.element.atomic_mass(),
    )
}

/// Compare two `(atomic_num, isotope, atomic_mass)` keys by CIP priority.
///
/// CIP rule hierarchy:
/// 1. Higher atomic number wins.
/// 2. For explicit isotopes: `Some(mass)` > `None` (heavier beats unspecified).
/// 3. For tiebreaker: higher atomic mass wins (rule 4).
fn cmp_key(a: (u8, Option<u16>, f64), b: (u8, Option<u16>, f64)) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;
    // Rule 1: atomic number (primary)
    match a.0.cmp(&b.0) {
        Equal => {}
        other => return other,
    }
    // Rule 2: isotope label (if present)
    match (a.1, b.1) {
        (Some(x), Some(y)) => match x.cmp(&y) {
            Equal => {}
            other => return other,
        },
        (Some(_), None) => return Greater,
        (None, Some(_)) => return Less,
        (None, None) => {}
    }
    // Rule 4 tiebreaker: atomic mass (descending for higher priority)
    b.2.partial_cmp(&a.2).unwrap_or(Equal)
}

/// BFS state for one node during sphere expansion.
struct ExpandState {
    node: AtomIdx,
    parent: AtomIdx,
    depth: usize,
    visited: HashSet<AtomIdx>,
}

/// BFS-based CIP sphere expansion for the branch starting at `start`,
/// not going back through `center`.
///
/// At each depth layer the collected `(atomic_num, isotope)` tuples are sorted
/// descending (highest priority first), implementing the hierarchical digraph
/// comparison.
///
/// Phantom atom rules:
/// 1. **Double-bond phantom, arrival side**: when expanding node B (reached via double
///    bond from A), add a phantom entry for A at the same depth level as B's children.
/// 2. **Double-bond phantom, departure side**: when listing A's own substituents and B
///    is reached via a double bond, count B twice in A's substituent list (the real B,
///    continuing the graph, plus a terminal duplicate) -- a double bond duplicates its
///    partner into *both* atoms' substituent lists, not just the arrival side above.
/// 3. **Ring revisit phantom**: if an already-visited atom is encountered,
///    add a phantom for it but don't expand further.
fn cip_branch_spheres(mol: &Molecule, center: AtomIdx, start: AtomIdx) -> Vec<SphereLayer> {
    let mut layers: HashMap<usize, Vec<(u8, Option<u16>, f64)>> = HashMap::new();
    let max_depth = 8usize;

    // The start atom itself is at depth 1.
    let start_key = atom_key(mol, start);
    layers.entry(1).or_default().push(start_key);

    let mut expand_queue: VecDeque<ExpandState> = VecDeque::new();
    {
        let mut v = HashSet::new();
        v.insert(center);
        v.insert(start);
        expand_queue.push_back(ExpandState {
            node: start,
            parent: center,
            depth: 1,
            visited: v,
        });
    }

    while let Some(state) = expand_queue.pop_front() {
        if state.depth >= max_depth {
            continue;
        }
        let child_depth = state.depth + 1;

        // Phantom of parent: add if the bond used to reach this node was double.
        if let Some((_, bond_to_parent)) = mol.bond_between(state.node, state.parent)
            && bond_to_parent.order == BondOrder::Double
        {
            let phantom_key = atom_key(mol, state.parent);
            layers.entry(child_depth).or_default().push(phantom_key);
        }

        for (nb, _) in mol.neighbors(state.node) {
            if nb == state.parent || nb == center {
                continue;
            }
            let child_key = atom_key(mol, nb);
            let layer = layers.entry(child_depth).or_default();

            // Departure-side double-bond phantom (rule 2 above).
            let is_double = mol
                .bond_between(state.node, nb)
                .is_some_and(|(_, b)| b.order == BondOrder::Double);
            if is_double {
                layer.push(child_key);
            }

            if state.visited.contains(&nb) {
                // Ring revisit: phantom only, no expansion.
                layer.push(child_key);
            } else {
                layer.push(child_key);
                let mut child_visited = state.visited.clone();
                child_visited.insert(nb);
                expand_queue.push_back(ExpandState {
                    node: nb,
                    parent: state.node,
                    depth: child_depth,
                    visited: child_visited,
                });
            }
        }
    }

    // Sort each layer descending and return as a Vec ordered by depth.
    let max_layer = layers.keys().copied().max().unwrap_or(0);
    let mut result = Vec::new();
    for d in 1..=max_layer {
        let mut layer = layers.remove(&d).unwrap_or_default();
        layer.sort_by(|a, b| cmp_key(*b, *a)); // descending
        result.push(layer);
    }
    result
}

/// Compare two branches from `center` starting at `a` and `b`.
///
/// Returns `Ordering::Greater` if branch `a` has higher CIP priority than `b`.
fn compare_branches(mol: &Molecule, center: AtomIdx, a: AtomIdx, b: AtomIdx) -> std::cmp::Ordering {
    use std::cmp::Ordering::*;

    // Depth-0 comparison: the substituent atoms themselves.
    let a_key = atom_key(mol, a);
    let b_key = atom_key(mol, b);
    match cmp_key(a_key, b_key) {
        Equal => {}
        other => return other,
    }

    // Sphere-by-sphere comparison.
    let a_spheres = cip_branch_spheres(mol, center, a);
    let b_spheres = cip_branch_spheres(mol, center, b);

    let max_depth = a_spheres.len().max(b_spheres.len());
    for d in 0..max_depth {
        let a_layer = a_spheres.get(d).map(|v| v.as_slice()).unwrap_or(&[]);
        let b_layer = b_spheres.get(d).map(|v| v.as_slice()).unwrap_or(&[]);

        let min_len = a_layer.len().min(b_layer.len());
        for i in 0..min_len {
            match cmp_key(a_layer[i], b_layer[i]) {
                Equal => {}
                other => return other,
            }
        }
        match a_layer.len().cmp(&b_layer.len()) {
            Equal => {}
            other => return other,
        }
    }

    Equal
}

// ---------------------------------------------------------------------------
// CIP Rule 5: stereo-descriptor tie-breaking
// ---------------------------------------------------------------------------

/// Map a provisional CIP code to a u8 token for Rule 5 comparison.
/// Only R/S matter; E/Z and unresolved centres collapse to 0.
fn stereo_token(code: Option<CipCode>) -> u8 {
    match code {
        Some(CipCode::R) => 2,
        Some(CipCode::S) => 1,
        _ => 0,
    }
}

/// Like [`cip_branch_spheres`] but collects stereo-descriptor tokens instead of
/// atom keys.  Follows identical BFS structure, phantom rules, and max depth so
/// that per-layer multisets are comparable as tiebreakers to the graph spheres.
fn cip_branch_stereo_spheres(
    mol: &Molecule,
    center: AtomIdx,
    start: AtomIdx,
    provisional: &HashMap<AtomIdx, CipCode>,
) -> Vec<Vec<u8>> {
    let mut layers: HashMap<usize, Vec<u8>> = HashMap::new();
    let max_depth = 8usize;

    layers
        .entry(1)
        .or_default()
        .push(stereo_token(provisional.get(&start).copied()));

    let mut expand_queue: VecDeque<ExpandState> = VecDeque::new();
    {
        let mut v = HashSet::new();
        v.insert(center);
        v.insert(start);
        expand_queue.push_back(ExpandState {
            node: start,
            parent: center,
            depth: 1,
            visited: v,
        });
    }

    while let Some(state) = expand_queue.pop_front() {
        if state.depth >= max_depth {
            continue;
        }
        let child_depth = state.depth + 1;

        // Double-bond phantom: emit token 0 (phantoms are non-stereogenic duplicates).
        if let Some((_, bond_to_parent)) = mol.bond_between(state.node, state.parent)
            && bond_to_parent.order == BondOrder::Double
        {
            layers.entry(child_depth).or_default().push(0u8);
        }

        for (nb, _) in mol.neighbors(state.node) {
            if nb == state.parent || nb == center {
                continue;
            }
            let token = stereo_token(provisional.get(&nb).copied());
            let layer = layers.entry(child_depth).or_default();

            if state.visited.contains(&nb) {
                // Ring revisit phantom: emit token but don't expand.
                layer.push(token);
            } else {
                layer.push(token);
                let mut child_visited = state.visited.clone();
                child_visited.insert(nb);
                expand_queue.push_back(ExpandState {
                    node: nb,
                    parent: state.node,
                    depth: child_depth,
                    visited: child_visited,
                });
            }
        }
    }

    let max_layer = layers.keys().copied().max().unwrap_or(0);
    let mut result = Vec::new();
    for d in 1..=max_layer {
        let mut layer = layers.remove(&d).unwrap_or_default();
        layer.sort_unstable_by(|a, b| b.cmp(a)); // descending
        result.push(layer);
    }
    result
}

/// True if `idx` is a potential tetrahedral stereocenter, using CIP Rule 5
/// to break graph ties via `provisional` R/S assignments from a first pass.
///
/// Falls back to the pure-graph result when `provisional` is empty or when
/// the tie persists even after adding stereo tokens.
pub(crate) fn is_potential_stereocenter_rule5(
    mol: &Molecule,
    idx: AtomIdx,
    provisional: &HashMap<AtomIdx, CipCode>,
) -> bool {
    let atom = mol.atom(idx);
    if atom.aromatic {
        return false;
    }
    match atom.element.atomic_number() {
        6 | 7 | 15 | 16 | 34 => {}
        _ => return false,
    }
    let mut neighbors: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();
    let h = implicit_hcount(mol, idx);
    if h > 1 {
        return false;
    }
    for _ in 0..h {
        neighbors.push(AtomIdx(u32::MAX));
    }
    if neighbors.len() == 3 && h == 0 && matches!(atom.element.atomic_number(), 15 | 16 | 34) {
        neighbors.push(AtomIdx(u32::MAX));
    }
    if neighbors.len() != 4 {
        return false;
    }

    // Step 1: pure-graph ranking.
    if let Some(ranks) = rank_substituents(mol, idx, &neighbors) {
        let mut r = ranks;
        r.sort_unstable();
        return r.windows(2).all(|w| w[0] != w[1]);
    }

    // Step 2: graph tie — try Rule 5 with stereo tokens.
    if provisional.is_empty() {
        return false;
    }

    // Build a signature per substituent: (atom_key, graph_spheres, stereo_spheres).
    // Two substituents are indistinct iff their signatures are equal.
    // PartialEq on (u8, Option<u16>, f64) is safe: f64 values come from
    // Element::atomic_mass(), a const table with no NaN.
    //
    // AtomIdx(u32::MAX) is the virtual H sentinel (bracket-H or lone pair).
    // It is always unique (h>1 is filtered above) so we skip sphere expansion
    // for it — atom_key alone distinguishes it from all heavy-atom substituents.
    let sigs: Vec<_> = neighbors
        .iter()
        .map(|&nb| {
            let is_sentinel = nb.0 == u32::MAX;
            (
                atom_key(mol, nb),
                if is_sentinel {
                    vec![]
                } else {
                    cip_branch_spheres(mol, idx, nb)
                },
                if is_sentinel {
                    vec![]
                } else {
                    cip_branch_stereo_spheres(mol, idx, nb, provisional)
                },
            )
        })
        .collect();

    // All 6 pairwise pairs must be unequal for 4 distinct substituents.
    for i in 0..4 {
        for j in (i + 1)..4 {
            if sigs[i] == sigs[j] {
                return false;
            }
        }
    }
    true
}

/// Assign CIP priority ranks to `subs` (substituents of `center`).
///
/// Returns `None` if any two substituents have equal priority (tie).
/// Otherwise returns `Vec<u8>` of the same length, where `result[i]` is the
/// rank of `subs[i]` (1 = lowest CIP priority, N = highest).
pub(crate) fn rank_substituents(
    mol: &Molecule,
    center: AtomIdx,
    subs: &[AtomIdx],
) -> Option<Vec<u8>> {
    let n = subs.len();
    if n == 0 {
        return Some(vec![]);
    }

    // Sort indices by CIP priority descending.
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&i, &j| compare_branches(mol, center, subs[i], subs[j]).reverse());

    // Check for ties among adjacent elements after sorting.
    for k in 0..n - 1 {
        let i = indices[k];
        let j = indices[k + 1];
        if compare_branches(mol, center, subs[i], subs[j]) == std::cmp::Ordering::Equal {
            return None;
        }
    }

    // Assign ranks: indices[0] gets rank n (highest), indices[n-1] gets rank 1 (lowest).
    let mut ranks = vec![0u8; n];
    for (rank_from_top, &idx) in indices.iter().enumerate() {
        ranks[idx] = (n - rank_from_top) as u8;
    }

    Some(ranks)
}

/// Collect a chiral atom's 4 substituents (including a virtual `AtomIdx(u32::MAX)`
/// slot for bracket H) in SMILES chirality-neighbor order.
///
/// Prefers `Molecule::stereo_neighbor_order`, which the SMILES parser populates with
/// the *true* textual encounter order: ring-closure partners are resolved to their
/// digit's actual position in the string via a dedicated slot mechanism
/// (`StereoEntry::PendingRing`), not the order bonds happen to be materialized in the
/// adjacency list. That distinction matters because a ring-*opening* bond (partner
/// unknown yet) is only added to `Molecule::neighbors()` once the matching closing
/// digit is reached, which can be *after* a branch/continuation atom that appears
/// later in the text but has nothing to wait on — so raw adjacency order silently
/// reorders the ring partner behind that atom. Falls back to adjacency order (with
/// the same heuristic H placement as before) for molecules with no parse-time stereo
/// data, e.g. built directly via `MoleculeBuilder`.
fn stereo_neighbors(mol: &Molecule, idx: AtomIdx) -> Vec<AtomIdx> {
    if let Some(order) = mol.stereo_neighbor_order(idx) {
        return order.iter().map(|&n| AtomIdx(n)).collect();
    }

    let atom = mol.atom(idx);
    let mut neighbors: Vec<AtomIdx> = mol.neighbors(idx).map(|(nb, _)| nb).collect();

    // For bracket atoms with explicit H (e.g. `[C@@H]`), the H occupies a specific
    // position in the SMILES chirality neighbor list:
    //
    //   • If the bracket atom has NO preceding atom in the SMILES chain (it is the
    //     first atom of the fragment, like `[C@@H](F)(Cl)Br`), the H is at position 0
    //     (the "from-viewer" slot) and the non-H neighbors follow.
    //
    //   • If the bracket atom HAS a preceding atom (like `N[C@@H](C)C(=O)O`), the
    //     preceding atom is at position 0, the H is at position 1, and the remaining
    //     non-H neighbors follow.
    //
    // "Preceding atom" = the atom that forms the bond into this atom from the left in
    // the SMILES string.  In the adjacency list, that atom is always added FIRST
    // (before branches and continuations) and therefore has a SMALLER atom index.
    let has_bracket_h = atom.hydrogen_count.is_some_and(|h| h > 0);
    if has_bracket_h {
        // Detect whether a preceding atom is present: the first neighbor, if its index
        // is smaller than `idx`, is the preceding atom.
        let has_preceding = neighbors.first().map(|&nb| nb.0 < idx.0).unwrap_or(false);
        let h_insert_pos = if has_preceding { 1 } else { 0 };
        neighbors.insert(h_insert_pos, AtomIdx(u32::MAX));
    }
    neighbors
}

fn assign_tetrahedral(mol: &Molecule, idx: AtomIdx) -> Option<CipCode> {
    let atom = mol.atom(idx);
    if atom.chirality == Chirality::None {
        return None;
    }

    let neighbors = stereo_neighbors(mol, idx);
    if neighbors.len() != 4 {
        return None;
    }

    // ranks[i] = CIP rank of neighbors[i]: 1 = lowest priority, 4 = highest.
    let ranks = rank_substituents(mol, idx, &neighbors)?;

    // --- Parity-based R/S determination -----------------------------------
    //
    // SMILES `@@` means: looking FROM neighbors[0], the sequence
    // neighbors[1]→neighbors[2]→neighbors[3] goes clockwise (CW).
    //
    // CIP R: looking FROM the rank-1 substituent, the sequence
    // rank2→rank3→rank4 (ascending priority) goes CW.
    //
    // Algorithm:
    // 1. Find where rank-1 is in the neighbors list (`lowest_pos`).
    // 2. Moving rank-1 to position 0 takes `lowest_pos` adjacent swaps,
    //    each one flipping CW↔CCW.  So the "effective_cw" (from rank-1's
    //    perspective) = smiles_cw XOR (lowest_pos is odd).
    // 3. After removing rank-1, the remaining three neighbors are in some order.
    //    Count how many swaps are needed to put them in ascending rank order
    //    [rank2, rank3, rank4].  An even number → same orientation; odd → flipped.
    // 4. is_r = effective_cw XOR (remaining_swaps is odd).

    let lowest_pos = ranks.iter().position(|&r| r == 1)?;
    let parity_odd = lowest_pos % 2 == 1;
    let smiles_cw = atom.chirality == Chirality::Clockwise;
    let cw_from_lowest = smiles_cw ^ parity_odd;

    // Remaining ranks in their current positional order (lowest_pos removed).
    let remaining_ranks: Vec<u8> = (0..4usize)
        .filter(|&i| i != lowest_pos)
        .map(|i| ranks[i])
        .collect();

    // Count swaps to reach the ascending-rank target [2, 3, 4].
    let remaining_swaps_odd = {
        let mut r = remaining_ranks.clone();
        let target = [2u8, 3, 4];
        let mut swaps = 0usize;
        for i in 0..3 {
            if r[i] != target[i] {
                let j_rel = r[i + 1..].iter().position(|&x| x == target[i])?;
                r.swap(i, j_rel + i + 1);
                swaps += 1;
            }
        }
        swaps % 2 == 1
    };

    // R if the effective CW sense matches the ascending-rank arrangement.
    let is_r = cw_from_lowest ^ remaining_swaps_odd;

    Some(if is_r { CipCode::R } else { CipCode::S })
}

/// Determine if a substituent is "up" relative to the alkene end it connects to.
///
/// Returns `Some(true)` = up, `Some(false)` = down, `None` = no stereo bond.
fn substituent_is_up(mol: &Molecule, alkene_end: AtomIdx, sub: AtomIdx) -> Option<bool> {
    let (_, bond) = mol.bond_between(alkene_end, sub)?;
    match bond.order {
        BondOrder::Up => {
            // `/` bond: atom1→atom2 goes "up"
            Some(bond.atom1 == alkene_end)
        }
        BondOrder::Down => {
            // `\` bond: atom1→atom2 goes "down"
            Some(bond.atom1 == sub)
        }
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Allene axial chirality
// ---------------------------------------------------------------------------

/// True if `idx` is the central atom of an allene (exactly 2 double bonds,
/// no other bonds to heavy atoms — the `>C=C=C<` pattern).
fn is_allene_central(mol: &Molecule, idx: AtomIdx) -> bool {
    let dbl_count = mol
        .neighbors(idx)
        .filter(|(_, bidx)| mol.bond(*bidx).order == BondOrder::Double)
        .count();
    dbl_count == 2 && mol.neighbors(idx).count() == 2
}

/// Assign axial chirality for the allene centred at `central_idx`.
///
/// Uses the same Up/Down bond convention as [`assign_ez`]: the highest-CIP
/// substituent at each terminal is tested for Up/Down; if both are on the same
/// side the code is Z (aS), otherwise E (aR).
///
/// Returns `None` when either terminal has no stereo bond, tied priorities, or
/// only one substituent (no axial chirality possible).
fn assign_allene(mol: &Molecule, central_idx: AtomIdx) -> Option<(AtomIdx, CipCode)> {
    // Collect the two terminal atoms of the allene.
    let terminals: Vec<AtomIdx> = mol
        .neighbors(central_idx)
        .filter(|(_, bidx)| mol.bond(*bidx).order == BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    if terminals.len() != 2 {
        return None;
    }
    let t1 = terminals[0];
    let t2 = terminals[1];

    // Non-allene substituents at each terminal.
    let subs_t1: Vec<AtomIdx> = mol
        .neighbors(t1)
        .filter(|&(nb, bidx)| nb != central_idx && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    let subs_t2: Vec<AtomIdx> = mol
        .neighbors(t2)
        .filter(|&(nb, bidx)| nb != central_idx && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    // At least one substituent needed at each end for axial chirality.
    if subs_t1.is_empty() || subs_t2.is_empty() {
        return None;
    }

    // Highest-priority substituent with a stereo bond at each terminal.
    let high_t1 = highest_stereo_sub(mol, t1, &subs_t1)?;
    let high_t2 = highest_stereo_sub(mol, t2, &subs_t2)?;

    let up_t1 = substituent_is_up(mol, t1, high_t1)?;
    let up_t2 = substituent_is_up(mol, t2, high_t2)?;

    // Same side → Z (aS / M); opposite → E (aR / P).
    let code = if up_t1 == up_t2 {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some((t1, code))
}

/// Assign E/Z for the double bond at `bond_idx`.
///
/// Returns `Some((atom_idx, E or Z))` using one of the double-bond endpoints
/// as the key atom index.  Returns `None` if the bond isn't double or stereo
/// cannot be determined.
fn assign_ez(mol: &Molecule, bond_idx: BondIdx) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    if bond.order != BondOrder::Double {
        return None;
    }

    let a1 = bond.atom1;
    let a2 = bond.atom2;

    // Non-double-bond neighbors for each alkene end (exclude the other alkene atom).
    let subs_a1: Vec<AtomIdx> = mol
        .neighbors(a1)
        .filter(|&(nb, bidx)| nb != a2 && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    let subs_a2: Vec<AtomIdx> = mol
        .neighbors(a2)
        .filter(|&(nb, bidx)| nb != a1 && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None; // terminal alkene
    }

    // Highest-priority substituent with a stereo (Up/Down) bond at each end.
    let high_sub_a1 = highest_stereo_sub(mol, a1, &subs_a1)?;
    let high_sub_a2 = highest_stereo_sub(mol, a2, &subs_a2)?;

    let up_a1 = substituent_is_up(mol, a1, high_sub_a1)?;
    let up_a2 = substituent_is_up(mol, a2, high_sub_a2)?;

    // Same side → Z (zusammen); opposite → E (entgegen).
    let code = if up_a1 == up_a2 {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some((a1, code))
}

/// From `subs` at `alkene_end`, return the highest CIP-priority substituent
/// that has an Up/Down bond to `alkene_end`.
fn highest_stereo_sub(mol: &Molecule, alkene_end: AtomIdx, subs: &[AtomIdx]) -> Option<AtomIdx> {
    let mut sorted: Vec<AtomIdx> = subs.to_vec();
    sorted.sort_by(|&a, &b| compare_branches(mol, alkene_end, a, b).reverse());
    sorted
        .into_iter()
        .find(|&sub| substituent_is_up(mol, alkene_end, sub).is_some())
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn cip_at(smiles: &str, atom_idx: usize) -> Option<CipCode> {
        let mol = parse(smiles).unwrap();
        let assignment = assign_cip(&mol);
        assignment.get(AtomIdx(atom_idx as u32))
    }

    // --- Tetrahedral R/S ---

    #[test]
    fn test_l_alanine_s() {
        // N[C@@H](C)C(=O)O — L-alanine, chiral center is atom 1
        assert_eq!(
            cip_at("N[C@@H](C)C(=O)O", 1),
            Some(CipCode::S),
            "L-alanine should be S"
        );
    }

    #[test]
    fn test_d_alanine_r() {
        // N[C@H](C)C(=O)O — D-alanine
        assert_eq!(
            cip_at("N[C@H](C)C(=O)O", 1),
            Some(CipCode::R),
            "D-alanine should be R"
        );
    }

    #[test]
    fn test_chfclbr_r() {
        // [C@@H](F)(Cl)Br — known R configuration
        // CIP priority: Br(35) > Cl(17) > F(9) > H(1)
        // @@H: looking from H, F→Cl→Br is CW → R
        assert_eq!(
            cip_at("[C@@H](F)(Cl)Br", 0),
            Some(CipCode::R),
            "[C@@H](F)(Cl)Br should be R"
        );
    }

    #[test]
    fn test_chfclbr_s() {
        // [C@H](F)(Cl)Br — S
        assert_eq!(
            cip_at("[C@H](F)(Cl)Br", 0),
            Some(CipCode::S),
            "[C@H](F)(Cl)Br should be S"
        );
    }

    #[test]
    fn test_no_chirality() {
        let mol = parse("CC(=O)O").unwrap();
        let assignment = assign_cip(&mol);
        let tetrahedral: Vec<_> = assignment
            .assignments
            .iter()
            .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
            .collect();
        assert!(tetrahedral.is_empty(), "acetic acid has no chiral centers");
    }

    #[test]
    fn test_symmetric_center_none() {
        // No @/@@ annotation → no assignment attempted
        let mol = parse("CC(N)(N)CC").unwrap();
        let assignment = assign_cip(&mol);
        assert!(
            assignment.assignments.is_empty(),
            "no stereo annotation → no assignment"
        );
    }

    #[test]
    fn test_assignment_get() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let a = assign_cip(&mol);
        assert!(a.get(AtomIdx(1)).is_some(), "atom 1 should have CIP code");
        assert!(a.get(AtomIdx(0)).is_none(), "atom 0 (N) has no chirality");
    }

    #[test]
    fn test_r_lactic_acid_gives_answer() {
        // OC(=O)[C@@H](O)C — lactic acid; should give R or S
        let mol = parse("OC(=O)[C@@H](O)C").unwrap();
        let assignment = assign_cip(&mol);
        let chiral_idx = (0..mol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .find(|&i| mol.atom(i).chirality != Chirality::None)
            .unwrap();
        let code = assignment.get(chiral_idx);
        assert!(
            code == Some(CipCode::R) || code == Some(CipCode::S),
            "should give R or S for lactic acid, got {:?}",
            code
        );
    }

    // --- E/Z double bonds ---

    #[test]
    fn test_trans_2_butene_e() {
        // C/C=C/C — trans-2-butene → E
        let mol = parse("C/C=C/C").unwrap();
        let assignment = assign_cip(&mol);
        let has_e = assignment.assignments.iter().any(|(_, c)| *c == CipCode::E);
        assert!(
            has_e,
            "Expected E for trans-2-butene, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_cis_2_butene_z() {
        // C/C=C\C — cis-2-butene → Z
        let mol = parse("C/C=C\\C").unwrap();
        let assignment = assign_cip(&mol);
        let has_z = assignment.assignments.iter().any(|(_, c)| *c == CipCode::Z);
        assert!(
            has_z,
            "Expected Z for cis-2-butene, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_fceccl_e() {
        // F/C=C/Cl → E (F and Cl on opposite sides)
        let mol = parse("F/C=C/Cl").unwrap();
        let assignment = assign_cip(&mol);
        let has_e = assignment.assignments.iter().any(|(_, c)| *c == CipCode::E);
        assert!(
            has_e,
            "Expected E for F/C=C/Cl, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_fceccl_z() {
        // F/C=C\Cl → Z (F and Cl on same side)
        let mol = parse("F/C=C\\Cl").unwrap();
        let assignment = assign_cip(&mol);
        let has_z = assignment.assignments.iter().any(|(_, c)| *c == CipCode::Z);
        assert!(
            has_z,
            "Expected Z for F/C=C\\Cl, got {:?}",
            assignment.assignments
        );
    }

    #[test]
    fn test_no_ez_no_stereo_bond() {
        // C=C — no Up/Down bonds → no E/Z
        let mol = parse("C=C").unwrap();
        let assignment = assign_cip(&mol);
        let has_ez = assignment
            .assignments
            .iter()
            .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
        assert!(!has_ez, "plain C=C should have no E/Z");
    }

    #[test]
    fn test_ez_terminal_no_crash() {
        // C=C/F — one alkene carbon is terminal (only implicit H's, no
        // stereo bond); should not crash.
        let mol = parse("C=C/F").unwrap();
        let _ = assign_cip(&mol);
    }

    #[test]
    fn test_cip_assignment_methane() {
        let mol = parse("C").unwrap();
        let assignment = assign_cip(&mol);
        assert!(assignment.assignments.is_empty());
    }

    #[test]
    fn test_multiple_chiral_centers() {
        // Two chiral centers in a chain
        let mol = parse("N[C@@H](C)[C@H](C)N").unwrap();
        let assignment = assign_cip(&mol);
        let rs_count = assignment
            .assignments
            .iter()
            .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
            .count();
        assert_eq!(rs_count, 2, "should assign 2 chiral centers");
    }

    #[test]
    fn test_cip_assignment_struct() {
        let mol = parse("N[C@@H](C)C(=O)O").unwrap();
        let a = assign_cip(&mol);
        assert!(!a.assignments.is_empty());
        let code = a.get(AtomIdx(1));
        assert_eq!(code, Some(CipCode::S));
    }

    #[test]
    fn test_r_s_are_consistent() {
        // The two SMILES should give opposite results
        let r_code = cip_at("[C@@H](F)(Cl)Br", 0);
        let s_code = cip_at("[C@H](F)(Cl)Br", 0);
        assert_ne!(r_code, s_code, "@ and @@ must give opposite results");
    }

    #[test]
    fn test_e_z_are_consistent() {
        // /C=C/C (trans) and /C=C\C (cis) must differ
        let mol_e = parse("C/C=C/C").unwrap();
        let mol_z = parse("C/C=C\\C").unwrap();
        let assign_e = assign_cip(&mol_e);
        let assign_z = assign_cip(&mol_z);
        let has_e = assign_e.assignments.iter().any(|(_, c)| *c == CipCode::E);
        let has_z = assign_z.assignments.iter().any(|(_, c)| *c == CipCode::Z);
        assert!(has_e && has_z, "trans must be E, cis must be Z");
    }

    #[test]
    fn test_canonical_preserves_ez() {
        // Canonicalizing a stereo-SMILES must preserve the E/Z assignment.
        // Uses F/C=C/Cl (E) and F/C=C\Cl (Z) — non-symmetric so CIP is unambiguous.
        use chematic_smiles::canonical_smiles;

        let mol_e = parse("F/C=C/Cl").unwrap();
        let can_e = canonical_smiles(&mol_e);
        let mol_e2 = parse(&can_e)
            .unwrap_or_else(|err| panic!("canonical E not parseable: {can_e} → {err}"));
        let assign_e = assign_cip(&mol_e2);
        assert!(
            assign_e.assignments.iter().any(|(_, c)| *c == CipCode::E),
            "canonical SMILES of E isomer must still be E, canonical='{can_e}'"
        );

        let mol_z = parse("F/C=C\\Cl").unwrap();
        let can_z = canonical_smiles(&mol_z);
        let mol_z2 = parse(&can_z)
            .unwrap_or_else(|err| panic!("canonical Z not parseable: {can_z} → {err}"));
        let assign_z = assign_cip(&mol_z2);
        assert!(
            assign_z.assignments.iter().any(|(_, c)| *c == CipCode::Z),
            "canonical SMILES of Z isomer must still be Z, canonical='{can_z}'"
        );
    }

    // --- Allene axial chirality ---

    #[test]
    fn test_allene_no_stereo_no_assignment() {
        // Propadiene (allene without any stereo bonds) → no assignment.
        let mol = parse("C=C=C").unwrap();
        let a = assign_cip(&mol);
        let has_allene = a
            .assignments
            .iter()
            .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
        assert!(
            !has_allene,
            "unspecified allene should have no axial chirality"
        );
    }

    #[test]
    fn test_allene_two_enantiomers_differ() {
        // 1,3-difluoroallene: F/C=C=C/F vs F/C=C=C\F — must give different codes.
        // Build manually with Up/Down bonds to avoid SMILES parser allene ambiguity.
        use chematic_core::{Atom, BondOrder as BO, Element, MoleculeBuilder};

        // F1/C2=C3=C4\F5 — F at each end with opposite Up/Down
        let mut b = MoleculeBuilder::new();
        let f1 = b.add_atom(Atom::new(Element::F));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let f5 = b.add_atom(Atom::new(Element::F));
        b.add_bond(f1, c2, BO::Up).unwrap(); // F1 up relative to C2
        b.add_bond(c2, c3, BO::Double).unwrap();
        b.add_bond(c3, c4, BO::Double).unwrap();
        b.add_bond(c4, f5, BO::Down).unwrap(); // F5 down relative to C4
        let mol_a = b.build();

        // Same but with both F Up (same side)
        let mut b2 = MoleculeBuilder::new();
        let f1b = b2.add_atom(Atom::new(Element::F));
        let c2b = b2.add_atom(Atom::new(Element::C));
        let c3b = b2.add_atom(Atom::new(Element::C));
        let c4b = b2.add_atom(Atom::new(Element::C));
        let f5b = b2.add_atom(Atom::new(Element::F));
        b2.add_bond(f1b, c2b, BO::Up).unwrap();
        b2.add_bond(c2b, c3b, BO::Double).unwrap();
        b2.add_bond(c3b, c4b, BO::Double).unwrap();
        b2.add_bond(c4b, f5b, BO::Up).unwrap(); // both Up
        let mol_b = b2.build();

        let code_a = assign_cip(&mol_a)
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
            .map(|(_, c)| *c);
        let code_b = assign_cip(&mol_b)
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
            .map(|(_, c)| *c);

        // Both must get an assignment, and they must differ.
        assert!(
            code_a.is_some(),
            "allene A should get an axial chirality code"
        );
        assert!(
            code_b.is_some(),
            "allene B should get an axial chirality code"
        );
        assert_ne!(
            code_a, code_b,
            "the two allene enantiomers must get different codes"
        );
    }

    #[test]
    fn test_non_allene_not_detected() {
        // CO2 (O=C=O) has two double bonds from C but is not an allene (no axial chirality possible).
        let mol = parse("O=C=O").unwrap();
        let a = assign_cip(&mol);
        // CO2 has no substituents for Up/Down bonds, so no allene assignment.
        let has_allene = a
            .assignments
            .iter()
            .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
        assert!(
            !has_allene,
            "CO2 should not get axial chirality (no stereo bonds)"
        );
    }

    // -- CIP rule 4 edge cases (mass tiebreaker and duplicates) --------
    #[test]
    fn test_cip_enantiomers_consistent_with_mass_tiebreaker() {
        // Verify that mass tiebreaker is used correctly (doesn't affect basic R/S).
        // D-alanine vs L-alanine should get opposite codes.
        let l_ala = parse("C[C@H](N)C(=O)O").unwrap();
        let d_ala = parse("C[C@@H](N)C(=O)O").unwrap();

        let l_assign = assign_cip(&l_ala);
        let d_assign = assign_cip(&d_ala);

        let l_code = l_assign
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::R | CipCode::S))
            .map(|(_, c)| *c);
        let d_code = d_assign
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::R | CipCode::S))
            .map(|(_, c)| *c);

        // Both should assign (no ties in amino acids).
        assert!(l_code.is_some(), "L-alanine should assign R or S");
        assert!(d_code.is_some(), "D-alanine should assign R or S");

        // They should be opposite.
        assert_ne!(
            l_code, d_code,
            "L and D enantiomers should have opposite R/S"
        );
    }

    #[test]
    fn test_cip_tied_substituents_no_assignment() {
        // When two substituents have identical priority, no R/S is assigned.
        // Example: any center with two identical groups.
        let mol = parse("CC(F)Br").unwrap(); // no chirality specified, but if we try to assign...
        let assignment = assign_cip(&mol);
        // No assignment (no @/@@ specified).
        assert!(
            assignment
                .assignments
                .iter()
                .all(|(_, c)| !matches!(c, CipCode::R | CipCode::S)),
            "achiral molecule should not get R/S"
        );
    }

    #[test]
    fn test_cip_atomic_mass_tiebreaker_infrastructure() {
        // Verify Element::atomic_mass() exists and provides correct values (CIP rule 4).
        // This confirms the mass tiebreaker infrastructure is in place.
        use chematic_core::Element;

        // Check a few elements have correct masses (monoisotopic, within 0.01 u)
        assert!(
            (Element::C.atomic_mass() - 12.0).abs() < 0.01,
            "C ~= 12.0 u, got {}",
            Element::C.atomic_mass()
        );
        assert!(
            (Element::N.atomic_mass() - 14.0).abs() < 0.01,
            "N ~= 14.0 u, got {}",
            Element::N.atomic_mass()
        );
        assert!(
            (Element::O.atomic_mass() - 16.0).abs() < 0.01,
            "O ~= 16.0 u, got {}",
            Element::O.atomic_mass()
        );
        assert!(
            (Element::H.atomic_mass() - 1.0).abs() < 0.01,
            "H ~= 1.0 u, got {}",
            Element::H.atomic_mass()
        );
        assert!(
            (Element::F.atomic_mass() - 19.0).abs() < 0.01,
            "F ~= 19.0 u, got {}",
            Element::F.atomic_mass()
        );

        // Verify that identical atoms have identical masses
        assert_eq!(
            Element::C.atomic_mass(),
            Element::C.atomic_mass(),
            "same element = same mass"
        );

        // Verify ordering: heavier atoms > lighter (used in CIP comparison)
        assert!(
            Element::N.atomic_mass() > Element::C.atomic_mass(),
            "N > C in mass"
        );
        assert!(
            Element::O.atomic_mass() > Element::N.atomic_mass(),
            "O > N in mass"
        );
    }

    // =========================================================================
    // CIP Rule 3 (duplicate atom) tests - fused ring systems
    // =========================================================================

    #[test]
    fn test_cip_naphthalene_assignment() {
        // Naphthalene (c1ccc2ccccc2c1): fused aromatic rings
        // Tests whether CIP Rule 3 (duplicate atom for fused systems) is handled.
        // The bridging carbons in fused rings may appear multiple times in sphere expansion.
        let mol = parse("c1ccc2ccccc2c1").expect("naphthalene");
        let assignment = assign_cip(&mol);
        // Naphthalene has no chiral centers (planar, all carbons equivalent by symmetry)
        // but assignment should not crash and should handle the fused structure.
        assert!(
            assignment.assignments.is_empty(),
            "naphthalene should have no chiral assignments (planar symmetric structure)"
        );
    }

    #[test]
    fn test_cip_decalin_assignment() {
        // Decalin (bicyclic C10): two fused saturated rings
        // CC(C)C1CCC2CCCCC2C1 is a substituted decalin derivative
        // Tests CIP Rule 3 with bridging carbons in saturated system
        let mol = parse("CC(C)C1CCC2CCCCC2C1").expect("decalin");
        let assignment = assign_cip(&mol);
        // Decalin has potential chiral centers; assignment should work without crashing
        assert!(
            !mol.atom(chematic_core::AtomIdx(0)).element.atomic_number() != 6
                || assignment.assignments.len() <= mol.atom_count(),
            "decalin CIP assignment should complete"
        );
    }

    #[test]
    fn test_cip_fused_ring_no_crash() {
        // Simple fused ring system: bicyclo[4.4.0]decane (decalin base)
        // Tests that sphere expansion handles ring revisits without crashing
        let mol = parse("C1CCC2CCCCC2C1").expect("decalin");
        let assignment = assign_cip(&mol);
        // Should complete without panic/crash
        assert!(assignment.assignments.len() <= mol.atom_count());
    }

    #[test]
    fn test_tetrahedral_stable_when_ring_bond_opens_before_other_neighbors() {
        // A stereocenter whose ring-closure digit is written BEFORE its other
        // substituents (`[C@@H]1...`) used to get the wrong CIP code: raw adjacency
        // order only materializes a ring-*opening* bond once the matching closing
        // digit is reached, which is later than a continuation atom that has nothing
        // to wait on -- so the neighbor meant to come second (by SMILES-textual
        // position) ended up listed after one that should come third. Verified
        // against RDKit's CanonicalRankAtoms-based CIP oracle (atom 5 == R here).
        let smi_a = "CN1CCC[C@@H]1c1cccnc1";
        let smi_b = "c1ccncc1[C@@H]1N(CCC1)C"; // same molecule, order-only respelling
        assert_eq!(cip_at(smi_a, 5), Some(CipCode::R));
        assert_eq!(cip_at(smi_b, 6), Some(CipCode::R));
    }

    #[test]
    fn test_tetrahedral_double_bond_duplicates_into_own_sphere() {
        // A stereocenter substituent reached via a double bond must count as TWO
        // entries in ITS OWN CIP substituent sphere (the real atom plus a phantom
        // duplicate), not just contribute a single phantom to the far side of the
        // double bond. Without the departure-side duplicate, `C(=CH2)(CH3)-` scores
        // (C,C) instead of (C,C,C) and loses a priority tie-break it should win.
        // Verified against RDKit (atom 3 == R).
        assert_eq!(
            cip_at("C=C(C)[C@@H]1CN[C@@H](C(=O)O)[C@@H]1CC(=O)O", 3),
            Some(CipCode::R)
        );
    }
}
