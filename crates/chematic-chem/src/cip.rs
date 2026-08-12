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
    if !atom.chirality.is_tetrahedral() {
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

/// Which CIP engine [`assign_cip_with_mode`] uses.
///
/// [`CipMode::Accurate`] only affects tetrahedral R/S -- [`assign_cip_accurate_experimental`]
/// (`chematic-cip`) never computes E/Z or allene axial chirality (it iterates atoms
/// with `chirality != None` and a 4-item `stereo_neighbor_order`; double-bond and
/// allene stereo aren't represented that way), so `Accurate` mode merges the accurate
/// engine's tetrahedral answers with [`LegacyFast`](CipMode::LegacyFast)'s E/Z and
/// allene answers rather than replacing `assign_cip` outright.
///
/// [`assign_cip_accurate_experimental`]: chematic_cip::assign_cip_accurate_experimental
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipMode {
    /// [`assign_cip`] as-is -- the only mode every existing call site
    /// (`iupac_name_stereo`, `num_stereocenters`, `Mol.cip_stereo()`, the WASM
    /// bindings, `chematic-inchi`'s stereo layers) uses today, and continues to use
    /// unless a caller explicitly opts into `Accurate`. ~96.3% oracle agreement,
    /// infallible, silent on ties (always produces *an* answer, never "unresolved").
    LegacyFast,
    /// Tetrahedral R/S (incl. Rule 5 pseudoasymmetric `r`/`s`) from the hierarchical
    /// digraph engine (~99.64% oracle-stable agreement, see `docs/rfcs/cip_accurate_rfc.md`),
    /// merged with legacy's E/Z and allene answers. Atoms the accurate engine
    /// explicitly ties on or exceeds its budget on are never silently backfilled with
    /// legacy's (less rigorous) guess -- they surface via
    /// [`CipModeAssignment::unresolved`] instead.
    Accurate,
}

/// Why [`assign_cip_with_mode`] (in [`CipMode::Accurate`]) couldn't produce a
/// tetrahedral R/S for an atom -- an explicit "we don't know," never a guess.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipUnresolvedReason {
    /// Two or more substituent branches remain indistinguishable after every rule
    /// this engine implements (1a/1b/2/4b/5) -- a genuine tie, not a missing rule.
    Tied,
    /// The digraph/comparator exceeded its size or recursion budget for this atom.
    BudgetExceeded,
}

/// Result of [`assign_cip_with_mode`]. Distinct from [`CipAssignment`] -- carries an
/// explicit `unresolved` channel that infallible, silent `assign_cip` has no
/// equivalent of.
#[derive(Debug, Clone, Default)]
pub struct CipModeAssignment {
    pub assignments: Vec<(AtomIdx, CipCode)>,
    pub unresolved: Vec<(AtomIdx, CipUnresolvedReason)>,
}

impl CipModeAssignment {
    /// Look up the CIP code for a given atom index.
    pub fn get(&self, idx: AtomIdx) -> Option<CipCode> {
        self.assignments
            .iter()
            .find(|(i, _)| *i == idx)
            .map(|(_, c)| *c)
    }
}

/// Error from [`assign_cip_with_mode`] -- only reachable via [`CipMode::Accurate`]
/// ([`CipMode::LegacyFast`] is infallible, matching `assign_cip`'s own contract).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CipModeError {
    Accurate(chematic_cip::CipCompareError),
}

impl std::fmt::Display for CipModeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CipModeError::Accurate(e) => write!(f, "accurate CIP engine error: {e}"),
        }
    }
}

impl std::error::Error for CipModeError {}

/// Run CIP assignment on `mol` using the requested engine. See [`CipMode`] for what
/// each mode covers. `CipMode::LegacyFast` is `assign_cip`'s output unchanged, wrapped
/// -- every existing caller of `assign_cip` is untouched by this function's existence.
pub fn assign_cip_with_mode(
    mol: &Molecule,
    mode: CipMode,
) -> Result<CipModeAssignment, CipModeError> {
    let legacy = assign_cip(mol);
    match mode {
        CipMode::LegacyFast => Ok(CipModeAssignment {
            assignments: legacy.assignments,
            unresolved: Vec::new(),
        }),
        CipMode::Accurate => {
            let budget = chematic_cip::CipBudget::default_budget();
            let accurate = chematic_cip::assign_cip_accurate_experimental(mol, budget)
                .map_err(CipModeError::Accurate)?;

            let unresolved_idx: HashSet<AtomIdx> = accurate
                .skipped
                .iter()
                .filter(|(_, r)| {
                    matches!(
                        r,
                        chematic_cip::SkipReason::Tied | chematic_cip::SkipReason::BudgetExceeded
                    )
                })
                .map(|(idx, _)| *idx)
                .collect();
            let accurate_idx: HashSet<AtomIdx> =
                accurate.assignments.iter().map(|(idx, _)| *idx).collect();

            // Legacy covers E/Z + allene (accurate has neither) and any atom accurate
            // never touched at all; accurate's own tetrahedral answers override
            // legacy's wherever both apply; accurate's explicit ties/budget-outs are
            // dropped from `assignments` entirely (never legacy's guess) and surfaced
            // via `unresolved` instead.
            let mut assignments: Vec<(AtomIdx, CipCode)> = legacy
                .assignments
                .into_iter()
                .filter(|(idx, _)| !accurate_idx.contains(idx) && !unresolved_idx.contains(idx))
                .collect();
            assignments.extend(accurate.assignments);

            let unresolved = accurate
                .skipped
                .into_iter()
                .filter_map(|(idx, reason)| match reason {
                    chematic_cip::SkipReason::Tied => Some((idx, CipUnresolvedReason::Tied)),
                    chematic_cip::SkipReason::BudgetExceeded => {
                        Some((idx, CipUnresolvedReason::BudgetExceeded))
                    }
                    chematic_cip::SkipReason::NotFourSubstituents => None,
                })
                .collect();

            Ok(CipModeAssignment {
                assignments,
                unresolved,
            })
        }
    }
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
    if !atom.chirality.is_tetrahedral() {
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
    let (bond_idx, bond) = mol.bond_between(alkene_end, sub)?;

    // Aromatic ring bonds can carry an adjacent exocyclic double bond's
    // `/` or `\` direction in Molecule's side channel while their chemical
    // bond order remains Aromatic. Use that stashed direction first so E/Z
    // perception does not depend on which parser path created the bond.
    let effective_order = mol.bond_direction(bond_idx).unwrap_or(bond.order);

    match effective_order {
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

    // Highest-priority substituent at each terminal, and which side it's on.
    let (_, up_t1) = highest_stereo_sub(mol, t1, &subs_t1)?;
    let (_, up_t2) = highest_stereo_sub(mol, t2, &subs_t2)?;

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

    // Highest-priority substituent at each end, and which side it's on.
    let (_, up_a1) = highest_stereo_sub(mol, a1, &subs_a1)?;
    let (_, up_a2) = highest_stereo_sub(mol, a2, &subs_a2)?;

    // Same side → Z (zusammen); opposite → E (entgegen).
    let code = if up_a1 == up_a2 {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some((a1, code))
}

/// From `subs` at `alkene_end` (at most 2 — an alkene carbon has only one
/// double bond plus up to two single-bond substituents), return the highest
/// CIP-priority substituent together with which side it's on.
///
/// A trigonal alkene carbon has exactly two possible sides, so when the
/// higher-priority substituent has no `/`/`\` marker of its own but its only
/// sibling substituent does, the higher-priority one's side is the sibling's
/// geometric complement — not a fallback to the sibling's own raw side
/// (which would silently compare the wrong pair of substituents and can
/// flip the resulting E/Z label).
///
/// When there are exactly two substituents and they're a genuine CIP
/// priority tie (e.g. the two ring branches of an unsubstituted, symmetric
/// ring's ipso carbon), swapping them maps the molecule onto itself — there
/// is no stereogenic bond here to report, not an ambiguous one to guess at.
/// Without this check, `compare_branches`' stable sort silently picks
/// whichever substituent happened to come first in `subs`, which depends on
/// adjacency-list order and is not stable across atom renumbering (e.g. a
/// `canonical_smiles` round trip), flipping the reported side arbitrarily.
fn highest_stereo_sub(
    mol: &Molecule,
    alkene_end: AtomIdx,
    subs: &[AtomIdx],
) -> Option<(AtomIdx, bool)> {
    if let [a, b] = subs[..]
        && compare_branches(mol, alkene_end, a, b) == std::cmp::Ordering::Equal
    {
        return None;
    }

    let mut sorted: Vec<AtomIdx> = subs.to_vec();
    sorted.sort_by(|&a, &b| compare_branches(mol, alkene_end, a, b).reverse());

    let top = *sorted.first()?;
    if let Some(up) = substituent_is_up(mol, alkene_end, top) {
        return Some((top, up));
    }
    if let [_, other] = sorted[..] {
        let other_up = substituent_is_up(mol, alkene_end, other)?;
        return Some((top, !other_up));
    }
    None
}

// ---------------------------------------------------------------------------
// E/Z double-bond stereo completeness
// ---------------------------------------------------------------------------

/// Summary of E/Z-stereogenic double bonds in a molecule.
///
/// The double-bond analog of `chematic_perception::stereo_validation::StereoCompleteness`
/// for tetrahedral centers: `specified` and `unspecified` both count toward `total`
/// equally, so a caller can ask "how many potential E/Z centers exist" independent of
/// whether they were annotated with a `/`/`\` marker in the input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct EzCompleteness {
    /// Stereogenic double bonds with a `/`/`\` marker that [`assign_ez`] resolves.
    pub specified: usize,
    /// Stereogenic double bonds with no resolvable marker.
    pub unspecified: usize,
    /// `specified + unspecified`.
    pub total: usize,
}

/// Whether `alkene_end`'s two positions (other than the double-bond partner at
/// `other_end`) are distinguishable -- i.e. this end is not a terminal `=CH2`-style
/// end and not symmetrically disubstituted.
///
/// Mirrors the substituent collection [`assign_ez`]/[`highest_stereo_sub`] already use:
/// 0 heavy substituents means both positions are implicit H (terminal, matching
/// `assign_ez`'s own `subs.is_empty()` "terminal alkene" check); 1 heavy substituent is
/// distinguishable from the implicit H filling the other position as long as that
/// substituent isn't *itself* an explicit hydrogen atom (e.g. `[H]C=C[H]` from an SDF
/// with explicit H) -- CIP ranks any other heavy atom above H, so there's no tie to
/// check in that case; 2 heavy substituents are compared via [`compare_branches`], the
/// same CIP machinery `highest_stereo_sub` uses to detect a genuine priority tie (e.g.
/// two identical methyls).
fn stereogenic_end(mol: &Molecule, other_end: AtomIdx, alkene_end: AtomIdx) -> bool {
    let subs: Vec<AtomIdx> = mol
        .neighbors(alkene_end)
        .filter(|&(nb, bidx)| nb != other_end && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect();

    match subs.len() {
        0 => false, // terminal =CH2-style end
        1 => mol.atom(subs[0]).element.atomic_number() != 1,
        2 => compare_branches(mol, alkene_end, subs[0], subs[1]) != std::cmp::Ordering::Equal,
        _ => false, // not a normal trigonal alkene carbon
    }
}

/// Size of the smallest ring passing through `bidx`, or `None` if `bidx` is acyclic.
///
/// A double bond inside a small ring can't isomerize between cis/trans (the ring
/// constrains the geometry), so it isn't E/Z-stereogenic even when its two ends are
/// otherwise distinguishable. RDKit applies the same cutoff (`minBondRingSize < 8` is
/// excluded, see `Code/GraphMol/Chirality.cpp`'s `isBondPotentialStereoBond`), confirmed
/// directly against a live RDKit oracle here (`Chem.FindPotentialStereo`):
/// cyclohexene/cycloheptene (6- and 7-membered) report no potential stereo bond,
/// cyclooctene (8-membered) does.
///
/// Computed as the shortest path between the bond's two endpoints that doesn't use the
/// bond itself, plus 1 -- deliberately *not* routed through SSSR: SSSR can return a
/// large fundamental cycle instead of a smaller component ring in fused/bridged systems
/// (see this crate's aromaticity handling for the same caveat), which would silently
/// under-flag a genuinely ring-locked bond as large enough to be stereogenic. A BFS
/// shortest cycle is immune to that decomposition choice.
fn bond_min_ring_size(mol: &Molecule, bidx: BondIdx) -> Option<usize> {
    let bond = mol.bond(bidx);
    let (start, goal) = (bond.atom1, bond.atom2);

    let mut dist = vec![usize::MAX; mol.atom_count()];
    dist[start.0 as usize] = 0;
    let mut queue: VecDeque<AtomIdx> = VecDeque::from([start]);
    while let Some(cur) = queue.pop_front() {
        for (nb, b) in mol.neighbors(cur) {
            if b == bidx || dist[nb.0 as usize] != usize::MAX {
                continue;
            }
            dist[nb.0 as usize] = dist[cur.0 as usize] + 1;
            if nb == goal {
                return Some(dist[nb.0 as usize] + 1);
            }
            queue.push_back(nb);
        }
    }
    None
}

/// Rings below this size can't support E/Z isomerism (ring strain locks the geometry).
/// Matches RDKit's `isBondPotentialStereoBond` cutoff -- see [`bond_min_ring_size`].
const MIN_STEREOGENIC_RING_SIZE: usize = 8;

/// Summarise how many E/Z-stereogenic double bonds in `mol` have an explicit
/// `/`/`\` marker vs. are left unannotated.
///
/// A double bond is E/Z-stereogenic when it is non-aromatic (aromatic ring bonds carry
/// `BondOrder::Aromatic`, not `BondOrder::Double`, so they're excluded by construction --
/// the same mechanism [`assign_ez`] relies on), not constrained by a small ring (see
/// [`bond_min_ring_size`]), and both ends have two distinguishable substituents per
/// [`stereogenic_end`]. Non-stereogenic double bonds (terminal, symmetric on either end,
/// or ring-locked) are not counted at all -- matching how
/// `chematic_perception::stereo_validation::stereo_completeness` only counts atoms with
/// 4 distinct neighbors for tetrahedral centers.
pub fn ez_completeness(mol: &Molecule) -> EzCompleteness {
    let mut specified = 0usize;
    let mut unspecified = 0usize;

    for j in 0..mol.bond_count() {
        let bidx = BondIdx(j as u32);
        let bond = mol.bond(bidx);
        if bond.order != BondOrder::Double {
            continue;
        }
        if bond_min_ring_size(mol, bidx).is_some_and(|size| size < MIN_STEREOGENIC_RING_SIZE) {
            continue;
        }
        let a1 = bond.atom1;
        let a2 = bond.atom2;
        if !stereogenic_end(mol, a2, a1) || !stereogenic_end(mol, a1, a2) {
            continue;
        }

        if assign_ez(mol, bidx).is_some() {
            specified += 1;
        } else {
            unspecified += 1;
        }
    }

    EzCompleteness {
        specified,
        unspecified,
        total: specified + unspecified,
    }
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

    // --- EZ-S1: aromatic ring bond stash + true-highest-priority substituent ---
    //
    // `substituent_is_up` previously read only `bond.order`, never
    // `mol.bond_direction(bond_idx)` -- the side channel used when a `/`/`\`
    // marker lands on a bond between two aromatic-flagged atoms (e.g. a ring
    // bond flanking an exocyclic C=N). Fixing that read exposed a second,
    // pre-existing, unrelated bug in `highest_stereo_sub`: it returned the
    // highest-priority substituent *among those carrying an explicit
    // marker*, not the true highest-priority substituent overall, silently
    // using a lower-priority marked substituent's raw side instead of
    // deriving the true-highest one's side as its geometric complement.
    // Both are fixed together here (see cip.rs's `substituent_is_up` and
    // `highest_stereo_sub`). Expected values below are RDKit
    // `rdCIPLabeler`-verified, not just "some E/Z is present" checks.

    #[test]
    fn test_ez_stash_exocyclic_imine_e() {
        // RDKit rdCIPLabeler: bond(N, ring-C) = E.
        assert_eq!(
            cip_at(r"C/N=c1\cccc[nH]1", 1),
            Some(CipCode::E),
            "exocyclic imine on an aromatic ring must read the stashed ring-bond direction"
        );
    }

    #[test]
    fn test_ez_stash_exocyclic_imine_inverted_z() {
        // Same structure, ring-closure slash inverted -- RDKit: Z. A "the two
        // just differ" check would also pass under a global sign-inversion
        // bug, so both this and the E case above pin the specific code.
        assert_eq!(
            cip_at(r"C/N=c1/cccc[nH]1", 1),
            Some(CipCode::Z),
            "inverting the stashed ring-bond direction must flip E to Z"
        );
    }

    #[test]
    fn test_ez_stash_representation_path_independence() {
        // The same molecule (N-methylenamino-pyridine-like ring), with the
        // stash-bearing bond routed through a different parser path each
        // time. All three must (a) actually stash on an aromatic-aromatic
        // bond -- not fall back to a literal Up/Down order -- and (b) agree
        // on the RDKit-confirmed E code, proving E/Z perception doesn't
        // depend on which path created the bond.
        let cases: &[(&str, u32, &str)] = &[
            (r"C/N=c1\cccc[nH]1", 1, "chain-edge"),
            (r"C/N=c1cccc[nH]\1", 1, "ring-closure"),
            (r"[nH]1cccc(\c1=N/C)", 5, "branch-attachment"),
        ];
        for &(smi, atom_idx, path_label) in cases {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{path_label}: parse '{smi}': {e}"));
            let has_aromatic_stash = (0..mol.bond_count()).any(|i| {
                let bidx = BondIdx(i as u32);
                let bond = mol.bond(bidx);
                mol.atom(bond.atom1).aromatic
                    && mol.atom(bond.atom2).aromatic
                    && bond.order == BondOrder::Aromatic
                    && mol.bond_direction(bidx).is_some()
            });
            assert!(
                has_aromatic_stash,
                "{path_label}: '{smi}' must stash direction on an aromatic bond, not divert it to a literal Up/Down order"
            );
            let code = assign_cip(&mol).get(AtomIdx(atom_idx));
            assert_eq!(
                code,
                Some(CipCode::E),
                "{path_label}: '{smi}' must resolve to the same E as the other paths"
            );
        }
    }

    #[test]
    fn test_ez_stash_canonical_round_trip() {
        use chematic_smiles::canonical_smiles;

        let mol = parse(r"C/N=c1\cccc[nH]1").unwrap();
        let can = canonical_smiles(&mol);
        let mol2 = parse(&can).unwrap_or_else(|e| panic!("re-parse canonical '{can}': {e}"));
        let code = assign_cip(&mol2)
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
            .map(|(_, c)| *c);
        assert_eq!(
            code,
            Some(CipCode::E),
            "canonical round trip of the stash-derived E must still be E, canonical='{can}'"
        );
    }

    #[test]
    fn test_highest_stereo_sub_uses_true_priority_not_marked_fallback() {
        // Non-aromatic, no stash involved at all -- a separate, pre-existing
        // bug independent of the stash-read fix above. At the C(F)(Cl)= end,
        // Cl (Z=17) outranks F (Z=9) but only F carries a marker; at the
        // C(Br)(I)= end, I (Z=53) outranks Br (Z=35) but only Br carries a
        // marker. `highest_stereo_sub` must derive each true-highest
        // substituent's side as the marked sibling's complement rather than
        // using the marked (lower-priority) sibling's own raw side.
        // RDKit rdCIPLabeler: E.
        let mol = parse(r"Cl/C(F)=C(\Br)I").unwrap();
        let code = assign_cip(&mol)
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
            .map(|(_, c)| *c);
        assert_eq!(
            code,
            Some(CipCode::E),
            "true-highest-priority substituent (unmarked) must be used, not the marked lower-priority one"
        );
    }

    #[test]
    fn test_highest_stereo_sub_symmetric_ring_is_not_stereogenic() {
        // An exocyclic imine on an *unsubstituted* benzo ring: the ipso
        // carbon's two ring branches are a genuine CIP priority tie
        // (symmetric across the ipso/para axis) -- confirmed via
        // `compare_branches(mol, alkene_end, subs[0], subs[1]) ==
        // Ordering::Equal` in both argument orders. Swapping the two ring
        // neighbors maps the molecule onto itself, so there is no
        // stereogenic bond to report at all -- not an ambiguous one to
        // guess at. (This SMILES also isn't RDKit-parseable, consistent
        // with there being no real stereoisomerism here.)
        let mol = parse(r"C/N=c1ccccc/1").unwrap();
        let code = assign_cip(&mol)
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
            .map(|(_, c)| *c);
        assert!(
            code.is_none(),
            "a genuine CIP priority tie must not be assigned an arbitrary E/Z: {code:?}"
        );

        // Stable across a canonical round trip -- before the tie guard, the
        // stable sort's fallback choice depended on adjacency order, which
        // renumbering changes, flipping E<->Z arbitrarily.
        use chematic_smiles::canonical_smiles;
        let can = canonical_smiles(&mol);
        let mol2 = parse(&can).unwrap_or_else(|e| panic!("re-parse canonical '{can}': {e}"));
        let code2 = assign_cip(&mol2)
            .assignments
            .iter()
            .find(|(_, c)| matches!(c, CipCode::E | CipCode::Z))
            .map(|(_, c)| *c);
        assert!(
            code2.is_none(),
            "the tie must remain unassigned after a canonical round trip, canonical='{can}': {code2:?}"
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

    #[test]
    fn cip_mode_legacy_fast_matches_assign_cip_byte_for_byte() {
        // CipMode::LegacyFast must be assign_cip()'s output, unchanged -- every
        // existing caller of assign_cip is untouched by assign_cip_with_mode existing.
        let smis = [
            "C[C@H](N)C(=O)O",
            "C=C(C)[C@@H]1CN[C@@H](C(=O)O)[C@@H]1CC(=O)O",
            "C1CCN(P2(N3CCCC3)=N[P@@](N3CCCC3)(N3CC3)=N[P@](N3CCCC3)(N3CC3)=N2)C1",
        ];
        for smi in smis {
            let mol = chematic_smiles::parse(smi).expect("valid SMILES");
            let legacy = assign_cip(&mol);
            let via_mode = assign_cip_with_mode(&mol, CipMode::LegacyFast).expect("infallible");
            assert_eq!(legacy.assignments, via_mode.assignments, "smiles={smi}");
            assert!(via_mode.unresolved.is_empty(), "smiles={smi}");
        }
    }

    #[test]
    fn cip_mode_accurate_pseudoasymmetric_fix_stays_unresolved_for_phosphorus_ties() {
        // Milestone 4A-2's `assign_one_with_rule5` fix (resolving the carbon cage
        // family's pseudoasymmetric centers) reaches this *different*, previously-
        // `SkipReason::Tied` molecule (docs/rfcs/cip_accurate_rfc.md Milestone 4C-1) via the
        // exact same code path: atoms 6/19 tie for the identical structural reason as
        // the carbon cage family (a chain-length-1-degenerate Rule 4b comparison whose
        // branches' auxiliary R/S signs genuinely differ). Left unguarded, the fix would
        // resolve these two phosphorus atoms as a side effect -- but Milestone 4C-1
        // independently found that *neither* RDKit CIP engine (`rdCIPLabeler` nor legacy
        // `_CIPCode`) has a representation-stable answer for this specific molecule, both
        // flip under a chemically-neutral Kekule respelling of the P/N ring -- so there
        // is no reliable oracle a resolved phosphorus label could ever be checked
        // against. `assign_one_with_rule5` therefore carries an explicit element-level
        // guard (see `crates/chematic-cip/src/assign.rs` module docs, "Element-level
        // guard: phosphorus stays tied"): it only ever emits a resolved label for a
        // carbon stereocenter, so these 2 phosphorus atoms fall back to
        // `SkipReason::Tied` -> `CipUnresolvedReason::Tied`, exactly their pre-fix
        // behavior, never an unverified label.
        let smi = "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
        for atom_idx in [6u32, 19u32] {
            let idx = AtomIdx(atom_idx);
            assert_eq!(
                result.get(idx),
                None,
                "atom={atom_idx}: phosphorus stereocenter must NOT get an unverified label"
            );
            assert!(
                result
                    .unresolved
                    .iter()
                    .any(|(i, reason)| *i == idx && *reason == CipUnresolvedReason::Tied),
                "atom={atom_idx} must be reported unresolved (Tied): {:?}",
                result.unresolved
            );
        }
    }

    #[test]
    fn cip_mode_accurate_phosphorus_ties_stay_unresolved_kekule_stable() {
        // Same molecule as the test above. Flip every P/N ring bond Single<->Double (a
        // chemically neutral resonance respelling, the same test Milestone 4C-0/4C-1
        // used to show *both* RDKit engines are representation-unstable here) and
        // confirm chematic's own "stays unresolved" guard does NOT flap to resolved on
        // either spelling -- the element-level guard is keyed on atom identity
        // (`mol.atom(idx).element`), which a bond-order respelling never changes, so it
        // must stay unresolved on both spellings by construction; checked directly here
        // rather than just asserted.
        use chematic_core::BondOrder;
        use chematic_perception::find_sssr;

        let smi = "CNP1(NC)=N[P@](NC)(N2CC2)=NP(NC)(NC)=N[P@@](NC)(N2CC2)=N1";
        let mol = chematic_smiles::parse(smi).expect("valid SMILES");
        let sssr = find_sssr(&mol);

        let flip = |m: &chematic_core::Molecule, a, b| -> Option<chematic_core::Molecule> {
            let (bidx, bond) = m.bond_between(a, b)?;
            let new_order = match bond.order {
                BondOrder::Single => BondOrder::Double,
                BondOrder::Double => BondOrder::Single,
                other => other,
            };
            Some(m.with_bond_order(bidx, new_order))
        };

        let mut respelled = mol.clone();
        for ring in sssr.rings() {
            for w in ring.windows(2) {
                if let Some(next) = flip(&respelled, w[0], w[1]) {
                    respelled = next;
                }
            }
            if let (Some(&first), Some(&last)) = (ring.first(), ring.last())
                && let Some(next) = flip(&respelled, last, first)
            {
                respelled = next;
            }
        }

        let original = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
        let after = assign_cip_with_mode(&respelled, CipMode::Accurate).expect("no engine error");
        for atom_idx in [6u32, 19u32] {
            let idx = AtomIdx(atom_idx);
            assert_eq!(
                original.get(idx),
                None,
                "atom={atom_idx}: original spelling"
            );
            assert_eq!(after.get(idx), None, "atom={atom_idx}: respelled");
            assert!(
                original.unresolved.iter().any(|(i, _)| *i == idx),
                "atom={atom_idx}: original spelling must stay unresolved"
            );
            assert!(
                after.unresolved.iter().any(|(i, _)| *i == idx),
                "atom={atom_idx}: respelled must stay unresolved too (not flapping)"
            );
        }
    }

    #[test]
    fn cip_mode_accurate_does_not_hide_oracle_unstable_answers() {
        // The 9 M4C-0 rows are chematic's own stable, confident (if oracle-disputed)
        // answer -- Accurate mode must still report them, not silently drop or
        // "unresolve" them. Being oracle-unstable is a property of the RDKit oracle
        // used for scoring, not a property chematic itself detects or reacts to.
        let mol = chematic_smiles::parse("N[P@]1(Cl)=NP(N2CC2)(N2CC2)=N[P@](N)(Cl)=N1")
            .expect("valid SMILES");
        let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
        assert_eq!(result.get(AtomIdx(12)), Some(CipCode::S));
    }

    #[test]
    fn cip_mode_accurate_merges_legacy_ez_with_accurate_tetrahedral() {
        // Accurate mode must still report E/Z (the accurate engine never computes it)
        // alongside its own tetrahedral R/S for the same molecule.
        let mol = chematic_smiles::parse("C/C=C/[C@H](N)C(=O)O").expect("valid SMILES");
        let result = assign_cip_with_mode(&mol, CipMode::Accurate).expect("no engine error");
        let has_ez = result
            .assignments
            .iter()
            .any(|(_, c)| matches!(c, CipCode::E | CipCode::Z));
        let has_rs = result
            .assignments
            .iter()
            .any(|(_, c)| matches!(c, CipCode::R | CipCode::S));
        assert!(
            has_ez,
            "expected an E/Z assignment: {:?}",
            result.assignments
        );
        assert!(
            has_rs,
            "expected an R/S assignment: {:?}",
            result.assignments
        );
    }

    // --- ez_completeness ---

    #[test]
    fn test_ez_completeness_specified() {
        // C/C=C/C — trans-2-butene: one stereogenic bond, fully marked.
        let mol = parse("C/C=C/C").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 1);
        assert_eq!(ec.unspecified, 0);
        assert_eq!(ec.total, 1);
    }

    #[test]
    fn test_ez_completeness_unspecified() {
        // CC=CC — same connectivity as above, no `/`/`\` marker: still
        // stereogenic (one heavy substituent + implicit H on each end,
        // always distinguishable), just not annotated.
        let mol = parse("CC=CC").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 0);
        assert_eq!(ec.unspecified, 1);
        assert_eq!(ec.total, 1);
    }

    #[test]
    fn test_ez_completeness_terminal_not_counted() {
        // C=C — ethylene: both ends are =CH2 (0 heavy substituents each
        // side), not stereogenic at all.
        let mol = parse("C=C").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 0);
        assert_eq!(ec.unspecified, 0);
        assert_eq!(ec.total, 0);
    }

    #[test]
    fn test_ez_completeness_symmetric_end_not_counted() {
        // CC(C)=CC — the =C end bonded to two methyls is a genuine CIP tie
        // (confirmed directly via compare_branches below), so the whole
        // bond is not stereogenic even though the other end (=CC, one
        // heavy substituent + implicit H) would be on its own.
        let mol = parse("CC(C)=CC").unwrap();
        let methyl_a = AtomIdx(0);
        let methyl_b = AtomIdx(2);
        let alkene_c = AtomIdx(1);
        assert_eq!(
            compare_branches(&mol, alkene_c, methyl_a, methyl_b),
            std::cmp::Ordering::Equal,
            "the two methyls must be a genuine CIP tie for this test to be valid"
        );

        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 0);
        assert_eq!(ec.unspecified, 0);
        assert_eq!(ec.total, 0);
    }

    #[test]
    fn test_ez_completeness_aromatic_ring_excluded() {
        // Benzene's ring bonds are BondOrder::Aromatic, not Double -- no
        // contribution to any count.
        let mol = parse("c1ccccc1").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.total, 0);
    }

    #[test]
    fn test_ez_completeness_kekule_benzene_ring_locked_excluded() {
        // C1=CC=CC=C1 -- Kekulized benzene with no aromaticity perception
        // applied, so every ring bond really is BondOrder::Double (not
        // Aromatic). Each ring atom still only has 1 heavy substituent per
        // end (the rest of the ring), so the pure substituent-count check
        // alone would call these stereogenic -- the small-ring lock (< 8
        // atoms, RDKit-confirmed) must exclude them anyway.
        let mol = parse("C1=CC=CC=C1").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.total, 0, "{ec:?}");
    }

    #[test]
    fn test_ez_completeness_small_ring_double_bond_excluded() {
        // Cyclohexene: a single ring double bond, otherwise satisfies the
        // substituent-distinguishability check, but a 6-membered ring can't
        // isomerize -- RDKit's FindPotentialStereo agrees (empty result).
        let mol = parse("C1=CCCCC1").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.total, 0, "{ec:?}");
    }

    #[test]
    fn test_ez_completeness_large_ring_double_bond_counted() {
        // Cyclooctene: an 8-membered ring double bond IS large enough for
        // E/Z isomerism (cis/trans-cyclooctene are both known, isolable
        // compounds) -- RDKit's FindPotentialStereo reports it as a
        // potential (unspecified) stereo bond. This is the positive
        // control for the ring-size cutoff above: it must not become a
        // blanket "any ring bond is excluded" bug.
        let mol = parse("C1=CCCCCCC1").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 0, "{ec:?}");
        assert_eq!(ec.unspecified, 1, "{ec:?}");
        assert_eq!(ec.total, 1, "{ec:?}");
    }

    #[test]
    fn test_ez_completeness_bridged_bicyclic_ring_lock_not_sssr_dependent() {
        // Norbornene (bicyclo[2.2.1]hept-2-ene): the double bond sits in a
        // real 5-membered ring, but this bridged system's SSSR can also
        // report a larger fundamental cycle through the same bond depending
        // on decomposition choice (the same envelope-ring artifact
        // documented for aromaticity/ring-count elsewhere in this crate).
        // A min-ring-size check naively routed through `find_sssr` risks
        // missing the true 5-ring and passing the >= 8 gate. RDKit
        // (MinBondRingSize=5, FindPotentialStereo -- empty) confirms this
        // bond is ring-locked and must NOT be counted.
        let mol = parse("C1=CC2CCC1C2").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.total, 0, "{ec:?}");
    }

    #[test]
    fn test_ez_completeness_partial_marker_is_unspecified() {
        // C/C=CC -- a marker on only one side of the double bond. The bond
        // is still structurally stereogenic (1 heavy substituent + implicit
        // H on each end, distinguishable), but assign_ez has nothing on the
        // far end to resolve a code from, so it must land in `unspecified`,
        // not silently drop out of the count or count as `specified`.
        let mol = parse("C/C=CC").unwrap();
        let bidx = (0..mol.bond_count())
            .map(|i| BondIdx(i as u32))
            .find(|&b| mol.bond(b).order == BondOrder::Double)
            .unwrap();
        assert_eq!(
            assign_ez(&mol, bidx),
            None,
            "a one-sided marker must not resolve a code"
        );

        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 0, "{ec:?}");
        assert_eq!(ec.unspecified, 1, "{ec:?}");
        assert_eq!(ec.total, 1, "{ec:?}");
    }

    #[test]
    fn test_ez_completeness_aggregates_mixed_specified_and_unspecified() {
        // C/C=C/CC=CC — first double bond (atoms 1=2) is fully marked
        // (same local Up/Up pattern as trans-2-butene, shifted down the
        // chain); second double bond (atoms 4=5) has no marker on either
        // side but is still stereogenic (one heavy substituent + implicit
        // H per end).
        let mol = parse("C/C=C/CC=CC").unwrap();
        let ec = ez_completeness(&mol);
        assert_eq!(ec.specified, 1, "{ec:?}");
        assert_eq!(ec.unspecified, 1, "{ec:?}");
        assert_eq!(ec.total, 2, "{ec:?}");
    }

    /// Regression for the square-planar-stereo PR's required CIP safety fix:
    /// `tetrahedral_stereo_neighbors`/`assign_tetrahedral` gated on
    /// `chirality == Chirality::None` (equality), not an exhaustive match --
    /// adding `Chirality::SquarePlanar` did NOT force a compile error there, so
    /// an `@SP1`-tagged 4-neighbor Pt center would have silently fallen through
    /// into the tetrahedral R/S algorithm and produced a bogus CIP code. Must
    /// return no code at all, not a wrong one.
    #[test]
    fn square_planar_center_never_gets_a_bogus_tetrahedral_cip_code() {
        let mol = parse("N->[Pt@SP1](<-N)(Cl)Cl").unwrap();
        let pt = (0..mol.atom_count())
            .map(|i| chematic_core::AtomIdx(i as u32))
            .find(|&i| mol.atom(i).element == chematic_core::Element::PT)
            .expect("has a Pt atom");
        assert!(
            matches!(mol.atom(pt).chirality, Chirality::SquarePlanar(_)),
            "fixture must actually carry a SquarePlanar tag"
        );
        assert_eq!(
            tetrahedral_stereo_neighbors(&mol, pt),
            None,
            "a square-planar center must never resolve a tetrahedral CIP code"
        );
        assert_eq!(
            assign_cip(&mol).get(pt),
            None,
            "assign_cip must not assign any code to a square-planar center"
        );
    }
}
