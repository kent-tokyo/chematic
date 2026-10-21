//! CIP (Cahn–Ingold–Prelog) stereochemistry assignment.
//!
//! Implements R/S (tetrahedral) and E/Z (double bond) assignment for
//! molecules parsed from SMILES with chirality annotations.

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, CipCode, Element, Molecule, implicit_hcount,
};

#[path = "cip_aromatic_ez.rs"]
mod aromatic_ez;

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

// Only two explicitly present substituents require digraph ranking. Most
// molecules have no marked E/Z bond, and a single heavy substituent can be
// read directly from its direction without preparing Kekule/MANCUDE data.
struct LazyEzRanker<'a> {
    mol: &'a Molecule,
    budget: chematic_cip::CipBudget,
    value: std::cell::OnceCell<chematic_cip::SubstituentRanker>,
}

impl<'a> LazyEzRanker<'a> {
    fn new(mol: &'a Molecule, budget: chematic_cip::CipBudget) -> Self {
        Self {
            mol,
            budget,
            value: std::cell::OnceCell::new(),
        }
    }

    fn get(&self) -> &chematic_cip::SubstituentRanker {
        self.value
            .get_or_init(|| chematic_cip::SubstituentRanker::new(self.mol, self.budget))
    }
}

/// E/Z labels keyed by the double bond they describe.
///
/// [`assign_cip`] (and the E/Z part of [`assign_cip_with_mode`]) report each
/// E/Z label under the double bond's first atom (`bond.atom1`), in bond-index
/// order. This returns the same labels, in the same order, keyed by the bond
/// itself, so callers need not reconstruct which bond an atom-keyed E/Z entry
/// refers to.
pub fn assign_ez_bonds(mol: &Molecule) -> Vec<(BondIdx, CipCode)> {
    assign_ez_bonds_with_mode(mol, CipMode::LegacyFast)
}

/// [`assign_ez_bonds`] for the E/Z labels of [`assign_cip_with_mode`] in
/// `mode` ([`CipMode::Accurate`] ranks substituents with the accurate engine).
pub fn assign_ez_bonds_with_mode(mol: &Molecule, mode: CipMode) -> Vec<(BondIdx, CipCode)> {
    match mode {
        CipMode::LegacyFast => (0..mol.bond_count())
            .filter_map(|j| {
                let bidx = BondIdx(j as u32);
                assign_ez(mol, bidx).map(|(_, code)| (bidx, code))
            })
            .collect(),
        CipMode::Accurate => {
            let ranker = LazyEzRanker::new(mol, chematic_cip::CipBudget::default_budget());
            (0..mol.bond_count())
                .filter_map(|j| {
                    let bidx = BondIdx(j as u32);
                    assign_ez_accurate(mol, bidx, &ranker).map(|(_, code)| (bidx, code))
                })
                .collect()
        }
    }
}

/// Which CIP engine [`assign_cip_with_mode`] uses.
///
/// [`assign_cip_accurate_experimental`] (`chematic-cip`) labels tetrahedral centres
/// only (it iterates atoms with `chirality != None` and a 4-item
/// `stereo_neighbor_order`); `Accurate` mode labels double-bond E/Z by ranking each
/// end's substituents with the same engine (`SubstituentRanker`, #634) and keeps
/// [`LegacyFast`](CipMode::LegacyFast)'s allene answers.
///
/// [`assign_cip_accurate_experimental`]: chematic_cip::assign_cip_accurate_experimental
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipMode {
    /// [`assign_cip`] as-is -- the only mode every existing call site
    /// (`iupac_name_stereo`, `num_stereocenters`, `Mol.cip_stereo()`, the WASM
    /// bindings, `chematic-inchi`'s stereo layers) uses today, and continues to use
    /// unless a caller explicitly opts into `Accurate`. Aromatic-substituted
    /// E/Z bonds use the accurate substituent ranker when it resolves both
    /// ends; this fixes MANCUDE priorities without changing the legacy
    /// tetrahedral engine. Infallible and silent on unresolved tetrahedral ties.
    LegacyFast,
    /// Tetrahedral R/S (incl. Rule 5 pseudoasymmetric `r`/`s`) from the hierarchical
    /// digraph engine (~99.64% oracle-stable agreement, see `docs/rfcs/cip_accurate_rfc.md`),
    /// E/Z ranked by the same engine, and legacy's allene answers. Atoms the accurate engine
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
    /// No chemically stable label is available. A wildcard/dummy atom can
    /// leave a direct ligand's priority unknown or reverse a more distant
    /// center's R/S label when its element is specified.
    OracleUnstable,
    /// A stereo-tagged centre with three explicit ligands whose fourth ligand
    /// would be a lone pair (a bridgehead amine, a sulfoxide); lone pairs are
    /// not modelled as CIP ligands, so no label is emitted.
    LonePairCenter,
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

/// Whether the accurate (`CipMode::Accurate`) label of tetrahedral centre
/// `atom` can change between neutral Kekulé spellings of the same molecule
/// (a phosphorus on an unsaturated ring, such as a cyclophosphazene P). RDKit's
/// CIPLabeler flips those labels the same way; identity checks should not
/// compare them across inputs.
pub fn cip_label_depends_on_kekule_spelling(mol: &Molecule, atom: AtomIdx) -> bool {
    chematic_cip::label_depends_on_kekule_spelling(mol, atom)
}

/// Whether [`CipMode::Accurate`]'s label of tetrahedral centre `atom` is
/// RDKit's reading of the SMILES spelling it was parsed from rather than a
/// property of the molecule: a ring centre with three single-bonded ligands
/// and a lone pair (bridgehead amine, cyclic phosphine or sulfonium), whose
/// `@`/`@@` RDKit inverts when one ring-closure digit is written on it
/// (`chematic_cip::label_follows_rdkit_smiles_reading`).
pub fn cip_label_depends_on_smiles_spelling(mol: &Molecule, atom: AtomIdx) -> bool {
    chematic_cip::label_follows_rdkit_smiles_reading(mol, atom)
}

/// Run CIP assignment on `mol` using the requested engine. See [`CipMode`] for what
/// each mode covers. `CipMode::LegacyFast` is `assign_cip`'s output unchanged, wrapped
/// -- every existing caller of `assign_cip` is untouched by this function's existence.
/// Legacy covers E/Z + allene (accurate has neither) and any atom accurate
/// never touched at all; accurate's own tetrahedral answers override
/// legacy's wherever both apply; accurate's explicit ties/budget-outs are
/// dropped from `assignments` entirely (never legacy's guess) and surfaced
/// via `unresolved` instead.
///
/// Double-bond E/Z uses the accurate engine's substituent ranking too
/// (#634): the legacy sphere comparison adds no duplicate atoms for
/// aromatic bonds, so an aryl substituent could lose to tert-butyl.
/// Legacy E/Z entries are replaced one-for-one (same atom1 keying,
/// same bond order); an end the accurate ranking cannot order drops
/// its label rather than keeping legacy's guess.
pub fn assign_cip_with_mode(
    mol: &Molecule,
    mode: CipMode,
) -> Result<CipModeAssignment, CipModeError> {
    match mode {
        CipMode::LegacyFast => Ok(CipModeAssignment {
            assignments: assign_cip(mol).assignments,
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
                        chematic_cip::SkipReason::Tied
                            | chematic_cip::SkipReason::BudgetExceeded
                            | chematic_cip::SkipReason::OracleUnstable
                            | chematic_cip::SkipReason::LonePairCenter
                    )
                })
                .map(|(idx, _)| *idx)
                .collect();
            let accurate_idx: HashSet<AtomIdx> =
                accurate.assignments.iter().map(|(idx, _)| *idx).collect();

            let legacy_ez: Vec<(AtomIdx, CipCode)> = (0..mol.bond_count())
                .filter_map(|j| assign_ez(mol, BondIdx(j as u32)))
                .collect();
            // An accurate answer or explicit abstention supersedes the legacy
            // tetrahedral result. Do not rank those same ligands twice. Keep
            // every legacy fallback and its original tetrahedral/EZ/allene order.
            let mut legacy_assignments: Vec<_> = (0..mol.atom_count())
                .map(|i| AtomIdx(i as u32))
                .filter(|idx| {
                    mol.atom(*idx).chirality.is_tetrahedral()
                        && !accurate_idx.contains(idx)
                        && !unresolved_idx.contains(idx)
                })
                .filter_map(|idx| assign_tetrahedral(mol, idx).map(|code| (idx, code)))
                .collect();
            legacy_assignments.extend(legacy_ez.iter().copied());
            legacy_assignments.extend((0..mol.atom_count()).filter_map(|i| {
                let idx = AtomIdx(i as u32);
                is_allene_central(mol, idx)
                    .then(|| assign_allene(mol, idx))
                    .flatten()
            }));
            let ranker = LazyEzRanker::new(mol, budget);
            let accurate_ez: Vec<(AtomIdx, CipCode)> = (0..mol.bond_count())
                .filter_map(|j| assign_ez_accurate(mol, BondIdx(j as u32), &ranker))
                .collect();
            let mut legacy_ez_left = legacy_ez.len();
            let mut assignments: Vec<(AtomIdx, CipCode)> = Vec::new();
            for (idx, code) in legacy_assignments {
                // assign_cip emits tetrahedral, then double-bond E/Z (in bond
                // order), then allene entries; the E/Z block is exactly
                // `legacy_ez`.
                let is_ez_block = legacy_ez_left > 0
                    && matches!(code, CipCode::E | CipCode::Z)
                    && legacy_ez[legacy_ez.len() - legacy_ez_left] == (idx, code);
                if is_ez_block {
                    legacy_ez_left -= 1;
                    if legacy_ez_left == 0 {
                        assignments.extend(accurate_ez.iter().copied().filter(|(i, _)| {
                            !accurate_idx.contains(i) && !unresolved_idx.contains(i)
                        }));
                    }
                    continue;
                }
                if !accurate_idx.contains(&idx) && !unresolved_idx.contains(&idx) {
                    assignments.push((idx, code));
                }
            }
            if legacy_ez.is_empty() {
                // No legacy E/Z block to replace; append any accurate labels.
                assignments.extend(
                    accurate_ez
                        .iter()
                        .copied()
                        .filter(|(i, _)| !accurate_idx.contains(i) && !unresolved_idx.contains(i)),
                );
            }
            assignments.extend(accurate.assignments);

            let unresolved = accurate
                .skipped
                .into_iter()
                .filter_map(|(idx, reason)| match reason {
                    chematic_cip::SkipReason::Tied => Some((idx, CipUnresolvedReason::Tied)),
                    chematic_cip::SkipReason::BudgetExceeded => {
                        Some((idx, CipUnresolvedReason::BudgetExceeded))
                    }
                    chematic_cip::SkipReason::OracleUnstable => {
                        Some((idx, CipUnresolvedReason::OracleUnstable))
                    }
                    chematic_cip::SkipReason::LonePairCenter => {
                        Some((idx, CipUnresolvedReason::LonePairCenter))
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
///
/// Eager reference implementation; production code uses the lazy
/// [`SphereGen`]. Kept as the differential-test oracle.
#[cfg_attr(not(test), allow(dead_code))]
fn cip_branch_spheres(mol: &Molecule, center: AtomIdx, start: AtomIdx) -> Vec<SphereLayer> {
    let mut layers: HashMap<usize, Vec<(u8, Option<u16>, f64)>> = HashMap::new();
    // Distinguish remote substituents on long ring/chain paths.  Eight
    // layers silently classified some otherwise distinct branches as tied;
    // keep the bound finite for pathological graphs while covering the
    // longest structures in the descriptor corpus.
    let max_depth = 16usize;

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

/// Lazy, layer-at-a-time equivalent of [`cip_branch_spheres`].
///
/// The sphere expansion is a FIFO walk whose states are dequeued in
/// nondecreasing depth, so layer `d + 1` is complete once every depth-`d`
/// state has been processed. Producing layers on demand lets a comparison
/// stop at the first differing sphere instead of enumerating every path to
/// depth 16 for both branches. Each produced layer (after sorting) is
/// identical to the corresponding entry of `cip_branch_spheres`, and the
/// generator is exhausted exactly where that vector ends. The visited set is
/// the path itself (at most 17 atoms), kept as a small vector.
struct SphereGen<'a> {
    mol: &'a Molecule,
    center: AtomIdx,
    queue: VecDeque<(AtomIdx, AtomIdx, usize, Vec<AtomIdx>)>,
    first: Option<(u8, Option<u16>, f64)>,
    depth: usize,
}

impl<'a> SphereGen<'a> {
    const MAX_DEPTH: usize = 16;

    fn new(mol: &'a Molecule, center: AtomIdx, start: AtomIdx) -> Self {
        let mut queue = VecDeque::new();
        // `u32::MAX` is the virtual hydrogen/lone-pair sentinel. It has an
        // atom key, but it is terminal and must never enter the molecular
        // adjacency walk. A real explicit H is terminal for the same reason,
        // so both representations compare identically after their first
        // sphere.
        if start.0 != u32::MAX {
            queue.push_back((start, center, 1usize, vec![center, start]));
        }
        Self {
            mol,
            center,
            queue,
            first: Some(atom_key(mol, start)),
            depth: 1,
        }
    }

    /// Next sphere layer (sorted descending), or `None` once exhausted.
    fn next_layer(&mut self) -> Option<SphereLayer> {
        if let Some(first) = self.first.take() {
            return Some(vec![first]);
        }
        let d = self.depth;
        let mut layer: SphereLayer = Vec::new();
        let mol = self.mol;
        while self.queue.front().is_some_and(|s| s.2 == d) {
            let (node, parent, depth, visited) = self.queue.pop_front().unwrap();
            if depth >= Self::MAX_DEPTH {
                continue;
            }
            let child_depth = depth + 1;
            if let Some((_, bond_to_parent)) = mol.bond_between(node, parent)
                && bond_to_parent.order == BondOrder::Double
            {
                layer.push(atom_key(mol, parent));
            }
            for (nb, _) in mol.neighbors(node) {
                if nb == parent || nb == self.center {
                    continue;
                }
                let child_key = atom_key(mol, nb);
                let is_double = mol
                    .bond_between(node, nb)
                    .is_some_and(|(_, b)| b.order == BondOrder::Double);
                if is_double {
                    layer.push(child_key);
                }
                layer.push(child_key);
                if !visited.contains(&nb) {
                    let mut child_visited = visited.clone();
                    child_visited.push(nb);
                    self.queue.push_back((nb, node, child_depth, child_visited));
                }
            }
        }
        self.depth += 1;
        if layer.is_empty() {
            // No entries at this depth means no states at this depth either,
            // so every deeper layer is empty too.
            self.queue.clear();
            return None;
        }
        layer.sort_by(|a, b| cmp_key(*b, *a)); // descending
        Some(layer)
    }
}

/// Lazy equality of two branches' sphere sequences (same result as
/// `cip_branch_spheres(a) == cip_branch_spheres(b)`).
fn sphere_sequences_equal(
    mol: &Molecule,
    center_a: AtomIdx,
    start_a: AtomIdx,
    center_b: AtomIdx,
    start_b: AtomIdx,
) -> bool {
    let mut ga = SphereGen::new(mol, center_a, start_a);
    let mut gb = SphereGen::new(mol, center_b, start_b);
    loop {
        match (ga.next_layer(), gb.next_layer()) {
            (None, None) => return true,
            (Some(la), Some(lb)) if la == lb => {}
            _ => return false,
        }
    }
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

    // Sphere-by-sphere comparison, generated lazily (see `SphereGen`).
    let mut a_gen = SphereGen::new(mol, center, a);
    let mut b_gen = SphereGen::new(mol, center, b);
    loop {
        let a_spheres = a_gen.next_layer();
        let b_spheres = b_gen.next_layer();
        if a_spheres.is_none() && b_spheres.is_none() {
            break;
        }
        let a_layer = a_spheres.as_deref().unwrap_or(&[]);
        let b_layer = b_spheres.as_deref().unwrap_or(&[]);

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

/// True when the bond from `center` to `neighbor` has an alternate path,
/// i.e. participates in at least one cycle. The bond itself is excluded from
/// the BFS so the direct edge cannot make the test trivially succeed.
fn bond_is_in_ring(mol: &Molecule, center: AtomIdx, neighbor: AtomIdx) -> bool {
    let Some((bond, _)) = mol.bond_between(center, neighbor) else {
        return false;
    };
    let mut visited = HashSet::from([center]);
    let mut queue = VecDeque::from([center]);
    while let Some(current) = queue.pop_front() {
        for (next, edge) in mol.neighbors(current) {
            if edge == bond {
                continue;
            }
            if next == neighbor {
                return true;
            }
            if visited.insert(next) {
                queue.push_back(next);
            }
        }
    }
    false
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
    let max_depth = 16usize;

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
    // Aromatic atoms are normally planar and cannot be tetrahedral centers.
    // RDKit's potential-center model is an exception for aromatic sulfur with
    // an exocyclic substituent (for example `[s+]([O-])`): its aromatic ring
    // neighbors plus the exocyclic branch form the four-site center counted by
    // CalcNumAtomStereoCenters. Keep the exception narrow; aromatic C/N/P and
    // other heteroatoms remain excluded.
    if atom.aromatic && atom.element.atomic_number() != 16 {
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
    // RDKit's ring-constrained neutral sp3-N rule is narrower than merely
    // saying "three distinct branches": each of the three neighbors must be
    // a ring degree-2 atom, so a fused bridgehead topology is required. This
    // excludes the common fused-ring tertiary amines whose one neighbor is a
    // degree-3 junction and whose inversion remains freely accessible.
    if atom.element.atomic_number() == 7
        && !atom.aromatic
        && atom.charge == 0
        && h == 0
        && neighbors.len() == 3
        && neighbors.iter().all(|&neighbor| {
            mol.atom(neighbor).element.atomic_number() != 1
                && mol.neighbors(neighbor).count() == 2
                && bond_is_in_ring(mol, idx, neighbor)
        })
        && !neighbors.iter().any(|&neighbor| {
            mol.neighbors(neighbor).any(|(other, bond_idx)| {
                mol.bond(bond_idx).order == BondOrder::Double
                    && matches!(mol.atom(other).element.atomic_number(), 8 | 16)
            })
        })
        && accurate_graph_branches_are_distinct(mol, idx, &neighbors)
    {
        return true;
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

    // The legacy pooled-sphere comparator can tie distinct carbon branches in
    // fused/bridged systems because it discards branch provenance. Re-check
    // only that unresolved carbon case with the provenance-preserving Rules
    // 1a/1b/2 comparator. This is a classification fallback, not a change to
    // the public legacy CIP labeler and not a Rule 5 guess.
    if atom.element == Element::C && accurate_graph_branches_are_distinct(mol, idx, &neighbors) {
        return true;
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
    //
    // Signatures are compared lazily — atom key, then graph spheres (layer by
    // layer, see `SphereGen`), then stereo spheres (computed only when
    // needed) — which is equivalent to comparing the eager
    // `(atom_key, cip_branch_spheres, cip_branch_stereo_spheres)` tuples.
    let keys: Vec<_> = neighbors.iter().map(|&nb| atom_key(mol, nb)).collect();
    let stereo: Vec<std::cell::OnceCell<_>> = (0..neighbors.len())
        .map(|_| std::cell::OnceCell::new())
        .collect();
    let stereo_of = |k: usize| {
        stereo[k].get_or_init(|| {
            let nb = neighbors[k];
            if nb.0 == u32::MAX {
                vec![]
            } else {
                cip_branch_stereo_spheres(mol, idx, nb, provisional)
            }
        })
    };
    let spheres_equal = |a: AtomIdx, b: AtomIdx| match (a.0 == u32::MAX, b.0 == u32::MAX) {
        (true, true) => true,
        // A sentinel has no spheres; a heavy branch always has layer 1.
        (true, false) | (false, true) => false,
        (false, false) => sphere_sequences_equal(mol, idx, a, idx, b),
    };

    // All 6 pairwise pairs must be unequal for 4 distinct substituents.
    for i in 0..4 {
        for j in (i + 1)..4 {
            if keys[i] == keys[j]
                && spheres_equal(neighbors[i], neighbors[j])
                && stereo_of(i) == stereo_of(j)
            {
                return false;
            }
        }
    }
    true
}

fn accurate_graph_branches_are_distinct(
    mol: &Molecule,
    center: AtomIdx,
    neighbors: &[AtomIdx],
) -> bool {
    let mut graph =
        match chematic_cip::CipDigraph::new(mol, center, chematic_cip::CipBudget::default_budget())
        {
            Ok(graph) => graph,
            Err(_) => return false,
        };
    let root_children = match graph.expand_children(graph.root()) {
        Ok(children) => children,
        Err(_) => return false,
    };
    if root_children.len() != neighbors.len() {
        return false;
    }
    let mut context = chematic_cip::CompareContext::new();
    match chematic_cip::rank_children(&mut graph, &root_children, &mut context) {
        Ok(groups) => groups.iter().all(|group| group.len() == 1),
        Err(_) => false,
    }
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
/// Whether `bond` lies on a ring of fewer than eight atoms.
fn double_bond_in_small_ring(mol: &Molecule, bond: BondIdx) -> bool {
    let (start, goal) = (mol.bond(bond).atom1, mol.bond(bond).atom2);
    let mut seen = vec![false; mol.atom_count()];
    seen[start.0 as usize] = true;
    let mut frontier = vec![start];
    for _ in 0..6 {
        let mut next = Vec::new();
        for &a in &frontier {
            for (nb, bi) in mol.neighbors(a) {
                if bi == bond {
                    continue;
                }
                if nb == goal {
                    return true;
                }
                if !seen[nb.0 as usize] {
                    seen[nb.0 as usize] = true;
                    next.push(nb);
                }
            }
        }
        frontier = next;
    }
    false
}

/// Whether `end` (double-bonded to `other`) has a substituent carrying a
/// `/`/`\` marker: without one at both ends a double bond has no E/Z to
/// read, so the ring and CIP checks can be skipped.
fn end_has_marked_substituent(mol: &Molecule, end: AtomIdx, other: AtomIdx) -> bool {
    mol.neighbors(end).any(|(nb, bidx)| {
        nb != other
            && mol.bond(bidx).order != BondOrder::Double
            && substituent_is_up(mol, end, nb).is_some()
    })
}

fn assign_ez(mol: &Molecule, bond_idx: BondIdx) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    // A ring of fewer than eight atoms holds its double bond cis: no E/Z
    // (RDKit drops such stereo too).
    if bond.order != BondOrder::Double
        || !end_has_marked_substituent(mol, bond.atom1, bond.atom2)
        || !end_has_marked_substituent(mol, bond.atom2, bond.atom1)
        || double_bond_in_small_ring(mol, bond_idx)
    {
        return None;
    }

    if let Some(result) = aromatic_ez::try_assign(mol, bond_idx) {
        return Some(result);
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
    // Without a `/`/`\` marker on either substituent there is no E/Z to
    // read, whatever the ranking: skip the CIP comparison (a carbonyl, an
    // unmarked alkene).
    if subs
        .iter()
        .all(|&sub| substituent_is_up(mol, alkene_end, sub).is_none())
    {
        return None;
    }
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

/// E/Z for the double bond at `bond_idx`, ranking each end's substituents with
/// the accurate engine's comparator ([`chematic_cip::SubstituentRanker`]:
/// hierarchical digraph, Rules 1a/2, MANCUDE duplicates for aromatic rings)
/// instead of the legacy sphere comparison. Same geometry reading and keying
/// as [`assign_ez`]. `None` when the bond is not a specified stereo double
/// bond or an end's two substituents tie (or exceed the comparison budget) —
/// never a guess.
fn assign_ez_accurate(
    mol: &Molecule,
    bond_idx: BondIdx,
    ranker: &LazyEzRanker<'_>,
) -> Option<(AtomIdx, CipCode)> {
    let bond = mol.bond(bond_idx);
    // A ring of fewer than eight atoms holds its double bond cis: no E/Z
    // (RDKit drops such stereo too).
    if bond.order != BondOrder::Double
        || !end_has_marked_substituent(mol, bond.atom1, bond.atom2)
        || !end_has_marked_substituent(mol, bond.atom2, bond.atom1)
        || double_bond_in_small_ring(mol, bond_idx)
    {
        return None;
    }
    let a1 = bond.atom1;
    let a2 = bond.atom2;
    let subs = |end: AtomIdx, other: AtomIdx| -> Vec<AtomIdx> {
        mol.neighbors(end)
            .filter(|&(nb, bidx)| nb != other && mol.bond(bidx).order != BondOrder::Double)
            .map(|(nb, _)| nb)
            .collect()
    };
    let subs_a1 = subs(a1, a2);
    let subs_a2 = subs(a2, a1);
    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None;
    }
    let (_, up_a1) = highest_stereo_sub_accurate(mol, a1, &subs_a1, ranker)?;
    let (_, up_a2) = highest_stereo_sub_accurate(mol, a2, &subs_a2, ranker)?;
    let code = if up_a1 == up_a2 {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some((a1, code))
}

/// [`highest_stereo_sub`] with the accurate comparator (at most two
/// substituents on a double-bond end).
fn highest_stereo_sub_accurate(
    mol: &Molecule,
    alkene_end: AtomIdx,
    subs: &[AtomIdx],
    ranker: &LazyEzRanker<'_>,
) -> Option<(AtomIdx, bool)> {
    // No marker on either substituent: nothing to read (see
    // `highest_stereo_sub`).
    if subs
        .iter()
        .all(|&sub| substituent_is_up(mol, alkene_end, sub).is_none())
    {
        return None;
    }
    let (top, other) = match subs {
        [only] => (*only, None),
        [a, b] => match ranker.get().compare(mol, alkene_end, *a, *b).ok()?? {
            std::cmp::Ordering::Greater => (*a, Some(*b)),
            std::cmp::Ordering::Less => (*b, Some(*a)),
            std::cmp::Ordering::Equal => return None,
        },
        _ => return None,
    };
    if let Some(up) = substituent_is_up(mol, alkene_end, top) {
        return Some((top, up));
    }
    let other_up = substituent_is_up(mol, alkene_end, other?)?;
    Some((top, !other_up))
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

/// RDKit's canonical atom ranks (`Chem.CanonicalRankAtoms(mol,
/// breakTies=True)`) including double-bond stereo: the E/Z labels of the
/// legacy engine stand for RDKit's `STEREOE`/`STEREOZ` bond stereo (see
/// [`chematic_perception::rdkit_canonical_atom_ranks_with_bond_stereo`]).
pub fn rdkit_canonical_atom_ranks(mol: &Molecule) -> Vec<u32> {
    let ez: Vec<(BondIdx, bool)> = assign_ez_bonds_with_mode(mol, CipMode::LegacyFast)
        .into_iter()
        .filter_map(|(b, code)| match code {
            CipCode::E => Some((b, true)),
            CipCode::Z => Some((b, false)),
            _ => None,
        })
        .collect();
    chematic_perception::rdkit_canonical_atom_ranks_with_bond_stereo(mol, &ez)
}

#[cfg(test)]
#[path = "cip_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "cip_square_planar_tests.rs"]
mod square_planar_tests;

#[cfg(test)]
#[path = "cip_preparation_tests.rs"]
mod preparation_tests;
