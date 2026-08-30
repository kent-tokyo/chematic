//! Maximum Common Substructure (MCS) search using the McGregor connected-growth algorithm.
//!
//! # Overview
//! Given a set of molecules, `find_mcs` returns a `QueryMolecule` that represents the
//! largest substructure common to all input molecules.  The algorithm is a branch-and-bound
//! search that grows the common subgraph one atom at a time, pruning subtrees whose upper
//! bound on reachable size is no better than the current best.

use std::collections::HashMap;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::{RingSet, find_sssr};

use crate::clock::Instant;
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

// ---------------------------------------------------------------------------
// Public API types
// ---------------------------------------------------------------------------

/// Atom-level comparison mode for MCS search.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum AtomCompare {
    /// Atoms match if they share the same atomic number, regardless of aromaticity
    /// (matches RDKit's `rdFMCS.AtomCompare.CompareElements` exactly -- confirmed
    /// against a live oracle, e.g. `[#6]` matches both aromatic and non-aromatic
    /// carbon; RDKit never encodes aromaticity as a per-atom constraint, only via
    /// bond-type queries). This is the default.
    #[default]
    Elements,
    /// Any heavy atom (atomic number > 1) matches any other heavy atom, regardless of
    /// element.  Useful for scaffold hopping across heterocycle series (C↔N, C↔O, etc.).
    AnyHeavyAtom,
    /// Any atom matches any other atom, including hydrogen.
    Any,
}

/// Bond-level comparison mode for MCS search.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub enum BondCompare {
    /// Bonds match only if they have the same normalized order (Up/Down → Single).
    /// This is the default.
    #[default]
    OrderOrAromatic,
    /// Any bond order matches any other bond order.
    Any,
}

/// Configuration options for MCS search.
#[derive(Debug, Clone)]
pub struct McsConfig {
    /// If `true`, bond orders must match between molecules. If `false`, any bond order is
    /// accepted (atom-only MCS).
    pub match_bonds: bool,
    /// Minimum number of atoms required for the result to be returned.  If the best
    /// common substructure has fewer atoms than `min_atoms`, an empty `QueryMolecule` is
    /// returned instead.
    pub min_atoms: usize,
    /// Optional time limit in milliseconds.  The search is aborted when the deadline is
    /// reached and the best result found so far is returned.
    pub timeout_ms: Option<u64>,
    /// If `true`, ring atoms can only be matched to ring atoms, and non-ring atoms only
    /// to non-ring atoms.  Prevents chemically meaningless matches across ring boundaries.
    pub ring_matches_ring_only: bool,
    /// If `true`, the MCS result must not contain a partial ring.  Any ring that is only
    /// partially present in the MCS is removed entirely (iterative post-processing).
    pub complete_rings_only: bool,
    /// Atom-level comparison strategy.  Defaults to [`AtomCompare::Elements`].
    pub atom_compare: AtomCompare,
    /// Bond-level comparison strategy.  Defaults to [`BondCompare::OrderOrAromatic`].
    pub bond_compare: BondCompare,
    /// If `true`, stereochemistry (chirality) must match between atoms. If `false`, chirality
    /// is ignored (default). Prevents matching of enantiomers in MCS.
    pub match_chiral_tag: bool,
    /// If `true`, formal charge must match between atoms. If `false` (default), charge is
    /// ignored -- e.g. a carboxylate `[O-]` matches a neutral `O`. Mirrors `match_chiral_tag`.
    pub match_charge: bool,
    /// If `true`, isotope label must match between atoms (including "no label" only matching
    /// "no label"). If `false` (default), isotope is ignored -- e.g. `[13CH4]` matches `C`.
    /// Mirrors `match_chiral_tag`.
    pub match_isotope: bool,
    /// If `true`, when multiple MCS of the same atom count exist, prefer the one with more bonds.
    /// Defaults to `true` (matching RDKit's default behavior).
    pub maximize_bonds: bool,
}

impl Default for McsConfig {
    fn default() -> Self {
        Self {
            match_bonds: true,
            min_atoms: 1,
            timeout_ms: None,
            ring_matches_ring_only: false,
            complete_rings_only: false,
            atom_compare: AtomCompare::Elements,
            bond_compare: BondCompare::OrderOrAromatic,
            match_chiral_tag: false,
            match_charge: false,
            match_isotope: false,
            maximize_bonds: true,
        }
    }
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Find the Maximum Common Substructure of `mols` using default configuration.
///
/// Returns an empty `QueryMolecule` if `mols` is empty.
pub fn find_mcs(mols: &[&Molecule]) -> QueryMolecule {
    find_mcs_with_config(mols, &McsConfig::default())
}

/// Find the Maximum Common Substructure of `mols` using `config`.
///
/// If `config.timeout_ms` is set and the deadline is reached, returns the
/// best result found so far -- silently, with no way to tell it apart from
/// an exhaustive result. Use [`find_mcs_with_config_checked`] /
/// [`McsOutcome`] when that distinction matters (mirrors
/// `match_vf2.rs`'s own `find_matches_with_rings_and_config` /
/// `_checked` split for the same reason).
pub fn find_mcs_with_config(mols: &[&Molecule], config: &McsConfig) -> QueryMolecule {
    find_mcs_with_config_checked(mols, config).into_query()
}

/// Outcome of an MCS search that may have been cut off by
/// `McsConfig::timeout_ms`.
///
/// The whole point of this type is that `TimedOut` must never be silently
/// treated the same as `Exhaustive`: a cut-off search's result is the best
/// found before the deadline, not proven to be the true maximum common
/// substructure.
#[derive(Debug, Clone)]
pub enum McsOutcome {
    /// The search explored the whole candidate space and confirmed `result`
    /// is the true (config-constrained) maximum common substructure.
    ///
    /// This claim depends on `grow`'s branch-and-bound actually being
    /// complete over connected subgraphs of `mols[0]` -- true since the
    /// include/exclude branch fix (every frontier atom is eventually either
    /// included or explicitly excluded, not just the first one tried at
    /// each node). Before that fix this variant was a real over-claim: a
    /// frontier atom incompatible with some other molecule could silently
    /// kill a branch that would have reached a strictly larger common
    /// substructure by excluding it instead.
    Exhaustive(QueryMolecule),
    /// `timeout_ms` was reached before the search finished; `result` is the
    /// best candidate found so far, not guaranteed optimal.
    TimedOut(QueryMolecule),
}

impl McsOutcome {
    /// Discard the exhaustive/timed-out distinction and return the result
    /// either way -- what [`find_mcs`]/[`find_mcs_with_config`] have always
    /// returned.
    pub fn into_query(self) -> QueryMolecule {
        match self {
            McsOutcome::Exhaustive(q) | McsOutcome::TimedOut(q) => q,
        }
    }

    /// `true` if the search was cut off by `timeout_ms` before finishing.
    pub fn was_timed_out(&self) -> bool {
        matches!(self, McsOutcome::TimedOut(_))
    }
}

/// Like [`find_mcs_with_config`], but reports whether `config.timeout_ms`
/// was reached before the search finished exhaustively -- see [`McsOutcome`].
pub fn find_mcs_with_config_checked(mols: &[&Molecule], config: &McsConfig) -> McsOutcome {
    if mols.is_empty() {
        return McsOutcome::Exhaustive(QueryMolecule::new());
    }
    if mols.len() == 1 {
        return McsOutcome::Exhaustive(molecule_to_query(mols[0]));
    }

    let deadline = config
        .timeout_ms
        .map(|ms| Instant::now() + std::time::Duration::from_millis(ms));

    let ring_sets: Vec<RingSet> = if config.ring_matches_ring_only || config.complete_rings_only {
        mols.iter().map(|m| find_sssr(m)).collect()
    } else {
        Vec::new()
    };

    let n0 = mols[0].atom_count();

    let mut state = McsState {
        mols,
        config,
        best: PartialMapping::empty(mols.len()),
        deadline,
        timed_out: false,
        ring_sets,
        excluded: vec![false; n0],
    };

    'outer: for a0 in 0..n0 {
        if state.timed_out {
            break;
        }
        if let Some(d) = state.deadline
            && Instant::now() >= d
        {
            state.timed_out = true;
            break;
        }

        // Gather compatible atoms from each other molecule for this seed atom.
        let seed_candidates =
            collect_seed_candidates(mols, AtomIdx(a0 as u32), config, &state.ring_sets);

        // Iterate the Cartesian product of candidate lists.
        for seed in CartesianProduct::new(&seed_candidates) {
            if state.timed_out {
                break 'outer;
            }
            if let Some(d) = state.deadline
                && Instant::now() >= d
            {
                state.timed_out = true;
                break 'outer;
            }

            // Prune: upper bound for a single-atom seed is at most min(atom counts by element).
            let ub = upper_bound_for_seeds(mols, AtomIdx(a0 as u32), &seed);
            if ub <= state.best.size {
                continue;
            }

            let mut mapping = PartialMapping::from_seed(mols, AtomIdx(a0 as u32), &seed);
            // Exclusion is a per-seed-tree decision; every push/pop in `grow`
            // unwinds back to all-`false`, but reset explicitly so a fresh
            // seed never inherits state from a prior seed's search tree.
            state.excluded.iter_mut().for_each(|b| *b = false);
            grow(&mut state, &mut mapping);
        }
    }

    // Apply complete_rings_only post-processing: remove partially-covered rings.
    if config.complete_rings_only && !state.ring_sets.is_empty() {
        prune_partial_rings(mols, &mut state.best, &state.ring_sets);
    }

    // Apply min_atoms filter.
    if state.best.size < config.min_atoms {
        return if state.timed_out {
            McsOutcome::TimedOut(QueryMolecule::new())
        } else {
            McsOutcome::Exhaustive(QueryMolecule::new())
        };
    }

    let result = build_query(mols[0], &state.best, config);
    if state.timed_out {
        McsOutcome::TimedOut(result)
    } else {
        McsOutcome::Exhaustive(result)
    }
}

// ---------------------------------------------------------------------------
// Internal data structures
// ---------------------------------------------------------------------------

/// A partial mapping between all molecules' atom indices and query atom indices.
#[derive(Clone)]
struct PartialMapping {
    /// `mol_map[mol_idx][atom_idx]` = `Some(query_idx)` if the atom is mapped.
    mol_map: Vec<Vec<Option<usize>>>,
    /// `query_to_mol[query_idx][mol_idx]` = `AtomIdx`.
    query_to_mol: Vec<Vec<AtomIdx>>,
    /// Number of atoms currently in the mapping.
    size: usize,
    /// Number of bonds in the common subgraph (for bond-match mode).
    bond_count: usize,
}

impl PartialMapping {
    /// Create an empty mapping for `n_mols` molecules (with `n_atoms[i]` atoms each).
    fn empty(n_mols: usize) -> Self {
        // We start with zero-length mol maps; they'll be resized lazily or we just
        // keep them as empty here and they're properly initialised in `from_seed`.
        Self {
            mol_map: vec![Vec::new(); n_mols],
            query_to_mol: Vec::new(),
            size: 0,
            bond_count: 0,
        }
    }

    /// Create a single-atom seed mapping: query atom 0 maps to `a0` in mol[0] and
    /// `seed[i]` in mol[i+1].
    fn from_seed(mols: &[&Molecule], a0: AtomIdx, seed: &[AtomIdx]) -> Self {
        let mut mol_map: Vec<Vec<Option<usize>>> =
            mols.iter().map(|m| vec![None; m.atom_count()]).collect();
        let mut query_to_mol: Vec<Vec<AtomIdx>> = Vec::new();

        // Map the seed tuple to query atom 0.
        let mut row = vec![a0];
        row.extend_from_slice(seed);
        for (mi, &ai) in row.iter().enumerate() {
            mol_map[mi][ai.0 as usize] = Some(0);
        }
        query_to_mol.push(row);

        Self {
            mol_map,
            query_to_mol,
            size: 1,
            bond_count: 0,
        }
    }

    /// Check whether `atom_idx` in molecule `mol_idx` is already mapped.
    fn is_mapped(&self, mol_idx: usize, atom_idx: AtomIdx) -> bool {
        let v = &self.mol_map[mol_idx];
        let idx = atom_idx.0 as usize;
        idx < v.len() && v[idx].is_some()
    }

    /// Look up the query index for `atom_idx` in molecule `mol_idx`.
    fn query_idx_of(&self, mol_idx: usize, atom_idx: AtomIdx) -> Option<usize> {
        let v = &self.mol_map[mol_idx];
        let idx = atom_idx.0 as usize;
        if idx < v.len() { v[idx] } else { None }
    }

    /// Extend the mapping with a new atom tuple (one atom per molecule).
    /// Returns the new query index.
    fn extend(&mut self, atoms: &[AtomIdx], extra_bonds: usize) -> usize {
        let q = self.query_to_mol.len();
        for (mi, &ai) in atoms.iter().enumerate() {
            self.mol_map[mi][ai.0 as usize] = Some(q);
        }
        self.query_to_mol.push(atoms.to_vec());
        self.size += 1;
        self.bond_count += extra_bonds;
        q
    }

    /// Retract the last-added atom tuple (must be the one with the highest query index).
    fn retract(&mut self, atoms: &[AtomIdx], extra_bonds: usize) {
        for (mi, &ai) in atoms.iter().enumerate() {
            self.mol_map[mi][ai.0 as usize] = None;
        }
        self.query_to_mol.pop();
        self.size -= 1;
        self.bond_count -= extra_bonds;
    }
}

// ---------------------------------------------------------------------------
// Search state
// ---------------------------------------------------------------------------

struct McsState<'a> {
    mols: &'a [&'a Molecule],
    config: &'a McsConfig,
    best: PartialMapping,
    deadline: Option<Instant>,
    timed_out: bool,
    ring_sets: Vec<RingSet>,
    /// `excluded[a.0]` = `true` if mol[0] atom `a` has been decided OUT of
    /// the current search-tree branch (see the include/exclude split in
    /// `grow`). Reset to all-`false` before growing from each new seed --
    /// exclusion is a per-branch decision, not global across seeds.
    excluded: Vec<bool>,
}

// ---------------------------------------------------------------------------
// Core search: grow
// ---------------------------------------------------------------------------

/// Canonical key for tie-breaking equal-size(-and-bond-count) mappings: the
/// mols[0]-side atom indices in query-discovery order. Comparing this key
/// (rather than "whichever was found first") makes the winner a pure
/// function of the mapping's content, so it doesn't depend on which order
/// branch-and-bound happens to visit tied candidates in.
///
/// Does NOT achieve full spelling-invariance: mols[0]'s AtomIdx numbering
/// IS its parse order, so a differently-spelled-but-equivalent mols[0] can
/// still produce a different winner on a genuine (non-automorphic,
/// element-different) tie -- confirmed empirically (see
/// scratch_element_different_tie_probe / commit history), same residual
/// tier as chematic-chem/standardize.rs's #9/#10 and the reverted #4 iupac
/// tie-break this round. This key still removes real internal
/// nondeterminism (DFS-visit-order dependence for a FIXED mols[0]
/// labeling) and is kept for that reason, not because it solves the
/// general problem.
fn mcs_tiebreak_key(mapping: &PartialMapping) -> Vec<u32> {
    mapping.query_to_mol.iter().map(|row| row[0].0).collect()
}

fn grow(state: &mut McsState<'_>, mapping: &mut PartialMapping) {
    // Update best if this mapping is larger, or same size but more bonds (if enabled),
    // or a same-size-and-bond-count tie broken by a canonical key (not "first found in
    // DFS order", which depends on mols[0]'s atom insertion order -- the same bug shape
    // as scaffold.rs Rule 8 / is_fused_ring found elsewhere in this project).
    let is_better = mapping.size > state.best.size
        || (mapping.size == state.best.size
            && ((state.config.maximize_bonds && mapping.bond_count > state.best.bond_count)
                || ((!state.config.maximize_bonds
                    || mapping.bond_count == state.best.bond_count)
                    && mcs_tiebreak_key(mapping) < mcs_tiebreak_key(&state.best))));
    if is_better {
        state.best = mapping.clone();
    }

    // Check timeout.
    if let Some(d) = state.deadline
        && Instant::now() >= d
    {
        state.timed_out = true;
        return;
    }
    if state.timed_out {
        return;
    }

    // Upper bound pruning: compute max additional atoms reachable from current state.
    let additional_ub = upper_bound_additional(state.mols, mapping, &state.excluded);
    if mapping.size + additional_ub <= state.best.size {
        return;
    }

    // Find frontier candidates: unmapped, non-excluded neighbors of already-mapped
    // atoms in mol[0]. For each such neighbor `n0`, decide (in a moment) whether to
    // include it in the mapping or exclude it from this branch entirely.
    let mol0 = state.mols[0];

    let mut frontier: Vec<AtomIdx> = Vec::new();
    for row in &mapping.query_to_mol {
        let a0 = row[0];
        for (nb, _) in mol0.neighbors(a0) {
            if !mapping.is_mapped(0, nb)
                && !state.excluded[nb.0 as usize]
                && !frontier.contains(&nb)
            {
                frontier.push(nb);
            }
        }
    }

    if frontier.is_empty() {
        return;
    }

    // Only the first frontier atom is branched on per call (the "connected growth"
    // property), but BOTH outcomes for it are explored -- include it now, or exclude
    // it and let a later frontier atom be tried in its place. Without the exclude
    // branch, a `n0` that happens to be incompatible with every other molecule kills
    // the whole subtree even when skipping it (and mapping a different frontier atom
    // instead) would reach a strictly larger common substructure -- see
    // scratch_repro_frontier0_only_incompleteness.
    let n0 = frontier[0];

    // Branch 1: include n0. For each other molecule, find compatible unmapped atoms
    // that are:
    // 1. Atom-compatible with n0
    // 2. Adjacent to every mol[i]-atom that corresponds to a mapped neighbor of n0 in mol[0]
    // 3. Bond-compatible with each such adjacency
    let candidates_per_mol =
        build_frontier_candidates(state.mols, mapping, n0, state.config, &state.ring_sets);

    // Cartesian product across molecules 1..n.
    for tuple in CartesianProduct::new(&candidates_per_mol) {
        if state.timed_out {
            return;
        }

        // Build the full atom tuple: [n0, tuple[0], tuple[1], ...]
        let mut atoms: Vec<AtomIdx> = Vec::with_capacity(state.mols.len());
        atoms.push(n0);
        atoms.extend_from_slice(&tuple);

        // Count new bonds (to already-mapped atoms in mol[0]) for bookkeeping.
        let extra_bonds = count_new_bonds(mol0, mapping, n0);

        mapping.extend(&atoms, extra_bonds);
        grow(state, mapping);
        mapping.retract(&atoms, extra_bonds);
    }

    // Branch 2: exclude n0 from this branch and recurse so the next frontier atom
    // gets a turn. Reset on the way back out -- exclusion only applies within this
    // subtree, not to sibling branches elsewhere in the search.
    if !state.timed_out {
        state.excluded[n0.0 as usize] = true;
        grow(state, mapping);
        state.excluded[n0.0 as usize] = false;
    }
}

// ---------------------------------------------------------------------------
// Frontier candidate generation
// ---------------------------------------------------------------------------

/// For each molecule i >= 1, collect the atoms that can be paired with `n0` from mol[0]
/// to extend the mapping.
///
/// A candidate atom `ni` in mol[i] is valid if:
/// - It is not yet mapped.
/// - It is atom-compatible with `n0`.
/// - For every already-mapped neighbor `m0` of `n0` in mol[0], atom `ni` is bonded to
///   `mapping.query_to_mol[q_m][i]` in mol[i] with a compatible bond order (when
///   `match_bonds=true`).
fn build_frontier_candidates(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    n0: AtomIdx,
    config: &McsConfig,
    ring_sets: &[RingSet],
) -> Vec<Vec<AtomIdx>> {
    let mol0 = mols[0];
    let atom0 = mol0.atom(n0);
    let n0_in_ring = !ring_sets.is_empty() && ring_sets[0].contains_atom(n0);

    // Collect mapped neighbors of n0 in mol[0] along with the bond order.
    let mut mapped_neighbors: Vec<(usize, AtomIdx, BondOrder)> = Vec::new();
    for (nb, bidx) in mol0.neighbors(n0) {
        if let Some(q) = mapping.query_idx_of(0, nb) {
            let bond = mol0.bond(bidx);
            mapped_neighbors.push((q, nb, bond.order));
        }
    }

    let mut result: Vec<Vec<AtomIdx>> = Vec::with_capacity(mols.len() - 1);

    for mi in 1..mols.len() {
        let mol_i = mols[mi];
        let mut candidates: Vec<AtomIdx> = Vec::new();

        'atom: for (ai, atom_i) in mol_i.atoms() {
            // Must not be mapped.
            if mapping.is_mapped(mi, ai) {
                continue;
            }
            // Must be atom-compatible.
            if !atoms_compatible(
                atom0,
                atom_i,
                &config.atom_compare,
                config.match_chiral_tag,
                config.match_charge,
                config.match_isotope,
            ) {
                continue;
            }
            // ring_matches_ring_only: ring atoms must match ring atoms only.
            if config.ring_matches_ring_only {
                let ai_in_ring = ring_sets[mi].contains_atom(ai);
                if n0_in_ring != ai_in_ring {
                    continue;
                }
            }
            // Must be bonded to every corresponding mapped neighbor with compatible bond.
            for &(q, _m0, bond_order_0) in &mapped_neighbors {
                let m_i = mapping.query_to_mol[q][mi];
                match mol_i.bond_between(ai, m_i) {
                    None => continue 'atom,
                    Some((_bidx, bond_entry)) => {
                        if config.match_bonds
                            && !bonds_compatible(
                                bond_order_0,
                                bond_entry.order,
                                &config.bond_compare,
                            )
                        {
                            continue 'atom;
                        }
                    }
                }
            }
            candidates.push(ai);
        }

        result.push(candidates);
    }

    result
}

// ---------------------------------------------------------------------------
// Seed collection
// ---------------------------------------------------------------------------

/// For each molecule i >= 1, collect atoms compatible with `a0` in mol[0].
fn collect_seed_candidates(
    mols: &[&Molecule],
    a0: AtomIdx,
    config: &McsConfig,
    ring_sets: &[RingSet],
) -> Vec<Vec<AtomIdx>> {
    let atom0 = mols[0].atom(a0);
    let a0_in_ring = !ring_sets.is_empty() && ring_sets[0].contains_atom(a0);
    let mut result = Vec::with_capacity(mols.len() - 1);
    for mi in 1..mols.len() {
        let cands: Vec<AtomIdx> = mols[mi]
            .atoms()
            .filter(|(ai, a)| {
                if !atoms_compatible(
                    atom0,
                    a,
                    &config.atom_compare,
                    config.match_chiral_tag,
                    config.match_charge,
                    config.match_isotope,
                ) {
                    return false;
                }
                if config.ring_matches_ring_only {
                    let ai_in_ring = ring_sets[mi].contains_atom(*ai);
                    if a0_in_ring != ai_in_ring {
                        return false;
                    }
                }
                true
            })
            .map(|(idx, _)| idx)
            .collect();
        result.push(cands);
    }
    result
}

// ---------------------------------------------------------------------------
// Upper bound heuristics
// ---------------------------------------------------------------------------

/// Upper bound on the total MCS size for a given seed (used at the seed level).
fn upper_bound_for_seeds(mols: &[&Molecule], a0: AtomIdx, seed: &[AtomIdx]) -> usize {
    // Count atoms of each element in each molecule.
    // Minimum across molecules gives max atoms of that element that can be matched.
    let mut counts: HashMap<u8, Vec<usize>> = HashMap::new();
    for (mi, mol) in mols.iter().enumerate() {
        let exclude = if mi == 0 { a0 } else { seed[mi - 1] };
        for (ai, atom) in mol.atoms() {
            if ai == exclude {
                continue; // seed atom already counts; we'll add 1 at the end
            }
            let en = atom.element.atomic_number();
            counts.entry(en).or_insert_with(|| vec![0; mols.len()])[mi] += 1;
        }
    }
    let ub: usize = counts.values().map(|v| *v.iter().min().unwrap_or(&0)).sum();
    ub + 1 // +1 for the seed atom itself
}

/// Upper bound on *additional* atoms reachable from the current mapping.
///
/// For each element, count unmapped atoms in each molecule (mol[0] atoms already
/// excluded from this branch don't count -- they can never be reached); the maximum
/// additional atoms of that element = min over all molecules.
fn upper_bound_additional(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    excluded: &[bool],
) -> usize {
    let mut counts: HashMap<u8, Vec<usize>> = HashMap::new();

    for (mi, mol) in mols.iter().enumerate() {
        for (ai, atom) in mol.atoms() {
            if mapping.is_mapped(mi, ai) {
                continue;
            }
            if mi == 0 && excluded[ai.0 as usize] {
                continue;
            }
            let en = atom.element.atomic_number();
            let entry = counts.entry(en).or_insert_with(|| vec![0; mols.len()]);
            entry[mi] += 1;
        }
    }

    counts.values().map(|v| *v.iter().min().unwrap_or(&0)).sum()
}

// ---------------------------------------------------------------------------
// Compatibility helpers
// ---------------------------------------------------------------------------

/// Atom compatibility according to `compare` mode.
fn atoms_compatible(
    a: &chematic_core::Atom,
    b: &chematic_core::Atom,
    compare: &AtomCompare,
    match_chiral: bool,
    match_charge: bool,
    match_isotope: bool,
) -> bool {
    let base = match compare {
        AtomCompare::Elements => a.element.atomic_number() == b.element.atomic_number(),
        AtomCompare::AnyHeavyAtom => a.element.atomic_number() > 1 && b.element.atomic_number() > 1,
        AtomCompare::Any => true,
    };
    if !base {
        return false;
    }
    if match_chiral && a.chirality != b.chirality {
        return false;
    }
    if match_charge && a.charge != b.charge {
        return false;
    }
    if match_isotope && a.isotope != b.isotope {
        return false;
    }
    true
}

/// Bond compatibility according to `compare` mode.
fn bonds_compatible(a: BondOrder, b: BondOrder, compare: &BondCompare) -> bool {
    match compare {
        BondCompare::OrderOrAromatic => normalize_bond(a) == normalize_bond(b),
        BondCompare::Any => true,
    }
}

/// Normalize Up/Down stereo bonds to Single for comparison purposes.
fn normalize_bond(o: BondOrder) -> BondOrder {
    match o {
        BondOrder::Up | BondOrder::Down => BondOrder::Single,
        other => other,
    }
}

// ---------------------------------------------------------------------------
// Bond counting helper
// ---------------------------------------------------------------------------

/// Count the number of already-mapped neighbors of `n0` in mol[0].
fn count_new_bonds(mol0: &Molecule, mapping: &PartialMapping, n0: AtomIdx) -> usize {
    mol0.neighbors(n0)
        .filter(|(nb, _)| mapping.is_mapped(0, *nb))
        .count()
}

// ---------------------------------------------------------------------------
// complete_rings_only post-processing
// ---------------------------------------------------------------------------

/// Remove atoms from `mapping` that would leave a ring partially covered.
///
/// A ring is "partially covered" if at least one but not all of its atoms appear in the
/// mapping.  Removing those atoms can expose further partial rings, so the pruning
/// iterates until convergence.
fn prune_partial_rings(mols: &[&Molecule], mapping: &mut PartialMapping, ring_sets: &[RingSet]) {
    loop {
        // Collect query indices that must be removed.
        let mut to_remove: Vec<usize> = Vec::new();

        // Only check mol[0]'s rings. The MCS is represented as a subgraph of mol[0],
        // so ring completeness is determined by mol[0]'s ring structure, not by how
        // the mapping happens to cover rings in the other (possibly larger) molecules.
        let ring_set = &ring_sets[0];
        for ring in ring_set.rings() {
            // Determine which ring atoms are currently mapped in mol[0].
            let mapped_in_ring: Vec<AtomIdx> = ring
                .iter()
                .copied()
                .filter(|&ai| mapping.is_mapped(0, ai))
                .collect();

            // Partial coverage: some but not all ring atoms are mapped.
            if !mapped_in_ring.is_empty() && mapped_in_ring.len() < ring.len() {
                for ai in mapped_in_ring {
                    if let Some(qi) = mapping.query_idx_of(0, ai)
                        && !to_remove.contains(&qi)
                    {
                        to_remove.push(qi);
                    }
                }
            }
        }

        if to_remove.is_empty() {
            break;
        }

        // Remove by rebuilding the mapping without the marked query indices.
        // Build a new query_to_mol (renumbered) and matching mol_maps.
        let n_mols = mols.len();
        let mut new_query_to_mol: Vec<Vec<AtomIdx>> = Vec::new();
        let mut new_mol_map: Vec<Vec<Option<usize>>> =
            mols.iter().map(|m| vec![None; m.atom_count()]).collect();

        for (old_qi, row) in mapping.query_to_mol.iter().enumerate() {
            if to_remove.contains(&old_qi) {
                continue;
            }
            let new_qi = new_query_to_mol.len();
            new_query_to_mol.push(row.clone());
            for mi in 0..n_mols {
                let ai = row[mi];
                new_mol_map[mi][ai.0 as usize] = Some(new_qi);
            }
        }

        mapping.query_to_mol = new_query_to_mol;
        mapping.mol_map = new_mol_map;
        mapping.size = mapping.query_to_mol.len();
    }
}

// ---------------------------------------------------------------------------
// Build output QueryMolecule
// ---------------------------------------------------------------------------

fn build_query(mol0: &Molecule, mapping: &PartialMapping, config: &McsConfig) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    // Add one query atom per mapped position. Aromaticity is deliberately NOT
    // encoded here as a per-atom constraint (matches RDKit's own `CompareElements`
    // representation, confirmed via live oracle): it's already fully captured by
    // the aromatic bond queries added below, and hard-coding mol0's own atom.aromatic
    // would make the query fail to re-match a different input molecule whose
    // corresponding atom disagrees on aromaticity (a real, confirmed regression when
    // AtomCompare::Elements permits cross-aromaticity matches).
    for row in &mapping.query_to_mol {
        let a0 = row[0];
        let atom = mol0.atom(a0);
        let an = atom.element.atomic_number();

        let query = AtomQuery::Primitive(AtomPrimitive::AtomicNum(an));
        qmol.add_atom(query);
    }

    // Add bonds between mapped atoms (based on mol[0] topology).
    for (q1, row1) in mapping.query_to_mol.iter().enumerate() {
        let a1 = row1[0];
        for (q2, row2) in mapping.query_to_mol.iter().enumerate() {
            if q2 <= q1 {
                continue;
            }
            let a2 = row2[0];
            if let Some((_bidx, bond_entry)) = mol0.bond_between(a1, a2) {
                let bq = bond_order_to_query(bond_entry.order, config.match_bonds);
                qmol.add_bond(q1, q2, bq);
            }
        }
    }

    qmol
}

/// Convert a target-molecule `BondOrder` into the corresponding `BondQuery`
/// used in a `QueryMolecule`.
///
/// When `match_bonds` is false, every bond becomes `Any` (atom-only MCS).
/// Up/Down/Dative stereo-like bonds collapse to Single. Query bond variants
/// become their closest SMARTS query equivalent.
fn bond_order_to_query(order: BondOrder, match_bonds: bool) -> BondQuery {
    if !match_bonds {
        return BondQuery::Primitive(BondPrimitive::Any);
    }
    match normalize_bond(order) {
        BondOrder::Single | BondOrder::Dative => BondQuery::Primitive(BondPrimitive::Single),
        BondOrder::Double => BondQuery::Primitive(BondPrimitive::Double),
        BondOrder::Triple => BondQuery::Primitive(BondPrimitive::Triple),
        BondOrder::Aromatic => BondQuery::Primitive(BondPrimitive::Aromatic),
        BondOrder::Quadruple | BondOrder::Zero | BondOrder::QueryAny => {
            BondQuery::Primitive(BondPrimitive::Any)
        }
        BondOrder::QuerySingleOrDouble => BondQuery::Or(
            Box::new(BondQuery::Primitive(BondPrimitive::Single)),
            Box::new(BondQuery::Primitive(BondPrimitive::Double)),
        ),
        BondOrder::QuerySingleOrAromatic => BondQuery::Or(
            Box::new(BondQuery::Primitive(BondPrimitive::Single)),
            Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
        ),
        BondOrder::QueryDoubleOrAromatic => BondQuery::Or(
            Box::new(BondQuery::Primitive(BondPrimitive::Double)),
            Box::new(BondQuery::Primitive(BondPrimitive::Aromatic)),
        ),
        BondOrder::Up | BondOrder::Down => BondQuery::Primitive(BondPrimitive::Single),
    }
}

// ---------------------------------------------------------------------------
// Convert a single molecule to a QueryMolecule (for the single-mol case)
// ---------------------------------------------------------------------------

fn molecule_to_query(mol: &Molecule) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for (_, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let query = AtomQuery::Primitive(AtomPrimitive::AtomicNum(an));
        qmol.add_atom(query);
    }

    for (_, bond) in mol.bonds() {
        let bq = bond_order_to_query(bond.order, true);
        qmol.add_bond(bond.atom1.0 as usize, bond.atom2.0 as usize, bq);
    }

    qmol
}

// ---------------------------------------------------------------------------
// Cartesian product iterator
// ---------------------------------------------------------------------------

/// Iterator over the Cartesian product of a slice of `Vec<T>` slices.
///
/// `CartesianProduct::new(&[vec![a,b], vec![x,y]])` yields `[a,x]`, `[a,y]`, `[b,x]`, `[b,y]`.
struct CartesianProduct<'a, T> {
    lists: &'a [Vec<T>],
    indices: Vec<usize>,
    done: bool,
}

impl<'a, T: Copy> CartesianProduct<'a, T> {
    fn new(lists: &'a [Vec<T>]) -> Self {
        let done = lists.iter().any(|l| l.is_empty());
        let indices = vec![0; lists.len()];
        Self {
            lists,
            indices,
            done,
        }
    }
}

impl<'a, T: Copy> Iterator for CartesianProduct<'a, T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Vec<T>> {
        if self.done {
            return None;
        }

        // Emit current combination.
        let item: Vec<T> = self
            .lists
            .iter()
            .zip(self.indices.iter())
            .map(|(l, &i)| l[i])
            .collect();

        // Advance indices (rightmost first).
        let mut carry = true;
        for k in (0..self.lists.len()).rev() {
            if carry {
                self.indices[k] += 1;
                if self.indices[k] < self.lists[k].len() {
                    carry = false;
                } else {
                    self.indices[k] = 0;
                }
            }
        }
        if carry {
            self.done = true;
        }

        Some(item)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_empty_input() {
        let result = find_mcs(&[]);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn test_single_molecule() {
        let benzene = parse("c1ccccc1").unwrap();
        let result = find_mcs(&[&benzene]);
        assert_eq!(result.atom_count(), 6);
    }

    #[test]
    fn test_ethane_propane() {
        // CC and CCC → MCS should be CC (2 atoms)
        let a = parse("CC").unwrap();
        let b = parse("CCC").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert!(
            result.atom_count() >= 2,
            "Expected >= 2 atoms, got {}",
            result.atom_count()
        );
    }

    #[test]
    fn test_methane_ethanol() {
        // C and CCO → MCS should be C (1 atom)
        let a = parse("C").unwrap();
        let b = parse("CCO").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert!(result.atom_count() >= 1);
    }

    #[test]
    fn test_benzene_toluene() {
        // c1ccccc1 and Cc1ccccc1 → MCS should be the benzene ring (6 atoms)
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert_eq!(
            result.atom_count(),
            6,
            "Expected benzene ring (6 atoms), got {}",
            result.atom_count()
        );
    }

    #[test]
    fn test_no_common_atoms() {
        // Only N atoms vs only O atoms → no MCS
        let a = parse("N").unwrap(); // ammonia
        let b = parse("O").unwrap(); // water
        let result = find_mcs(&[&a, &b]);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn test_identical_molecules() {
        let a = parse("CC").unwrap();
        let b = parse("CC").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert_eq!(result.atom_count(), 2);
    }

    #[test]
    fn test_match_bonds_false() {
        // With match_bonds=false, single and double bond C-C should match
        let config = McsConfig {
            match_bonds: false,
            ..McsConfig::default()
        };
        let a = parse("CC").unwrap();
        let b = parse("C=C").unwrap();
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert!(result.atom_count() >= 2);
    }

    #[test]
    fn test_min_atoms_filter() {
        // CC and CCC have MCS of size 2, but if min_atoms=5, result should be empty
        let config = McsConfig {
            min_atoms: 5,
            ..McsConfig::default()
        };
        let a = parse("CC").unwrap();
        let b = parse("CCC").unwrap();
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn test_timeout_does_not_panic() {
        let config = McsConfig {
            timeout_ms: Some(1),
            ..McsConfig::default()
        }; // 1 ms
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        // Should return without panic (may return partial result)
        let _ = find_mcs_with_config(&[&a, &b], &config);
    }

    #[test]
    fn test_checked_reports_exhaustive_without_timeout() {
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let unchecked_result = find_mcs_with_config(&[&a, &b], &McsConfig::default());
        let checked = find_mcs_with_config_checked(&[&a, &b], &McsConfig::default());
        assert!(
            !checked.was_timed_out(),
            "no timeout_ms set: search must report Exhaustive"
        );
        assert!(matches!(checked, McsOutcome::Exhaustive(_)));
        assert_eq!(
            checked.into_query().atom_count(),
            unchecked_result.atom_count(),
            "find_mcs_with_config_checked's result must match find_mcs_with_config's"
        );
    }

    #[test]
    fn test_checked_reports_timed_out_with_zero_deadline() {
        // A 0ms deadline is reached on the very first check, guaranteeing a
        // cut-off search (unlike test_timeout_does_not_panic's 1ms, which
        // isn't guaranteed to actually trigger on a small molecule pair --
        // that test only checks "no panic", this one checks the outcome tag).
        let config = McsConfig {
            timeout_ms: Some(0),
            ..McsConfig::default()
        };
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let checked = find_mcs_with_config_checked(&[&a, &b], &config);
        assert!(
            checked.was_timed_out(),
            "0ms deadline must be reported as TimedOut, not silently treated as Exhaustive"
        );
        assert!(matches!(checked, McsOutcome::TimedOut(_)));
    }

    #[test]
    fn test_result_matches_all_inputs() {
        // The MCS result should be findable as a substructure in all input molecules
        use crate::find_matches;
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let result = find_mcs(&[&a, &b]);
        if result.atom_count() > 0 {
            assert!(!find_matches(&result, &a).is_empty());
            assert!(!find_matches(&result, &b).is_empty());
        }
    }

    #[test]
    fn test_aspirin_benzoic_acid() {
        // CC(=O)Oc1ccccc1C(=O)O and OC(=O)c1ccccc1
        // Common: benzoic acid core → >= 9 atoms
        let a = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let b = parse("OC(=O)c1ccccc1").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert!(
            result.atom_count() >= 7,
            "Expected >= 7 atoms for aspirin/benzoic acid MCS, got {}",
            result.atom_count()
        );
    }

    // -----------------------------------------------------------------------
    // ring_matches_ring_only tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_ring_matches_ring_only_blocks_ring_to_nonring() {
        // cyclohexane (all ring C's) vs hexane (all chain C's)
        // Without constraint: MCS = 5 atoms (C-C-C-C-C path, both have it)
        // With ring_matches_ring_only: ring C (cyclohexane) can't match non-ring C (hexane) → MCS = 0
        let a = parse("C1CCCCC1").unwrap();
        let b = parse("CCCCCC").unwrap();

        // Confirm default finds a path match.
        let default_result = find_mcs(&[&a, &b]);
        assert!(
            default_result.atom_count() > 0,
            "Default MCS should be non-zero"
        );

        let config = McsConfig {
            ring_matches_ring_only: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert_eq!(
            result.atom_count(),
            0,
            "ring_matches_ring_only: ring C must not match chain C, expected 0 got {}",
            result.atom_count()
        );
    }

    #[test]
    fn test_ring_matches_ring_only_preserves_ring_only_mcs() {
        // benzene vs toluene: ring atoms should still match each other
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let config = McsConfig {
            ring_matches_ring_only: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert_eq!(
            result.atom_count(),
            6,
            "benzene ring should still be the MCS"
        );
    }

    // -----------------------------------------------------------------------
    // complete_rings_only tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_complete_rings_only_quinoline_series() {
        // With complete_rings_only: partial rings are removed from the MCS result.
        // In the quinoline series the exocyclic C would create a partial ring if it
        // happened to be included; enabling complete_rings_only ensures the quinoline
        // scaffold (10 atoms) is returned intact.
        let a = parse("c1ccc2nc(CC)ccc2c1").unwrap();
        let b = parse("c1ccc2nc(CO)ccc2c1").unwrap();
        let c = parse("c1ccc2nc(CN)ccc2c1").unwrap();

        let config = McsConfig {
            complete_rings_only: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&a, &b, &c], &config);
        // Result must be >= 10 atoms (quinoline) and no partial ring
        assert!(
            result.atom_count() >= 10,
            "Expected at least quinoline scaffold (10) with complete_rings_only, got {}",
            result.atom_count()
        );
    }

    #[test]
    fn test_complete_rings_only_does_not_break_full_ring_mcs() {
        // benzene vs toluene: full benzene ring should survive complete_rings_only
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let config = McsConfig {
            complete_rings_only: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert_eq!(
            result.atom_count(),
            6,
            "Full benzene ring should be preserved"
        );
    }

    // -----------------------------------------------------------------------
    // Diagnostic: verify complete_rings_only ACTUALLY fires pruning
    // -----------------------------------------------------------------------

    #[test]
    fn test_complete_rings_only_pruning_fires() {
        // Benzene vs naphthalene: raw MCS without constraint is 6 ring atoms (one full ring).
        // With complete_rings_only, fused-ring bridgeheads are in two rings each; verify
        // the full benzene ring is preserved (6 atoms), not over-pruned.
        let a = parse("c1ccccc1").unwrap();
        let b = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene

        let default_result = find_mcs(&[&a, &b]);
        let config = McsConfig {
            complete_rings_only: true,
            ..McsConfig::default()
        };
        let pruned_result = find_mcs_with_config(&[&a, &b], &config);

        eprintln!(
            "benzene vs naphthalene — raw MCS: {}, complete_rings_only: {}",
            default_result.atom_count(),
            pruned_result.atom_count()
        );
        // Both should be 6 (one complete benzene ring)
        assert_eq!(
            default_result.atom_count(),
            6,
            "Default MCS should be benzene ring (6 atoms)"
        );
        assert_eq!(
            pruned_result.atom_count(),
            6,
            "complete_rings_only should preserve the full 6-ring (no partial ring)"
        );
    }

    #[test]
    fn test_complete_rings_only_prunes_partial_fused_ring() {
        // C1Cc2ccccc21 is benzocyclobutene (4-ring fused to 6-ring), not indane.
        // Toluene: Cc1ccccc1.
        // Raw MCS = 7 atoms: the 6 aromatic C's + aliphatic atom1 from benzocyclobutene
        // matched to the methyl C of toluene (both non-aromatic, atoms_compatible passes).
        // Atom1 of benzocyclobutene is in the 4-membered ring — partial ring in mol[0].
        // complete_rings_only fires: remove atom1 + bridgehead atoms (shared with 4-ring) →
        // the 6-ring becomes partially covered → cascading removal → result = 0.
        let benzocyclobutene = parse("C1Cc2ccccc21").unwrap();
        let toluene = parse("Cc1ccccc1").unwrap();

        let default_result = find_mcs(&[&benzocyclobutene, &toluene]);
        let config = McsConfig {
            complete_rings_only: true,
            ..McsConfig::default()
        };
        let pruned_result = find_mcs_with_config(&[&benzocyclobutene, &toluene], &config);

        // Raw MCS includes a partial 4-ring; cascading pruning removes everything.
        assert!(
            default_result.atom_count() > 0,
            "default MCS should be non-empty"
        );
        assert_eq!(
            pruned_result.atom_count(),
            0,
            "complete_rings_only must prune to 0 when partial ring cascades"
        );
    }

    #[test]
    fn test_complete_rings_only_quinoline_raw_vs_pruned() {
        // Verify quinoline series: check whether complete_rings_only actually changes
        // the atom count vs default (determines if pruning fires at all for this case)
        let a = parse("c1ccc2nc(CC)ccc2c1").unwrap();
        let b = parse("c1ccc2nc(CO)ccc2c1").unwrap();
        let c = parse("c1ccc2nc(CN)ccc2c1").unwrap();

        let default_result = find_mcs(&[&a, &b, &c]);
        let config = McsConfig {
            complete_rings_only: true,
            ..McsConfig::default()
        };
        let pruned_result = find_mcs_with_config(&[&a, &b, &c], &config);

        eprintln!(
            "quinoline series raw MCS: {}, with complete_rings_only: {}",
            default_result.atom_count(),
            pruned_result.atom_count()
        );
        // Document whether complete_rings_only triggered a change
        // (if equal, pruning was not triggered for this series)
    }

    #[test]
    fn test_ring_matches_ring_only_and_complete_rings_only_combined() {
        // Both flags together: quinoline series should return the full quinoline scaffold
        let a = parse("c1ccc2nc(CC)ccc2c1").unwrap();
        let b = parse("c1ccc2nc(CO)ccc2c1").unwrap();
        let config = McsConfig {
            ring_matches_ring_only: true,
            complete_rings_only: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert!(
            result.atom_count() >= 10,
            "Combined flags should yield at least the quinoline scaffold, got {}",
            result.atom_count()
        );
    }

    // -----------------------------------------------------------------------
    // match_chiral_tag tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_match_chiral_tag_false_matches_enantiomers() {
        // R-alanine vs S-alanine: with match_chiral_tag=false (default),
        // should match all atoms (ignoring chirality)
        let r_ala = parse("[C@H](C)(N)C(=O)O").unwrap();
        let s_ala = parse("[C@@H](C)(N)C(=O)O").unwrap();
        let default_result = find_mcs(&[&r_ala, &s_ala]);
        assert_eq!(
            default_result.atom_count(),
            r_ala.atom_count(),
            "Default (match_chiral_tag=false) should match all atoms of enantiomers"
        );
    }

    #[test]
    fn test_match_chiral_tag_true_blocks_stereocenters() {
        // R-alanine vs S-alanine: with match_chiral_tag=true,
        // the stereocenter C should not match (different chirality).
        // However, non-chiral atoms (N, O, carboxyl C) can still match.
        // So MCS should be smaller than without the constraint.
        let r_ala = parse("[C@H](C)(N)C(=O)O").unwrap();
        let s_ala = parse("[C@@H](C)(N)C(=O)O").unwrap();

        let default_result = find_mcs(&[&r_ala, &s_ala]);
        let config = McsConfig {
            match_chiral_tag: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&r_ala, &s_ala], &config);

        assert!(
            result.atom_count() < default_result.atom_count(),
            "match_chiral_tag=true should find smaller MCS than default: {} vs {}",
            result.atom_count(),
            default_result.atom_count()
        );
    }

    #[test]
    fn test_match_chiral_tag_true_matches_same_stereoisomer() {
        // R-alanine vs R-alanine: with match_chiral_tag=true,
        // should still match all atoms (same stereochemistry)
        let r_ala1 = parse("[C@H](C)(N)C(=O)O").unwrap();
        let r_ala2 = parse("[C@H](C)(N)C(=O)O").unwrap();
        let config = McsConfig {
            match_chiral_tag: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&r_ala1, &r_ala2], &config);
        assert_eq!(
            result.atom_count(),
            r_ala1.atom_count(),
            "match_chiral_tag=true should still match molecules with same stereochemistry"
        );
    }

    // -----------------------------------------------------------------------
    // match_charge / match_isotope tests (mirror match_chiral_tag exactly)
    // -----------------------------------------------------------------------

    #[test]
    fn test_match_charge_false_matches_charge_states() {
        // Acetate vs neutral acetic acid: with match_charge=false (default),
        // the carboxylate/carboxylic-acid oxygen should still match.
        let acetate = parse("CC(=O)[O-]").unwrap();
        let acetic_acid = parse("CC(=O)O").unwrap();
        let default_result = find_mcs(&[&acetate, &acetic_acid]);
        assert_eq!(
            default_result.atom_count(),
            acetate.atom_count(),
            "Default (match_charge=false) should match all atoms across charge states"
        );
    }

    #[test]
    fn test_match_charge_true_blocks_differently_charged_atoms() {
        // Acetate vs neutral acetic acid: with match_charge=true, the [O-]
        // atom must not match the neutral O, so MCS should be smaller.
        let acetate = parse("CC(=O)[O-]").unwrap();
        let acetic_acid = parse("CC(=O)O").unwrap();
        let default_result = find_mcs(&[&acetate, &acetic_acid]);
        let config = McsConfig {
            match_charge: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&acetate, &acetic_acid], &config);
        assert!(
            result.atom_count() < default_result.atom_count(),
            "match_charge=true should find smaller MCS than default: {} vs {}",
            result.atom_count(),
            default_result.atom_count()
        );
    }

    #[test]
    fn test_match_charge_true_matches_same_charge_state() {
        // Acetate vs acetate: with match_charge=true, should still match
        // all atoms (identical charge state).
        let acetate1 = parse("CC(=O)[O-]").unwrap();
        let acetate2 = parse("CC(=O)[O-]").unwrap();
        let config = McsConfig {
            match_charge: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&acetate1, &acetate2], &config);
        assert_eq!(
            result.atom_count(),
            acetate1.atom_count(),
            "match_charge=true should still match molecules with the same charge state"
        );
    }

    #[test]
    fn test_match_isotope_false_matches_isotope_labels() {
        // 13C-labeled methanol vs unlabeled methanol: with match_isotope=false
        // (default), the labeled carbon should still match.
        let labeled = parse("[13CH3]O").unwrap();
        let unlabeled = parse("CO").unwrap();
        let default_result = find_mcs(&[&labeled, &unlabeled]);
        assert_eq!(
            default_result.atom_count(),
            labeled.atom_count(),
            "Default (match_isotope=false) should match all atoms across isotope labels"
        );
    }

    #[test]
    fn test_match_isotope_true_blocks_differently_labeled_atoms() {
        // 13C-labeled methanol vs unlabeled methanol: with match_isotope=true,
        // the labeled carbon must not match the unlabeled one, so MCS should
        // be smaller.
        let labeled = parse("[13CH3]O").unwrap();
        let unlabeled = parse("CO").unwrap();
        let default_result = find_mcs(&[&labeled, &unlabeled]);
        let config = McsConfig {
            match_isotope: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&labeled, &unlabeled], &config);
        assert!(
            result.atom_count() < default_result.atom_count(),
            "match_isotope=true should find smaller MCS than default: {} vs {}",
            result.atom_count(),
            default_result.atom_count()
        );
    }

    #[test]
    fn test_match_isotope_true_matches_same_isotope_label() {
        // Two identically 13C-labeled methanol molecules: with
        // match_isotope=true, should still match all atoms.
        let labeled1 = parse("[13CH3]O").unwrap();
        let labeled2 = parse("[13CH3]O").unwrap();
        let config = McsConfig {
            match_isotope: true,
            ..McsConfig::default()
        };
        let result = find_mcs_with_config(&[&labeled1, &labeled2], &config);
        assert_eq!(
            result.atom_count(),
            labeled1.atom_count(),
            "match_isotope=true should still match molecules with the same isotope label"
        );
    }

    #[test]
    fn frontier_exclude_branch_finds_larger_mcs_than_include_only() {
        // Regression test for a real completeness bug: `grow` used to try
        // only `frontier[0]` at each node with no way to exclude it and try
        // a different frontier atom instead.
        //
        // mol0: O-C(-N)(-N) — hub C's neighbor list (parse/bond-insertion
        // order) is [O, N(branch), N]. mol1: N-C(-N) — hub' with 2 N leaves,
        // no O at all. True MCS is the full 3-atom C(N)(N) (mol1 in its
        // entirety): every seed that reaches "hub mapped, one N-leaf still
        // unmapped" puts unmatchable O ahead of the remaining N leaf in
        // frontier order. Without an exclude branch, growth died at size 2
        // (stuck on O) regardless of which atom was used as the seed --
        // fixed by branching include/exclude on `frontier[0]` instead of
        // only ever including it.
        let mol0 = parse("OC(N)N").unwrap();
        let mol1 = parse("NC(N)").unwrap();
        let result = find_mcs(&[&mol0, &mol1]);
        assert_eq!(
            result.atom_count(),
            3,
            "expected true MCS (C(N)(N), 3 atoms = all of mol1), got {}",
            result.atom_count()
        );
    }

    #[test]
    fn test_element_different_tie_is_spelling_dependent_known_gap() {
        // Two disconnected 3-atom fragments differing only in the
        // heteroatom (N vs O): both are size-3/bond-count-2 candidates of
        // equal quality, a genuine (non-automorphic, element-different) tie
        // -- unlike most ties elsewhere in this crate's tests, the choice
        // here IS externally observable (AtomicNum(7) vs AtomicNum(8) in
        // the returned QueryMolecule), because AtomCompare::Elements records
        // atomic number per atom.
        //
        // Spot-checked per user request after the #4 (iupac) tie-break
        // fix turned out to be a net-negative regression: confirmed via a
        // temporary revert that mcs_tiebreak_key's choice here is IDENTICAL
        // to the pre-#8-fix (first-found-in-DFS-wins) behavior -- this is
        // not a regression introduced by #8, it's a pre-existing gap that
        // #8's tie-break doesn't resolve (see mcs_tiebreak_key's doc
        // comment: no full spelling-invariance is claimed for it).
        let n_first = find_mcs(&[&parse("NCC.OCC").unwrap(), &parse("NCC.OCC").unwrap()]);
        let o_first = find_mcs(&[&parse("OCC.NCC").unwrap(), &parse("OCC.NCC").unwrap()]);

        assert_eq!(n_first.atom_count(), 3);
        assert_eq!(o_first.atom_count(), 3);
        assert_eq!(
            n_first.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(7))
        );
        assert_eq!(
            o_first.atoms[0].query,
            AtomQuery::Primitive(AtomPrimitive::AtomicNum(8))
        );
    }
}
