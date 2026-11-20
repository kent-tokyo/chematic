//! Maximum Common Substructure (MCS) search using the McGregor connected-growth algorithm.
//!
//! # Overview
//! Given a set of molecules, `find_mcs` returns a `QueryMolecule` that represents the
//! largest substructure common to all input molecules.  The algorithm is a branch-and-bound
//! search that grows the common subgraph one atom at a time, pruning subtrees whose upper
//! bound on reachable size is no better than the current best.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::{RingSet, find_sssr};

use crate::clock::Instant;
use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

#[path = "mcs_bounds.rs"]
mod bounds;
use bounds::{AtomUpperBound, BondUpperBound};

#[path = "mcs_frontier.rs"]
mod frontier;
use frontier::{extend_frontier_atom, mapped_rows_share_common_bond};

#[path = "mcs_stereo.rs"]
mod stereo;
use stereo::mapped_chirality_compatible;

#[cfg(test)]
#[path = "mcs_comparator_tests.rs"]
mod comparator_tests;

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
    /// Minimum number of atoms for a candidate when comparing two or more molecules.
    /// If none reaches `min_atoms`, return an empty `QueryMolecule`.
    pub min_atoms: usize,
    /// Optional time limit in milliseconds.  The search is aborted when the deadline is
    /// reached and the best result found so far is returned.
    pub timeout_ms: Option<u64>,
    /// Maximum number of search nodes. `None` disables this bound.
    pub max_search_nodes: Option<usize>,
    /// Maximum atoms accepted from any input molecule.
    pub max_atoms_per_molecule: usize,
    /// Maximum bonds accepted from any input molecule.
    pub max_bonds_per_molecule: usize,
    /// If `true`, ring atoms can only be matched to ring atoms, and non-ring atoms only
    /// to non-ring atoms.  Prevents chemically meaningless matches across ring boundaries.
    pub ring_matches_ring_only: bool,
    /// If `true`, a candidate is accepted only when every ring it touches in the
    /// first molecule is complete in the common subgraph.
    pub complete_rings_only: bool,
    /// Atom-level comparison strategy.  Defaults to [`AtomCompare::Elements`].
    pub atom_compare: AtomCompare,
    /// Bond-level comparison strategy.  Defaults to [`BondCompare::OrderOrAromatic`].
    pub bond_compare: BondCompare,
    /// If `true`, fully mapped tetrahedral centers must have the same stereochemistry
    /// after accounting for ligand order. If `false`, chirality is ignored (default).
    pub match_chiral_tag: bool,
    /// If `true`, formal charge must match between atoms. If `false` (default), charge is
    /// ignored -- e.g. a carboxylate `[O-]` matches a neutral `O`. Mirrors `match_chiral_tag`.
    pub match_charge: bool,
    /// If `true`, isotope label must match between atoms (including "no label" only matching
    /// "no label"). If `false` (default), isotope is ignored -- e.g. `[13CH4]` matches `C`.
    /// Mirrors `match_chiral_tag`.
    pub match_isotope: bool,
    /// If `true`, when multiple MCS of the same atom count exist, prefer the one with more bonds.
    /// Defaults to `true` for the atom-first objective. Select
    /// [`McsObjective::BondCount`] to make bonds the primary objective.
    pub maximize_bonds: bool,
}

impl Default for McsConfig {
    fn default() -> Self {
        Self {
            match_bonds: true,
            min_atoms: 1,
            timeout_ms: None,
            max_search_nodes: Some(1_000_000),
            max_atoms_per_molecule: 100_000,
            max_bonds_per_molecule: 200_000,
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

/// What the MCS search maximizes before resolving ties.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum McsObjective {
    /// Maximize atoms, using `McsConfig::maximize_bonds` for tied atom counts.
    AtomCount,
    /// Maximize bonds first, then prefer more atoms.
    BondCount,
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
    /// A configured input or search-node limit was reached; `result` is the
    /// best candidate found so far and is not proven to be maximal.
    ResourceLimited(QueryMolecule),
}

impl McsOutcome {
    /// Discard the exhaustive/timed-out distinction and return the result
    /// either way -- what [`find_mcs`]/[`find_mcs_with_config`] have always
    /// returned.
    pub fn into_query(self) -> QueryMolecule {
        match self {
            McsOutcome::Exhaustive(q)
            | McsOutcome::TimedOut(q)
            | McsOutcome::ResourceLimited(q) => q,
        }
    }

    /// `true` if the search was cut off by `timeout_ms` before finishing.
    pub fn was_timed_out(&self) -> bool {
        matches!(self, McsOutcome::TimedOut(_))
    }

    /// `true` if an input-size or search-node limit stopped the search.
    pub fn was_resource_limited(&self) -> bool {
        matches!(self, McsOutcome::ResourceLimited(_))
    }
}

/// Like [`find_mcs_with_config`], but reports whether `config.timeout_ms`
/// was reached before the search finished exhaustively -- see [`McsOutcome`].
pub fn find_mcs_with_config_checked(mols: &[&Molecule], config: &McsConfig) -> McsOutcome {
    find_mcs_with_objective_checked(mols, config, McsObjective::AtomCount)
}

/// Find an MCS using RDKit's `CompareAnyHeavyAtom` atom policy while retaining
/// every other option in `config`. Heavy atoms match each other, explicit
/// hydrogen matches explicit hydrogen, and wildcard matches wildcard.
/// The policy overrides `config.atom_compare`; [`find_mcs_with_config_checked`]
/// retains the original heavy-only behavior.
pub fn find_mcs_with_rdkit_any_heavy_checked(mols: &[&Molecule], config: &McsConfig) -> McsOutcome {
    find_mcs_with_rdkit_any_heavy_objective_checked(mols, config, McsObjective::AtomCount)
}

/// Maximum common substructure with an explicit primary objective.
pub fn find_mcs_with_objective_checked(
    mols: &[&Molecule],
    config: &McsConfig,
    objective: McsObjective,
) -> McsOutcome {
    find_mcs_impl(mols, config, false, objective)
}

/// Explicit objective under the RDKit any-heavy-atom matching policy.
pub fn find_mcs_with_rdkit_any_heavy_objective_checked(
    mols: &[&Molecule],
    config: &McsConfig,
    objective: McsObjective,
) -> McsOutcome {
    find_mcs_impl(mols, config, true, objective)
}

fn find_mcs_impl(
    mols: &[&Molecule],
    config: &McsConfig,
    rdkit_heavy: bool,
    objective: McsObjective,
) -> McsOutcome {
    if mols.iter().any(|mol| {
        mol.atom_count() > config.max_atoms_per_molecule
            || mol.bond_count() > config.max_bonds_per_molecule
    }) {
        return McsOutcome::ResourceLimited(QueryMolecule::new());
    }
    if mols.is_empty() {
        return McsOutcome::Exhaustive(QueryMolecule::new());
    }
    if mols.len() == 1 {
        return if config.max_search_nodes == Some(0) {
            McsOutcome::ResourceLimited(QueryMolecule::new())
        } else {
            McsOutcome::Exhaustive(molecule_to_query(mols[0], config))
        };
    }

    let n0 = mols[0].atom_count();
    let mut state = McsState::new(mols, config, rdkit_heavy, objective);

    'outer: for a0 in 0..n0 {
        if state.timed_out || state.resource_limited || state.proven_optimal {
            break;
        }
        if let Some(d) = state.deadline
            && Instant::now() >= d
        {
            state.timed_out = true;
            break;
        }

        // The edge-based partition used by bond-first search can reach every
        // connected mapping from its smallest-index vertex.
        if objective == McsObjective::BondCount && a0 > 0 {
            state.exclude_prior_seed(AtomIdx((a0 - 1) as u32));
        }

        // Gather compatible atoms from each other molecule for this seed atom.
        let seed_candidates = collect_seed_candidates(
            mols,
            AtomIdx(a0 as u32),
            config,
            &state.ring_sets,
            rdkit_heavy,
        );

        // Iterate the Cartesian product of candidate lists.
        for seed in CartesianProduct::new(&seed_candidates) {
            if state.timed_out || state.resource_limited || state.proven_optimal {
                break 'outer;
            }
            if let Some(d) = state.deadline
                && Instant::now() >= d
            {
                state.timed_out = true;
                break 'outer;
            }

            if !state.bounds_can_improve() {
                continue;
            }

            let mut mapping = PartialMapping::from_seed(mols, AtomIdx(a0 as u32), &seed);
            grow(&mut state, &mut mapping);
        }
    }

    let result = if state.best.size < config.min_atoms {
        QueryMolecule::new()
    } else {
        build_query(mols, &state.best, config)
    };
    if state.resource_limited {
        McsOutcome::ResourceLimited(result)
    } else if state.timed_out {
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
    rdkit_heavy: bool,
    objective: McsObjective,
    best: PartialMapping,
    deadline: Option<Instant>,
    timed_out: bool,
    resource_limited: bool,
    search_nodes: usize,
    /// Upper bounds on the size of any mapping under the selected comparators.
    absolute_atom_upper_bound: usize,
    absolute_bond_upper_bound: usize,
    /// Set once a valid mapping reaches both absolute bounds. At that point
    /// no later branch can improve either optimization objective, so symmetric
    /// seed permutations can be skipped without weakening completeness.
    proven_optimal: bool,
    ring_sets: Vec<RingSet>,
    atom_bound: AtomUpperBound,
    bond_bound: Option<BondUpperBound>,
}

impl<'a> McsState<'a> {
    fn new(
        mols: &'a [&'a Molecule],
        config: &'a McsConfig,
        rdkit_heavy: bool,
        objective: McsObjective,
    ) -> Self {
        let deadline = config
            .timeout_ms
            .map(|ms| Instant::now() + std::time::Duration::from_millis(ms));

        let ring_sets: Vec<RingSet> = if config.ring_matches_ring_only || config.complete_rings_only
        {
            mols.iter().map(|m| find_sssr(m)).collect()
        } else {
            Vec::new()
        };

        let atom_bound = AtomUpperBound::new(mols, &config.atom_compare, rdkit_heavy);
        let bond_bound = (objective == McsObjective::BondCount)
            .then(|| BondUpperBound::new(mols, config, rdkit_heavy));
        let (absolute_atom_upper_bound, absolute_bond_upper_bound) = match &bond_bound {
            Some(bonds) => (atom_bound.total().min(bonds.total() + 1), bonds.total()),
            None => (
                mols.iter().map(|mol| mol.atom_count()).min().unwrap_or(0),
                mols.iter().map(|mol| mol.bond_count()).min().unwrap_or(0),
            ),
        };

        Self {
            mols,
            config,
            rdkit_heavy,
            objective,
            best: PartialMapping::empty(mols.len()),
            deadline,
            timed_out: false,
            resource_limited: false,
            search_nodes: 0,
            absolute_atom_upper_bound,
            absolute_bond_upper_bound,
            proven_optimal: false,
            ring_sets,
            atom_bound,
            bond_bound,
        }
    }

    fn source_bond_allowed(&self, bond: BondIdx) -> bool {
        self.bond_bound
            .as_ref()
            .is_none_or(|bound| bound.is_allowed(bond))
    }

    fn exclude_prior_seed(&mut self, atom: AtomIdx) {
        self.atom_bound.set_excluded(atom, true);
        if let Some(bonds) = &mut self.bond_bound {
            for (_, bond) in self.mols[0].neighbors(atom) {
                if bonds.is_allowed(bond) {
                    bonds.set_excluded(bond, true);
                }
            }
        }
    }

    fn score(&self, mapping: &PartialMapping) -> (usize, usize) {
        match self.objective {
            McsObjective::AtomCount => (
                mapping.size,
                if self.config.maximize_bonds {
                    mapping.bond_count
                } else {
                    0
                },
            ),
            McsObjective::BondCount => (mapping.bond_count, mapping.size),
        }
    }

    fn bounds_can_improve(&self) -> bool {
        let mut atoms = self.atom_bound.total();
        if let Some(bond_bound) = &self.bond_bound {
            let bonds = bond_bound.total();
            // A connected result has at most one more vertex than edges.
            atoms = atoms.min(bonds + 1);
            if bonds != self.best.bond_count {
                return atoms >= self.config.min_atoms && bonds > self.best.bond_count;
            }
        }
        atoms >= self.config.min_atoms && atoms >= self.best.size
    }
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
fn mcs_tiebreak_key(mapping: &PartialMapping) -> impl Iterator<Item = u32> + '_ {
    mapping.query_to_mol.iter().map(|row| row[0].0)
}

/// Reject a fully mapped stereocenter with incompatible ligand parity, and
/// keep the best candidate that satisfies all configured constraints.
fn consider_mapping(state: &mut McsState<'_>, mapping: &PartialMapping) -> bool {
    // A raw @/@@ tag is relative to the SMILES neighbor order. Only a mapped
    // ligand permutation can establish whether two fully represented centers
    // have the same stereochemistry. A partial center remains eligible.
    if state.config.match_chiral_tag && !mapped_chirality_compatible(state.mols, mapping) {
        return false;
    }
    if mapping.size < state.config.min_atoms {
        return true;
    }
    let score = state.score(mapping);
    let current_best = state.score(&state.best);
    let is_better = score > current_best
        || (score == current_best && mcs_tiebreak_key(mapping).lt(mcs_tiebreak_key(&state.best)));
    // Pruning only the unconstrained winner after the search can discard a
    // smaller, valid candidate found earlier (for example a common chain next
    // to two incompatible rings). Keep the best *valid* mapping instead.
    if is_better
        && (!state.config.complete_rings_only
            || mapping_has_complete_rings(state.mols, mapping, &state.ring_sets, state.config))
    {
        state.best = mapping.clone();
        state.proven_optimal = state.best.size == state.absolute_atom_upper_bound
            && ((state.objective == McsObjective::AtomCount && !state.config.maximize_bonds)
                || state.best.bond_count == state.absolute_bond_upper_bound);
    }
    true
}

fn grow(state: &mut McsState<'_>, mapping: &mut PartialMapping) {
    let mut excluded_atoms = Vec::new();
    let mut excluded_bonds = Vec::new();
    while prepare_growth_node(state, mapping) {
        let Some(frontier) = next_frontier_atom(state, mapping) else {
            break;
        };
        extend_frontier_atom(state, mapping, frontier);
        if state.timed_out || state.resource_limited || state.proven_optimal {
            break;
        }
        match state.objective {
            McsObjective::AtomCount => {
                state.atom_bound.set_excluded(frontier, true);
                excluded_atoms.push(frontier);
            }
            McsObjective::BondCount => {
                exclude_frontier_bonds(state, mapping, frontier, &mut excluded_bonds);
            }
        }
    }
    for atom in excluded_atoms.into_iter().rev() {
        state.atom_bound.set_excluded(atom, false);
    }
    if let Some(bonds) = &mut state.bond_bound {
        for bond in excluded_bonds {
            bonds.set_excluded(bond, false);
        }
    }
}

fn prepare_growth_node(state: &mut McsState<'_>, mapping: &PartialMapping) -> bool {
    if state.proven_optimal {
        return false;
    }
    if let Some(limit) = state.config.max_search_nodes {
        if state.search_nodes >= limit {
            state.resource_limited = true;
            return false;
        }
        state.search_nodes += 1;
    }
    if !consider_mapping(state, mapping) {
        return false;
    }
    if let Some(deadline) = state.deadline
        && Instant::now() >= deadline
    {
        state.timed_out = true;
        return false;
    }
    !state.timed_out
        && !state.resource_limited
        && !state.proven_optimal
        && state.bounds_can_improve()
}

fn next_frontier_atom(state: &McsState<'_>, mapping: &PartialMapping) -> Option<AtomIdx> {
    let source = state.mols[0];
    mapping.query_to_mol.iter().find_map(|row| {
        source
            .neighbors(row[0])
            .find(|&(atom, bond)| {
                !mapping.is_mapped(0, atom)
                    && !state.atom_bound.is_excluded(atom)
                    && state.source_bond_allowed(bond)
            })
            .map(|(atom, _)| atom)
    })
}

fn exclude_frontier_bonds(
    state: &mut McsState<'_>,
    mapping: &PartialMapping,
    frontier: AtomIdx,
    excluded: &mut Vec<BondIdx>,
) {
    // The include branch covered every assignment in which any current
    // edge from the mapping to this vertex is common. Exclude those edges.
    // Later vertices can still connect the frontier by other edges.
    if let Some(bonds) = &mut state.bond_bound {
        for (neighbor, bond) in state.mols[0].neighbors(frontier) {
            if mapping.is_mapped(0, neighbor) && bonds.is_allowed(bond) {
                bonds.set_excluded(bond, true);
                excluded.push(bond);
            }
        }
    }
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
    rdkit_heavy: bool,
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
                    rdkit_heavy,
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
// Compatibility helpers
// ---------------------------------------------------------------------------

fn rdkit_heavy_atom_class(atom: &chematic_core::Atom) -> u8 {
    match crate::match_vf2::effective_atomic_number(atom) {
        0 => 0,
        1 => 1,
        _ => 2,
    }
}

/// Atom compatibility according to `compare` mode.
fn atoms_compatible(
    a: &chematic_core::Atom,
    b: &chematic_core::Atom,
    compare: &AtomCompare,
    rdkit_heavy: bool,
    match_chiral: bool,
    match_charge: bool,
    match_isotope: bool,
) -> bool {
    let class = bounds::atom_class(a, compare, rdkit_heavy);
    if class.is_none() || class != bounds::atom_class(b, compare, rdkit_heavy) {
        return false;
    }
    if match_chiral
        && !(a.chirality.is_tetrahedral() && b.chirality.is_tetrahedral())
        && a.chirality != b.chirality
    {
        return false;
    }
    if match_charge && a.charge != b.charge {
        return false;
    }
    if match_isotope && a.isotope.unwrap_or(0) != b.isotope.unwrap_or(0) {
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
// complete_rings_only candidate validation
// ---------------------------------------------------------------------------

/// A candidate touching a ring must contain all its atoms and common bonds.
/// Rings are checked on mol[0], which defines the query subgraph.
fn mapping_has_complete_rings(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    ring_sets: &[RingSet],
    config: &McsConfig,
) -> bool {
    ring_sets[0].rings().iter().all(|ring| {
        let mapped = ring.iter().filter(|&&ai| mapping.is_mapped(0, ai)).count();
        mapped == 0
            || (mapped == ring.len()
                && ring.iter().enumerate().all(|(index, &left)| {
                    let right = ring[(index + 1) % ring.len()];
                    let Some(q1) = mapping.query_idx_of(0, left) else {
                        return false;
                    };
                    let Some(q2) = mapping.query_idx_of(0, right) else {
                        return false;
                    };
                    mapped_rows_share_common_bond(mols, mapping, q1, q2, config)
                }))
    })
}

// ---------------------------------------------------------------------------
// Build output QueryMolecule
// ---------------------------------------------------------------------------

fn constrain_atom(
    mut query: AtomQuery,
    atom: &chematic_core::Atom,
    config: &McsConfig,
) -> AtomQuery {
    if config.match_charge {
        query = AtomQuery::And(
            Box::new(query),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Charge(atom.charge))),
        );
    }
    if config.match_isotope {
        query = AtomQuery::And(
            Box::new(query),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Isotope(
                atom.isotope.unwrap_or(0),
            ))),
        );
    }
    query
}

fn build_query(mols: &[&Molecule], mapping: &PartialMapping, config: &McsConfig) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();
    let mol0 = mols[0];

    // Add one query atom per mapped position. Aromaticity is deliberately NOT
    // encoded here as a per-atom constraint (matches RDKit's own `CompareElements`
    // representation, confirmed via live oracle): it's already fully captured by
    // the aromatic bond queries added below, and hard-coding mol0's own atom.aromatic
    // would make the query fail to re-match a different input molecule whose
    // corresponding atom disagrees on aromaticity (a real, confirmed regression when
    // AtomCompare::Elements permits cross-aromaticity matches).
    for row in &mapping.query_to_mol {
        // RDKit's FMCS query includes every element observed at this mapped
        // position when the atom comparator accepts unlike elements. Keeping
        // only mol0's element would produce an MCS that cannot match another
        // input even though the search mapped it successfully.
        let mut numbers = Vec::new();
        for (&mol, &atom_idx) in mols.iter().zip(row) {
            let number = crate::match_vf2::effective_atomic_number(mol.atom(atom_idx));
            if !numbers.contains(&number) {
                numbers.push(number);
            }
        }
        let query = numbers
            .into_iter()
            .map(|number| AtomQuery::Primitive(AtomPrimitive::AtomicNum(number)))
            .reduce(|left, right| AtomQuery::Or(Box::new(left), Box::new(right)))
            .expect("an MCS mapping includes its first input");
        qmol.add_atom(constrain_atom(query, mol0.atom(row[0]), config));
    }

    // Add only bonds present in every input.  Extra edges in any one input
    // are intentionally omitted: MCS is non-induced subgraph matching.
    for (q1, row1) in mapping.query_to_mol.iter().enumerate() {
        let a1 = row1[0];
        for (q2, row2) in mapping.query_to_mol.iter().enumerate() {
            if q2 <= q1 {
                continue;
            }
            let a2 = row2[0];
            if mol0.bond_between(a1, a2).is_some()
                && mapped_rows_share_common_bond(mols, mapping, q1, q2, config)
            {
                // CompareAny and match_bonds=false may map different bond
                // orders. Emit the observed alternatives, not mol0's order
                // or an unrestricted bond that also matches unseen orders.
                let mut alternatives = Vec::new();
                for (mol, (&left, &right)) in mols.iter().zip(row1.iter().zip(row2)) {
                    let (_, bond) = mol
                        .bond_between(left, right)
                        .expect("every MCS query bond is present in each input");
                    let query = bond_order_to_query(bond.order);
                    if !alternatives.contains(&query) {
                        alternatives.push(query);
                    }
                }
                let bq = alternatives
                    .into_iter()
                    .reduce(|left, right| BondQuery::Or(Box::new(left), Box::new(right)))
                    .expect("an MCS query bond includes its first input");
                qmol.add_bond(q1, q2, bq);
            }
        }
    }

    qmol
}

/// Convert a target-molecule `BondOrder` into the corresponding `BondQuery`
/// used in a `QueryMolecule`.
///
/// Up/Down/Dative stereo-like bonds collapse to Single. Query bond variants
/// become their closest SMARTS query equivalent. `build_query` combines the
/// orders observed across inputs when bond matching was relaxed.
fn bond_order_to_query(order: BondOrder) -> BondQuery {
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

fn molecule_to_query(mol: &Molecule, config: &McsConfig) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for (_, atom) in mol.atoms() {
        let an = crate::match_vf2::effective_atomic_number(atom);
        let query = AtomQuery::Primitive(AtomPrimitive::AtomicNum(an));
        qmol.add_atom(constrain_atom(query, atom, config));
    }

    for (_, bond) in mol.bonds() {
        let bq = bond_order_to_query(bond.order);
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
mod tests;
