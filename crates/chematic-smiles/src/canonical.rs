//! Canonical SMILES generation via the Morgan (extended connectivity) algorithm.
//!
//! A canonical SMILES is a unique string representation of a molecule:
//! two molecules that are identical (same graph, same atom properties) will
//! always produce the same canonical SMILES string.
//!
//! Algorithm:
//! 1. Assign initial invariants to each atom (atomic number, degree, charge, …).
//! 2. Iteratively update ranks using Morgan-style neighbor aggregation until
//!    the number of distinct ranks stabilises.
//! 3. Use the resulting ranks to impose a canonical DFS traversal order.
//!    Critically, both the ring-closure discovery DFS and the write DFS
//!    use the *same* canonical traversal order so the output is stable.
//! 4. Tie-breaking when two atoms have equal Morgan rank is resolved by
//!    atomic_number → isotope → charge → aromaticity → degree.
//!
//! Reference: Weininger, D. (1988) J. Chem. Inf. Comput. Sci. 28, 31-36.

use std::collections::{HashMap, HashSet};

use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL, implicit_hcount,
    remap_square_planar_tag, remap_tetrahedral_parity, valence_inferred_hcount,
};

use crate::writer::{
    bond_token_from, emit_bracket_hydrogens, square_planar_token, suppress_standalone_wedge,
};

/// Return the atom indices sorted into canonical (Morgan-rank) order.
///
/// The returned `Vec<usize>` lists atom positions (0-based) in the order they
/// would be encountered during a canonical DFS write.  Atoms with higher
/// Morgan rank appear earlier.  This is the same ordering `canonical_smiles`
/// uses internally: raw `morgan_ranks` ties are resolved via the same
/// individualize-refine + lexicographically-smallest-string selection, not
/// left as an input-order-dependent plateau.
///
/// Useful for normalizing atom-indexed property arrays to a canonical order.
pub fn canonical_atom_order(mol: &Molecule) -> Vec<usize> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }
    let (ranks, _) = winning_individualized_ranks(mol);
    let mut order: Vec<usize> = (0..n).collect();
    // Sort descending by rank (highest rank first, as in canonical DFS).
    order.sort_unstable_by(|&a, &b| ranks[b].cmp(&ranks[a]));
    order
}

/// Resolve `morgan_ranks` ties via individualize-refine and return the fully
/// discrete per-atom ranks of whichever branch produces the
/// lexicographically smallest canonical SMILES, *and* that winning string --
/// shared by `canonical_smiles` and `canonical_atom_order` so both use the
/// identical tie-break, instead of `canonical_atom_order` silently falling
/// back to raw (tie-break-free) `morgan_ranks`.
///
/// Returning the already-written winning string (not just its ranks) lets
/// [`canonical_smiles`] reuse it directly instead of calling
/// `CanonicalWriter::write_all` a second time on the same ranks -- every
/// branch considered here is written exactly once, tied or not, so a caller
/// that also wrote the winner again was paying for a fully redundant
/// traversal on every single call (perf; issue found bisecting the
/// `run_reactants`/`apply_retro` regression between chematic 0.4.25 and
/// 0.4.30 -- see docs/rfcs/reaction_transform_perf.md).
///
/// As of this PR, the actual search is `canonical_search::
/// winning_individualized_ranks_with_limits` (automorphism-orbit-pruned,
/// streaming, called here with `CanonicalizationLimits::unbounded()`),
/// which is provably equivalent to (a strict subset of the same candidate
/// (ranks, string) pairs considered by) the legacy exhaustive enumeration
/// below -- see `docs/rfcs/canonical_automorphism_pruning.md`. The legacy
/// exhaustive path (`legacy_winning_individualized_ranks`) is kept as a
/// last-resort fallback for the astronomically-rare case of an internal bug
/// in the new engine (`CanonicalizationError::InvalidInternalMapping`): it
/// degrades to "slow but correct" rather than panicking or returning wrong
/// output. `CanonicalizationError::SearchBudgetExceeded` cannot occur here
/// since `unbounded()` never checks either budget.
fn winning_individualized_ranks(mol: &Molecule) -> (Vec<u64>, String) {
    match crate::canonical_search::winning_individualized_ranks_with_limits(
        mol,
        &crate::canonical_search::CanonicalizationLimits::unbounded(),
    ) {
        Ok(result) => result,
        Err(_) => legacy_winning_individualized_ranks(mol),
    }
}

/// The original (pre-orbit-pruning) individualize-refine search: exhaustive,
/// capped at `MAX_INDIVIDUALIZE_BRANCHES`. Kept as (a) the fallback of last
/// resort for `winning_individualized_ranks` and (b) test-only exhaustive
/// oracle support (`canonical_smiles_exhaustive_oracle`, `usize::MAX`
/// budget). No longer the default hot path.
fn legacy_winning_individualized_ranks(mol: &Molecule) -> (Vec<u64>, String) {
    let plateaued = morgan_ranks(mol);
    let mut budget = MAX_INDIVIDUALIZE_BRANCHES;
    let branches = enumerate_discrete_ranks(mol, plateaued, &mut budget);
    branches
        .into_iter()
        .map(|ranks| {
            let s = CanonicalWriter::new(mol, &ranks).write_all();
            (ranks, s)
        })
        .min_by(|(_, a), (_, b)| a.cmp(b))
        .unwrap_or_default()
}

/// Benchmark-only: the pre-orbit-pruning exhaustive search (same code as
/// [`legacy_winning_individualized_ranks`], capped at
/// `MAX_INDIVIDUALIZE_BRANCHES` exactly like before this PR), exposed so
/// `examples/canonical_orbit_perf.rs` can measure a genuine before/after
/// comparison against the exact prior algorithm from within the same crate
/// version, instead of requiring a separate git-worktree checkout. Not part
/// of the stable public API (`#[doc(hidden)]`); never call this from
/// production code -- use [`canonical_smiles`].
#[doc(hidden)]
pub fn legacy_canonical_smiles_for_benchmark(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }
    legacy_winning_individualized_ranks(mol).1
}

/// Benchmark-only: the number of individualize-refine branches (leaves) the
/// pre-orbit-pruning exhaustive search would explore for `mol`, capped at
/// `MAX_INDIVIDUALIZE_BRANCHES` exactly like before this PR. Paired with
/// [`legacy_canonical_smiles_for_benchmark`] so `examples/
/// canonical_orbit_perf.rs` can report a precise leaf-count reduction
/// (old branch count vs the new engine's `leaves_written` instrumentation
/// counter), not just wall-clock timing. `#[doc(hidden)]`, benchmark-only.
#[doc(hidden)]
pub fn legacy_branch_count_for_benchmark(mol: &Molecule) -> usize {
    if mol.atom_count() == 0 {
        return 0;
    }
    let plateaued = morgan_ranks(mol);
    let mut budget = MAX_INDIVIDUALIZE_BRANCHES;
    enumerate_discrete_ranks(mol, plateaued, &mut budget).len()
}

/// Test-only, unbounded, orbit-pruning-free exhaustive individualize-refine
/// oracle (section 12): the legacy enumeration with an effectively
/// unlimited budget, used to cross-check the orbit-pruned engine's output.
/// Never called "the" canonicalization implementation -- this is a
/// deliberately slow ground truth for tests only.
#[cfg(test)]
pub(crate) fn canonical_smiles_exhaustive_oracle(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }
    let plateaued = morgan_ranks(mol);
    let mut budget = usize::MAX;
    let branches = enumerate_discrete_ranks(mol, plateaued, &mut budget);
    branches
        .into_iter()
        .map(|ranks| CanonicalWriter::new(mol, &ranks).write_all())
        .min()
        .unwrap_or_default()
}

/// Same ground truth as [`canonical_smiles_exhaustive_oracle`], but also
/// returns the winning rank vector -- found necessary during independent
/// Round-1 correctness review (PR #193): the existing oracle checks only
/// ever compared the winning canonical *string*, never the *rank vector*
/// `canonical_search::winning_individualized_ranks_with_limits` also
/// returns and the public `canonical_atom_order` API consumes. Two branches
/// within one automorphism orbit can legitimately share a minimal string
/// via different rank vectors, so string equality alone does not prove
/// rank-vector equality. Compare the orbit-pruned engine's rank vector
/// against *this* (the unbounded exhaustive oracle's), not against
/// `legacy_winning_individualized_ranks`'s -- the legacy engine's own
/// winning leaf can legitimately differ too (same reason, independently
/// flagged by Round 2), so oracle-vs-oracle is the correct ground truth.
#[cfg(test)]
pub(crate) fn canonical_smiles_exhaustive_oracle_with_ranks(mol: &Molecule) -> (Vec<u64>, String) {
    if mol.atom_count() == 0 {
        return (Vec::new(), String::new());
    }
    let plateaued = morgan_ranks(mol);
    let mut budget = usize::MAX;
    let branches = enumerate_discrete_ranks(mol, plateaued, &mut budget);
    branches
        .into_iter()
        .map(|ranks| {
            let s = CanonicalWriter::new(mol, &ranks).write_all();
            (ranks, s)
        })
        .min_by(|(_, a), (_, b)| a.cmp(b))
        .unwrap_or_default()
}

/// Return `true` if atoms `a` and `b` are topologically equivalent (symmetric).
///
/// Two atoms are considered equivalent when they have the same Morgan rank —
/// meaning no graph-based feature (element, charge, degree, neighbour
/// environment, …) can distinguish them.
///
/// # Example
/// All six carbons of benzene are equivalent; the two carbons of ethane are
/// equivalent; the two oxygens of acetic acid are **not** (different degree).
/// Assign a symmetry class number to every atom.
///
/// Atoms with the same class number are topologically equivalent (symmetric).
/// Class numbers are consecutive integers starting at 0, ordered by increasing
/// Morgan rank (lowest rank = class 0).
///
/// # Example
/// Benzene returns `[0,0,0,0,0,0]` (all 6 carbons equivalent).
/// Toluene returns `[0,1,1,1,1,1,2]` (methyl-C, ring-Cs, ipso-C).
pub fn equivalent_atom_classes(mol: &Molecule) -> Vec<usize> {
    let ranks = morgan_ranks(mol);
    // Sort unique rank values to assign stable class numbers.
    let mut unique: Vec<u64> = ranks.clone();
    unique.sort_unstable();
    unique.dedup();
    ranks
        .iter()
        .map(|r| unique.partition_point(|&u| u < *r))
        .collect()
}

pub fn are_atoms_equivalent(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let ranks = morgan_ranks(mol);
    let ia = a.0 as usize;
    let ib = b.0 as usize;
    if ia >= ranks.len() || ib >= ranks.len() {
        return false;
    }
    ranks[ia] == ranks[ib]
}

/// Return the canonical SMILES for a molecule.
///
/// For molecules with no atoms, returns an empty string.
/// Disconnected fragments (multiple components) are joined with `.`.
///
/// Atom ordering is fully discretized before writing (individualize-refine,
/// see `enumerate_discrete_ranks`): when the plain Morgan refinement in
/// [`morgan_ranks`] plateaus with genuine (non-automorphism) ties still
/// present, every possible resolution is tried and the lexicographically
/// smallest resulting string is returned. This makes the output invariant to
/// which input atom ordering/spelling the molecule was parsed from, not just
/// idempotent under repeated self-canonicalization.
pub fn canonical_smiles(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }

    let (_, winning_string) = winning_individualized_ranks(mol);
    winning_string
}

/// Compute Morgan (extended connectivity) ranks for all atoms.
///
/// Returns a vector of normalised ordinal ranks (0-based, gap-free)
/// indexed by atom position (same order as `mol.atoms()`). This is pure
/// neighbor-hash refinement to a fixpoint -- it does NOT individualize
/// remaining ties, so atoms in the same non-trivial automorphism orbit (or
/// in a refinement cell that merely *contains* an orbit) keep equal ranks.
/// That is the correct, useful notion of "rank" for topological-symmetry
/// queries (see [`equivalent_atom_classes`], [`are_atoms_equivalent`]).
///
/// [`canonical_smiles`] does NOT use this directly for atom ordering when
/// ties remain -- see `enumerate_discrete_ranks` for the individualize-refine
/// step that resolves ties before writing.
pub fn morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();
    let initial: Vec<u64> = (0..n)
        .map(|i| initial_invariant(mol, AtomIdx(i as u32)))
        .collect();
    refine_ranks(mol, initial)
}

/// Refine `ranks` (any starting coloring, not necessarily the initial
/// invariant) via neighbor-hash iteration until the number of distinct
/// classes stops increasing. Used both for the plain (tie-preserving) ranks
/// in [`morgan_ranks`] and, with a perturbed starting coloring, as the
/// "refine" half of individualize-refine in `enumerate_discrete_ranks`.
pub(crate) fn refine_ranks(mol: &Molecule, mut ranks: Vec<u64>) -> Vec<u64> {
    let n = ranks.len();
    let max_iter = n + 2;
    for _ in 0..max_iter {
        let old_distinct = count_distinct(&ranks);

        let new_ranks: Vec<u64> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                // Include bond order in the neighbor contribution so that atoms
                // bonded via different bond types (e.g. O= vs O-H in acetic acid)
                // receive distinct Morgan ranks even when neighbor atom ranks are
                // otherwise identical.
                let mut neighbor_contributions: Vec<u64> = mol
                    .neighbors(idx)
                    .map(|(nb, bidx)| {
                        let bond_val = bond_order_value(mol.bond(bidx).order);
                        fnv_hash_sequence(ranks[nb.0 as usize], &[bond_val])
                    })
                    .collect();
                neighbor_contributions.sort_unstable();
                fnv_hash_sequence(ranks[i], &neighbor_contributions)
            })
            .collect();

        let new_distinct = count_distinct(&new_ranks);
        ranks = new_ranks;

        if new_distinct <= old_distinct {
            break;
        }
    }

    normalize_ranks(&ranks)
}

/// Safety cap on the number of discrete rank assignments
/// `enumerate_discrete_ranks` will explore. It exists to guarantee
/// termination on pathologically symmetric inputs (fullerene fragments,
/// deep dendrimers) where the principled fix is automorphism-aware branch
/// pruning (nauty-style), not attempted here. Once exhausted, the
/// remaining ties in that branch fall back to `canonical_cmp`'s finite
/// tie-break chain (deterministic, but not guaranteed order-independent).
///
/// This IS hit in practice, corrected 2026-07-12 after a claim of "never
/// hit" here went unverified: measured on 5,000 real ChEMBL-derived
/// molecules, 3 (0.06%) exceeded this cap, needing up to 168,219 branches
/// (16.8x). All three are real drug-synthesis intermediates with multiple
/// Boc/pivaloyl tert-butyl protecting groups, each an independent 3-way
/// symmetric orbit that multiplies combinatorially across the molecule.
/// For all three, truncation at 10,000 was confirmed (against an
/// unbounded run, and separately against 32 independent re-spellings) to
/// still find the correct lexicographically-smallest winner -- the
/// exhausted cells in these cases are true automorphism orbits (every
/// individualization within the cell writes the same string, so the
/// blowup is redundant duplicates, not competing candidates). This is not
/// a guarantee for the general case: the failure mode this cap is meant
/// to bound is a cell that merely *contains* an orbit (genuinely different
/// candidates truncated away), which no real molecule in this corpus
/// happened to exercise. Raising the constant is not a principled fix (the
/// observed distribution has a cliff -- p99.9 is ~4,922, but the 3
/// offenders need 74k-168k -- so no fixed multiple closes the gap); the
/// real fix is the orbit-aware pruning mentioned above.
const MAX_INDIVIDUALIZE_BRANCHES: usize = 10_000;

/// Group atom indices by their current (gap-free ordinal) rank value.
/// `by_rank[r]` lists, in ascending atom-index order, every atom whose rank
/// equals `r`. Shared by the legacy exhaustive `enumerate_discrete_ranks`
/// and the orbit-pruned `canonical_search::search_canonical`, so both agree
/// on cell membership and (within a cell) traversal order.
pub(crate) fn group_by_rank(ranks: &[u64]) -> Vec<Vec<usize>> {
    let mut by_rank: Vec<Vec<usize>> = Vec::new();
    for (i, &r) in ranks.iter().enumerate() {
        let r = r as usize;
        if by_rank.len() <= r {
            by_rank.resize(r + 1, Vec::new());
        }
        by_rank[r].push(i);
    }
    by_rank
}

/// Individualize atom `atom_idx` within its current rank class: insert a new
/// rank strictly between its class and the next-higher class, so a
/// subsequent refinement pass can propagate the distinction through the rest
/// of the graph. `ranks` must be gap-free ordinals (as produced by
/// `refine_ranks`/`normalize_ranks`).
pub(crate) fn individualize(ranks: &[u64], atom_idx: usize) -> Vec<u64> {
    let v = ranks[atom_idx];
    ranks
        .iter()
        .enumerate()
        .map(|(i, &r)| {
            if i == atom_idx {
                v + 1
            } else if r > v {
                r + 1
            } else {
                r
            }
        })
        .collect()
}

/// Enumerate every discrete (all-singleton) rank assignment reachable from
/// `ranks` (already refined to a fixpoint) via individualize-refine
/// branching.
///
/// Refinement cells are always unions of automorphism orbits (a standard
/// 1-WL / equitable-partition fact: refinement can never split an orbit).
/// So: if a cell IS an orbit, every choice of which atom to individualize
/// yields an automorphic result -- the resulting SMILES strings are
/// identical, so exploring all of them is correct but redundant. If a cell
/// properly CONTAINS an orbit, no order-independent rule can select a single
/// representative (if one existed, refinement would already have used it as
/// an invariant and the cell would not be tied) -- the only order-independent
/// resolution is to try every atom in the cell and let the caller take the
/// lexicographically smallest resulting string. Cell SELECTION (which
/// non-singleton class to branch on next) is itself order-independent:
/// always the lowest-ranked non-singleton cell.
fn enumerate_discrete_ranks(mol: &Molecule, ranks: Vec<u64>, budget: &mut usize) -> Vec<Vec<u64>> {
    let by_rank = group_by_rank(&ranks);

    // `by_rank` is indexed by ordinal rank value, so the first multi-member
    // entry found is the lowest-ranked non-singleton cell.
    let Some(members) = by_rank.iter().find(|m| m.len() > 1) else {
        return vec![ranks];
    };

    let mut results = Vec::new();
    for &atom_idx in members {
        if *budget == 0 {
            break;
        }
        *budget -= 1;
        let individualized = individualize(&ranks, atom_idx);
        let re_refined = refine_ranks(mol, individualized);
        results.extend(enumerate_discrete_ranks(mol, re_refined, budget));
    }
    results
}

/// Initial per-atom invariant packed into a u64.
fn initial_invariant(mol: &Molecule, idx: AtomIdx) -> u64 {
    let atom = mol.atom(idx);

    if atom.wildcard {
        return 0;
    }

    let an = atom.element.atomic_number() as u64;
    let degree = mol.degree(idx) as u64;
    let charge = (atom.charge as i64 + 128) as u64;
    let iso = atom.isotope.unwrap_or(0) as u64;
    let arom = atom.aromatic as u64;
    // The *effective* H count (explicit bracket value if present, else the
    // valence-inferred implicit count -- the same unification
    // `emit_bracket_hydrogens` already applies when writing) must seed the
    // invariant, not raw `atom.hydrogen_count`. An explicit bracket atom
    // whose stored H count merely repeats what valence inference would have
    // produced anyway (e.g. `[Cl]` vs organic-subset `Cl`, both H=0) is the
    // same chemical species and must get the same invariant; using
    // `hydrogen_count.is_some()` as part of the seed instead let bracket
    // "spelling" alone change Morgan ranks and therefore `canonical_smiles`'s
    // output for two representations of one molecule (see issue #205).
    let h_flag = implicit_hcount(mol, idx) as u64;

    (an << 56) | (degree << 48) | (charge << 40) | (iso << 24) | (h_flag << 16) | (arom << 8)
}

/// Map a BondOrder to a stable integer for use in Morgan rank hashing.
fn bond_order_value(order: BondOrder) -> u64 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Quadruple => 5,
        _ => 0,
    }
}

fn fnv_hash_sequence(base: u64, values: &[u64]) -> u64 {
    const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    let mut h = FNV_OFFSET ^ base.wrapping_mul(FNV_PRIME);
    for &v in values {
        h ^= v;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

fn count_distinct(ranks: &[u64]) -> usize {
    let mut seen: Vec<u64> = ranks.to_vec();
    seen.sort_unstable();
    seen.dedup();
    seen.len()
}

fn normalize_ranks(ranks: &[u64]) -> Vec<u64> {
    let mut sorted: Vec<(u64, usize)> = ranks
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    sorted.sort_unstable_by_key(|&(v, _)| v);

    let mut result = vec![0u64; ranks.len()];
    let mut current_rank: u64 = 0;
    let mut prev_val = sorted[0].0;

    for (val, idx) in sorted {
        if val != prev_val {
            current_rank += 1;
            prev_val = val;
        }
        result[idx] = current_rank;
    }

    result
}

pub(crate) struct CanonicalWriter<'a> {
    mol: &'a Molecule,
    ranks: &'a [u64],
    written: Vec<bool>,
    ring_bonds: HashSet<BondIdx>,
    /// (ring_num, is_open_side, ring_partner_atom, physical_bond). Deliberately
    /// does NOT store a precomputed `BondOrder`: which endpoint of `bidx` is
    /// "open" vs "close" is fixed at discovery time (`dfs_mark`), but the
    /// actual `/`/`\` token can only be correctly computed once written --
    /// see `normalize_ez`'s doc comment for why baking in a write-direction
    /// re-orientation this early corrupts E/Z groups spanning bonds visited
    /// in different directions (issue #390).
    atom_ring_nums: HashMap<AtomIdx, Vec<(u32, bool, AtomIdx, BondIdx)>>,
    next_ring: u32,
    out: String,
    /// Union-find groups of directional (`/`/`\`) bonds that jointly encode
    /// one connected E/Z system — flipping every member preserves geometry,
    /// flipping a subset does not. Keyed/rooted by `BondIdx`.
    ez_group: HashMap<BondIdx, BondIdx>,
    /// Groups whose first-encountered bond (in write order) came out `Down`;
    /// every remaining bond in the group is flipped so the first directional
    /// bond of each system is always `/`, regardless of input spelling.
    ez_flip: HashMap<BondIdx, bool>,
    /// Resolved "which bond carries the `/`/`\` marker" override, computed
    /// once up front by [`Self::resolve_ez_markers`] and consulted by every
    /// site that would otherwise read the bond's raw parse-time direction
    /// (literal `Up`/`Down` order, or the aromatic-bond-direction stash).
    /// Bonds not present here fall back to that raw reading unchanged — this
    /// map only ever contains entries for the substituent bonds of a
    /// tri-/tetra-substituted stereo alkene end, where the *choice* of which
    /// substituent carries the mark is otherwise input-order-dependent (see
    /// `resolve_ez_markers` for why). A demoted (no-longer-marked) bond is
    /// stored here too, pinned to its own plain (non-directional) order, so
    /// a stray parse-time marker on it can't leak through.
    ez_marker: HashMap<BondIdx, BondOrder>,
    /// Test-only instrumentation: every alkene end for which
    /// `resolve_ez_marker_for_end` hit the shared-candidate-bond abstain
    /// guard specifically (as opposed to "not ambiguous" or "no direction
    /// info at all"). Lets regression tests assert the guard actually
    /// fired for a given fixture by reading production's own record of
    /// which branch it took, rather than re-deriving the topology
    /// condition in the test and hoping it matches. `cfg(test)`-gated —
    /// zero cost/size impact on release builds.
    #[cfg(test)]
    ez_shared_bond_abstains: Vec<AtomIdx>,
}

impl<'a> CanonicalWriter<'a> {
    pub(crate) fn new(mol: &'a Molecule, ranks: &'a [u64]) -> Self {
        let n = mol.atom_count();
        Self {
            mol,
            ranks,
            written: vec![false; n],
            ring_bonds: HashSet::new(),
            atom_ring_nums: HashMap::new(),
            next_ring: 1,
            out: String::new(),
            ez_group: HashMap::new(),
            ez_flip: HashMap::new(),
            ez_marker: HashMap::new(),
            #[cfg(test)]
            ez_shared_bond_abstains: Vec::new(),
        }
    }

    /// The raw parse-time direction of `bidx`, if any: a literal `Up`/`Down`
    /// bond order, or (when both endpoints are aromatic) the stashed
    /// direction of an adjacent exocyclic double bond. This is the
    /// *unresolved* reading — [`Self::effective_order`] is what emission
    /// code should use instead, since it additionally applies
    /// `resolve_ez_markers`'s carrier choice on top of this.
    ///
    /// Delegates to `crate::writer::raw_bond_direction` -- the exact same
    /// rule the plain (non-canonical) writer uses, kept as one shared
    /// implementation on purpose so the two writers can never silently
    /// disagree on what a bond's stored/stashed direction means.
    fn raw_input_direction(&self, bidx: BondIdx) -> Option<BondOrder> {
        crate::writer::raw_bond_direction(self.mol, bidx)
    }

    /// Strip directionality from a bond order, leaving its "plain" chemical
    /// order untouched: `Up`/`Down` becomes `Single` (the literal-marker
    /// case), anything else (notably `Aromatic`, the stash case) is
    /// returned as-is. Used to demote a substituent bond that must no
    /// longer carry the `/`/`\` marker in the output.
    fn plain_order(order: BondOrder) -> BondOrder {
        if matches!(order, BondOrder::Up | BondOrder::Down) {
            BondOrder::Single
        } else {
            order
        }
    }

    /// The effective direction to use when emitting `bidx`: `resolve_ez_markers`'s
    /// resolved choice if this bond participates in a tri-/tetra-substituted
    /// stereo alkene end, otherwise the unresolved raw parse-time direction
    /// (literal order or aromatic-bond-direction stash), otherwise the
    /// bond's own real chemical order. Every write-time site that needs "the
    /// order to show for this bond" (ring-closure emission and tree-edge
    /// emission alike) must go through this, not read `bond_direction`/
    /// `order` directly, or the two sites can disagree on which bond a
    /// moved E/Z marker landed on.
    fn effective_order(&self, bidx: BondIdx) -> BondOrder {
        if let Some(&resolved) = self.ez_marker.get(&bidx) {
            return resolved;
        }
        self.raw_input_direction(bidx)
            .unwrap_or(self.mol.bond(bidx).order)
    }

    /// Return `true` if `dir` encodes "up" (`/`, read atom1→atom2) as seen
    /// from `alkene_end`'s side of `bond` — the same atom1/atom2-relative
    /// convention `chematic_chem::cip::substituent_is_up` uses to read a
    /// marker back out, so a value written here round-trips identically.
    fn direction_is_up(dir: BondOrder, bond_atom1: AtomIdx, alkene_end: AtomIdx) -> bool {
        match dir {
            BondOrder::Up => bond_atom1 == alkene_end,
            BondOrder::Down => bond_atom1 != alkene_end,
            _ => false, // never called with a non-directional `dir`
        }
    }

    /// The inverse of `direction_is_up`: which `Up`/`Down` order to store on
    /// `bond` so that, read from `alkene_end`'s side, it encodes `want_up`.
    fn direction_for_up(bond_atom1: AtomIdx, alkene_end: AtomIdx, want_up: bool) -> BondOrder {
        if (bond_atom1 == alkene_end) == want_up {
            BondOrder::Up
        } else {
            BondOrder::Down
        }
    }

    /// Maximum coupled-component size [`Self::resolve_component_jointly`]
    /// will attempt to jointly resolve via bounded enumeration
    /// (`2^size` candidate assignments). Measured, not guessed, and fully
    /// reproducible from a corpus committed to this repo (no external
    /// download required):
    ///
    /// ```text
    /// cargo run -p chematic-smiles --release --example ez_shared_carrier_component_audit -- \
    ///     scripts/descriptor_census_corpus.smi
    /// ```
    ///
    /// Expected output (re-run to confirm, not assumed): 18/18 pinned
    /// fixtures each produce exactly one size-2 component; the 5,000-line
    /// committed corpus produces 31 coupling components total, still every
    /// one of them size exactly 2, 0 cycles. (An independent run against
    /// the larger, non-committed ChEMBL corpus used for this crate's own
    /// full-corpus residual measurements — see
    /// `docs/rfcs/canonical_smiles_residual_rfc.md` — found the same: every
    /// component observed has been size exactly 2.) Every node has at most
    /// 2 candidate substituent bonds, so every component is a simple path
    /// or cycle, never a general graph. This cap gives an 8x margin over
    /// the measured maximum. Worst case is bounded at `2^16` = 65,536
    /// candidate-assignment evaluations, each O(size) with small constants
    /// (no independent wall-clock measurement is claimed here — every
    /// corpus scanned so far never exceeds 4 evaluations, i.e. component
    /// size 2, so the `2^16` bound is a safety margin, not an observed
    /// cost). A component that exceeds the cap abstains as a whole --
    /// never partially applied.
    const MAX_JOINT_COMPONENT_SIZE: usize = 16;

    /// Resolve, for every tri-/tetra-substituted stereo alkene end, which ONE
    /// of its ≥2 non-double-bond substituent bonds should carry the `/`/`\`
    /// marker in the output — deterministically, from the molecule's own
    /// (already fully individualized) canonical ranks, rather than just
    /// inheriting whichever bond the original parse happened to mark.
    ///
    /// **The bug this fixes** (canonical-residual RFC, Root cause 1): SMILES
    /// only requires ONE of an alkene carbon's ≥2 substituent bonds to carry
    /// an explicit marker (the geometry of the other is implied — at a
    /// trigonal alkene carbon with two substituents, they are always on
    /// opposite sides). *Which* substituent gets marked is a free choice at
    /// write time, but the writer previously just checked "does this
    /// specific bond already carry a direction", inherited straight from
    /// parse time — so two RDKit-valid respellings of the identical molecule
    /// that mark different substituents produced two different (but
    /// chemically identical) canonical outputs. Resolving the carrier here,
    /// once, from topology + rank alone, makes the choice depend only on the
    /// molecule itself, not on which substituent the input happened to mark.
    ///
    /// Two independently stereogenic double bonds can additionally share one
    /// physical candidate substituent bond between them (issue #149 —
    /// typically a ring-closure bond connecting two adjacent stereocenters
    /// each bearing their own exocyclic stereo double bond). Resolving one
    /// end's carrier in isolation from such a coupled partner risks moving
    /// or demoting a mark the *other* system's own resolution simultaneously
    /// relies on for the same physical bond. [`Self::coupling_components`]
    /// partitions every ambiguous end into maximal connected clusters joined
    /// by exactly this kind of shared-bond edge; [`Self::resolve_component_
    /// jointly`] solves each cluster as one unit (a lone, uncoupled end is
    /// simply a cluster of size 1, resolved by the same code path at
    /// negligible extra cost — see that method's own doc comment).
    ///
    /// Runs before `build_ez_groups`/ring discovery/DFS write, so it depends
    /// only on molecule topology, never on write order.
    fn resolve_ez_markers(&mut self) {
        if self.ranks.is_empty() {
            return;
        }
        let stereo_alkene_ends = Self::compute_stereo_alkene_ends(self.mol);
        for component in Self::coupling_components(self.mol, &stereo_alkene_ends) {
            self.resolve_component_jointly(&component);
        }
    }

    fn end_has_substituent(mol: &Molecule, end: AtomIdx) -> bool {
        mol.neighbors(end)
            .any(|(_, b)| mol.bond(b).order != BondOrder::Double)
    }

    /// Non-double-bond neighbors of `end` — an alkene carbon's up-to-two
    /// sigma substituents. Filtered by bond order, not by comparing against
    /// a specific double-bond `BondIdx`, so an allene/cumulene terminus
    /// correctly excludes *every* double bond at `end`, matching
    /// `chematic_chem::cip::assign_ez`'s own substituent-collection
    /// convention.
    fn substituents(mol: &Molecule, end: AtomIdx) -> Vec<(AtomIdx, BondIdx)> {
        mol.neighbors(end)
            .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
            .collect()
    }

    /// Whether the double bond between `a` and `b` is endocyclic in a ring
    /// smaller than 8 atoms — i.e. some SSSR ring contains both atoms and
    /// has fewer than 8 members. Such a bond's real-world geometry is fixed
    /// by the ring itself, not a free stereochemical choice (issue #149's
    /// ring-constrained residual, `docs/rfcs/ez_ring_constrained_residual_
    /// audit.md`): RDKit independently agrees the bond is not a real
    /// stereocenter (`STEREONONE`/absent from `FindPotentialStereo`) with
    /// 0/1,387 row-level disagreements on a 5,000-molecule corpus at exactly
    /// this threshold. Deliberately checks ring membership of the BOND (both
    /// endpoints in the same ring), not merely of the atom — a ring carbon
    /// with a genuinely free *exocyclic* double bond (its own C=X pointing
    /// outward) must not be excluded, and the audit found the naive
    /// atom-membership version over-excludes by ~9 percentage points of a
    /// real corpus for exactly this reason.
    fn double_bond_endocyclic_in_small_ring(
        rings: &[Vec<AtomIdx>],
        a: AtomIdx,
        b: AtomIdx,
    ) -> bool {
        rings
            .iter()
            .any(|r| r.len() < 8 && r.contains(&a) && r.contains(&b))
    }

    /// Every atom that is a stereogenic double bond's terminus with exactly
    /// 2 candidate substituents — the ambiguity precondition a marker-carrier
    /// choice exists for at all (0 or 1 substituent has nothing to choose
    /// between; a real double bond's end never has more than 2). A double
    /// bond only counts when BOTH ends have at least one substituent
    /// (matching `chematic_chem::cip::assign_ez`'s own
    /// `subs_a1.is_empty() || subs_a2.is_empty()` guard — e.g. a
    /// ketone/aldehyde `C=O` never does, the O side has none), so the carbon
    /// side of such a bond is never treated as ambiguous purely because it
    /// happens to have 2 ring-bond substituents of its own.
    ///
    /// Also excludes both ends of a double bond that is itself endocyclic in
    /// a ring smaller than 8 atoms (see
    /// [`Self::double_bond_endocyclic_in_small_ring`]) — such a bond has no
    /// real stereochemical freedom, so treating it as ambiguous only
    /// destabilizes marker-carrier selection for a genuinely stereogenic
    /// bond it happens to share a candidate bond with (issue #149).
    ///
    /// This runs on every `canonical_smiles()` call via `resolve_ez_markers`,
    /// so `chematic_perception::find_sssr` (Horton: O(V·E) candidates + GF(2)
    /// elimination) is only ever computed when there is at least one
    /// double-bond candidate that could need the ring check — most
    /// molecules (no double bonds, or only non-candidate ones like a
    /// ketone's `C=O`) skip it entirely.
    fn compute_stereo_alkene_ends(mol: &Molecule) -> HashSet<AtomIdx> {
        let candidates: Vec<_> = mol
            .bonds()
            .filter(|(_, bond)| bond.order == BondOrder::Double)
            .filter(|(_, bond)| {
                Self::end_has_substituent(mol, bond.atom1)
                    && Self::end_has_substituent(mol, bond.atom2)
            })
            .collect();
        if candidates.is_empty() {
            return HashSet::new();
        }

        let rings = chematic_perception::find_sssr(mol);
        let rings = rings.rings();
        let mut ends = HashSet::new();
        for (_, bond) in candidates {
            if Self::double_bond_endocyclic_in_small_ring(rings, bond.atom1, bond.atom2) {
                continue;
            }
            for end in [bond.atom1, bond.atom2] {
                if Self::substituents(mol, end).len() == 2 {
                    ends.insert(end);
                }
            }
        }
        ends
    }

    /// Partition every stereo-alkene-end atom into maximal connected
    /// components of the shared-candidate-bond coupling graph: nodes are the
    /// ends themselves, edges connect two ends directly bonded via one of
    /// their own (≤2) candidate substituent bonds — the only way one
    /// physical bond can be a marker-carrier candidate for two different
    /// double bonds' own resolution at once (a candidate bond is by
    /// definition incident to its owning end, so it can only ever be
    /// "shared" with the atom on its *other* side). Since every node has ≤2
    /// candidate bonds, every component is a simple path or cycle, never a
    /// general graph (measured directly — see [`Self::MAX_JOINT_COMPONENT_
    /// SIZE`]'s doc comment). A molecule with no coupling at all yields one
    /// singleton component per ambiguous end (the overwhelming common
    /// case). Component *membership* is a pure function of molecule
    /// topology; the order components (and atoms within them) are returned
    /// in has no effect on the result — [`Self::resolve_component_jointly`]
    /// re-sorts by canonical rank internally before making any decision.
    ///
    /// Traversal order (which node BFS starts from, and the order its
    /// neighbors are pushed) is made deterministic by sorting on raw
    /// `AtomIdx` — not because that order is canonical (it isn't), but
    /// because `ends` is a `HashSet` and iterating it directly would leak
    /// process-random hash-seed order into `component`'s element order.
    /// [`Self::resolve_component_jointly`]'s own rank-based sort washes
    /// that out for every rank-*distinct* member regardless, but ties in
    /// `self.ranks` are only possible for genuinely automorphic ends, and
    /// leaving *that* residual order-dependence to an unseeded HashSet
    /// would be process-random input to an otherwise-deterministic
    /// decision — so it is removed here at the source instead.
    fn coupling_components(mol: &Molecule, ends: &HashSet<AtomIdx>) -> Vec<Vec<AtomIdx>> {
        let mut starts: Vec<AtomIdx> = ends.iter().copied().collect();
        starts.sort_by_key(|a| a.0);
        let mut visited: HashSet<AtomIdx> = HashSet::new();
        let mut components = Vec::new();
        for start in starts {
            if visited.contains(&start) {
                continue;
            }
            let mut queue = vec![start];
            let mut component = Vec::new();
            visited.insert(start);
            while let Some(cur) = queue.pop() {
                component.push(cur);
                let mut nbs: Vec<AtomIdx> = Self::substituents(mol, cur)
                    .into_iter()
                    .filter(|(nb, _)| ends.contains(nb))
                    .map(|(nb, _)| nb)
                    .collect();
                nbs.sort_by_key(|a| a.0);
                for nb in nbs {
                    if visited.insert(nb) {
                        queue.push(nb);
                    }
                }
            }
            components.push(component);
        }
        components
    }

    /// Jointly resolve one coupling component (see
    /// [`Self::coupling_components`]) — a cluster of ambiguous alkene ends,
    /// possibly of size 1 (no coupling at all, the common case: this method
    /// is the *only* code path that assigns marker carriers, unified rather
    /// than kept as a separate fast path, since a singleton component can
    /// never conflict with anything and this method's enumeration collapses
    /// to exactly the old single-end algorithm's own choice at negligible
    /// extra cost — 2 evaluations instead of a direct pick).
    ///
    /// Enumerates every combination of "which of its 2 candidate bonds does
    /// end `i` use as its own marker carrier" (bounded by
    /// [`Self::MAX_JOINT_COMPONENT_SIZE`]), keeps only combinations where
    /// every physical bond in the component receives a single, mutually
    /// consistent literal value from every end that votes on it (see
    /// [`Self::evaluate_choice`]), and — among those — picks the one
    /// deviating from each end's own canonical-rank-preferred choice the
    /// *fewest* times (primary key: total deviation count, ascending), with
    /// a **molecule-intrinsic** tie-break among equal-count combinations:
    /// ends sorted by ascending `self.ranks` (never `AtomIdx`/`BondIdx`/
    /// traversal order), comparing each combination's per-end "did I
    /// deviate" bit vector in that order, lexicographically smallest wins
    /// (deviations are pushed onto the higher-ranked end first).
    ///
    /// That two-part key is a bijection of `choice` (distinct combinations
    /// always produce distinct keys), so it alone can never report a tie —
    /// a real one still exists when two *different* minimal-count
    /// combinations differ only by which of two `self.ranks`-equal (hence
    /// provably automorphic) ends carries a given deviation: no signal this
    /// molecule offers distinguishes them, so [`Self::coupling_components`]'s
    /// deterministic-but-arbitrary traversal order must not be allowed to
    /// pick a "winner" by accident. Detected explicitly below rather than
    /// left to fall out of the key comparison. If no consistent combination
    /// exists at all, or a genuine tie is detected, the whole component is
    /// left exactly as the input spelled it: always safe, never guessed.
    fn resolve_component_jointly(&mut self, component: &[AtomIdx]) {
        let mut ordered: Vec<AtomIdx> = component.to_vec();
        ordered.sort_by_key(|&a| self.ranks[a.0 as usize]);
        let k = ordered.len();

        if k > Self::MAX_JOINT_COMPONENT_SIZE {
            #[cfg(test)]
            self.ez_shared_bond_abstains.extend_from_slice(&ordered);
            return; // too large to jointly resolve -- leave untouched, never partially applied
        }

        let subs: Vec<[(AtomIdx, BondIdx); 2]> = ordered
            .iter()
            .map(|&end| {
                let s = Self::substituents(self.mol, end);
                debug_assert_eq!(s.len(), 2, "component member must be a 2-substituent end");
                [s[0], s[1]]
            })
            .collect();

        // Each end's own rank-preferred candidate index. A tie (both
        // candidates equal rank) can only mean the two substituents are
        // automorphic (per `self.ranks`'s own fully-discrete-for-genuinely-
        // distinct-atoms invariant, established at `winning_individualized_
        // ranks`) — either index then writes the same canonical string.
        let pref: Vec<usize> = subs
            .iter()
            .map(|s| {
                if self.ranks[s[1].0.0 as usize] < self.ranks[s[0].0.0 as usize] {
                    1
                } else {
                    0
                }
            })
            .collect();

        // Each end's own geometry fact — computed ONCE per end, via a FIXED
        // rank-based reference substituent (never "whichever candidate
        // happens to carry a raw mark in the CURRENT input text") — see
        // [`Self::reference_up`]. This is what makes the whole enumeration
        // below idempotent: two different valid spellings of the same real
        // molecule always mark *some* substituent at each end, but not
        // necessarily the same one, and an earlier version of this method
        // read `raw_input_direction` per-candidate inside the enumeration
        // itself, so the set of non-conflicting combinations (and hence the
        // winner) could differ between "the original input" and "this
        // method's own previous output re-parsed" — an idempotence bug, not
        // just a permutation-invariance one. `reference_up` fixes this by
        // reading the SAME fact regardless of which of the two candidates
        // the input happened to mark.
        let ref_up: Vec<Option<bool>> = ordered
            .iter()
            .zip(subs.iter())
            .zip(pref.iter())
            .map(|((&end, s), &p)| self.reference_up(end, s, p))
            .collect();

        // Every valid (non-conflicting) global assignment, as its deviation
        // bit vector (relative to `pref`) paired with the bond-value plan
        // it produces.
        let mut valid: Vec<(Vec<u8>, HashMap<BondIdx, BondOrder>)> = Vec::new();
        for mask in 0u32..(1u32 << k) {
            let choice: Vec<usize> = (0..k).map(|i| ((mask >> i) & 1) as usize).collect();
            if let Some(plan) = self.evaluate_choice(&ordered, &subs, &pref, &ref_up, &choice) {
                let bits: Vec<u8> = (0..k).map(|i| u8::from(choice[i] != pref[i])).collect();
                valid.push((bits, plan));
            }
        }

        let Some(min_count) = valid
            .iter()
            .map(|(bits, _)| bits.iter().map(|&b| u32::from(b)).sum::<u32>())
            .min()
        else {
            #[cfg(test)]
            self.ez_shared_bond_abstains.extend_from_slice(&ordered);
            return; // no consistent assignment exists at all
        };
        let mut winners: Vec<&(Vec<u8>, HashMap<BondIdx, BondOrder>)> = valid
            .iter()
            .filter(|(bits, _)| bits.iter().map(|&b| u32::from(b)).sum::<u32>() == min_count)
            .collect();
        winners.sort_by(|a, b| a.0.cmp(&b.0));
        let (winner_bits, winner_plan) = winners[0];

        // Genuine tie: swapping the deviation between two ends this
        // molecule's own canonical rank cannot tell apart (rank-equal, so
        // provably automorphic) also reaches an equally-minimal, equally
        // valid assignment. No available intrinsic signal picks a real
        // winner between them, so abstain rather than let traversal-order
        // arbitrarily decide.
        let is_tie = (0..k).any(|i| {
            (i + 1..k).any(|j| {
                self.ranks[ordered[i].0 as usize] == self.ranks[ordered[j].0 as usize]
                    && winner_bits[i] != winner_bits[j]
                    && winners.iter().any(|(bits, _)| {
                        bits[i] == winner_bits[j]
                            && bits[j] == winner_bits[i]
                            && (0..k).all(|p| p == i || p == j || bits[p] == winner_bits[p])
                    })
            })
        });

        if is_tie {
            #[cfg(test)]
            self.ez_shared_bond_abstains.extend_from_slice(&ordered);
            return;
        }

        for (bidx, order) in winner_plan.clone() {
            self.ez_marker.insert(bidx, order);
        }
    }

    /// One end's own geometry fact: whether its rank-preferred ("reference")
    /// candidate substituent (`subs[pref_idx]`) reads as "up", read via a
    /// FIXED reference — the reference bond's own raw direction if it has
    /// one, else its sibling's raw direction inverted (the trigonal-carbon
    /// sibling-complement identity: the two substituents of a stereogenic
    /// alkene end always sit on opposite sides of the double-bond axis).
    /// Returns `None` only when *neither* candidate carries any raw
    /// direction at all (no geometry info at this end).
    ///
    /// This is deliberately never "whichever candidate happens to carry a
    /// raw mark" — always the SAME (rank-based) candidate — so its result
    /// depends only on molecule topology and encoded geometry, never on
    /// which of the two candidates the CURRENT input text happened to mark.
    /// Mirrors the test-only `up_of_reference` oracle's own logic
    /// (kept as two copies, not shared, because that one intentionally
    /// stays test-only scaffolding — see its doc comment).
    fn reference_up(
        &self,
        alkene_end: AtomIdx,
        subs: &[(AtomIdx, BondIdx); 2],
        pref_idx: usize,
    ) -> Option<bool> {
        let reference = subs[pref_idx];
        let sibling = subs[1 - pref_idx];
        if let Some(dir) = self.raw_input_direction(reference.1) {
            Some(Self::direction_is_up(
                dir,
                self.mol.bond(reference.1).atom1,
                alkene_end,
            ))
        } else {
            let dir = self.raw_input_direction(sibling.1)?;
            Some(!Self::direction_is_up(
                dir,
                self.mol.bond(sibling.1).atom1,
                alkene_end,
            ))
        }
    }

    /// One candidate global assignment for a coupled component:
    /// `choice[i]` selects which of `subs[i]`'s two candidate bonds end
    /// `ordered[i]` uses as its own marker carrier. Returns the combined
    /// bond-value plan if every physical bond this component touches
    /// receives a single, consistent literal value from every end that
    /// votes on it (via [`Self::end_votes`]); `None` if any bond gets two
    /// disagreeing votes — the corruption [`Self::resolve_component_
    /// jointly`] exists to prevent.
    fn evaluate_choice(
        &self,
        ordered: &[AtomIdx],
        subs: &[[(AtomIdx, BondIdx); 2]],
        pref: &[usize],
        ref_up: &[Option<bool>],
        choice: &[usize],
    ) -> Option<HashMap<BondIdx, BondOrder>> {
        let mut votes: HashMap<BondIdx, BondOrder> = HashMap::new();
        for (i, &end) in ordered.iter().enumerate() {
            // A carrier election is not geometry-neutral when the losing
            // (demoted) candidate is itself load-bearing for a *different*
            // stereo double bond -- demoting it would silently strip that
            // other double bond's own explicit geometry, and electing the
            // winner would (via `build_ez_groups`) hand that other double
            // bond's side a spurious value it never had (issue #390's actual
            // root cause). See `is_load_bearing_elsewhere`.
            let other_bidx = subs[i][1 - choice[i]].1;
            if self.is_load_bearing_elsewhere(other_bidx, end) {
                return None;
            }
            for (bidx, value) in self.end_votes(end, &subs[i], pref[i], ref_up[i], choice[i]) {
                match votes.get(&bidx) {
                    None => {
                        votes.insert(bidx, value);
                    }
                    Some(&existing) if existing == value => {}
                    Some(_) => return None,
                }
            }
        }
        Some(votes)
    }

    /// True if `bidx` (one of `owning_end`'s own two candidate substituent
    /// bonds) carries a genuine raw input mark *and* also flanks a
    /// *different* stereo double bond's end that is itself *not* ambiguous
    /// (exactly one substituent, so `bidx` is its ONLY possible source of
    /// geometric information) -- i.e. electing `owning_end`'s *other*
    /// candidate as carrier would demote `bidx` to a plain bond, silently
    /// stripping that other, unrelated double bond of the one explicit mark
    /// it has no alternative way to recover.
    ///
    /// Deliberately excludes the case where the far end is itself ambiguous
    /// (2 substituents, i.e. also in [`Self::compute_stereo_alkene_ends`]):
    /// such an end has its own resolution machinery (possibly jointly, via
    /// [`Self::coupling_components`]/[`Self::resolve_component_jointly`])
    /// and is not solely dependent on this one candidate's raw mark --
    /// forbidding its demotion here would block genuinely coupled,
    /// resolvable systems (`EZ_SHARED_CARRIER_FULLY_RESOLVED`) for no
    /// reason, since a real conflict there is already caught by
    /// `evaluate_choice`'s own vote-consistency check.
    ///
    /// This is deliberately narrower than "does `bidx` have a raw mark":
    /// demoting a raw-marked candidate that flanks nothing else (the common,
    /// intended case `resolve_ez_markers` exists for -- picking whichever of
    /// two candidates the canonical rank prefers, regardless of which one
    /// the input happened to mark) must stay allowed. Only a candidate that
    /// is *also* someone else's sole, non-ambiguous load-bearing mark is
    /// protected.
    fn is_load_bearing_elsewhere(&self, bidx: BondIdx, owning_end: AtomIdx) -> bool {
        if self.raw_input_direction(bidx).is_none() {
            return false;
        }
        let bond = self.mol.bond(bidx);
        let Some(other_end) = bond.other(owning_end) else {
            return false;
        };
        if Self::substituents(self.mol, other_end).len() != 1 {
            return false; // other_end is itself ambiguous -- has its own resolution path
        }
        self.mol.neighbors(other_end).any(|(_, nb_bidx)| {
            nb_bidx != bidx
                && self.mol.bond(nb_bidx).order == BondOrder::Double
                && Self::end_has_substituent(self.mol, other_end)
                && self
                    .mol
                    .bond(nb_bidx)
                    .other(other_end)
                    .is_some_and(|far| Self::end_has_substituent(self.mol, far))
        })
    }

    /// One end's own votes on its two candidate bonds, given its
    /// precomputed [`Self::reference_up`] fact — generalizes the original
    /// single-end carrier/sibling branching to an arbitrary (not
    /// necessarily rank-preferred) `chosen` index, so a lone (uncoupled)
    /// end and a jointly-resolved coupled end share exactly the same
    /// per-end semantics. Always votes on *both* of the end's candidate
    /// bonds (the non-chosen one demoted to its plain order) — not
    /// conditionally, so the vote set depends only on `ref_up` (itself
    /// input-mark-placement-invariant), never on whether the non-chosen
    /// bond happened to carry a raw mark in this particular spelling.
    fn end_votes(
        &self,
        alkene_end: AtomIdx,
        subs: &[(AtomIdx, BondIdx); 2],
        pref_idx: usize,
        ref_up: Option<bool>,
        chosen_idx: usize,
    ) -> Vec<(BondIdx, BondOrder)> {
        let Some(ref_up) = ref_up else {
            return Vec::new(); // no direction info at this end at all
        };
        let chosen = subs[chosen_idx];
        let other = subs[1 - chosen_idx];
        let chosen_up = if chosen_idx == pref_idx {
            ref_up
        } else {
            !ref_up
        };
        let picked = Self::direction_for_up(self.mol.bond(chosen.1).atom1, alkene_end, chosen_up);
        vec![
            (chosen.1, picked),
            (other.1, Self::plain_order(self.mol.bond(other.1).order)),
        ]
    }

    /// Union all directional single bonds flanking each stereo double bond
    /// into one group per connected E/Z system (order-independent — depends
    /// only on molecule topology, not on canonical ranks or write order).
    fn build_ez_groups(&mut self) {
        fn find(group: &mut HashMap<BondIdx, BondIdx>, x: BondIdx) -> BondIdx {
            let parent = *group.get(&x).unwrap_or(&x);
            if parent == x {
                x
            } else {
                let root = find(group, parent);
                group.insert(x, root);
                root
            }
        }

        fn union(group: &mut HashMap<BondIdx, BondIdx>, a: BondIdx, b: BondIdx) {
            let ra = find(group, a);
            let rb = find(group, b);
            if ra != rb {
                group.insert(ra, rb);
            }
        }

        for bidx in 0..self.mol.bond_count() {
            let bidx = BondIdx(bidx as u32);
            let bond = self.mol.bond(bidx);
            if bond.order != BondOrder::Double {
                continue;
            }
            let mut side_bonds = Vec::new();
            for endpoint in [bond.atom1, bond.atom2] {
                for (_, nb_bidx) in self.mol.neighbors(endpoint) {
                    if nb_bidx == bidx {
                        continue;
                    }
                    // `effective_order` (not the raw bond order/stash) so a
                    // marker `resolve_ez_markers` moved onto a different
                    // substituent is grouped by its NEW carrier bond, not
                    // its old (now-demoted) one.
                    if matches!(
                        self.effective_order(nb_bidx),
                        BondOrder::Up | BondOrder::Down
                    ) {
                        side_bonds.push(nb_bidx);
                    }
                }
            }
            let Some(&first) = side_bonds.first() else {
                continue;
            };
            self.ez_group.entry(first).or_insert(first);
            for &b in &side_bonds[1..] {
                self.ez_group.entry(b).or_insert(b);
                union(&mut self.ez_group, first, b);
            }
        }
    }

    /// Normalize a directional bond so every member of its E/Z system is
    /// flipped consistently -- in mol-relative (`atom1`->`atom2`) terms -- to
    /// preserve geometry, with the group's shared sign pinned so the first
    /// occurrence in canonical write order prints as `Up` (`/`).
    ///
    /// This function does two genuinely different jobs, on two different
    /// frames, and conflating them is exactly how issue #390 happened:
    ///
    /// 1. **Propagation** (every call): flip [`Self::effective_order`] --
    ///    mol-relative, topology-fixed -- by the group's shared bit.
    ///    `build_ez_groups` unions bonds by molecule topology alone, so this
    ///    MUST run in that same topology-fixed frame. A bond's `atom1`/
    ///    `atom2` struct fields are fixed at parse/build time, but which
    ///    endpoint a *specific* DFS write visits first varies per occurrence
    ///    and per candidate canonical numbering; flipping an
    ///    already-write-oriented value mixes bonds visited "forward" with
    ///    bonds visited "backward" relative to their own storage, silently
    ///    inverting the very same/different relationship this function
    ///    exists to preserve whenever a group spans bonds with different
    ///    write directions (confirmed against RDKit: a coupled system with
    ///    one member's DFS direction reversed relative to another's produced
    ///    a genuinely wrong, not just differently spelled, configuration).
    /// 2. **Anchoring** (first call per group only): seed the group's shared
    ///    bit from whether *this specific occurrence*, re-oriented for
    ///    `from_atom` via [`Self::reorient_for_write`], would print as
    ///    `Down` -- because the point of the bit is to pin a *printed*
    ///    character, which only exists in write-perspective terms. Seeding
    ///    it from the mol-relative value instead pins an artifact of
    ///    whichever input text happened to be parsed (`atom1`/`atom2`
    ///    ordering is parse-order, not a canonical property): a different
    ///    spelling of the same molecule can swap `atom1`/`atom2`, seed a
    ///    different bit, and sign-flip the entire group -- geometry still
    ///    internally consistent (propagation was already fixed), but no
    ///    longer canonical or idempotent. `from_atom` is used **only** here,
    ///    never mixed into propagation.
    ///
    /// Every call site MUST separately apply [`Self::reorient_for_write`] to
    /// this function's *return value* for its own occurrence -- passing
    /// `from_atom` here does not do that; it only seeds the anchor.
    fn normalize_ez(&mut self, bidx: BondIdx, from_atom: AtomIdx) -> BondOrder {
        let order = self.effective_order(bidx);
        if !matches!(order, BondOrder::Up | BondOrder::Down) {
            return order;
        }
        let root = {
            let mut x = *self.ez_group.get(&bidx).unwrap_or(&bidx);
            while let Some(&p) = self.ez_group.get(&x) {
                if p == x {
                    break;
                }
                x = p;
            }
            x
        };
        let flip = if let Some(&f) = self.ez_flip.get(&root) {
            f
        } else {
            let bond_atom1 = self.mol.bond(bidx).atom1;
            let printed = Self::reorient_for_write(bond_atom1, from_atom, order);
            let seed = printed == BondOrder::Down;
            self.ez_flip.insert(root, seed);
            seed
        };
        if flip {
            match order {
                BondOrder::Up => BondOrder::Down,
                BondOrder::Down => BondOrder::Up,
                other => other,
            }
        } else {
            order
        }
    }

    /// Re-orient a mol-relative (`atom1`->`atom2`) directional order for
    /// reading "from `from_atom` toward the bond's other endpoint" -- the
    /// per-occurrence step every [`Self::normalize_ez`] caller must apply to
    /// its return value. Mirrors `crate::writer::direction_from` exactly
    /// (kept as a second copy since this one only ever runs after
    /// `normalize_ez`, which has no equivalent in the plain writer).
    fn reorient_for_write(bond_atom1: AtomIdx, from_atom: AtomIdx, order: BondOrder) -> BondOrder {
        match order {
            BondOrder::Up => {
                if bond_atom1 == from_atom {
                    BondOrder::Up
                } else {
                    BondOrder::Down
                }
            }
            BondOrder::Down => {
                if bond_atom1 == from_atom {
                    BondOrder::Down
                } else {
                    BondOrder::Up
                }
            }
            other => other,
        }
    }

    pub(crate) fn write_all(mut self) -> String {
        // Phase -1: pick, for every tri-/tetra-substituted stereo alkene
        // end, which substituent bond canonically carries the `/`/`\`
        // marker (topology + rank only — independent of write order, and of
        // which substituent the original parse happened to mark).
        self.resolve_ez_markers();

        // Phase 0: group directional bonds into connected E/Z systems
        // (topology-only, independent of canonical order).
        self.build_ez_groups();

        // Phase 1: discover ring-closure back-edges using the SAME canonical DFS
        // order that the writer will use. This ensures ring-closure numbers are
        // stable across re-parses.
        self.find_ring_closures();

        // Phase 2: canonical DFS serialization.
        let starts = self.canonical_atom_list();
        let mut first = true;
        for start in starts {
            if self.written[start.0 as usize] {
                continue;
            }
            if !first {
                self.out.push('.');
            }
            first = false;
            self.write_chain(start, None, None);
        }

        self.out
    }

    /// Return all atoms sorted in canonical order: highest rank first, ties
    /// broken by chemical properties invariant across re-parses.
    fn canonical_atom_list(&self) -> Vec<AtomIdx> {
        let mut atoms: Vec<AtomIdx> = (0..self.mol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .collect();
        atoms.sort_by(|&a, &b| self.canonical_cmp(b, a)); // descending
        atoms
    }

    /// Canonical ordering comparator (ascending; negate for descending).
    /// Tie-breaking uses chemical properties only (not atom indices),
    /// so the order is invariant between runs on chemically identical molecules.
    fn canonical_cmp(&self, a: AtomIdx, b: AtomIdx) -> std::cmp::Ordering {
        let ra = self.ranks[a.0 as usize];
        let rb = self.ranks[b.0 as usize];
        if ra != rb {
            return ra.cmp(&rb);
        }

        let atom_a = self.mol.atom(a);
        let atom_b = self.mol.atom(b);

        // Break ties with: atomic_number → isotope → charge → aromatic → degree
        atom_a
            .element
            .atomic_number()
            .cmp(&atom_b.element.atomic_number())
            .then_with(|| {
                atom_a
                    .isotope
                    .unwrap_or(0)
                    .cmp(&atom_b.isotope.unwrap_or(0))
            })
            .then_with(|| atom_a.charge.cmp(&atom_b.charge))
            .then_with(|| (atom_a.aromatic as u8).cmp(&(atom_b.aromatic as u8)))
            .then_with(|| self.mol.degree(a).cmp(&self.mol.degree(b)))
    }

    /// Discover back-edges by running the same canonical DFS as the writer.
    /// Using identical traversal order ensures ring-closure numbers are stable.
    fn find_ring_closures(&mut self) {
        let n = self.mol.atom_count();
        let mut visited = vec![false; n];
        let mut in_stack = vec![false; n];

        // Iterate in canonical order (same as write_all).
        let starts = self.canonical_atom_list();
        for start in starts {
            if !visited[start.0 as usize] {
                self.dfs_mark(start, None, &mut visited, &mut in_stack);
            }
        }
    }

    fn dfs_mark(
        &mut self,
        atom: AtomIdx,
        from_bond: Option<BondIdx>,
        visited: &mut Vec<bool>,
        in_stack: &mut Vec<bool>,
    ) {
        visited[atom.0 as usize] = true;
        in_stack[atom.0 as usize] = true;

        let mut neighbors: Vec<(AtomIdx, BondIdx)> = self.mol.neighbors(atom).collect();
        self.sort_neighbors_canonical(&mut neighbors);

        for (neighbor, bidx) in neighbors {
            if Some(bidx) == from_bond {
                continue;
            }
            if self.ring_bonds.contains(&bidx) {
                continue;
            }

            if !visited[neighbor.0 as usize] {
                self.dfs_mark(neighbor, Some(bidx), visited, in_stack);
            } else if in_stack[neighbor.0 as usize] {
                self.ring_bonds.insert(bidx);
                let rn = self.next_ring;
                self.next_ring += 1;
                // Discovery only decides ring numbering and which endpoint is
                // "open" (carries the marker, written when `neighbor` is
                // reached) vs "close" (suppressed, written when `atom` is
                // reached) -- NOT the actual `/`/`\` token. That is computed
                // later, in `write_chain`, from `normalize_ez`'s mol-relative
                // result re-oriented for whichever atom is actually being
                // written at consumption time (see `normalize_ez`'s doc
                // comment and `atom_ring_nums`'s field comment).
                self.atom_ring_nums
                    .entry(neighbor)
                    .or_default()
                    .push((rn, true, atom, bidx)); // open side; partner = close atom
                self.atom_ring_nums
                    .entry(atom)
                    .or_default()
                    .push((rn, false, neighbor, bidx)); // close side; partner = open atom
            }
        }

        in_stack[atom.0 as usize] = false;
    }

    /// `incoming_bond` is the already-oriented token for the edge leading to
    /// `atom` (`None` for the root or an implicit bond), not a `BondOrder`:
    /// a dative arrow's direction depends on which endpoint is written
    /// first, and `BondOrder::smiles_char` would truncate `"->"` to `'-'`
    /// anyway (issue #194). See `crate::writer::bond_token_from`.
    fn write_chain(
        &mut self,
        atom: AtomIdx,
        from_atom: Option<AtomIdx>,
        incoming_bond: Option<&'static str>,
    ) {
        self.written[atom.0 as usize] = true;

        if let Some(token) = incoming_bond {
            self.out.push_str(token);
        }

        // Compute parity-corrected chirality before ring data is consumed.
        let corrected_chirality = self.corrected_chirality(atom, from_atom);
        self.emit_atom(atom, corrected_chirality);

        // Ring-closure digits.
        if let Some(rings) = self.atom_ring_nums.remove(&atom) {
            for (rn, is_open, partner, bidx) in rings {
                // The open side carries the marker (normalize_ez's
                // mol-relative result, re-oriented for `atom` -- the
                // endpoint actually being written right now); the close
                // side is always suppressed to its plain, non-directional
                // order to avoid printing conflicting ring-closure chars at
                // both ends of the same back-edge. Mirrors dfs_mark's
                // former order_at_open/order_at_close split exactly, just
                // computed now instead of at discovery time (see
                // `normalize_ez`'s doc comment for why).
                let bond_order = if is_open {
                    let normalized = self.normalize_ez(bidx, atom);
                    Self::reorient_for_write(self.mol.bond(bidx).atom1, atom, normalized)
                } else {
                    match self.effective_order(bidx) {
                        BondOrder::Up | BondOrder::Down => {
                            Self::plain_order(self.mol.bond(bidx).order)
                        }
                        other => other,
                    }
                };
                let bond_order = suppress_standalone_wedge(self.mol, bidx, bond_order);
                // Whether a ring-closure digit needs an explicit bond-order
                // prefix depends on BOTH endpoints' aromaticity, exactly like
                // a tree-edge's own `implicit` computation below: a bare
                // digit between two aromatic atoms is read back as an
                // *aromatic* bond by the parser, so a genuinely Single
                // ring-closure bond between two aromatic atoms (e.g. a
                // non-aromatic fusion bond joining two separately-aromatic
                // ring systems, `c1-2`) must carry the `-` marker or it
                // silently becomes aromatic on re-parse -- issue #395. The
                // previous check only inspected `atom`'s own aromaticity,
                // never the ring-closure partner's.
                let atom_arom = self.mol.atom(atom).aromatic;
                let partner_arom = self.mol.atom(partner).aromatic;
                let implicit = match bond_order {
                    BondOrder::Single => !(atom_arom && partner_arom),
                    BondOrder::Aromatic => atom_arom && partner_arom,
                    _ => false,
                };
                if !implicit {
                    // Oriented from the atom being written right now, so a
                    // dative ring closure prints `->` at its donor end and
                    // `<-` at its acceptor end (the same bond read from
                    // opposite directions).
                    self.out
                        .push_str(bond_token_from(self.mol, bidx, bond_order, atom));
                }
                // SMILES ring-closure numbers are limited to 1–99.
                // Molecules needing ≥ 100 simultaneous open ring closures are
                // exotic beyond any known organic chemistry; skip extras rather
                // than panic from `char::from_digit` overflow.
                if rn > 99 {
                    continue;
                }
                if rn >= 10 {
                    self.out.push('%');
                    self.out.push(char::from_digit(rn / 10, 10).unwrap());
                    self.out.push(char::from_digit(rn % 10, 10).unwrap());
                } else {
                    self.out.push(char::from_digit(rn, 10).unwrap());
                }
            }
        }

        // Tree-edge children, sorted canonically.
        let mut children: Vec<(AtomIdx, BondIdx)> = self
            .mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .collect();

        // Sort children by canonical rank (ascending → highest rank = main chain).
        children.sort_by(|&(a, ..), &(b, ..)| self.canonical_cmp(a, b));

        let n = children.len();
        for (i, (child, bidx)) in children.into_iter().enumerate() {
            // Normalized here (not before sorting) so the flip decision is
            // made in true left-to-right write order: this atom's earlier
            // (lower-rank) children have already fully recursed by the time
            // a later sibling's direction is decided. `effective_order`
            // applies `resolve_ez_markers`'s resolved carrier choice and the
            // aromatic-bond-direction stash, both ahead of the bond's own
            // literal order; `normalize_ez` (mol-relative) must run before
            // re-orienting for `atom`'s write direction, never after -- see
            // its doc comment.
            let normalized = self.normalize_ez(bidx, atom);
            let bond_order = Self::reorient_for_write(self.mol.bond(bidx).atom1, atom, normalized);
            let bond_order = suppress_standalone_wedge(self.mol, bidx, bond_order);
            let is_last = i == n - 1;
            let parent_arom = self.mol.atom(atom).aromatic;
            let child_arom = self.mol.atom(child).aromatic;
            let implicit = match bond_order {
                BondOrder::Single => !(parent_arom && child_arom),
                BondOrder::Aromatic => parent_arom && child_arom,
                _ => false,
            };
            let written_bond = if implicit {
                None
            } else {
                Some(bond_token_from(self.mol, bidx, bond_order, atom))
            };

            if !is_last {
                self.out.push('(');
                self.write_chain(child, Some(atom), written_bond);
                self.out.push(')');
            } else {
                self.write_chain(child, Some(atom), written_bond);
            }
        }
    }

    /// Sort a neighbor list in canonical order (for consistent DFS traversal).
    fn sort_neighbors_canonical(&self, neighbors: &mut [(AtomIdx, BondIdx)]) {
        neighbors.sort_by(|&(a, _), &(b, _)| self.canonical_cmp(b, a)); // descending
    }

    fn emit_atom(&mut self, idx: AtomIdx, chirality: Chirality) {
        let atom = self.mol.atom(idx);

        if atom.wildcard {
            self.out.push_str("[*]");
            return;
        }

        // An explicit `hydrogen_count` only forces bracket notation here if it
        // carries genuine disambiguating information -- i.e. differs from
        // what organic-subset valence inference would produce anyway if the
        // atom were unbracketed. An atom whose stored H count merely repeats
        // that inferred value (e.g. `[Cl]` vs organic-subset `Cl`, both H=0)
        // is the same chemical species and must canonicalize to the same
        // (unbracketed) spelling regardless of which notation it was parsed
        // from -- unlike `writer.rs`'s plain (non-canonical) `emit_atom`,
        // which intentionally preserves original notation and is unaffected
        // by this distinction.
        let h_is_disambiguating = atom
            .hydrogen_count
            .is_some_and(|h| h != valence_inferred_hcount(self.mol, idx));

        let needs_bracket = atom.isotope.is_some()
            || atom.charge != 0
            || h_is_disambiguating
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some()
            || chirality != Chirality::None;

        if needs_bracket {
            self.out.push('[');
            if let Some(iso) = atom.isotope {
                self.out.push_str(&iso.to_string());
            }
            let sym = if atom.aromatic {
                atom.element.symbol().to_lowercase()
            } else {
                atom.element.symbol().to_string()
            };
            self.out.push_str(&sym);

            match chirality {
                Chirality::CounterClockwise => self.out.push('@'),
                Chirality::Clockwise => self.out.push_str("@@"),
                Chirality::None => {}
                Chirality::SquarePlanar(p) => self.out.push_str(square_planar_token(p)),
            }

            emit_bracket_hydrogens(&mut self.out, self.mol, idx);

            match atom.charge {
                0 => {}
                1 => self.out.push('+'),
                -1 => self.out.push('-'),
                c if c > 0 => self.out.push_str(&format!("+{c}")),
                c => self.out.push_str(&c.to_string()),
            }

            if let Some(m) = atom.atom_map {
                self.out.push(':');
                self.out.push_str(&m.to_string());
            }

            self.out.push(']');
        } else if atom.aromatic {
            self.out.push_str(&atom.element.symbol().to_lowercase());
        } else {
            self.out.push_str(atom.element.symbol());
        }
    }

    /// Compute the parity-corrected chirality for `atom` when it is written
    /// with `from_atom` as the predecessor in the canonical DFS.
    ///
    /// Tetrahedral: returns the stored chirality unchanged when no stereo
    /// neighbor order is recorded (e.g. programmatically constructed
    /// molecules) -- with only 2 possible states, "unchanged, no better
    /// information" is a safe no-op. Square-planar: drops to
    /// [`Chirality::None`] instead whenever the neighbor order can't be
    /// verified (no recorded order, size mismatch, duplicate/foreign ids) --
    /// for a 3-state tag, passing the *original* tag through against a
    /// *reordered* neighbor list can silently describe a different,
    /// plausible-but-wrong stereoisomer, exactly the failure category this
    /// mechanism exists to eliminate. See
    /// [`chematic_core::remap_square_planar_tag`] (the generalized
    /// stereo-geometry module, `docs/rfcs/generalized_stereo_geometry_rfc.md`).
    fn corrected_chirality(&self, atom: AtomIdx, from_atom: Option<AtomIdx>) -> Chirality {
        let stored = self.mol.atom(atom).chirality;
        if stored == Chirality::None {
            return Chirality::None;
        }

        let Some(original) = self.mol.stereo_neighbor_order(atom) else {
            return if stored.is_tetrahedral() {
                stored
            } else {
                Chirality::None
            };
        };

        let atom_data = self.mol.atom(atom);
        let has_h = atom_data.hydrogen_count.is_some_and(|h| h > 0);

        // Build canonical neighbor sequence in SMILES output order:
        // 1. from_atom   (or H_SENTINEL if root and has bracket H)
        // 2. bracket H   (only when from_atom is Some and has_h)
        // 3. ring-closure partners in ring-number order
        // 4. children in ascending canonical rank (branches first, main chain last)
        let mut canonical: Vec<u32> = Vec::with_capacity(original.len());

        match from_atom {
            Some(prev) => {
                canonical.push(prev.0);
                if has_h {
                    canonical.push(STEREO_H_SENTINEL);
                }
            }
            None => {
                if has_h {
                    canonical.push(STEREO_H_SENTINEL);
                }
            }
        }

        if let Some(rings) = self.atom_ring_nums.get(&atom) {
            for &(_, _, partner, _) in rings {
                canonical.push(partner.0);
            }
        }

        let mut children: Vec<AtomIdx> = self
            .mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .map(|(nb, _)| nb)
            .collect();
        children.sort_by(|&a, &b| self.canonical_cmp(a, b)); // ascending rank
        for child in children {
            canonical.push(child.0);
        }

        if canonical.len() != original.len() {
            return if stored.is_tetrahedral() {
                stored // size mismatch → tetrahedral fallback (unchanged)
            } else {
                Chirality::None // square-planar: unverifiable → drop, don't guess
            };
        }

        // Square-planar centers are always genuinely 4-coordinate, so
        // `original`/`canonical` are always 4-element there. Tetrahedral
        // centers are 4-element in the common case, but an allene *end*
        // carbon (sp2, one real double-bond partner standing in for the
        // 4th tetrahedral-like position) legitimately has only 3 entries
        // (e.g. `F[C@@H]=[C]=[C@H]Cl`'s F-bearing atom: [F, implicit-H
        // sentinel, =C partner]) -- `StereoGeometry::Tetrahedral` is fixed
        // at 4 slots and doesn't model that case, so it keeps using the
        // length-generic `permutation_is_odd` fallback unchanged for any
        // non-4 length (never reachable for square-planar in practice,
        // since that geometry has no cumulated-bond analog, but handled
        // the same way there for symmetry/safety).
        let original_arr = <[u32; 4]>::try_from(original).ok();
        let canonical_arr = <[u32; 4]>::try_from(canonical.as_slice()).ok();

        match stored {
            Chirality::CounterClockwise | Chirality::Clockwise => {
                // `unwrap_or(false)` -- i.e. "no parity flip, pass the
                // stored tag through unchanged" -- is the deliberate,
                // documented fallback for ANY `remap_tetrahedral_parity`
                // error here (`DuplicateSlotId` or `MismatchedLigandSet`),
                // not an oversight that swallows a newly-possible error
                // class. This is the exact same "unchanged, no better
                // information" fallback `corrected_chirality` already uses
                // a few lines above for "no recorded order at all" and
                // "size mismatch" -- consistent treatment for every flavor
                // of "this tetrahedral parity computation couldn't be
                // trusted," matching the 2-state-tag philosophy documented
                // on this function: passing the ORIGINAL tag through
                // against an unverified/malformed order is still a valid
                // state for a 2-state tag, just not provably correct (see
                // `docs/rfcs/square_planar_stereo_rfc.md`'s "fail-closed on
                // data-integrity problems" section for why this reasoning
                // does NOT extend to square-planar's 3-state tag, which
                // drops to `Chirality::None` instead a few lines below).
                // `MismatchedLigandSet` specifically is unreachable from any
                // real parsed molecule (this function's `original`/
                // `canonical` are always built from the same molecule's
                // actual bonded-neighbor id set), so this fallback is inert
                // in practice, not a compromise made for a real production
                // case. `debug_assert!`/`unreachable!` were deliberately
                // NOT used here instead: this exact PR's own `chematic-3d`
                // fix replaced a `debug_assert!` that was a silent
                // release-mode no-op, and reintroducing that anti-pattern
                // here -- for an error class that, unlike that bug, this
                // function itself already reports as a typed `Err` rather
                // than silently miscomputing -- would be inconsistent with
                // the lesson that fix encodes.
                let is_odd = match (original_arr, canonical_arr) {
                    (Some(o), Some(c)) => remap_tetrahedral_parity(o, c).unwrap_or(false),
                    _ => permutation_is_odd(original, &canonical),
                };
                if is_odd {
                    if stored == Chirality::CounterClockwise {
                        Chirality::Clockwise
                    } else {
                        Chirality::CounterClockwise
                    }
                } else {
                    stored
                }
            }
            Chirality::SquarePlanar(tag) => match (original_arr, canonical_arr) {
                (Some(o), Some(c)) => remap_square_planar_tag(tag, o, c)
                    .map(Chirality::SquarePlanar)
                    .unwrap_or(Chirality::None),
                _ => Chirality::None,
            },
            Chirality::None => Chirality::None,
        }
    }
}

/// Return `true` if the permutation mapping `original` order to `canonical` order
/// has odd parity (i.e. requires an odd number of transpositions).
///
/// Both slices must contain the same multiset of values (no duplicates —
/// a repeated value would collide in the by-value index below). Generic
/// over `T` (any length) so it backs two different callers: `corrected_chirality`'s
/// non-4-element fallback (an allene *end* carbon's 3-element
/// `stereo_neighbor_order` -- see that function's comment; the common
/// 4-element tetrahedral case goes through `chematic_core::remap_tetrahedral_parity`
/// instead, the generalized stereo-geometry module,
/// `docs/rfcs/generalized_stereo_geometry_rfc.md`) and the test-only
/// rank-based tetrahedral fingerprint (`u64` canonical ranks, arbitrary
/// length).
fn permutation_is_odd<T: Eq + std::hash::Hash + Copy>(original: &[T], canonical: &[T]) -> bool {
    let n = original.len();
    let mut pos: HashMap<T, usize> = HashMap::with_capacity(n);
    for (i, &v) in original.iter().enumerate() {
        pos.insert(v, i);
    }
    // perm[i] = position in `original` of the element at `canonical[i]`
    let perm: Vec<usize> = canonical
        .iter()
        .map(|v| *pos.get(v).unwrap_or(&0))
        .collect();

    // Count cycles in the permutation; parity = (n - #cycles) % 2
    let mut visited = vec![false; n];
    let mut num_cycles = 0usize;
    for start in 0..n {
        if !visited[start] {
            num_cycles += 1;
            let mut j = start;
            while !visited[j] {
                visited[j] = true;
                j = perm[j];
            }
        }
    }
    (n - num_cycles) % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    /// Positive control (section 3/19): directly verify -- by calling the
    /// legacy exhaustive enumeration itself, not by trusting a doc comment
    /// -- that budget exhaustion in the *old* design silently returns a
    /// partial/wrong result rather than failing closed. With `budget = 0`,
    /// `enumerate_discrete_ranks`'s very first non-singleton cell hits the
    /// `*budget == 0 { break; }` guard on its first loop iteration and
    /// returns an empty `Vec`; the caller's `.min_by(..).unwrap_or_default()`
    /// then silently produces `("", vec![])` -- an EMPTY canonical SMILES,
    /// not an error. This is exactly the anti-pattern the new
    /// `canonical_smiles_with_limits`/`CanonicalizationError::
    /// SearchBudgetExceeded` API (in `canonical_search.rs`) exists to
    /// retire for callers who opt into a bounded search.
    #[test]
    fn old_design_zero_budget_silently_returns_empty_string_not_an_error() {
        let mol = parse("c1ccccc1").unwrap(); // benzene: 1 non-singleton cell of all 6 atoms
        let plateaued = morgan_ranks(&mol);
        let mut budget = 0usize;
        let branches = enumerate_discrete_ranks(&mol, plateaued, &mut budget);
        assert!(
            branches.is_empty(),
            "zero budget should yield zero explored branches from the old design"
        );

        // Reproduce exactly what `legacy_winning_individualized_ranks` (and,
        // before this PR, `winning_individualized_ranks`) does with that
        // empty branch list.
        let (ranks, s) = branches
            .into_iter()
            .map(|r: Vec<u64>| {
                let s = CanonicalWriter::new(&mol, &r).write_all();
                (r, s)
            })
            .min_by(|(_, a), (_, b)| a.cmp(b))
            .unwrap_or_default();
        assert_eq!(
            s, "",
            "old design: zero budget silently yields an EMPTY canonical string"
        );
        assert!(ranks.is_empty());

        // Confirm the NEW engine does NOT reproduce this failure mode: the
        // same molecule under an equivalently tiny bounded search fails
        // closed (a typed error), never an empty string.
        let bounded = crate::canonical_search::canonical_smiles_with_limits(
            &mol,
            &crate::canonical_search::CanonicalizationLimits {
                max_search_nodes: Some(0),
                max_automorphism_tests: None,
            },
        );
        assert!(
            matches!(
                bounded,
                Err(crate::canonical_search::CanonicalizationError::SearchBudgetExceeded { .. })
            ),
            "new engine must fail closed, got {bounded:?}"
        );
    }

    /// Build a copy of `mol` with atoms reordered by `perm` (perm[new_idx] = old_idx).
    /// Bonds are remapped to the new indices; stereo/direction metadata is
    /// intentionally dropped since this helper only exists to test whether the
    /// *skeleton* rank partition is invariant under atom relabeling.
    fn permute_molecule(mol: &Molecule, perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        for (_, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let _ = builder.add_bond(a, b, bond.order);
        }
        builder.build()
    }

    /// Like [`permute_molecule`], but preserves bond order (including
    /// literal `Up`/`Down` E/Z markers) and the `bond_direction` stash,
    /// remapped onto the new atom labeling -- needed to test canonical
    /// *convergence* across genuinely different atom-labelings of the same
    /// molecule (unlike `permute_molecule`, which deliberately drops this
    /// metadata since it only tests skeleton rank-partition invariance).
    /// Bonds are added in the same relative order as `mol.bonds()` yields
    /// them, so each new bond's index equals its original one -- the
    /// `bond_direction` stash (keyed by `BondIdx`) can therefore be
    /// reapplied directly after `build()`.
    fn relabel_molecule_preserving_ez(mol: &Molecule, perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        let mut direction_stash = Vec::new();
        for (bidx, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let new_bidx = builder
                .add_bond(a, b, bond.order)
                .expect("relabeling a valid molecule's own bonds cannot fail");
            if let Some(dir) = mol.bond_direction(bidx) {
                direction_stash.push((new_bidx, dir));
            }
        }
        // Tetrahedral chirality (`Atom.chirality`) is meaningless without
        // its accompanying `stereo_neighbor_order` (the original neighbor
        // listing @/@@ is interpreted against) -- remapped here the same
        // way `Molecule::remove_atom` remaps it on atom removal (see that
        // method's own doc comment), just with no removed atom to drop.
        for &old_idx in perm {
            let old_atom = AtomIdx(old_idx as u32);
            if let Some(order) = mol.stereo_neighbor_order(old_atom) {
                let new_order: Vec<u32> = order
                    .iter()
                    .map(|&v| {
                        if v == chematic_core::STEREO_H_SENTINEL {
                            v
                        } else {
                            old_to_new[v as usize]
                        }
                    })
                    .collect();
                builder.set_stereo_neighbor_order(AtomIdx(old_to_new[old_idx]), new_order);
            }
        }
        let mut relabeled = builder.build();
        for (bidx, dir) in direction_stash {
            relabeled.set_bond_direction(bidx, dir);
        }
        relabeled
    }

    /// Relabel a rank vector into "first-seen order" group ids, so partitions
    /// can be compared structurally (which atoms are grouped together)
    /// independent of the actual numeric rank values assigned.
    fn partition_key(ranks: &[u64]) -> Vec<usize> {
        let mut seen: Vec<u64> = Vec::new();
        ranks
            .iter()
            .map(|&r| match seen.iter().position(|&s| s == r) {
                Some(pos) => pos,
                None => {
                    seen.push(r);
                    seen.len() - 1
                }
            })
            .collect()
    }

    /// Sanity check demanded before any individualize-refine rewrite: the
    /// refinement-to-plateau partition (which atoms share a rank, not the raw
    /// numeric values) must be invariant under input atom permutation. If this
    /// fails, the root cause is a non-invariant initial invariant / refinement
    /// step itself (over-splitting orbits), not a missing individualize step
    /// (under-splitting ties) -- a completely different bug to fix.
    #[test]
    fn morgan_ranks_partition_is_permutation_invariant() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
        ];

        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let n = mol.atom_count();
            let part_orig = partition_key(&morgan_ranks(&mol));

            // A few deterministic, non-identity permutations (no RNG dependency).
            let perms: Vec<Vec<usize>> = vec![
                (0..n).rev().collect(),
                {
                    let mut p: Vec<usize> = (0..n).collect();
                    if n > 2 {
                        p.rotate_left(n / 3 + 1);
                    }
                    p
                },
                {
                    let mut p: Vec<usize> = (0..n).rev().collect();
                    if n > 3 {
                        p.swap(1, n - 2);
                        p.rotate_right(2);
                    }
                    p
                },
            ];

            for perm in perms {
                let permuted = permute_molecule(&mol, &perm);
                let part_perm = partition_key(&morgan_ranks(&permuted));

                // inverse: where did old atom `old_idx` land in the permuted molecule?
                let mut new_of_old = vec![0usize; n];
                for (new_idx, &old_idx) in perm.iter().enumerate() {
                    new_of_old[old_idx] = new_idx;
                }

                for i in 0..n {
                    for j in (i + 1)..n {
                        let same_orig = part_orig[i] == part_orig[j];
                        let same_perm = part_perm[new_of_old[i]] == part_perm[new_of_old[j]];
                        assert_eq!(
                            same_orig, same_perm,
                            "partition not permutation-invariant for '{smi}': \
                             atoms {i},{j} (perm {perm:?})"
                        );
                    }
                }
            }
        }
    }

    /// `canonical_atom_order` must be permutation-invariant: relabeling the
    /// same molecule's atoms (different parse order) must not change WHICH
    /// symmetry class of atom appears 1st, 2nd, 3rd, ... in the returned
    /// order. Unlike `canonical_smiles`, `canonical_atom_order` does not run
    /// individualize-refine -- it sorts raw `morgan_ranks` with no tie-break,
    /// so this is expected to fail on any molecule with a genuine
    /// (non-singleton) rank tie. This is a diagnostic probe for that gap, not
    /// an already-passing invariant.
    #[test]
    fn canonical_atom_order_permutation_invariance_probe() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
        ];

        let mut bad = 0;
        let mut total = 0;
        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let n = mol.atom_count();
            // Ground-truth atom-class labels, in the ORIGINAL molecule's index
            // space only -- never relabeled independently for the permuted
            // copy, or two arbitrary first-seen-order numberings would be
            // compared against each other and any mismatch would be a
            // methodology artifact, not a real instability (see the
            // project's canonicalization-tie-break-theory / measurement-
            // harness-controls lessons).
            let part_orig = partition_key(&morgan_ranks(&mol));
            let order_orig = canonical_atom_order(&mol);
            let profile_orig: Vec<usize> = order_orig.iter().map(|&i| part_orig[i]).collect();

            let perms: Vec<Vec<usize>> = vec![(0..n).rev().collect(), {
                let mut p: Vec<usize> = (0..n).collect();
                if n > 2 {
                    p.rotate_left(n / 3 + 1);
                }
                p
            }];

            for perm in perms {
                total += 1;
                // perm[new_idx] = old_idx (see permute_molecule's contract).
                let permuted = permute_molecule(&mol, &perm);
                let order_perm = canonical_atom_order(&permuted);
                // Map each returned NEW index back to the class of the
                // corresponding OLD atom, via `part_orig` -- the same
                // ground-truth labeling used for profile_orig.
                let profile_perm: Vec<usize> = order_perm
                    .iter()
                    .map(|&new_i| part_orig[perm[new_i]])
                    .collect();

                if profile_orig != profile_perm {
                    bad += 1;
                    eprintln!(
                        "canonical_atom_order NOT permutation-invariant for '{smi}' (perm {perm:?}): \
                         {profile_orig:?} != {profile_perm:?}"
                    );
                }
            }
        }
        eprintln!("canonical_atom_order instability: {bad}/{total} permutation trials");
        assert_eq!(
            bad, 0,
            "{bad}/{total} permutation trials were unstable -- see stderr"
        );
    }

    /// Direct probe (no permutation needed): does `canonical_atom_order`'s
    /// naive `morgan_ranks`-only sort ever disagree with the FULLY
    /// individualized/resolved rank order that `canonical_smiles` actually
    /// verified-correct output is built from? Comparing against
    /// same-partition-class labels (as the permutation probe above does) is
    /// blind to intra-class reordering among genuinely symmetric
    /// (automorphism-equivalent) atoms, where any order is harmless -- this
    /// test instead reconstructs the winning individualized branch directly,
    /// so it also catches disagreement WITHIN a Morgan-rank tie that is not
    /// a true automorphism (exactly the class of bug individualize-refine
    /// was built to fix for `canonical_smiles`, see Round 10-12 history).
    #[test]
    fn canonical_atom_order_matches_individualized_ranks_probe() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
            // Extra cases picked for symmetry that Weisfeiler-Leman-style
            // refinement is known to struggle with (fused/bridged systems).
            "C1CC2CCC1CC2",
            "C1CC2CC1CC2",
            "c1ccc(-c2ccccc2)cc1",
            "OC1CCC(O)CC1",
            "C12CC3CC(CC(C3)C1)C2",
        ];

        let mut needed_individualization = 0;
        let mut mismatched = 0;
        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let plateaued = morgan_ranks(&mol);
            let mut budget = MAX_INDIVIDUALIZE_BRANCHES;
            let branches = enumerate_discrete_ranks(&mol, plateaued, &mut budget);
            if branches.len() > 1 {
                needed_individualization += 1;
            }
            let winning_ranks = branches
                .into_iter()
                .min_by_key(|ranks| CanonicalWriter::new(&mol, ranks).write_all())
                .expect("at least one branch");

            let n = mol.atom_count();
            let mut winning_order: Vec<usize> = (0..n).collect();
            winning_order.sort_by(|&a, &b| winning_ranks[b].cmp(&winning_ranks[a]));

            let naive_order = canonical_atom_order(&mol);
            if naive_order != winning_order {
                mismatched += 1;
                eprintln!(
                    "canonical_atom_order disagrees with resolved canonical order for '{smi}': \
                     naive={naive_order:?} resolved={winning_order:?}"
                );
            }
        }
        eprintln!(
            "{needed_individualization}/{} molecules needed individualization; \
             {mismatched}/{} disagreed with canonical_atom_order",
            corpus.len(),
            corpus.len()
        );
        assert_eq!(
            mismatched, 0,
            "canonical_atom_order must match the individualized/resolved order -- see stderr"
        );
    }

    /// Canonical SMILES must be stable: applying it twice gives the same result.
    fn is_stable(smiles: &str) -> bool {
        let mol1 = parse(smiles).expect(smiles);
        let c1 = canonical_smiles(&mol1);
        assert!(
            !c1.is_empty(),
            "canonical_smiles returned empty for '{smiles}'"
        );
        let mol2 =
            parse(&c1).unwrap_or_else(|e| panic!("canonical SMILES '{c1}' is not parseable: {e}"));
        let c2 = canonical_smiles(&mol2);
        c1 == c2
    }

    /// Two SMILES representing the same molecule must give the same canonical form.
    fn same_canonical(a: &str, b: &str) -> bool {
        let mol_a = parse(a).expect(a);
        let mol_b = parse(b).expect(b);
        canonical_smiles(&mol_a) == canonical_smiles(&mol_b)
    }

    #[test]
    fn test_methane_stable() {
        assert!(is_stable("C"));
    }
    #[test]
    fn test_ethane_stable() {
        assert!(is_stable("CC"));
    }
    #[test]
    fn test_ethanol_stable() {
        assert!(is_stable("CCO"));
    }
    #[test]
    fn test_acetic_acid_stable() {
        assert!(is_stable("CC(=O)O"));
    }
    #[test]
    fn test_benzene_stable() {
        assert!(is_stable("c1ccccc1"));
    }
    #[test]
    fn test_pyridine_stable() {
        assert!(is_stable("c1ccncc1"));
    }
    #[test]
    fn test_naphthalene_stable() {
        assert!(is_stable("c1ccc2ccccc2c1"));
    }
    #[test]
    fn test_aspirin_stable() {
        assert!(is_stable("CC(=O)Oc1ccccc1C(=O)O"));
    }
    #[test]
    fn test_caffeine_stable() {
        assert!(is_stable("Cn1cnc2c1c(=O)n(c(=O)n2C)C"));
    }

    #[test]
    fn test_ethanol_same_from_different_starts() {
        assert!(same_canonical("CCO", "OCC"));
    }

    #[test]
    fn test_isobutane_same_canonical() {
        // CC(C)C and C(C)(C)C are the same molecule.
        assert!(same_canonical("CC(C)C", "C(C)(C)C"));
    }

    #[test]
    fn test_wildcard_roundtrip() {
        let mol = parse("[*]CC").unwrap();
        let c = canonical_smiles(&mol);
        assert!(!c.is_empty());
        let mol2 = parse(&c).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert!(is_stable("[*]CC"));
    }

    #[test]
    fn test_disconnected_stable() {
        assert!(is_stable("[Na+].[Cl-]"));
    }

    // E/Z stereo bond direction tests.
    #[test]
    fn test_ez_e_stable() {
        assert!(is_stable("C/C=C/C"));
    }
    #[test]
    fn test_ez_z_stable() {
        assert!(is_stable("C/C=C\\C"));
    }
    #[test]
    fn test_ez_fluoro_e_stable() {
        assert!(is_stable("F/C=C/Cl"));
    }
    #[test]
    fn test_ez_fluoro_z_stable() {
        assert!(is_stable("F/C=C\\Cl"));
    }
    #[test]
    fn test_ez_e_ne_z() {
        // E and Z isomers of 1-fluoro-2-chloroethylene must yield different canonical forms.
        let mol_e = parse("F/C=C/Cl").unwrap();
        let mol_z = parse("F/C=C\\Cl").unwrap();
        assert_ne!(canonical_smiles(&mol_e), canonical_smiles(&mol_z));
    }

    // ── Tetrahedral stereo parity tests ─────────────────────────────────────

    #[test]
    fn test_tetrahedral_stable_no_from_atom() {
        // Bracket-H form at start of fragment — no from-atom.
        assert!(is_stable("[C@@H](F)(Cl)Br"));
        assert!(is_stable("[C@H](F)(Cl)Br"));
    }

    #[test]
    fn test_tetrahedral_stable_with_from_atom() {
        // L-alanine: chiral atom has a from-atom (N).
        assert!(is_stable("N[C@@H](C)C(=O)O"));
        assert!(is_stable("N[C@H](C)C(=O)O"));
    }

    #[test]
    fn test_enantiomers_differ() {
        // R and S configurations must give distinct canonical SMILES.
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
        assert!(!same_canonical("[C@@H](F)(Cl)Br", "[C@H](F)(Cl)Br"));
    }

    #[test]
    fn test_tetrahedral_same_from_different_starts() {
        // L-alanine from N vs from methyl — odd permutation, parity correction required.
        // RDKit: N[C@@H](C)C(=O)O and C[C@H](N)C(=O)O both → C[C@H](N)C(=O)O.
        assert!(same_canonical("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"));
        // D-alanine must differ from L-alanine.
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
    }

    #[test]
    fn test_rdkit_agreement_alanine() {
        // Pairs where the Morgan ranks distinguish all atoms unambiguously.
        // N[C@@H](C)C(=O)O and C[C@H](N)C(=O)O: same L-alanine (RDKit agrees).
        assert!(same_canonical("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"));
        // Enantiomers must differ (RDKit: C[C@@H](N)C(=O)O for D-alanine).
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
        // Stability: L-alanine canonical is self-stable.
        assert!(is_stable("N[C@@H](C)C(=O)O"));
        assert!(is_stable("C[C@H](N)C(=O)O"));
    }

    #[test]
    fn test_tetrahedral_all_heavy_substituents_stable() {
        // Chiral centre with no bracket H (all four heavy substituents).
        assert!(is_stable("[C@](F)(Cl)(Br)I"));
        assert!(is_stable("[C@@](F)(Cl)(Br)I"));
    }

    #[test]
    fn test_tetrahedral_all_heavy_enantiomers_differ() {
        assert!(!same_canonical("[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I"));
    }

    #[test]
    fn test_ring_stereocentre_stable() {
        // Chiral atom inside a ring — tests ring-closure partner resolution.
        assert!(is_stable("[C@@H]1CCCC1F"));
        assert!(is_stable("[C@H]1CCCC1F"));
    }

    #[test]
    fn test_ring_stereocentre_enantiomers_differ() {
        assert!(!same_canonical("[C@@H]1CCCC1F", "[C@H]1CCCC1F"));
    }

    #[test]
    fn test_chirality_from_different_entry_points() {
        // Same chiral molecule, two SMILES with different traversal order.
        // F[C@@H](Cl)Br  ≡  Cl[C@H](F)Br  (same S-configuration, just written
        // from different entry atoms — verified by signed-tetrahedral-volume).
        // Their canonical SMILES must be identical.
        let c1 = canonical_smiles(&parse("F[C@@H](Cl)Br").unwrap());
        let c2 = canonical_smiles(&parse("Cl[C@H](F)Br").unwrap());
        assert_eq!(c1, c2, "same molecule from different starts should match");

        // Cross-check: the enantiomer gives a different canonical form.
        let c3 = canonical_smiles(&parse("F[C@H](Cl)Br").unwrap());
        assert_ne!(c1, c3, "enantiomers must differ");
    }

    // ── Bond-order canonicality tests (#14 fix) ──────────────────────────

    #[test]
    fn test_acetic_acid_canonical_same_from_different_starts() {
        // Bug #14: both oxygens in acetic acid had the same Morgan rank because
        // the refinement loop omitted bond orders.  After the fix, O= (double)
        // and O-H (single) get distinct ranks regardless of atom insertion order.
        assert!(same_canonical("CC(=O)O", "OC(C)=O"));
        assert!(same_canonical("CC(=O)O", "O=C(O)C"));
        assert!(same_canonical("CC(=O)O", "C(C)(=O)O"));
    }

    #[test]
    fn test_oxygens_in_acetic_acid_not_equivalent() {
        // The two oxygens (O= vs O-H) are chemically distinct and must receive
        // different Morgan symmetry classes.
        let mol = parse("CC(=O)O").unwrap();
        let classes = equivalent_atom_classes(&mol);
        let o_classes: Vec<usize> = mol
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 8)
            .map(|(i, _)| classes[i.0 as usize])
            .collect();
        assert_eq!(o_classes.len(), 2);
        assert_ne!(
            o_classes[0], o_classes[1],
            "O= and O-H must be in different symmetry classes"
        );
    }

    #[test]
    fn test_formic_acid_canonical_consistent() {
        // OC=O and O=CO — same formic acid, should canonicalize identically.
        assert!(same_canonical("OC=O", "O=CO"));
    }

    // ── RDKit PR #9066: conjugated E/Z round-trip ────────────────────────────

    #[test]
    fn conjugated_double_bond_ez_round_trip() {
        // RDKit PR #9066: removeRedundantBondDirSpecs() could strip bond directions
        // on conjugated double bonds, losing E/Z stereo.  Chematic does not apply
        // aggressive direction removal, but this test guards against regressions.
        for smi in &[
            r"F/C=C/C=C/Cl", // all-E conjugated diene
            r"F/C=C\C=C\Cl", // E then Z
            r"F/C=C/C=C\Cl", // E then inverted-Z
        ] {
            let mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e:?}"));
            let out = canonical_smiles(&mol);
            let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
            let out2 = canonical_smiles(&mol2);
            assert_eq!(
                out, out2,
                "conjugated E/Z must be stable after two rounds: {smi} → {out} → {out2}"
            );
        }
    }

    // ── Round 10: ring-closure directional bond flip ─────────────────────────
    //
    // A directional marker (`/`, `\`) is read "toward" the ring digit from
    // wherever it's written. At the ring-OPENING occurrence that's already
    // the open->close direction; at the CLOSING occurrence it's close->open
    // (the opposite traversal direction over the same physical bond) and must
    // be flipped before use (parser.rs `close_or_open_ring`). Before the fix,
    // the closing-side marker was stored raw/unflipped, which silently
    // produced a *different stereoisomer* whenever a random SMILES spelling
    // routed a conjugated system's connecting single bond through a
    // ring-closure digit instead of a plain adjacent chain bond -- confirmed
    // via a corpus-wide worst-of-10 sweep (RDKit-checked structural
    // correctness, not just self-stability/idempotency, which this bug class
    // passed trivially since it was deterministic-but-wrong on each input).

    #[test]
    fn ring_closure_direction_flip_real_world_repro() {
        // Real molecule found via corpus sweep. `variant` is an RDKit
        // doRandom=True re-spelling of the exact same molecule as `orig`,
        // routing the diene's connecting single bond through ring-closure
        // digit "1" instead of a plain chain bond. Before the parser fix,
        // chematic silently emitted a different (RDKit-confirmed
        // non-equivalent) stereoisomer for `variant`.
        let orig = r"CC1CCOC(=O)/C=C/C=C\C(=O)O[C@@H]2C[C@H]3O[C@@H]4C[C@@H](C)C(=O)C[C@]4(COC(=O)C1O)[C@]2(C)C31CO1";
        let variant = r"C1=C\C(=O)O[C@@H]2C[C@H]3O[C@H]4[C@@]([C@@]2(C32OC2)C)(CC(=O)[C@H](C)C4)COC(=O)C(O)C(C)CCOC(=O)/C=C/1";
        assert!(
            same_canonical(orig, variant),
            "ring-closure-routed diene must canonicalize identically to the \
             chain-form spelling of the same molecule"
        );
    }

    #[test]
    fn ring_closure_direction_minimal_ez_agreement() {
        // Minimal case isolating the same mechanism: a ring-closure bond
        // (distinct from the exocyclic C=C double bond itself) whose
        // directional markers are specified at BOTH the opening and closing
        // occurrences of the ring digit. Per the flip rule, opposite raw
        // symbols (one `/`, one `\`) describe one consistent bond and must
        // parse successfully; same-symbol at both ends is the conflicting
        // case (unchanged by this fix -- only Up/Down are flipped, so a
        // same-vs-different Double/Single conflict, e.g. "C=1CC-1", is
        // unaffected).
        let mol = parse(r"F/C=C/1CCCC\1").unwrap_or_else(|e| panic!("{e:?}"));
        let out = canonical_smiles(&mol);
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
        assert_eq!(
            canonical_smiles(&mol2),
            out,
            "ring-closure E/Z with opposite-symbol agreement must round-trip stably"
        );

        // Same-symbol at both ends of a ring-closure directional bond is now
        // (correctly) the conflicting combination.
        assert!(matches!(
            parse(r"F/C=C/1CCCC/1"),
            Err(crate::error::SmilesError::ConflictingRingBond { ring_num: 1, .. })
        ));
    }

    // ── Round 10: ring-digit reuse racing PendingRing resolution ─────────────
    //
    // A stereocenter that OPENS a ring whose partner closes INSIDE the
    // stereocenter's own branch subtree (e.g. `[C@]1(...[C@H]1...)`) has its
    // own stereo record still unfinalized at the moment of that first
    // closure -- the immediate-resolution fast path in `close_or_open_ring`
    // only patches already-finalized records, so this case falls through to
    // the end-of-parse fallback. Before the fix, that fallback resolved by
    // raw ring DIGIT via `ring_close_partners: HashMap<u8, AtomIdx>` -- if the
    // same digit was reused later for an unrelated ring (e.g. a trailing
    // phenyl `c1ccccc1`), the later reuse's closer silently overwrote the
    // earlier, still-pending resolution, corrupting the stereocenter's
    // neighbor order with a foreign atom index (confirmed: the wrong index
    // pointed at an aromatic carbon in the unrelated trailing ring, not
    // anywhere near the stereocenter). Fixed by keying resolution on a
    // per-occurrence slot id (`next_ring_slot`) that is never reused,
    // regardless of how many times the same ring digit is.

    #[test]
    fn ring_digit_reuse_inside_stereocenter_branch_real_world_repro() {
        // Real molecule found via corpus sweep. `variant` is an RDKit
        // doRandom=True re-spelling of the exact same molecule as `orig`,
        // where the stereocenter's ring-1 partner closes inside its own
        // branch AND ring digit 1 is reused later for a trailing phenyl.
        let orig = r"COc1ccc2c3c1OC1[C@H](O)[C@](CO)(CCCCCc4ccccc4)CC4C(C2)N(C)CCC341";
        let variant = r"C([C@@]1(CC2C34CCN(C2Cc2ccc(c(c24)OC3[C@@H]1O)OC)C)CO)CCCCc1ccccc1";
        assert!(
            same_canonical(orig, variant),
            "ring-digit reuse must not corrupt a stereocenter whose own ring \
             partner closes inside its branch"
        );
    }

    #[test]
    fn ring_digit_reuse_inside_stereocenter_branch_minimal() {
        // Minimal case matching the real repro's precondition exactly:
        // `[C@H]1` opens ring 1 and its partner closes INSIDE its own first
        // branch `(CC1)` -- i.e. before the parser ever advances to a new
        // *chain* atom for atom0, so atom0's stereo record is still
        // unfinalized at the moment of that closure (the immediate-resolution
        // fast path in `close_or_open_ring` cannot catch it; only the
        // end-of-parse fallback does). Ring digit 1 is then reused by an
        // unrelated, disconnected fragment. Before the fix, the fallback
        // resolved by raw digit and the benzene's ring closure silently
        // stole atom0's still-pending resolution.
        let smi = r"[C@H]1(CC1)Cl.c1ccccc1";
        let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
        let out = canonical_smiles(&mol);
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
        assert_eq!(
            canonical_smiles(&mol2),
            out,
            "stereocenter with in-branch ring closure + later digit reuse must be stable"
        );
    }

    // ── Allene cumulated double bond stereo ──────────────────────────────────

    #[test]
    fn allene_stereo_two_enantiomers_differ() {
        // F[C@@H]=[C]=[C@H]Cl and F[C@H]=[C]=[C@@H]Cl must produce different canonical SMILES.
        let mol_r = parse("F[C@@H]=[C]=[C@H]Cl").unwrap();
        let mol_s = parse("F[C@H]=[C]=[C@@H]Cl").unwrap();
        let smi_r = canonical_smiles(&mol_r);
        let smi_s = canonical_smiles(&mol_s);
        assert_ne!(
            smi_r, smi_s,
            "allene enantiomers must produce different canonical SMILES: {smi_r}"
        );
    }

    #[test]
    fn allene_stereo_round_trip_stable() {
        for smi in &["F[C@@H]=[C]=[C@H]Cl", "F[C@H]=[C]=[C@@H]Cl"] {
            let mol = parse(smi).unwrap();
            let out = canonical_smiles(&mol);
            let mol2 = parse(&out).unwrap();
            let out2 = canonical_smiles(&mol2);
            assert_eq!(
                out, out2,
                "allene stereo must be stable: {smi} -> {out} -> {out2}"
            );
        }
    }

    /// Pinned golden-value regression for the allene-end-carbon bug this PR
    /// found and fixed: an allene end carbon (sp2, one real double-bond
    /// partner standing in for the 4th tetrahedral-like position) has a
    /// **3-element** `stereo_neighbor_order`, not 4 -- routing it through
    /// `chematic_core::remap_tetrahedral_parity` (fixed at `[u32; 4]`) would
    /// silently fall through to the "unchanged" fallback for the wrong
    /// reason (array-conversion failure, not "no verifiable order"), which
    /// produces a DIFFERENT valid-looking tag than the length-generic
    /// `permutation_is_odd` correctly computes -- this was caught only by a
    /// byte-identical before/after diff during development
    /// (`allene_stereo_two_enantiomers_differ`/`allene_stereo_round_trip_stable`
    /// above kept passing throughout, since neither checks an exact golden
    /// value, only relative invariants -- see the RFC's §8/§12.1 for the
    /// full incident writeup). `corrected_chirality`'s fallback dispatch
    /// (`canonical.rs`, the `original_arr`/`canonical_arr` `<[u32; 4]>::try_from`
    /// pair) explicitly falls back to `permutation_is_odd` whenever either
    /// side isn't exactly 4 elements -- exercised here, not just described.
    /// Pinned exact values so a future regression fails a normal
    /// `cargo test`, not only a throwaway diff tool.
    #[test]
    fn allene_stereo_exact_canonical_value_is_pinned() {
        assert_eq!(
            canonical_smiles(&parse("F[C@@H]=[C]=[C@H]Cl").unwrap()),
            "C(=[C@H]Cl)=[C@H]F"
        );
        assert_eq!(
            canonical_smiles(&parse("F[C@H]=[C]=[C@@H]Cl").unwrap()),
            "C(=[C@@H]Cl)=[C@@H]F"
        );
    }

    /// Mirrors `square_planar_stereo.rs`'s own
    /// `untagged_four_coordinate_pt_is_never_auto_promoted` for the
    /// tetrahedral case: a plain, untagged 4-distinct-substituent carbon
    /// must stay `Chirality::None` through parsing and canonicalization --
    /// never auto-assigned a `@`/`@@` tag just because it happens to have 4
    /// distinct substituents that *could* form a stereocenter. Guaranteed
    /// structurally by `corrected_chirality`'s very first line (`if stored
    /// == Chirality::None { return Chirality::None; }`, unchanged by this
    /// PR): the new `stereo_geometry` module's `canonicalize_configuration`/
    /// `remap_*` functions are never even called for an atom with no
    /// declared chirality, so there is no code path in this PR that could
    /// invent one.
    #[test]
    fn untagged_tetrahedral_center_is_never_auto_promoted() {
        let mol = parse("C(F)(Cl)(Br)I").unwrap();
        for (idx, atom) in mol.atoms() {
            assert_eq!(
                atom.chirality,
                Chirality::None,
                "atom {idx:?} unexpectedly carries stereo in an untagged fixture"
            );
        }
        let smi = canonical_smiles(&mol);
        assert!(
            !smi.contains('@'),
            "untagged input must never gain a stereo tag on write: {smi}"
        );
    }

    // ── RDKit PR #8957: fused-ring stereo round-trip ────────────────────────

    #[test]
    fn ring_stereo_stable_in_fused_system() {
        // RDKit PR #8957: "modern stereo" perception inverted R/S in fused polycyclic
        // systems with multiple stereocenters.  The canonical SMILES form may use @
        // instead of @@ (both are valid encodings of the same stereoisomer depending
        // on traversal order), so the invariant is round-trip stability, not literal
        // @@ count.
        let smi = r"CC[C@@]1(C)C[C@@](CC)(c2ccccc2)CCO1";
        let mol = parse(smi).expect("fused ring stereo mol");
        let out = canonical_smiles(&mol);
        // Round-trip must be stable: canonical(canonical(x)) == canonical(x).
        let mol2 = parse(&out).expect("canonical re-parse");
        let out2 = canonical_smiles(&mol2);
        assert_eq!(
            out, out2,
            "fused ring stereo must be stable after canonical round-trip"
        );
        // The canonical SMILES must still contain at least 2 stereocenters (@/@@ count ≥ 2).
        let stereo_count = out.matches('@').count();
        assert!(
            stereo_count >= 2,
            "both stereocenters must be encoded (got {stereo_count}): {out}"
        );
    }

    // ── E/Z directional-bond canonical stability (issue: Sprint 8) ──────────
    //
    // The canonical writer emits `/`,`\` directional bonds with traversal-direction
    // correction but no separate "normalization" pass. These tests lock in that the
    // direction choice is already deterministic and idempotent for stable skeletons,
    // so a future writer change cannot silently regress E/Z output. (The residual
    // canonical_diff idempotency failures are large fused-polycyclic atom-ranking
    // non-convergence, not a `/`,`\` direction bug — see docs/rdkit_compat.md.)

    /// E/Z parity of the first stereo double bond: `Some(true)` = E (the two
    /// *reference* substituents are on opposite sides), `Some(false)` = Z,
    /// `None` = no specified geometry.
    ///
    /// At each end, the reference substituent is whichever of its (at most
    /// two) non-double-bond neighbors has the lower `morgan_ranks` value —
    /// a deterministic, structure-derived choice (not CIP priority, but
    /// consistent across any re-parse of the same molecule, which is all
    /// this self-consistency check needs) — falling back to the sibling's
    /// marker via the trigonal-carbon complement identity when the
    /// reference substituent's own bond isn't the one marked.
    ///
    /// A prior version of this helper picked "whichever substituent bond
    /// happens to be marked, in adjacency order" instead of a fixed
    /// structural reference. That is parse-order-dependent: at a
    /// trisubstituted end, adjacency order is just first-encountered-in-the-
    /// input order, not tied to which substituent is higher priority, so
    /// comparing (say) the low-priority substituent at one end against the
    /// high-priority one at the other can silently invert the apparent
    /// parity. It happened to still pass before `resolve_ez_markers`
    /// existed only because canonicalization never relocated a marker to a
    /// *different* substituent, so both sides of the before/after
    /// comparison picked the same (arbitrary) atom by construction. Once a
    /// marker can legitimately move to a different, equally-valid carrier,
    /// that construction no longer holds, and the fixed rank-based
    /// reference is what actually verifies "the encoded geometry didn't
    /// change" rather than "the same accidental artifact survived".
    fn double_bond_is_e(smiles: &str) -> Option<bool> {
        let mol = parse(smiles).unwrap();
        double_bond_is_e_mol(&mol)
    }

    fn raw_dir(mol: &Molecule, bidx: BondIdx) -> Option<BondOrder> {
        let order = mol.bond(bidx).order;
        if matches!(order, BondOrder::Up | BondOrder::Down) {
            return Some(order);
        }
        mol.bond_direction(bidx)
    }

    /// Whether `end`'s substituent side of a double bond points "up",
    /// reading whichever substituent actually carries a direction marker
    /// via the fixed, rank-based reference (lower-`ranks` substituent) with
    /// sibling-complement fallback -- the same reference `resolve_ez_marker_
    /// for_end` uses to pick a canonical carrier. Shared by
    /// `double_bond_is_e_mol` (single representative bond) and
    /// `geometry_fingerprint` (every stereogenic bond) so there is exactly
    /// one implementation of "which substituent to trust" for tests, not
    /// two that could silently disagree.
    fn up_of_reference(mol: &Molecule, ranks: &[u64], end: AtomIdx) -> Option<bool> {
        let subs: Vec<(AtomIdx, BondIdx)> = mol
            .neighbors(end)
            .filter(|&(_, b)| mol.bond(b).order != BondOrder::Double)
            .collect();
        match subs.len() {
            1 => {
                let (_, bidx) = subs[0];
                let dir = raw_dir(mol, bidx)?;
                Some(CanonicalWriter::direction_is_up(
                    dir,
                    mol.bond(bidx).atom1,
                    end,
                ))
            }
            2 => {
                let reference = *subs
                    .iter()
                    .min_by_key(|&&(a, _)| ranks[a.0 as usize])
                    .expect("subs has 2 elements");
                let sibling = if reference.0 == subs[0].0 {
                    subs[1]
                } else {
                    subs[0]
                };
                if let Some(dir) = raw_dir(mol, reference.1) {
                    Some(CanonicalWriter::direction_is_up(
                        dir,
                        mol.bond(reference.1).atom1,
                        end,
                    ))
                } else {
                    let dir = raw_dir(mol, sibling.1)?;
                    Some(!CanonicalWriter::direction_is_up(
                        dir,
                        mol.bond(sibling.1).atom1,
                        end,
                    ))
                }
            }
            _ => None,
        }
    }

    /// Like `double_bond_is_e`, but tries every double bond in the molecule
    /// (not just the first one found) and returns the first that yields a
    /// defined answer -- a molecule can have an earlier, non-stereogenic
    /// double bond (e.g. a ketone's `C=O`, whose oxygen side has no
    /// substituent at all) ahead of its actual stereo alkene in bond order.
    fn double_bond_is_e_mol(mol: &Molecule) -> Option<bool> {
        let ranks = morgan_ranks(mol);
        mol.bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .find_map(|(_, b)| {
                let ua = up_of_reference(mol, &ranks, b.atom1)?;
                let ub = up_of_reference(mol, &ranks, b.atom2)?;
                Some(ua != ub)
            })
    }

    /// A geometry fingerprint suitable for comparing two *different* parses
    /// of the same molecule (e.g. original input vs. a canonical-output
    /// reparse) -- one E/Z fact per stereogenic double bond (both ends have
    /// at least one substituent), via `up_of_reference`. Ordered by each
    /// bond's lower-ranked endpoint: `ranks` is a molecule-intrinsic, permutation-
    /// invariant key (unlike `BondIdx`/bond-parse-order, which differs
    /// between two independently-atom-ordered spellings of the same
    /// molecule), so this ordering lines up positionally across the two
    /// parses being compared as long as neither parse gained or lost a
    /// stereogenic double bond.
    fn geometry_fingerprint(mol: &Molecule) -> Vec<Option<bool>> {
        let ranks = morgan_ranks(mol);
        let mut doubles: Vec<(u64, BondIdx)> = mol
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Double)
            .filter(|(_, b)| {
                CanonicalWriter::end_has_substituent(mol, b.atom1)
                    && CanonicalWriter::end_has_substituent(mol, b.atom2)
            })
            .map(|(bidx, b)| {
                let key = ranks[b.atom1.0 as usize].min(ranks[b.atom2.0 as usize]);
                (key, bidx)
            })
            .collect();
        doubles.sort_by_key(|&(k, _)| k);

        doubles
            .into_iter()
            .map(|(_, bidx)| {
                let bond = mol.bond(bidx);
                let ua = up_of_reference(mol, &ranks, bond.atom1)?;
                let ub = up_of_reference(mol, &ranks, bond.atom2)?;
                Some(ua != ub)
            })
            .collect()
    }

    const EZ_STABLE_CORPUS: &[&str] = &[
        "C/C=C/C",     // (E)-2-butene
        "C/C=C\\C",    // (Z)-2-butene
        "F/C=C/F",     // (E)-1,2-difluoroethene
        "F/C=C\\F",    // (Z)
        "CC/C=C/CC",   // (E)-3-hexene
        "CC/C=C\\CC",  // (Z)-3-hexene
        "C/C=C/C=C/C", // (2E,4E)-hexadiene
        "Cl/C=C/Br",
        "C/C=C/c1ccccc1", // (E)-propenylbenzene
        "C/C(F)=C(\\F)C",
    ];

    #[test]
    fn ez_canonical_smiles_is_idempotent() {
        for s in EZ_STABLE_CORPUS {
            assert!(
                is_stable(s),
                "E/Z canonical SMILES must be idempotent for {s}"
            );
        }
    }

    #[test]
    fn ez_geometry_preserved_through_canonicalization() {
        for s in EZ_STABLE_CORPUS {
            let want = double_bond_is_e(s)
                .unwrap_or_else(|| panic!("input {s} must have specified geometry"));
            let canon = canonical_smiles(&parse(s).unwrap());
            let got = double_bond_is_e(&canon)
                .unwrap_or_else(|| panic!("canonical {canon} dropped geometry from {s}"));
            assert_eq!(got, want, "E/Z geometry changed: {s} -> {canon}");
        }
    }

    #[test]
    fn ez_e_and_z_differ_for_each_skeleton() {
        // Each E form must canonicalize differently from its Z form.
        for (e, z) in [
            ("C/C=C/C", "C/C=C\\C"),
            ("F/C=C/F", "F/C=C\\F"),
            ("CC/C=C/CC", "CC/C=C\\CC"),
        ] {
            assert_ne!(
                canonical_smiles(&parse(e).unwrap()),
                canonical_smiles(&parse(z).unwrap()),
                "E and Z must produce different canonical SMILES ({e} vs {z})"
            );
        }
    }

    // ── Fused-aromatic canonical idempotency (Sprint 9) ─────────────────────
    //
    // Lock in the fused aromatics that DO round-trip consistently. The residual
    // ~1.6% canonical idempotency failures on large fused polycyclics are caused
    // by aromaticity-perception round-trip inconsistency (a molecule vs the
    // re-parse of its own canonical SMILES can disagree on which bonds are
    // aromatic — e.g. 16 vs 17 on a fluorene-type linkage — which shifts Morgan
    // ranks). That is an aromaticity/parser-core issue, not a canonical-ranking
    // bug; see docs/rdkit_compat.md. These cases are stable and guarded here.

    #[test]
    fn fused_aromatic_canonical_is_idempotent() {
        for s in [
            "c1ccc2ccccc2c1",         // naphthalene
            "c1ccc2ncccc2c1",         // quinoline
            "c1ccc2c(c1)cc[nH]2",     // indole
            "c1ccc2cc3ccccc3cc2c1",   // anthracene
            "c1ccc2[nH]c3ccccc3c2c1", // carbazole
            "c1ccc2c(c1)oc1ccccc12",  // dibenzofuran
        ] {
            assert!(
                is_stable(s),
                "fused-aromatic canonical SMILES must be idempotent for {s}"
            );
        }
    }

    // ── Round 12: simple E/Z direction normalization ────────────────────────
    //
    // `/N=N/` and `\N=N\` are two equally valid SMILES spellings of the same
    // geometry (flipping every directional bond of one connected E/Z system
    // preserves meaning). Before this fix, the writer just propagated
    // whichever direction the parser happened to read, so two spellings of
    // the same molecule could canonicalize to two different strings. The fix
    // normalizes each connected E/Z system so its first directional bond (in
    // canonical write order) is always `/`.

    #[test]
    fn ez_simple_bond_direction_normalized_azo() {
        assert!(
            same_canonical("CN(C)/N=N/c1ccccc1", r"CN(C)\N=N\c1ccccc1",),
            "isolated E/Z double bond must canonicalize identically regardless \
             of which of the two equally-valid slash spellings was parsed"
        );
    }

    #[test]
    fn ez_simple_bond_direction_normalized_symmetric() {
        assert!(same_canonical("F/C=C/F", r"F\C=C\F"));
        assert!(same_canonical("C(/F)=C/F", r"C(\F)=C\F"));
    }

    // Bracket atoms with `hydrogen_count: None` must still get their implicit
    // hydrogens written by the canonical writer too — same bug/fixtures as
    // writer::tests::test_bracket_implicit_h_*.
    use chematic_core::{Atom, Element};

    #[test]
    fn canonical_bracket_implicit_h_ammonium_charge_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut n = Atom::new(Element::N);
        n.charge = 1;
        b.add_atom(n);
        assert_eq!(canonical_smiles(&b.build()), "[NH4+]");
    }

    #[test]
    fn canonical_bracket_implicit_h_isotope_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c = Atom::new(Element::C);
        c.isotope = Some(13);
        b.add_atom(c);
        assert_eq!(canonical_smiles(&b.build()), "[13CH4]");
    }

    #[test]
    fn canonical_bracket_implicit_h_atom_map_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(canonical_smiles(&b.build()), "C[CH3:7]");
    }

    #[test]
    fn canonical_bracket_implicit_h_isotope_and_atom_map() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.isotope = Some(13);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(canonical_smiles(&b.build()), "C[13CH3:7]");
    }

    #[test]
    fn canonical_bracket_implicit_h_hydroxide_charge_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut o = Atom::new(Element::O);
        o.charge = -1;
        b.add_atom(o);
        assert_eq!(canonical_smiles(&b.build()), "[OH-]");
    }

    // ── Standalone wedge/hash bond must not emit a meaningless SMILES
    // directional token in canonical_smiles either (same fix as
    // `writer::tests`, applied to the canonical writer's own emission
    // sites) ──────────────────────────────────────────────────────────────

    /// Minimal repro from docs/rfcs/stereo2d_reader_integration_rfc.md §3: a
    /// wedge bond with no adjacent double bond must not become `/`/`\`.
    #[test]
    fn canonical_standalone_solid_wedge_not_written_as_slash() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone solid wedge (no adjacent double bond) must not be \
             written as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);

        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// Same shape but with a hash wedge (`BondOrder::Down`).
    #[test]
    fn canonical_standalone_hash_wedge_not_written_as_backslash() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Down).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone hash wedge (no adjacent double bond) must not be \
             written as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Down);

        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// A wedge bond landing on a ring-closure edge must be suppressed
    /// identically in the canonical writer's ring-closure emission site.
    #[test]
    fn canonical_standalone_wedge_on_ring_closure_bond() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c1, c3, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a ring-closure bond must not be written as \
             a directional token: got '{out}'"
        );
        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.bond_count(), 3);
    }

    /// The fix must not disturb `Atom.chirality`/`Molecule::stereo_neighbor_order`
    /// while suppressing a standalone wedge bond's spurious token on the
    /// same stereocenter -- mirrors `writer::tests::
    /// test_standalone_wedge_does_not_disturb_stereocenter_chirality` for
    /// the canonical writer's own (rank-dependent) chirality-correction path.
    #[test]
    fn canonical_standalone_wedge_does_not_disturb_stereocenter_chirality() {
        let mut center = Atom::new(Element::C);
        center.chirality = Chirality::Clockwise;
        // Force bracket notation (h=0 still forces `[..]` without printing
        // an `H` token) so the chirality symbol is actually emitted.
        center.hydrogen_count = Some(0);
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(center);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let iodine = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, iodine, BondOrder::Up).unwrap();
        b.set_stereo_neighbor_order(c, vec![f.0, cl.0, br.0, iodine.0]);
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a stereocenter's substituent must still be \
             suppressed: got '{out}'"
        );
        assert!(
            out.contains('@'),
            "the stereocenter's own chirality symbol must still be emitted \
             alongside the (now-suppressed) wedge bond: got '{out}'"
        );

        assert_eq!(mol.atom(c).chirality, Chirality::Clockwise);
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
        assert_eq!(
            mol.stereo_neighbor_order(c),
            Some([f.0, cl.0, br.0, iodine.0].as_slice())
        );
    }

    /// Regression guard: genuine E/Z directional markers must still be
    /// emitted by canonical_smiles (already covered by `ez_e_stable` etc.
    /// above; this pins the specific double-bond-adjacency mechanism the
    /// standalone-wedge fix now gates on).
    #[test]
    fn canonical_genuine_ez_directional_bond_still_written() {
        let mol = parse("F/C=C/F").unwrap();
        let out = canonical_smiles(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine double-bond-adjacent directional marker must survive \
             canonical_smiles(): got '{out}'"
        );
    }

    /// The pre-existing "E/Z direction stashed on an aromatic bond" fixtures
    /// (`parser::tests::direction_stash_*`) all involve a real exocyclic
    /// double bond and must keep round-tripping unaffected by this fix --
    /// exercised end-to-end here via `canonical_smiles` directly (those
    /// tests live in parser.rs and already re-run on every `cargo test`).
    #[test]
    fn canonical_aromatic_stash_with_real_double_bond_still_emits_direction() {
        let mol = parse(r"N=c1\c(O)c(O)\c1=N").unwrap();
        let out = canonical_smiles(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine exocyclic-double-bond-adjacent aromatic stash must \
             still emit a directional token: got '{out}'"
        );
    }

    // ── E/Z marker-carrier normalization (fix/canonical-ez-carrier-
    // normalization) ─────────────────────────────────────────────────────
    //
    // At a tri-/tetra-substituted stereo alkene end, SMILES only requires
    // ONE of its ≥2 substituent bonds to carry the `/`/`\` marker -- the
    // other's position is implied (a trigonal alkene carbon has exactly two
    // sides). *Which* substituent gets the mark used to be whatever the
    // parser happened to read, so two RDKit-valid respellings of the same
    // molecule that mark different substituents produced two different
    // canonical outputs (docs/rfcs/canonical_smiles_residual_rfc.md, Root cause
    // 1). `resolve_ez_markers` picks the marker carrier deterministically
    // from canonical rank instead. Every case below is a real molecule from
    // the residual corpus (`validation/results/
    // canonical_residual_diagnosis_summary.json`'s `permutation_invariance_
    // failures_sample`), pinned as a regression guard, not a synthetic
    // approximation of the bug.

    /// Assert `a` and `b` -- two real, already-observed-divergent canonical
    /// outputs of the SAME molecule -- now canonicalize identically, AND
    /// that the E/Z geometry each encodes is unchanged by that
    /// canonicalization (checked via `double_bond_is_e`, not string
    /// comparison, per the "zero new semantic changes" requirement).
    #[track_caller]
    fn assert_ez_carrier_pair_resolved(a: &str, b: &str) {
        let want_a = double_bond_is_e(a);
        let want_b = double_bond_is_e(b);
        let canon_a = canonical_smiles(&parse(a).unwrap());
        let canon_b = canonical_smiles(&parse(b).unwrap());
        assert_eq!(
            canon_a, canon_b,
            "two divergent real-corpus spellings of the same molecule must \
             now canonicalize identically: '{a}' -> '{canon_a}' vs '{b}' -> '{canon_b}'"
        );
        assert_eq!(
            double_bond_is_e(&canon_a),
            want_a,
            "canonicalizing '{a}' -> '{canon_a}' changed its E/Z geometry"
        );
        assert_eq!(
            double_bond_is_e(&canon_b),
            want_b,
            "canonicalizing '{b}' -> '{canon_b}' changed its E/Z geometry"
        );
    }

    /// Simplest case: a disubstituted-vs-trisubstituted amidine carbon with
    /// exactly two non-H substituents (methylamino / 4-iodobenzylamino), no
    /// ring at all. The two spellings mark different substituents.
    #[test]
    fn ez_carrier_trisub_no_ring() {
        assert_ez_carrier_pair_resolved("Ic1ccc(cc1)CN/C(=N/C)NC", "Ic1ccc(cc1)CNC(=N/C)/NC");
    }

    /// Tetrasubstituted alkene: both ends have two real substituents (no
    /// implicit H on either alkene carbon), so both ends independently pick
    /// a canonical carrier.
    #[test]
    fn ez_carrier_tetrasub_no_ring() {
        assert_ez_carrier_pair_resolved(
            "COc1ccc(/C=C(\\C)C(=O)c2cc(OC)c(OC)c(OC)c2)cc1O",
            "COc1ccc(/C=C(C(=O)c2cc(OC)c(OC)c(OC)c2)\\C)cc1O",
        );
    }

    /// One candidate substituent is reached via a plain tree edge, the
    /// other via the SAME alkene end's ring-closure bond (an aliphatic
    /// cyclopentylidene hydrazone) -- covers "E/Z bond that is also part of
    /// a ring closure" directly, with no aromaticity involved at all.
    #[test]
    fn ez_carrier_ring_closure_candidate() {
        assert_ez_carrier_pair_resolved(
            "OC(=O)CCCCCCC/1CCCC1=N/NCCCCCC",
            "OC(=O)CCCCCCC1CCC/C1=N\\NCCCCCC",
        );
    }

    /// The aromatic-bond-direction stash (`resolve_aromatic_direction_stash`
    /// in parser.rs) combined with a ring-closure candidate: a four-membered
    /// ring carbon carries an exocyclic C=N imine, and the mark can validly
    /// land on either its plain ring tree-edge or its ring-closure bond
    /// (both aromatic-aromatic, so both route through the stash side
    /// channel, not a literal `Up`/`Down` bond order). This is the exact
    /// intersection the RFC's dominant root cause calls out: tri-substituted
    /// + aromatic stash + ring closure, all at once.
    ///
    /// This same ring also has a ketone (`c(=O)`) adjacent to the imine,
    /// sharing a ring bond with it -- a ketone's carbon has two ring-bond
    /// "substituents" of its own even though a ketone has no E/Z stereo at
    /// all (the oxygen side has none). `resolve_ez_markers` must recognize
    /// the ketone end is not stereogenic (both ends of a double bond need a
    /// substituent) and never touch that shared bond on the ketone's
    /// account, or it can move/erase the imine's own marker -- this pin
    /// exercises that guard too, not just the ring-closure/stash mechanism.
    #[test]
    fn ez_carrier_aromatic_stash_ring_closure() {
        assert_ez_carrier_pair_resolved(
            "c4(c(cncc4Cl)Cl)C(=O)Nc3ccc(cc3)C[C@H](/N=c/2c(=O)c(c2N1CCSCC1)O)C(=O)O",
            "c4(c(cncc4Cl)Cl)C(=O)Nc3ccc(cc3)C[C@H](/N=c2\\c(=O)c(c2N1CCSCC1)O)C(=O)O",
        );
    }

    /// A deterministic (seeded, reproducible -- never true randomness)
    /// Fisher-Yates permutation of `0..n`, for exercising many atom
    /// relabelings without external randomness.
    ///
    /// Deliberately NOT built on `chematic_smiles::random_smiles`: probing
    /// it during this investigation found that `random_smiles(&mol, 2)` on
    /// `CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N` (one of the
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures) produces a
    /// canonical key differing at a SECOND, unrelated ring stereocenter's
    /// `@`/`@@` tag from the original -- independently confirmed via RDKit
    /// (`Chem.MolToSmiles`/`compare_molecules(..., StandardInchiString)` ->
    /// `Distinct`, a genuinely different diastereomer, not just a
    /// respelling). Whether that is a real `random_smiles` bug or a
    /// legitimate edge case was NOT investigated further here -- it is
    /// out of scope for this PR -- but using it as this gate's relabeling
    /// source would risk silently conflating that separate, unconfirmed
    /// question with the E/Z-carrier claim this test exists to check.
    fn deterministic_permutation(n: usize, seed: u64) -> Vec<usize> {
        let mut state = seed ^ 0x9E37_79B9_7F4A_7C15;
        let mut next_u32 = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) as u32
        };
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = (next_u32() as usize) % (i + 1);
            perm.swap(i, j);
        }
        perm
    }

    /// Every alternate encoding of `mol` reachable by moving ONE raw
    /// `/`/`\` mark from its current candidate bond to the SIBLING
    /// candidate at the same stereo-alkene end, via the same sibling-
    /// complement identity `end_votes` itself relies on -- i.e. a
    /// different, but geometrically identical, valid SMILES-level
    /// spelling of the exact same real molecule.
    ///
    /// This is a DIFFERENT degree of freedom than atom-index relabeling
    /// (`relabel_molecule_preserving_ez`/`deterministic_permutation`),
    /// which always preserves which physical bond carries the mark.
    /// Measured directly: pure index relabeling alone (even 16+ deterministic
    /// permutations) never reproduces the divergence the corpus diagnosis
    /// found for several [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures,
    /// because that divergence is specifically about the solver's decision
    /// changing based on WHICH candidate bond happens to carry the input's
    /// raw mark -- this helper is what actually exercises that axis.
    ///
    /// Restricted to the plain (non-aromatic-stash) case: both candidates'
    /// real chemical order must be `Single` (one currently `Up`/`Down`, the
    /// other genuinely unmarked) -- sufficient for every fixture in this
    /// file (none of the 18 issue #149 fixtures route an E/Z-relevant
    /// candidate through the aromatic-bond-direction stash), not a claim
    /// this covers every possible molecule.
    ///
    /// Self-validating: a candidate bond can be shared between TWO
    /// different ends (the exact coupling this whole module resolves), so
    /// moving its mark away for ONE end's sake, considered in isolation,
    /// can silently strip the OTHER end's only source of geometry --
    /// measured directly on one of the
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures (the
    /// enol-ether double bond's own mark vanishing as a side effect of
    /// relocating the SEPARATE, merely-bond-adjacent coupled imine's own
    /// mark). Every candidate alternate is checked against
    /// `geometry_fingerprint` before being returned -- only alternates that
    /// provably encode the exact same real-world geometry as `mol` are
    /// ever handed to a caller, so this helper can never itself be the
    /// source of a "lost stereo" false positive in the tests that use it.
    fn alternate_ez_markings(mol: &Molecule) -> Vec<Molecule> {
        let baseline_geo = geometry_fingerprint(mol);
        let ends = CanonicalWriter::compute_stereo_alkene_ends(mol);
        let mut alternates = Vec::new();
        for &end in &ends {
            let subs = CanonicalWriter::substituents(mol, end);
            if subs.len() != 2 {
                continue;
            }
            for i in 0..2 {
                let chosen = subs[i];
                let other = subs[1 - i];
                let chosen_order = mol.bond(chosen.1).order;
                let other_order = mol.bond(other.1).order;
                let chosen_dir = match chosen_order {
                    BondOrder::Up | BondOrder::Down => chosen_order,
                    _ => continue,
                };
                if other_order != BondOrder::Single {
                    continue; // other already marked, or non-plain -- not this helper's case
                }

                let up =
                    CanonicalWriter::direction_is_up(chosen_dir, mol.bond(chosen.1).atom1, end);
                let new_other_dir =
                    CanonicalWriter::direction_for_up(mol.bond(other.1).atom1, end, !up);

                let alt = mol
                    .with_bond_order(chosen.1, BondOrder::Single)
                    .with_bond_order(other.1, new_other_dir);
                if geometry_fingerprint(&alt) != baseline_geo {
                    continue; // this end's own move stripped a coupled partner's geometry -- discard
                }
                alternates.push(alt);
            }
        }
        alternates
    }

    /// The full set of relabelings + alternate spellings this file's
    /// permutation-invariance gates check: the original, its reversed-
    /// atom-order relabeling, 16 deterministic Fisher-Yates relabelings
    /// (all preserving which bond carries each mark), PLUS every
    /// mark-relocated alternate from [`alternate_ez_markings`] (which
    /// instead moves a mark to its sibling candidate, atom labeling
    /// unchanged) -- the two axes together are what the corpus-measured
    /// residual divergence actually needs to be reproduced locally.
    fn ez_carrier_test_variants(mol: &Molecule) -> Vec<Molecule> {
        let n = mol.atom_count();
        let mut variants: Vec<Molecule> = vec![
            relabel_molecule_preserving_ez(mol, &(0..n).collect::<Vec<_>>()),
            relabel_molecule_preserving_ez(mol, &(0..n).rev().collect::<Vec<_>>()),
        ];
        for seed in 0..16u64 {
            variants.push(relabel_molecule_preserving_ez(
                mol,
                &deterministic_permutation(n, seed),
            ));
        }
        variants.extend(alternate_ez_markings(mol));
        variants
    }

    /// A permutation-invariant tetrahedral-stereo fact per stereocenter,
    /// ordered by the CENTER atom's own ascending canonical rank so it
    /// lines up positionally across two different atom-labelings of the
    /// same molecule -- the tetrahedral counterpart to
    /// `geometry_fingerprint` (E/Z). Reuses [`permutation_is_odd`], the
    /// same primitive `corrected_chirality` uses to compute the write-time
    /// parity correction, but against a FIXED rank-sorted reference order
    /// instead of "whichever order the writer happened to traverse in," so
    /// relabeling can never change the reported polarity. Implicit H
    /// (`STEREO_H_SENTINEL`) gets a fixed extreme rank purely as an
    /// internally-consistent tie-break -- this does not need to match real
    /// CIP H-priority, only be the same choice on both sides of a
    /// comparison. `None` for a stereocenter whose 4 neighbor ranks are
    /// not fully discriminated (would make the by-value permutation lookup
    /// ambiguous) -- rare (only possible for locally-automorphic
    /// substituents) and never silently miscompared.
    fn tetrahedral_fingerprint(mol: &Molecule) -> Vec<Option<bool>> {
        let ranks = morgan_ranks(mol);
        let mut centers: Vec<(u64, AtomIdx)> = mol
            .atoms()
            .filter(|(_, a)| a.chirality.is_tetrahedral())
            .map(|(idx, _)| (ranks[idx.0 as usize], idx))
            .collect();
        centers.sort_by_key(|&(r, _)| r);

        centers
            .into_iter()
            .map(|(_, idx)| {
                let stored = mol.atom(idx).chirality;
                let original = mol.stereo_neighbor_order(idx)?;
                let rank_of = |v: u32| -> u64 {
                    if v == STEREO_H_SENTINEL {
                        u64::MAX
                    } else {
                        ranks[v as usize]
                    }
                };
                let original_ranks: Vec<u64> = original.iter().map(|&v| rank_of(v)).collect();
                let mut sorted_ranks = original_ranks.clone();
                sorted_ranks.sort_unstable();
                let mut deduped = sorted_ranks.clone();
                deduped.dedup();
                if deduped.len() != sorted_ranks.len() {
                    return None; // neighbor ranks not fully discriminated
                }
                let is_odd = permutation_is_odd(&original_ranks, &sorted_ranks);
                Some(match (stored, is_odd) {
                    (Chirality::Clockwise, false) => true,
                    (Chirality::Clockwise, true) => false,
                    (Chirality::CounterClockwise, false) => false,
                    (Chirality::CounterClockwise, true) => true,
                    (Chirality::None, _) => unreachable!("filtered out above"),
                    (Chirality::SquarePlanar(_), _) => {
                        unreachable!("centers is filtered to is_tetrahedral() above")
                    }
                })
            })
            .collect()
    }

    /// For every shared-candidate coupling component in `mol` recorded as
    /// abstained (production's own `ez_shared_bond_abstains`), NONE of its
    /// ends' candidate bonds carry a resolved marker -- proving "abstain"
    /// really means untouched, never a partially-applied plan left behind
    /// by an early-return path that ran after some `ez_marker` inserts
    /// already happened.
    ///
    /// Deliberately NOT asserting the stronger "0 or all of an end's own
    /// 2 candidates are marked" for every end unconditionally: that is
    /// FALSE in a real, legitimate case measured directly on one of the
    /// [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures -- an end
    /// with zero raw geometry of its own (`reference_up` returns `None`,
    /// so `end_votes` contributes nothing) can still be a bystander whose
    /// shared bond gets marked purely by its *coupled partner's* own vote.
    /// That end's OTHER (private) candidate correctly stays unmarked. This
    /// is not corruption -- it is what `geometry_fingerprint`/
    /// `tetrahedral_fingerprint` (checked separately) confirm produces the
    /// right molecule either way -- so it must not be flagged here as if
    /// it were a bug.
    fn assert_no_partial_marker_application(mol: &Molecule) {
        let ends = CanonicalWriter::compute_stereo_alkene_ends(mol);
        let components = CanonicalWriter::coupling_components(mol, &ends);
        let (ranks, _) = winning_individualized_ranks(mol);
        let mut writer = CanonicalWriter::new(mol, &ranks);
        writer.resolve_ez_markers();

        for component in &components {
            let component_abstained = component
                .iter()
                .any(|end| writer.ez_shared_bond_abstains.contains(end));
            if !component_abstained {
                continue;
            }
            for &end in component {
                let subs = CanonicalWriter::substituents(mol, end);
                let marked = subs
                    .iter()
                    .filter(|(_, b)| writer.ez_marker.contains_key(b))
                    .count();
                assert_eq!(
                    marked, 0,
                    "end {end:?} belongs to a component recorded as \
                     abstained, but still has a marked candidate bond"
                );
            }
        }
    }

    /// 18 of these 19 are real-corpus molecules from the 282-molecule
    /// `has_ez_marker` diagnosis subset (issue #149) where two independently
    /// stereogenic double bonds share one physical candidate carrier bond --
    /// the joint component solver (`resolve_component_jointly`) resolves
    /// **all** of them fully: canonical output is invariant under every
    /// relabeling tested below, not just one. The 19th (see this list's own
    /// doc comment) is a separately-found real corpus molecule with the same
    /// shared-carrier shape, added once measured to resolve with the same
    /// rigor.
    ///
    /// Originally split 10/8: the last 8 (`CC1=C2CC[C@H](/C=N/N=C(N)N)...`
    /// through `CCO/C(O)=C(\C1=NCCN1)...`) were a genuine, still-open
    /// residual until `compute_stereo_alkene_ends` gained a ring-size gate
    /// (`double_bond_endocyclic_in_small_ring`, issue #149's Wave 2C audit,
    /// `docs/rfcs/ez_ring_constrained_residual_audit.md`) that excludes an
    /// end whose own double bond is endocyclic in a ring smaller than 8
    /// atoms -- exactly the shape all 8 shared (an endocyclic culprit with
    /// no real stereochemical freedom, coupled via a shared candidate bond
    /// to a genuinely stereogenic partner). Measured directly: all 8 now
    /// converge to one canonical output, and the original 10 (re-verified,
    /// not assumed -- 5 of them share the identical endocyclic shape per
    /// the audit's own finding) remain fully resolved. Merged into one list
    /// since the split no longer describes anything real.
    const EZ_SHARED_CARRIER_FULLY_RESOLVED: &[&str] = &[
        r"CCCCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/Cc2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
        r"CCC/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"O=C(Nc1ccc(C[C@H](/N=c2\c(O)c(O)\c2=N/c2ccccc2)C(=O)O)cc1)c1c(Cl)cncc1Cl",
        r"CC(C)(C)/N=c1\c(O)c(O)\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(I)c1",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1C(F)(F)F",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccccc1OC",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccc([N+](=O)[O-])cc1",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1cccc(C(F)(F)F)c1",
        r"CC1=C2CC[C@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
        r"CC1=C2CC[C@@H](/C=N/N=C(N)N)[C@@]2(C)CC/C1=N\N=C(N)N",
        r"COC(=O)/C=C/[C@H]1CCC2=C(C)/C(=N/N=C(N)N)CC[C@@]21C",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1ccc(I)cc1",
        r"CCO/C(O)=C1\C(C)=NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(C)c1",
        r"CCOC(=O)C1=C(C)N=C(C)/C(=C(/O)OCC)C1c1cccc(OC)c1",
        r"CCO/C(O)=C(\C1=NCCN1)c1nnc(N)s1",
        r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC",
    ];

    /// Proves all 19 [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] fixtures are
    /// genuinely, fully permutation-invariant -- not just under the one
    /// relabeling a weaker test might check. Per fixture: the original
    /// parse, its reversed-atom-order relabeling, and 16 deterministic
    /// (seeded) Fisher-Yates relabelings (see [`deterministic_permutation`]
    /// for why `random_smiles` is not used) -- 18 spellings total -- must
    /// all canonicalize to exactly ONE string. That string must be
    /// idempotent (`canonical(canonical(x)) == canonical(x)`), must
    /// re-parse without error, and a reparse of it must preserve both the
    /// E/Z ([`geometry_fingerprint`]) and tetrahedral
    /// ([`tetrahedral_fingerprint`]) stereo facts read from the original
    /// parse. Also asserts the joint solver never abstains for these 19
    /// (production's own `ez_shared_bond_abstains` record, not re-derived
    /// topology) and never applies a partial marker plan.
    ///
    /// The 19th fixture (two independently stereogenic exocyclic imines on
    /// the same four-membered ring, sharing the ring-closure bond each
    /// would otherwise use as a candidate carrier) was previously believed
    /// to be a genuine, documented residual that a fully general carrier
    /// choice could not resolve to one canonical string -- tracked by a
    /// separate, weaker test that only checked E/Z-preservation across two
    /// specific hand-picked respellings, not full permutation-invariance.
    /// Re-measured directly (not assumed) while auditing a stale-doc-comment
    /// flag in ROADMAP.md: it passes this test's full 18-spelling rigor,
    /// same as the other 18 -- the residual is resolved, and the split has
    /// been merged into this one list for the same reason the original 8/10
    /// split above was.
    #[test]
    fn ez_shared_carrier_fully_resolved_are_permutation_invariant() {
        for &s in EZ_SHARED_CARRIER_FULLY_RESOLVED {
            let mol = parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"));

            let (ranks, _) = winning_individualized_ranks(&mol);
            let mut writer = CanonicalWriter::new(&mol, &ranks);
            writer.resolve_ez_markers();
            assert!(
                writer.ez_shared_bond_abstains.is_empty(),
                "'{s}': expected the joint component solver to resolve \
                 without abstaining -- this fixture may need re-diagnosis"
            );
            assert_no_partial_marker_application(&mol);

            let variants = ez_carrier_test_variants(&mol);
            assert!(
                variants.len() >= 18,
                "sanity: at least 16 deterministic + 2 fixed relabelings"
            );

            let canonical_outputs: Vec<String> = variants.iter().map(canonical_smiles).collect();
            let unique: HashSet<&String> = canonical_outputs.iter().collect();
            assert_eq!(
                unique.len(),
                1,
                "'{s}': expected ONE canonical output across {} relabelings/\
                 markings (fully resolved), got {} distinct: {:?}",
                variants.len(),
                unique.len(),
                unique
            );
            let canon = canonical_outputs[0].clone();

            // Idempotence.
            let reparsed_once =
                parse(&canon).unwrap_or_else(|e| panic!("'{s}': reparse of '{canon}': {e}"));
            let canon_twice = canonical_smiles(&reparsed_once);
            assert_eq!(
                canon, canon_twice,
                "'{s}': canonical(canonical(x)) must equal canonical(x)"
            );

            // Zero corruption: E/Z and tetrahedral facts survive the
            // canonicalize -> reparse round trip.
            let before_geo = geometry_fingerprint(&mol);
            let after_geo = geometry_fingerprint(&reparsed_once);
            assert_eq!(
                before_geo, after_geo,
                "'{s}': E/Z geometry must be preserved by canonicalization"
            );
            let before_tet = tetrahedral_fingerprint(&mol);
            let after_tet = tetrahedral_fingerprint(&reparsed_once);
            assert_eq!(
                before_tet, after_tet,
                "'{s}': tetrahedral stereo must be preserved by canonicalization"
            );
            assert!(
                before_geo.iter().any(|f| f.is_some()),
                "test setup sanity: '{s}' must have at least one defined \
                 E/Z geometry fact"
            );
        }
    }

    /// Issue #390 witness: `O/N=C/C(C=N/O)=N\NC`. A prior version of
    /// `normalize_ez` silently changed this molecule's atom3=atom7 double
    /// bond from Z to E -- confirmed genuinely wrong (not just a different,
    /// equally valid spelling) via an independent RDKit `MolToInchi`/
    /// `GetStereo()` check, before any fix existed.
    ///
    /// Root cause, in two parts:
    /// 1. `resolve_ez_markers`'s carrier election for the ambiguous end
    ///    (2 candidates: the bond toward the atom1=atom2 double bond's own
    ///    side, raw-marked and load-bearing for atom1=atom2's geometry; the
    ///    bond toward the atom4=atom5 double bond's side, unmarked --
    ///    atom4=atom5 is itself genuinely undefined in this molecule, only
    ///    one of its two flanking bonds ever carries a mark, confirmed via
    ///    InChI's own `?` stereo descriptor for it) could elect the
    ///    unmarked candidate. That demotes the marked one -- silently
    ///    under-specifying atom1=atom2 -- while handing atom4=atom5 a
    ///    geometry it never had. Neither move is neutral. Fixed by
    ///    [`Self::is_load_bearing_elsewhere`]: an election must not demote a
    ///    raw-marked candidate that is some other, non-ambiguous double
    ///    bond's only geometric anchor.
    /// 2. Independently, `build_ez_groups`/`normalize_ez` chained the
    ///    elected carrier's group together with atom4=atom5's own group
    ///    (both touch the same connecting bond), and the group's sign was
    ///    seeded from a value that had already been re-oriented for one
    ///    specific DFS write direction -- which bond within a group gets
    ///    visited "forward" vs "backward" varies across candidate canonical
    ///    numberings for reasons unrelated to this bond's own geometry, so
    ///    the seeded sign could vary too. Fixed by splitting `normalize_ez`
    ///    into a mol-relative propagation step (always operates on
    ///    [`Self::effective_order`], never an already-write-oriented value)
    ///    and a write-perspective anchor-seeding step (uses the write atom
    ///    ONLY to decide the group's shared sign once, never to decide what
    ///    gets propagated).
    #[test]
    fn issue390_witness_geometry_preserved_and_stable() {
        let smi = "O/N=C/C(C=N/O)=N\\NC";
        let mol = parse(smi).unwrap();
        let canon = canonical_smiles(&mol);
        assert_eq!(
            canon, smi,
            "canonical form must match this already-canonically-written input exactly"
        );

        let reparsed = parse(&canon).unwrap();
        let canon_twice = canonical_smiles(&reparsed);
        assert_eq!(
            canon, canon_twice,
            "canonical(canonical(x)) must equal canonical(x)"
        );

        let before_geo = geometry_fingerprint(&mol);
        let after_geo = geometry_fingerprint(&reparsed);
        assert_eq!(
            before_geo, after_geo,
            "E/Z geometry must be preserved by canonicalization"
        );
        assert!(
            before_geo.iter().any(|f| f.is_some()),
            "test setup sanity: witness must have at least one defined E/Z fact"
        );
    }

    /// The issue #390 witness's ambiguous end has exactly one *safe*
    /// candidate: the bond toward atom1=atom2, which is load-bearing for
    /// that unrelated double bond's own geometry. Its sibling (toward
    /// atom4=atom5, itself undefined) is therefore NOT a valid
    /// geometry-preserving alternate carrier -- moving the mark there would
    /// silently strip atom1=atom2's geometry, exactly the corruption
    /// [`Self::is_load_bearing_elsewhere`] exists to forbid.
    /// [`alternate_ez_markings`] independently verifies geometry
    /// preservation before offering an alternate, so it returning none here
    /// is itself a direct, positive confirmation the fix is in effect --
    /// not a vacuous negative.
    #[test]
    fn issue390_witness_has_no_valid_alternate_marking() {
        let mol = parse("O/N=C/C(C=N/O)=N\\NC").unwrap();
        let alternates = alternate_ez_markings(&mol);
        assert!(
            alternates.is_empty(),
            "the only other candidate is load-bearing for atom1=atom2's own geometry; \
             moving the mark there must not be offered as a valid alternate"
        );
    }

    /// Mirror/stereoisomer distinctness: the issue #390 witness's E and Z
    /// forms must remain genuinely distinguishable after canonicalization,
    /// never collapse into the same string. Direct regression test for a
    /// mid-investigation defect (never shipped) where an earlier attempted
    /// fix made canonicalization informationally lossy for this coupled
    /// shape -- both the true-Z and true-E spellings canonicalized to
    /// identical output regardless of which geometry was actually input.
    #[test]
    fn issue390_mirror_stereoisomers_stay_distinct() {
        let z = canonical_smiles(&parse("O/N=C/C(C=N/O)=N\\NC").unwrap());
        let e = canonical_smiles(&parse("O/N=C/C(C=N/O)=N/NC").unwrap());
        assert_ne!(
            z, e,
            "genuinely different E/Z witnesses must not canonicalize to the same string"
        );
    }

    /// Atom-order permutation invariance for the issue #390 witness: every
    /// deterministic relabeling [`ez_carrier_test_variants`] produces
    /// (identity, reversed, 16 seeded Fisher-Yates shuffles, plus any
    /// geometry-preserving alternate markings) must canonicalize to the
    /// identical string, and that string must be idempotent and preserve
    /// the input's own E/Z geometry.
    #[test]
    fn issue390_witness_permutation_invariant() {
        let mol = parse("O/N=C/C(C=N/O)=N\\NC").unwrap();
        let variants = ez_carrier_test_variants(&mol);
        assert!(
            variants.len() >= 18,
            "sanity: at least 16 deterministic + 2 fixed relabelings"
        );

        let outputs: Vec<String> = variants.iter().map(canonical_smiles).collect();
        let unique: HashSet<&String> = outputs.iter().collect();
        assert_eq!(
            unique.len(),
            1,
            "expected ONE canonical output across {} relabelings/markings, got {}: {:?}",
            variants.len(),
            unique.len(),
            unique
        );

        let canon = outputs[0].clone();
        let reparsed = parse(&canon).unwrap_or_else(|e| panic!("reparse of '{canon}': {e}"));
        assert_eq!(
            canonical_smiles(&reparsed),
            canon,
            "canonical(canonical(x)) must equal canonical(x)"
        );

        let before_geo = geometry_fingerprint(&mol);
        let after_geo = geometry_fingerprint(&reparsed);
        assert_eq!(
            before_geo, after_geo,
            "E/Z geometry must be preserved across permutation + canonicalization"
        );
    }

    /// Direct empirical check of the invariant [`Self::resolve_component_
    /// jointly`]'s tie-abstain path relies on: `winning_individualized_
    /// ranks` is a strict total order over EVERY atom, with no ties, even
    /// for a genuinely symmetric molecule (swapping the two coupled ends
    /// is an actual graph automorphism here, not just a superficial
    /// resemblance) -- built via `MoleculeBuilder` directly: methyl-
    /// hydrazone-hydrazone-methyl, `P_A-A(=N_A-N_A2)-B(=N_B-N_B2)-P_B` with
    /// `A`/`B` directly bonded (the shared/coupled candidate) and each end
    /// independently wedge-marked so there is real E/Z content at stake.
    /// The extra terminal `N_A2`/`N_B2` atoms are required, not decorative
    /// -- `compute_stereo_alkene_ends` only counts a double bond as
    /// stereogenic when BOTH ends have >=1 substituent (matching
    /// `chematic_chem::cip::assign_ez`'s own guard), and a bare `C=N` imine
    /// nitrogen has none. If this ever starts reporting equal ranks for
    /// the two automorphic ends, the tie-abstain path documented on
    /// `resolve_component_jointly` becomes reachable in practice, not just
    /// defensively coded -- this test exists to catch that change, not to
    /// assert it can never happen.
    #[test]
    fn individualized_ranks_never_tie_even_for_symmetric_molecules() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        let mut b = MoleculeBuilder::new();
        let pa = b.add_atom(Atom::new(Element::C)); // P_A methyl
        let a = b.add_atom(Atom::new(Element::C)); // A (coupled alkene end)
        let na = b.add_atom(Atom::new(Element::N)); // N_A
        let na2 = b.add_atom(Atom::new(Element::N)); // N_A's own substituent
        let bb = b.add_atom(Atom::new(Element::C)); // B (coupled alkene end)
        let nb = b.add_atom(Atom::new(Element::N)); // N_B
        let nb2 = b.add_atom(Atom::new(Element::N)); // N_B's own substituent
        let pb = b.add_atom(Atom::new(Element::C)); // P_B methyl
        b.add_bond(pa, a, BondOrder::Up).unwrap();
        b.add_bond(a, na, BondOrder::Double).unwrap();
        b.add_bond(na, na2, BondOrder::Single).unwrap();
        b.add_bond(a, bb, BondOrder::Single).unwrap(); // shared/coupled candidate
        b.add_bond(bb, nb, BondOrder::Double).unwrap();
        b.add_bond(nb, nb2, BondOrder::Single).unwrap();
        b.add_bond(bb, pb, BondOrder::Up).unwrap();
        let mol = b.build();

        let (ranks, _) = winning_individualized_ranks(&mol);
        assert_ne!(
            ranks[a.0 as usize], ranks[bb.0 as usize],
            "individualized ranks are documented as a strict total order \
             even for automorphic atoms -- if this now ties, the \
             tie-abstain path this test's sibling exercises via a forced \
             rank array has become reachable through real input too"
        );

        // Sanity: confirm the two ends really are coupled (share the A-B
        // bond as a mutual candidate) and the solver resolves them (no
        // rank tie means no abstain), rather than this being a vacuous
        // check on a topology that never reaches the tie-break logic.
        let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
        let components = CanonicalWriter::coupling_components(&mol, &ends);
        assert!(
            components.iter().any(|c| c.len() == 2),
            "test setup sanity: A and B must form a genuine size-2 coupling \
             component"
        );
        let mut writer = CanonicalWriter::new(&mol, &ranks);
        writer.resolve_ez_markers();
        assert!(
            writer.ez_shared_bond_abstains.is_empty(),
            "test setup sanity: with distinct ranks, the solver should \
             resolve rather than abstain"
        );
    }

    /// White-box probe of the tie-abstain path, and a genuine (not assumed)
    /// finding about why it cannot be exercised through a size-2 component
    /// -- the ONLY size ever observed in real molecules (measured, see
    /// [`Self::MAX_JOINT_COMPONENT_SIZE`]'s doc comment).
    ///
    /// Built the SAME symmetric molecule as
    /// [`individualized_ranks_never_tie_even_for_symmetric_molecules`],
    /// then FORCED a rank tie between the two coupled ends via a
    /// hand-crafted `ranks` array (bypassing `winning_individualized_ranks`
    /// entirely, since that function never produces one on its own) --
    /// tried under two different mark configurations (`P_A`/`P_B` both
    /// `Up`, and `Up`/`Down`). In BOTH cases the solver still resolved
    /// uniquely, never abstaining, even with the forced tie. Traced why:
    /// for a size-2 component, a "mixed" combination (one end's chosen
    /// candidate is the shared bond, the other's is its own private one)
    /// ALWAYS conflicts by construction -- the end treating the shared
    /// bond as "chosen" votes a real direction on it, while the end
    /// treating it as "other" always votes to demote it to plain
    /// (`Self::plain_order`), and a direction can never equal "plain". So
    /// only two combinations are ever valid at size 2: BOTH ends use the
    /// shared bond, or BOTH use their own private one -- never a case
    /// where the two ends' deviation bits differ, which is the specific
    /// precondition `is_tie` checks for (`winner_bits[i] != winner_bits[j]`).
    /// A rank tie between the two ends is therefore insufficient on its
    /// own to create genuine ambiguity at size 2, REGARDLESS of which
    /// marks are present -- there is no "swap which one deviates" case to
    /// be ambiguous about, only "both agree" (resolved) or "both would
    /// conflict, so both fall back" (also resolved, uniquely).
    ///
    /// This does NOT prove `is_tie` can never fire -- a genuinely
    /// unreachable-in-practice size-3+ component (never observed in any
    /// corpus scanned) could still expose a real swap ambiguity between
    /// two of its rank-tied members while a third breaks the symmetry
    /// differently; that combination was not constructed here (todo,
    /// tracked as a documented gap rather than silently assumed safe).
    /// What IS established: the defensive `is_tie` check is unreachable
    /// for the only topology this crate's own corpus scans have ever
    /// produced, for two independent, now-verified reasons (ranks don't
    /// tie; and even a forced tie can't create a swap-style ambiguity at
    /// size 2) -- not merely "happens not to trigger today".
    #[test]
    fn size_two_component_forced_rank_tie_still_resolves_uniquely() {
        use chematic_core::{Atom, Element, MoleculeBuilder};
        fn build(pb_order: BondOrder) -> (Molecule, AtomIdx, AtomIdx) {
            let mut b = MoleculeBuilder::new();
            let pa = b.add_atom(Atom::new(Element::C));
            let a = b.add_atom(Atom::new(Element::C));
            let na = b.add_atom(Atom::new(Element::N));
            let na2 = b.add_atom(Atom::new(Element::N));
            let bb = b.add_atom(Atom::new(Element::C));
            let nb = b.add_atom(Atom::new(Element::N));
            let nb2 = b.add_atom(Atom::new(Element::N));
            let pb = b.add_atom(Atom::new(Element::C));
            b.add_bond(pa, a, BondOrder::Up).unwrap();
            b.add_bond(a, na, BondOrder::Double).unwrap();
            b.add_bond(na, na2, BondOrder::Single).unwrap();
            b.add_bond(a, bb, BondOrder::Single).unwrap();
            b.add_bond(bb, nb, BondOrder::Double).unwrap();
            b.add_bond(nb, nb2, BondOrder::Single).unwrap();
            b.add_bond(bb, pb, pb_order).unwrap();
            (b.build(), a, bb)
        }

        for pb_order in [BondOrder::Up, BondOrder::Down] {
            let (mol, a, bb) = build(pb_order);
            let (mut ranks, _) = winning_individualized_ranks(&mol);
            // Force a tie between the two coupled ends, preserving every
            // other atom's real rank -- the minimal, targeted way to
            // inject the one precondition `is_tie` checks for.
            ranks[a.0 as usize] = ranks[bb.0 as usize];

            let mut writer = CanonicalWriter::new(&mol, &ranks);
            writer.resolve_ez_markers();
            assert!(
                writer.ez_shared_bond_abstains.is_empty(),
                "pb_order={pb_order:?}: a forced rank tie at size 2 should \
                 still resolve uniquely (see this test's own doc comment \
                 for why no swap-style ambiguity is possible here) -- an \
                 abstain here would mean the structural argument above no \
                 longer holds and needs re-examining, not just re-asserting"
            );
            assert_no_partial_marker_application(&mol);
        }
    }

    /// [`Self::coupling_components`] deliberately sorts its `HashSet`
    /// iteration on raw `AtomIdx` purely to avoid process-random traversal
    /// order (see that method's own doc comment) -- never as a canonical
    /// tie-break. This proves that choice is sufficient: rebuilding the
    /// SAME set of stereo-alkene ends via several different insertion
    /// orders always yields the SAME partition into components (compared
    /// as sets of atoms, never as the returned `Vec`'s order), across
    /// every fixture in [`EZ_SHARED_CARRIER_FULLY_RESOLVED`] -- i.e.
    /// component MEMBERSHIP is a pure function of molecule structure, never
    /// of `HashSet` build/iteration order.
    #[test]
    fn coupling_components_membership_is_insertion_order_independent() {
        fn normalized_components(
            mol: &Molecule,
            ends: &HashSet<AtomIdx>,
        ) -> Vec<std::collections::BTreeSet<u32>> {
            let mut components: Vec<std::collections::BTreeSet<u32>> =
                CanonicalWriter::coupling_components(mol, ends)
                    .into_iter()
                    .map(|c| c.into_iter().map(|a| a.0).collect())
                    .collect();
            components.sort();
            components
        }

        for &s in EZ_SHARED_CARRIER_FULLY_RESOLVED {
            let mol = parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"));
            let ends = CanonicalWriter::compute_stereo_alkene_ends(&mol);
            let baseline = normalized_components(&mol, &ends);

            let mut as_vec: Vec<AtomIdx> = ends.iter().copied().collect();
            as_vec.sort_by_key(|a| a.0);
            let mut orders: Vec<Vec<AtomIdx>> = vec![as_vec.clone()];
            let mut reversed = as_vec.clone();
            reversed.reverse();
            orders.push(reversed);
            for seed in 0..8u64 {
                let perm = deterministic_permutation(as_vec.len(), seed);
                orders.push(perm.into_iter().map(|i| as_vec[i]).collect());
            }

            for order in orders {
                let rebuilt: HashSet<AtomIdx> = order.into_iter().collect();
                let variant = normalized_components(&mol, &rebuilt);
                assert_eq!(
                    baseline, variant,
                    "'{s}': coupling_components membership must not depend \
                     on HashSet insertion order"
                );
            }
        }
    }

    /// Positive control for [`tetrahedral_fingerprint`]: without this, the
    /// preservation checks above would pass vacuously if the fingerprint
    /// function were broken and always returned e.g. all-`None` or some
    /// other input-insensitive constant. Flipping one real `[C@H]` ->
    /// `[C@@H]` tag in one of the fully-resolved fixtures must change the
    /// fingerprint, proving the function actually reads the input.
    #[test]
    fn tetrahedral_fingerprint_is_sensitive_to_a_real_flip() {
        let original = "CCCCC/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let flipped = "CCCCC/N=c1\\c(O)c(O)\\c1=N/[C@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let fp_original = tetrahedral_fingerprint(&parse(original).unwrap());
        let fp_flipped = tetrahedral_fingerprint(&parse(flipped).unwrap());
        assert_ne!(
            fp_original, fp_flipped,
            "flipping a real [C@H]->[C@@H] tag must change the fingerprint"
        );
    }

    /// Positive control for `geometry_fingerprint`: without this, the
    /// "before == after" check above would pass vacuously if the fingerprint
    /// function were broken and always returned e.g. all-`None` or some
    /// other input-insensitive constant. Flipping one real mark in one of
    /// the residual fixtures (`/N=c1\...` -> `/N=c1/...`) must change the
    /// fingerprint, proving the function actually reads the input.
    #[test]
    fn geometry_fingerprint_is_sensitive_to_a_real_flip() {
        let original = "CCCCC/N=c1\\c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let flipped = "CCCCC/N=c1/c(O)c(O)\\c1=N/[C@@H](Cc1ccc(NC(=O)c2c(Cl)cncc2Cl)cc1)C(=O)O";
        let fp_original = geometry_fingerprint(&parse(original).unwrap());
        let fp_flipped = geometry_fingerprint(&parse(flipped).unwrap());
        assert_ne!(
            fp_original, fp_flipped,
            "geometry_fingerprint must be sensitive to a real E/Z flip, or \
             the no-corruption check above is vacuous"
        );
    }
}

// ── Explicit-vs-implicit hydrogen-count canonicalization invariance ────────
//
// Regression suite for issue #205 (kent-tokyo/renkin PR #65 finding): an
// atom's `hydrogen_count` field is `Some(n)` when it was written with
// bracket notation and `None` when written organic-subset, but these can
// represent the *same* chemical state (n happens to equal what valence
// inference would give anyway). Before this fix, `initial_invariant`
// treated `None` and `Some(n)` as different Morgan-rank seeds, and
// `emit_atom` treated any `Some(_)` as forcing bracket notation regardless
// of whether it was redundant -- so two representations of one molecule
// (e.g. a bracket `[Cl]` substituent vs. organic-subset `Cl`) could
// canonicalize to different strings. Fixed by routing both decisions
// through `implicit_hcount`/`valence_inferred_hcount` instead of the raw
// field. These tests prove the invariant mechanically (structural/data
// comparisons), not by assuming two SMILES "look like" the same molecule.
#[cfg(test)]
mod explicit_implicit_h_invariance {
    use super::*;
    use crate::parser::parse;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn atom_multiset(mol: &Molecule) -> Vec<(u8, i8, Option<u16>, bool, u8)> {
        // (atomic_number, charge, isotope, aromatic, effective_h_count) --
        // effective (not raw) H count, since that's the structural fact that
        // must be invariant; two representations legitimately differ in
        // whether hydrogen_count is None or Some(same value).
        let mut v: Vec<_> = mol
            .atoms()
            .map(|(idx, a)| {
                (
                    a.element.atomic_number(),
                    a.charge,
                    a.isotope,
                    a.aromatic,
                    chematic_core::implicit_hcount(mol, idx),
                )
            })
            .collect();
        v.sort();
        v
    }

    fn bond_multiset(mol: &Molecule) -> Vec<(u8, u8, u8)> {
        // (min_atomic_number, max_atomic_number, bond_order) at each edge,
        // order-independent identification of the bond multiset.
        let mut v: Vec<_> = mol
            .bonds()
            .map(|(_, b)| {
                let e1 = mol.atom(b.atom1).element.atomic_number();
                let e2 = mol.atom(b.atom2).element.atomic_number();
                (e1.min(e2), e1.max(e2), bond_order_value(b.order) as u8)
            })
            .collect();
        v.sort();
        v
    }

    /// Asserts `a` and `b` are the same chemical graph (same atom multiset,
    /// same bond multiset) *before* checking canonical string equality --
    /// this is what makes the string-equality assertions below a proof of
    /// canonicalization invariance rather than an assumption that two
    /// differently-written SMILES denote the same molecule.
    fn assert_same_graph_then_same_canonical(a: &Molecule, b: &Molecule, label: &str) {
        assert_eq!(
            a.atom_count(),
            b.atom_count(),
            "{label}: atom count must match for this to be a meaningful test"
        );
        assert_eq!(
            atom_multiset(a),
            atom_multiset(b),
            "{label}: atom multisets (element/charge/isotope/aromatic/effective-H) differ -- \
             these are not actually the same molecule, so canonical equality would prove nothing"
        );
        assert_eq!(
            bond_multiset(a),
            bond_multiset(b),
            "{label}: bond multisets differ -- not the same graph"
        );
        assert_eq!(
            canonical_smiles(a),
            canonical_smiles(b),
            "{label}: same graph must canonicalize identically"
        );
    }

    /// Isolated version of the kent-tokyo/renkin PR #65 fixture with zero
    /// reaction machinery involved -- pure SMILES parsing, bracket vs
    /// organic notation. `chematic-rxn`'s own test suite
    /// (`transform.rs::reaction_derived_matches_direct_parse_chlorobenzene`)
    /// covers the actual `run_reactants`-derived case, since chematic-smiles
    /// cannot depend on chematic-rxn.
    #[test]
    fn bracket_vs_organic_notation_isolated_no_reaction() {
        let cases = [
            ("[Cl]c1ccccc1", "Clc1ccccc1"), // aromatic, symmetric ring
            ("[Cl]CCC", "ClCCC"),           // aliphatic, asymmetric
            ("[OH]CC", "OCC"),              // hydroxyl group
            ("C[NH2]", "CN"),               // amine
            ("[CH3]C", "CC"),               // methyl written explicit
        ];
        for (bracket, organic) in cases {
            let a = parse(bracket).unwrap();
            let b = parse(organic).unwrap();
            assert_same_graph_then_same_canonical(&a, &b, &format!("{bracket} vs {organic}"));
        }
    }

    /// Randomized atom relabeling: build the same graph via `MoleculeBuilder`
    /// with atoms inserted in several different orders (including a
    /// deterministic pseudo-random shuffle), holding chemistry fixed. Proves
    /// insertion order alone -- independent of bracket/organic notation --
    /// cannot change `canonical_smiles`'s output.
    #[test]
    fn randomized_atom_relabeling_does_not_change_canonical_form() {
        // A 5-membered ring with a substituent breaking full symmetry:
        // cyclopentane with one CH2 replaced conceptually by using a
        // pendant Cl so every construction order is exercised meaningfully.
        fn ring_with_substituent(insertion_order: &[usize]) -> Molecule {
            // Logical atoms: 0..5 ring carbons, 5 = Cl substituent on ring atom 0.
            let atoms: Vec<Atom> = (0..5)
                .map(|_| Atom {
                    element: Element::C,
                    isotope: None,
                    charge: 0,
                    hydrogen_count: None,
                    aromatic: false,
                    chirality: Chirality::None,
                    wildcard: false,
                    atom_map: None,
                    cip_code: None,
                })
                .chain(std::iter::once(Atom {
                    element: Element::CL,
                    isotope: None,
                    charge: 0,
                    hydrogen_count: None,
                    aromatic: false,
                    chirality: Chirality::None,
                    wildcard: false,
                    atom_map: None,
                    cip_code: None,
                }))
                .collect();
            // Logical bonds: ring 0-1-2-3-4-0, plus 0-5 (Cl).
            let bonds = [
                (0usize, 1usize, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
                (4, 0, BondOrder::Single),
                (0, 5, BondOrder::Single),
            ];

            let mut logical_to_new = [0u32; 6];
            let mut b = MoleculeBuilder::new();
            for (new_idx, &logical_idx) in insertion_order.iter().enumerate() {
                logical_to_new[logical_idx] = new_idx as u32;
                b.add_atom(atoms[logical_idx].clone());
            }
            for (l1, l2, order) in bonds {
                let a1 = AtomIdx(logical_to_new[l1]);
                let a2 = AtomIdx(logical_to_new[l2]);
                let _ = b.add_bond(a1, a2, order);
            }
            b.build()
        }

        let orders: Vec<Vec<usize>> = vec![
            vec![0, 1, 2, 3, 4, 5],
            vec![5, 4, 3, 2, 1, 0],
            vec![2, 0, 4, 1, 5, 3],
            vec![5, 0, 1, 2, 3, 4],
            vec![3, 1, 5, 0, 4, 2],
        ];
        let reference = canonical_smiles(&ring_with_substituent(&orders[0]));
        for order in &orders[1..] {
            let mol = ring_with_substituent(order);
            assert_eq!(
                canonical_smiles(&mol),
                reference,
                "insertion order {order:?} must not change canonical_smiles"
            );
        }
    }

    /// Disconnected structures: a salt/mixture where one fragment uses
    /// bracket notation and the other uses organic-subset notation for
    /// what is, on the relevant atom, the same effective H count.
    #[test]
    fn disconnected_structures_bracket_organic_mix() {
        let a = parse("[Cl]CC.[Na+]").unwrap();
        let b = parse("ClCC.[Na+]").unwrap();
        assert_same_graph_then_same_canonical(&a, &b, "disconnected salt, bracket vs organic Cl");
    }

    /// A Kekulized ring (explicit alternating single/double bonds, no
    /// aromaticity perception applied), combined with a bracket/organic
    /// notation difference on its substituent, must still canonicalize
    /// identically. Kept separate from the aromatic-flagged case in
    /// `bracket_vs_organic_notation_isolated_no_reaction` -- raw `parse()`
    /// does not itself perceive aromaticity from Kekulized input (that is a
    /// separate, opt-in `chematic-perception` step), so an aromatic-flagged
    /// molecule and a Kekulized one are genuinely different `Atom.aromatic`
    /// states, not merely different spellings of the same data; comparing
    /// them here would test aromaticity perception, not this fix.
    #[test]
    fn kekulized_ring_bracket_vs_organic_substituent() {
        let bracket = parse("[Cl]C1=CC=CC=C1").unwrap();
        let organic = parse("ClC1=CC=CC=C1").unwrap();
        assert_same_graph_then_same_canonical(
            &bracket,
            &organic,
            "Kekulized ring, bracket vs organic Cl substituent",
        );
    }

    /// Isotope and charge must NOT be unified away -- only redundant
    /// explicit H-count information is. An explicit isotope or nonzero
    /// charge is always genuinely distinguishing (there is no "implicit
    /// isotope"/"implicit charge" to fall back to), so these must still
    /// force bracket notation and distinct canonical forms from their
    /// natural-abundance/neutral counterparts.
    #[test]
    fn isotope_and_charge_remain_distinguishing_not_unified() {
        let normal = parse("CCl").unwrap();
        let isotopic = parse("C[37Cl]").unwrap();
        assert_ne!(
            canonical_smiles(&normal),
            canonical_smiles(&isotopic),
            "an explicit isotope must remain distinguishing"
        );

        let neutral = parse("CC(=O)[O-]").unwrap();
        let anion_str = canonical_smiles(&neutral);
        assert!(
            anion_str.contains('-'),
            "explicit negative charge must remain visible: {anion_str}"
        );
    }

    /// A real stereocenter (explicit chirality) is unaffected by the
    /// explicit/implicit-H unification: canonicalizing a chiral molecule
    /// from two different atom orderings must still agree, and the
    /// canonical form must still be stable under repeated canonicalization
    /// (idempotence).
    #[test]
    fn stereocenter_canonicalization_idempotent_and_order_independent() {
        let a = parse("N[C@@H](C)C(=O)O").unwrap(); // L-alanine
        let b = parse("OC(=O)[C@H](C)N").unwrap(); // same molecule, other end first
        assert_same_graph_then_same_canonical(&a, &b, "L-alanine, two parse orders");

        let canon = canonical_smiles(&a);
        let reparsed = parse(&canon).unwrap();
        assert_eq!(
            canonical_smiles(&reparsed),
            canon,
            "canonical_smiles must be idempotent: canonicalize(parse(canonicalize(x))) == canonicalize(x)"
        );
    }

    /// General idempotence check across the bracket/organic pairs above:
    /// canonicalizing an already-canonical string must reproduce it
    /// byte-for-byte, regardless of which side of a bracket/organic pair it
    /// started from.
    #[test]
    fn canonicalization_is_idempotent_for_all_bracket_organic_pairs() {
        for (bracket, organic) in [
            ("[Cl]c1ccccc1", "Clc1ccccc1"),
            ("[Cl]CCC", "ClCCC"),
            ("[OH]CC", "OCC"),
        ] {
            for smi in [bracket, organic] {
                let mol = parse(smi).unwrap();
                let canon = canonical_smiles(&mol);
                let reparsed = parse(&canon).unwrap();
                assert_eq!(
                    canonical_smiles(&reparsed),
                    canon,
                    "not idempotent starting from {smi}: {canon}"
                );
            }
        }
    }

    // ── Round 11: ring-closure explicit-bond-order aromaticity check (#395) ──
    //
    // `write_chain`'s ring-closure marker decision checked only the
    // currently-written atom's own aromaticity, never its ring-closure
    // partner's -- unlike the equivalent tree-edge decision, which correctly
    // checks both endpoints. A bare ring-closure digit between two aromatic
    // atoms is read back as an *aromatic* bond, so a genuinely `Single`
    // ring-closure bond joining two atoms that each individually happen to
    // be aromatic (e.g. a non-aromatic fusion bond between two separately-
    // aromatic ring systems, `c1-2`) silently became aromatic on re-parse
    // whenever the writer omitted the `-` marker. Confirmed via corpus sweep
    // (130/10,000 molecules) and RDKit InChI cross-check.

    #[test]
    fn ring_closure_explicit_single_bond_between_aromatic_atoms_real_world_repro() {
        // Smallest real-corpus repro from issue #395: the final atom closes
        // both ring 1 (implicit aromatic, correct) and ring 2 (explicit
        // `-2`, a non-aromatic fusion bond) at once.
        let smi = "Oc1[nH]c(Br)nc2nnc(Br)c1-2";
        let mol = parse(smi).unwrap();
        let once = canonical_smiles(&mol);
        let reparsed = parse(&once).unwrap();
        let twice = canonical_smiles(&reparsed);
        assert_eq!(
            once, twice,
            "explicit-bond-order ring closure not idempotent: {smi} -> once={once} twice={twice}"
        );
    }

    #[test]
    fn ring_closure_explicit_bond_order_survives_reparse_for_various_orders() {
        // A minimal two-ring system where the fusion bond (ring digit `2`)
        // is deliberately non-aromatic, exercised at each explicit bond
        // order the writer can emit on a ring closure.
        for smi in [
            "c1ccc2c1-c1ccccc-21", // single fusion bond
            "c1ccc2c1=CC=CC2",     // double fusion bond into a non-aromatic ring
        ] {
            let mol = match parse(smi) {
                Ok(m) => m,
                Err(_) => continue, // not every hand-written combination is valid; skip malformed ones
            };
            let once = canonical_smiles(&mol);
            let reparsed = parse(&once).unwrap_or_else(|e| {
                panic!("re-parse of canonical output '{once}' (from {smi}) failed: {e}")
            });
            let twice = canonical_smiles(&reparsed);
            assert_eq!(
                once, twice,
                "not idempotent: {smi} -> once={once} twice={twice}"
            );
        }
    }
}
