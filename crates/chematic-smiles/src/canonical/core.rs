/// Return the atom indices sorted into canonical (Morgan-rank) order.
///
/// The returned `Vec<usize>` lists atom positions (0-based) in the order they
/// would be encountered during a canonical DFS write. Atoms with higher
/// Morgan rank appear earlier. This is the same ordering `canonical_smiles`
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
    // Memoized on `mol`: canonical SMILES, canonical atom order and every
    // InChI layer derive from this same (deterministic) search.
    let shared = mol.derived(chematic_core::DerivedSlot::CanonicalRanks, || {
        winning_individualized_ranks_uncached(mol)
    });
    (shared.0.clone(), shared.1.clone())
}

fn winning_individualized_ranks_uncached(mol: &Molecule) -> (Vec<u64>, String) {
    let mut best = winning_individualized_ranks_single(mol);
    for alternate in alternate_ez_carrier_spellings(mol) {
        let candidate = winning_individualized_ranks_single(&alternate);
        if candidate.1 < best.1 {
            best = candidate;
        }
    }
    best
}

/// Return geometrically equivalent spellings that move one plain E/Z marker
/// to the sibling substituent at the same alkene end.  SMILES permits either
/// substituent to carry the marker; choosing the lexicographically smallest
/// result across both spellings removes the last input-carrier dependency in
/// cyclic conjugated systems.
///
/// A marked bond is not moved when its far endpoint is the sole carrier for a
/// different stereogenic double bond.  That conservative guard prevents this
/// normalization from stripping unrelated E/Z information.  Coupled ends
/// remain handled by the writer's joint resolver.
fn alternate_ez_carrier_spellings(mol: &Molecule) -> Vec<Molecule> {
    let ends = CanonicalWriter::compute_stereo_alkene_ends(mol);
    // Multi-end coupled systems have additional global carrier constraints;
    // leave those to the writer's joint resolver until a complete equivalent
    // spelling proof exists for that larger state space. Without a stereo
    // alkene end there is nothing to respell (and no ranks to compute).
    if ends.is_empty() || ends.len() > 2 {
        return Vec::new();
    }
    let ranks = morgan_ranks(mol);
    let writer = CanonicalWriter::new(mol, &ranks);
    let mut alternates = Vec::new();
    for &end in &ends {
        let subs = CanonicalWriter::substituents(mol, end);
        if subs.len() != 2 {
            continue;
        }
        for i in 0..2 {
            let marked = subs[i];
            let sibling = subs[1 - i];
            let Some(direction) = writer.raw_input_direction(marked.1) else {
                continue;
            };
            // A bond leading to another stereo end is that end's candidate
            // too (the two are coupled, left to the writer's joint
            // resolver): a marker moved onto it also states that end's
            // geometry, and one moved off it can take that end's only
            // marker. Squaraine diimines in Kekulé form came out with one
            // imine inverted or unspecified (`C/N=C1/C(O)=C(O)\C1=N/C`,
            // #734).
            if !matches!(direction, BondOrder::Up | BondOrder::Down)
                || mol.bond(sibling.1).order != BondOrder::Single
                || ends.contains(&sibling.0)
                || ends.contains(&marked.0)
                || writer.is_load_bearing_elsewhere(marked.1, end)
            {
                continue;
            }
            let up = CanonicalWriter::direction_is_up(
                direction,
                writer.raw_direction_anchor(marked.1),
                end,
            );
            let sibling_direction =
                CanonicalWriter::direction_for_up(writer.raw_direction_anchor(sibling.1), end, !up);
            alternates.push(
                mol.with_bond_order(marked.1, BondOrder::Single)
                    .with_bond_order(sibling.1, sibling_direction),
            );
        }
    }
    alternates
}

fn winning_individualized_ranks_single(mol: &Molecule) -> (Vec<u64>, String) {
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
    let oriented = atom1_relative_stashes(mol);
    let mol = oriented.as_ref().unwrap_or(mol);
    let cleaned = without_lone_markers(mol);
    let mol = cleaned.as_ref().unwrap_or(mol);

    let (_, winning_string) = winning_individualized_ranks(mol);
    winning_string
}

/// `mol` with every stashed direction ([`Molecule::bond_direction`]) read
/// from the bond's atom1, or `None` when none is read from atom2. A stash on
/// a plain single bond becomes its literal `Up`/`Down` order; one on an
/// aromatic (or wedge) bond keeps its stash with atom1 as anchor. Both are
/// the same geometry; the writer's carrier resolution produces atom1-relative
/// markers and its emission reads them against the anchor, so a stash read
/// from atom2 came out mirrored (92 of 272 E/Z molecules with their markers
/// moved into atom2-anchored stashes changed geometry, #734).
fn atom1_relative_stashes(mol: &Molecule) -> Option<Molecule> {
    let flip = |d: BondOrder| match d {
        BondOrder::Up => BondOrder::Down,
        BondOrder::Down => BondOrder::Up,
        other => other,
    };
    let todo: Vec<(BondIdx, BondOrder, bool)> = mol
        .bonds()
        .filter_map(|(bidx, bond)| {
            let direction = mol.bond_direction(bidx)?;
            let anchor = mol.bond_direction_anchor(bidx).unwrap_or(bond.atom1);
            let direction = if anchor == bond.atom1 {
                direction
            } else {
                flip(direction)
            };
            let to_literal = bond.order == BondOrder::Single;
            (to_literal || anchor != bond.atom1).then_some((bidx, direction, to_literal))
        })
        .collect();
    if todo.is_empty() {
        return None;
    }
    let mut out = mol.clone();
    for (bidx, direction, to_literal) in todo {
        let atom1 = out.bond(bidx).atom1;
        out.clear_bond_direction(bidx);
        if to_literal {
            out.set_bond_order(bidx, direction);
        } else {
            out.set_bond_direction(bidx, direction);
            out.set_bond_direction_anchor(bidx, atom1);
        }
    }
    Some(out)
}

/// `mol` without `/`/`\\` markers that take part in no E/Z specification
/// (a marker flanking only double bonds whose other end has no marker, or
/// that sit in rings of fewer than eight atoms), or `None` when it has
/// none. RDKit drops such markers when parsing; removing them first makes
/// a spelling with and without them canonicalize alike. A literal
/// `Up`/`Down` order becomes `Single`.
fn without_lone_markers(mol: &Molecule) -> Option<Molecule> {
    let mut out: Option<Molecule> = None;
    loop {
        let current = out.as_ref().unwrap_or(mol);
        let lone: Vec<BondIdx> = current
            .bonds()
            .filter(|&(bidx, _)| crate::writer::raw_bond_direction(current, bidx).is_some())
            .filter(|&(bidx, _)| !crate::writer::flanks_marked_double_bond(current, bidx))
            .map(|(bidx, _)| bidx)
            .collect();
        if lone.is_empty() {
            return out;
        }
        let next = out.get_or_insert_with(|| mol.clone());
        for bidx in lone {
            next.clear_bond_direction(bidx);
            if matches!(next.bond(bidx).order, BondOrder::Up | BondOrder::Down) {
                next.set_bond_order(bidx, BondOrder::Single);
            }
        }
    }
}

/// [`canonical_smiles`] together with the order in which the atoms appear in
/// the string.
///
/// The string is exactly [`canonical_smiles`]'s. `order[k]` is the index in
/// `mol` of the `k`-th atom written, which is also the index that atom
/// receives when the string is parsed back with [`crate::parse`]. Every atom
/// is written exactly once, so `order` is a permutation of
/// `0..mol.atom_count()`. This is the canonical DFS visit order, which is not
/// the same permutation as [`canonical_atom_order`] (a rank sort).
pub fn canonical_smiles_with_atom_order(mol: &Molecule) -> (String, Vec<AtomIdx>) {
    if mol.atom_count() == 0 {
        return (String::new(), Vec::new());
    }
    let oriented = atom1_relative_stashes(mol);
    let mol = oriented.as_ref().unwrap_or(mol);
    let cleaned = without_lone_markers(mol);
    let mol = cleaned.as_ref().unwrap_or(mol);
    let (ranks, winning_string) = winning_individualized_ranks(mol);
    // The winner may have been written from an equivalent E/Z carrier
    // spelling of `mol` (same atoms and bonds, different directional
    // markers); the traversal only depends on the ranks and the graph.
    let (_, order) = CanonicalWriter::new(mol, &ranks).write_all_with_order();
    (winning_string, order)
}

/// Return a canonical SMILES only when its representation is self-stable.
///
/// Canonical E/Z carrier placement remains fail-closed for coupled systems
/// except the bounded aromatic-direction-stash planner: that planner searches
/// every eligible token slot, reparses each candidate, and selects a semantic
/// lexicographic minimum. This helper still reparses and requires idempotence
/// before admitting that narrow, proven path as a deduplication or cache key.
/// Enhanced stereo groups are refused because ordinary SMILES omits their
/// membership and relative/mixture semantics.
pub fn canonical_smiles_stable_key(mol: &Molecule) -> Option<String> {
    if !mol.stereo_groups().is_empty() {
        return None;
    }
    let candidate = canonical_smiles(mol);
    let reparsed = crate::parser::parse(&candidate).ok()?;
    if candidate != canonical_smiles(&reparsed) {
        return None;
    }

    let has_direction = |bidx: BondIdx| {
        matches!(reparsed.bond(bidx).order, BondOrder::Up | BondOrder::Down)
            || reparsed.bond_direction(bidx).is_some()
    };
    let ez_double_bonds = reparsed
        .bonds()
        .filter(|(_, bond)| bond.order == BondOrder::Double)
        .filter(|(double_idx, bond)| {
            [bond.atom1, bond.atom2].into_iter().all(|end| {
                reparsed
                    .neighbors(end)
                    .any(|(_, side_bond)| side_bond != *double_idx && has_direction(side_bond))
            })
        })
        .count();
    let has_aromatic_direction_stash = mol.bonds().any(|(bond, value)| {
        value.order == BondOrder::Aromatic && mol.bond_direction(bond).is_some()
    });
    (ez_double_bonds <= 1 || has_aromatic_direction_stash).then_some(candidate)
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
    let mut next_ranks = vec![0; n];
    let mut distinct_scratch = Vec::with_capacity(n);
    // Carry the previous distinct-count forward.  Recomputing it from the
    // unchanged `ranks` at the top of every iteration needlessly cloned and
    // sorted the same vector a second time.
    let mut old_distinct = count_distinct(&ranks, &mut distinct_scratch);
    for _ in 0..max_iter {
        for (i, next) in next_ranks.iter_mut().enumerate() {
            let idx = AtomIdx(i as u32);
            // Include bond order in the neighbor contribution so that atoms
            // bonded via different bond types (e.g. O= vs O-H in acetic acid)
            // receive distinct Morgan ranks even when neighbor atom ranks are
            // otherwise identical.
            // Molecular graphs overwhelmingly have degree <= 8. Keep
            // the per-atom refinement scratch space on the stack for
            // that common case; fall back to a Vec only for unusual
            // highly branched atoms.
            let mut inline = [0_u64; 8];
            let mut inline_len = 0usize;
            let mut overflow: Option<Vec<u64>> = None;
            for (nb, bidx) in mol.neighbors(idx) {
                let bond_val = bond_order_value(mol.bond(bidx).order);
                let contribution = fnv_hash_sequence(ranks[nb.0 as usize], &[bond_val]);
                if let Some(values) = overflow.as_mut() {
                    values.push(contribution);
                } else if inline_len < inline.len() {
                    inline[inline_len] = contribution;
                    inline_len += 1;
                } else {
                    let mut values = Vec::with_capacity(mol.degree(idx));
                    values.extend_from_slice(&inline);
                    values.push(contribution);
                    overflow = Some(values);
                }
            }
            *next = match overflow {
                Some(mut values) => {
                    values.sort_unstable();
                    fnv_hash_sequence(ranks[i], &values)
                }
                None => {
                    inline[..inline_len].sort_unstable();
                    fnv_hash_sequence(ranks[i], &inline[..inline_len])
                }
            };
        }

        let new_distinct = count_distinct(&next_ranks, &mut distinct_scratch);
        std::mem::swap(&mut ranks, &mut next_ranks);

        if new_distinct <= old_distinct {
            break;
        }
        old_distinct = new_distinct;
    }

    // The final distinct-count pass already sorted/deduplicated these exact
    // ranks. Reuse its dictionary instead of sorting an index vector again.
    for rank in &mut ranks {
        *rank = distinct_scratch
            .binary_search(rank)
            .expect("rank in dictionary") as u64;
    }
    ranks
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
/// `refine_ranks`).
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
        // Charge and isotope still tell dummy atoms apart (`[*-]`, `[13*]`).
        // A plain `[*]` stays 0.
        return (atom.isotope.unwrap_or(0) as u64) << 8 | (atom.charge as u8) as u64;
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
pub(crate) fn bond_order_value(order: BondOrder) -> u64 {
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

fn count_distinct(ranks: &[u64], scratch: &mut Vec<u64>) -> usize {
    scratch.clear();
    scratch.extend_from_slice(ranks);
    scratch.sort_unstable();
    scratch.dedup();
    scratch.len()
}

#[derive(Clone)]
pub(crate) struct CanonicalWriter<'a> {
    mol: &'a Molecule,
    ranks: &'a [u64],
    written: Vec<bool>,
    ring_bonds: Vec<bool>,
    /// (ring_num, is_open_side, ring_partner_atom, physical_bond). Deliberately
    /// does NOT store a precomputed `BondOrder`: which endpoint of `bidx` is
    /// "open" vs "close" is fixed at discovery time (`dfs_mark`), but the
    /// actual `/`/`\` token can only be correctly computed once written --
    /// see `normalize_ez`'s doc comment for why baking in a write-direction
    /// re-orientation this early corrupts E/Z groups spanning bonds visited
    /// in different directions (issue #390).
    atom_ring_nums: Vec<Vec<(u32, bool, AtomIdx, BondIdx)>>,
    /// Selected directional ring bonds whose marker is emitted at the closing
    /// occurrence. The bounded aromatic E/Z planner varies this together with
    /// carrier polarity because either legal digit occurrence can encode the
    /// same bond direction.
    ring_marker_on_close: HashSet<BondIdx>,
    /// Double bonds checked by [`Self::suppress_lone_marker`] for lying in
    /// a ring of fewer than eight atoms.
    small_ring_double: std::cell::RefCell<HashMap<BondIdx, bool>>,
    next_ring: u32,
    out: String,
    /// Atoms in the order `write_chain` emitted them.
    order: Vec<AtomIdx>,
    /// Union-find groups of directional (`/`/`\`) bonds that jointly encode
    /// one connected E/Z system — flipping every member preserves geometry,
    /// flipping a subset does not. Keyed/rooted by `BondIdx`.
    ez_group: HashMap<BondIdx, BondIdx>,
    /// Groups whose first-encountered bond (in write order) came out `Down`;
    /// every remaining bond in the group is flipped so the first directional
    /// bond of each system is always `/`, regardless of input spelling.
    ez_flip: HashMap<BondIdx, bool>,
    /// Optional component-level flip choices used when comparing equivalent
    /// aromatic-stash spellings.  A connected E/Z component may be inverted
    /// as a whole without changing its geometry; choosing the lexicographically
    /// smallest serialized candidate removes input-carrier dependence.
    forced_ez_flip: Option<HashMap<BondIdx, bool>>,
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
    /// Geometry facts captured before E/Z marker carrier resolution.  The
    /// component solver reads these facts instead of re-reading whichever
    /// raw marker location happened to occur in the input spelling.
    ez_geometry_facts: Vec<EzGeometryFact>,
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
    /// Test-only DFS edge preference used to determine whether the remaining
    /// issue #503 residual has a geometry-preserving common spelling when a
    /// whole coupled component is allowed to influence traversal.  Production
    /// still uses rank-only traversal until this bounded search has a proven
    /// acceptance rule.
    #[cfg(test)]
    traversal_bond_preference: Option<HashMap<BondIdx, bool>>,
    /// Complete directional-token assignment selected by the bounded aromatic
    /// E/Z planner. Unlike `ez_marker`, this deliberately bypasses local
    /// carrier election and component-global normalization: every eligible
    /// token slot is assigned explicitly, then the serialized candidate is
    /// reparsed and accepted only if it preserves rank-keyed E/Z facts.
    direct_ez_orders: Option<HashMap<BondIdx, BondOrder>>,
    /// Parent edge for the canonical DFS skeleton. Besides test diagnostics,
    /// allene parity needs the exact terminal-ligand order that serialization
    /// will use.
    canonical_tree_parent: Vec<Option<(AtomIdx, BondIdx)>>,
}

/// One input-representation-independent E/Z relation, retained before the
/// writer decides where `/` or `\\` tokens will appear.  `AtomIdx` and
/// `BondIdx` locate the fact within the current molecule only; any future
/// canonical choice must use the rank keys, never those parse-order indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct EzGeometryFact {
    double_bond: BondIdx,
    double_key: (u64, u64),
    ends: [EzGeometryEnd; 2],
    same_side: bool,
}

/// The fixed, rank-selected reference substituent at one double-bond end and
/// its encoded side of the alkene axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct EzGeometryEnd {
    atom: AtomIdx,
    reference_bond: BondIdx,
    reference_rank: u64,
    up: bool,
}

/// A parse-order-independent E/Z signature used to validate a candidate
/// serialization without invoking canonicalization again.  It deliberately
/// uses Morgan ranks only: full individualized canonical ranks would recurse
/// through the writer while the bounded planner is already serializing.
type EzSemanticSignature = Vec<((u64, u64), (u64, u64), bool)>;

/// Return whether every doubly-marked alkene end carries one marker on each
/// side of the double-bond axis.
///
/// SMILES permits either substituent bond at an alkene end to carry the
/// directional token.  When both substituent bonds are marked, however, the
/// two tokens must describe opposite sides.  Keeping two same-side markers
/// is internally self-consistent if a reader simply picks one of them, but
/// independent readers such as RDKit correctly treat it as a conflicting
/// specification and discard the affected E/Z assignment.  The complete
/// aromatic planner enumerates explicit token placements, so reject that
/// syntax before considering a candidate semantically admissible.
fn ez_direction_assignments_consistent(mol: &Molecule) -> bool {
    mol.bonds()
        .filter(|(_, bond)| bond.order == BondOrder::Double)
        .all(|(_, bond)| {
            [bond.atom1, bond.atom2].into_iter().all(|end| {
                let marked_sides: Vec<bool> = mol
                    .neighbors(end)
                    .filter(|(_, adjacent)| mol.bond(*adjacent).order != BondOrder::Double)
                    .filter_map(|(_, adjacent)| {
                        crate::writer::raw_bond_direction(mol, adjacent).map(|direction| {
                            let anchor = mol
                                .bond_direction_anchor(adjacent)
                                .unwrap_or(mol.bond(adjacent).atom1);
                            CanonicalWriter::direction_is_up(direction, anchor, end)
                        })
                    })
                    .collect();
                match marked_sides.as_slice() {
                    [] | [_] => true,
                    [first, second] => first != second,
                    _ => false,
                }
            })
        })
}

fn ez_semantic_signature(mol: &Molecule) -> EzSemanticSignature {
    let ranks = morgan_ranks(mol);
    let writer = CanonicalWriter::new(mol, &ranks);
    writer
        .extract_ez_geometry_facts()
        .into_iter()
        .map(|fact| {
            let mut references = [fact.ends[0].reference_rank, fact.ends[1].reference_rank];
            references.sort_unstable();
            (
                fact.double_key,
                (references[0], references[1]),
                fact.same_side,
            )
        })
        .collect()
}
