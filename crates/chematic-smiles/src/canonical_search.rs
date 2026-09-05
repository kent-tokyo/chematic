//! Streaming, automorphism-orbit-pruned canonical search.
//!
//! Replaces full branch enumeration (`crate::canonical::
//! enumerate_discrete_ranks`, still kept for the test-only exhaustive
//! oracle and as a last-resort fallback -- see below) with a DFS that:
//! 1. Never explores more than one representative per PROVEN automorphism
//!    orbit within a target cell (`canonical_automorphism::
//!    has_colored_automorphism_mapping` decides "proven").
//! 2. Never collects the full leaf set into memory -- it streams leaves
//!    into a single running incumbent (lexicographically smallest string).
//!
//! Every leaf's `ranks` vector is produced by the exact same
//! `crate::canonical::individualize` + `crate::canonical::refine_ranks`
//! primitives the legacy enumeration uses. Pruning only ever *skips*
//! sibling branches; it never perturbs the rank vector of a branch it does
//! explore. So any surviving branch is byte-for-byte identical to what
//! unpruned enumeration would have produced for that same individualization
//! sequence -- see `docs/rfcs/canonical_automorphism_pruning.md`, "why pruning
//! cannot change surviving output", for the full argument and its
//! `canonical_smiles_exhaustive_oracle` cross-check.

use chematic_core::{AtomIdx, Molecule};
use smallvec::{SmallVec, smallvec};

use crate::canonical::{CanonicalWriter, individualize, refine_ranks};
use crate::canonical_automorphism::has_colored_automorphism_mapping;
use crate::canonical_partition::{
    CanonicalColoredGraph, Partition, exact_refine, initial_partition,
};

mod stats {
    //! Feature-gated internal work counters, following the exact pattern
    //! established by `chematic-rxn`'s `perf_counters.rs`
    //! (`docs/rfcs/reaction_transform_perf.md`): process-global `AtomicUsize`s,
    //! zero cost (every bump compiles to an empty inline fn) unless the
    //! `canonical-search-instrumentation` feature is enabled.
    #[cfg(feature = "canonical-search-instrumentation")]
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[cfg(feature = "canonical-search-instrumentation")]
    static NODES_VISITED: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static LEAVES_WRITTEN: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static ORBIT_TESTS: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static ORBIT_UNIONS: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static CHILDREN_PRUNED: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static MAX_DEPTH: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static LARGEST_TARGET_CELL: AtomicUsize = AtomicUsize::new(0);
    #[cfg(feature = "canonical-search-instrumentation")]
    static BUDGET_EXHAUSTIONS: AtomicUsize = AtomicUsize::new(0);

    /// Snapshot of every counter. All-zero when the feature is disabled.
    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    pub struct CanonicalSearchStats {
        pub nodes_visited: usize,
        pub leaves_written: usize,
        pub orbit_tests: usize,
        pub orbit_unions: usize,
        pub children_pruned: usize,
        pub max_depth: usize,
        pub largest_target_cell: usize,
        pub budget_exhaustions: usize,
    }

    #[cfg(feature = "canonical-search-instrumentation")]
    pub fn snapshot() -> CanonicalSearchStats {
        CanonicalSearchStats {
            nodes_visited: NODES_VISITED.load(Ordering::Relaxed),
            leaves_written: LEAVES_WRITTEN.load(Ordering::Relaxed),
            orbit_tests: ORBIT_TESTS.load(Ordering::Relaxed),
            orbit_unions: ORBIT_UNIONS.load(Ordering::Relaxed),
            children_pruned: CHILDREN_PRUNED.load(Ordering::Relaxed),
            max_depth: MAX_DEPTH.load(Ordering::Relaxed),
            largest_target_cell: LARGEST_TARGET_CELL.load(Ordering::Relaxed),
            budget_exhaustions: BUDGET_EXHAUSTIONS.load(Ordering::Relaxed),
        }
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    pub fn snapshot() -> CanonicalSearchStats {
        CanonicalSearchStats::default()
    }

    #[cfg(feature = "canonical-search-instrumentation")]
    pub fn reset() {
        NODES_VISITED.store(0, Ordering::Relaxed);
        LEAVES_WRITTEN.store(0, Ordering::Relaxed);
        ORBIT_TESTS.store(0, Ordering::Relaxed);
        ORBIT_UNIONS.store(0, Ordering::Relaxed);
        CHILDREN_PRUNED.store(0, Ordering::Relaxed);
        MAX_DEPTH.store(0, Ordering::Relaxed);
        LARGEST_TARGET_CELL.store(0, Ordering::Relaxed);
        BUDGET_EXHAUSTIONS.store(0, Ordering::Relaxed);
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    pub fn reset() {}

    #[cfg(feature = "canonical-search-instrumentation")]
    pub(crate) fn record_node(depth: usize) {
        NODES_VISITED.fetch_add(1, Ordering::Relaxed);
        MAX_DEPTH.fetch_max(depth, Ordering::Relaxed);
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    #[inline(always)]
    pub(crate) fn record_node(_depth: usize) {}

    #[cfg(feature = "canonical-search-instrumentation")]
    pub(crate) fn record_leaf() {
        LEAVES_WRITTEN.fetch_add(1, Ordering::Relaxed);
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    #[inline(always)]
    pub(crate) fn record_leaf() {}

    #[cfg(feature = "canonical-search-instrumentation")]
    pub(crate) fn record_orbit_test(succeeded: bool) {
        ORBIT_TESTS.fetch_add(1, Ordering::Relaxed);
        if succeeded {
            ORBIT_UNIONS.fetch_add(1, Ordering::Relaxed);
        }
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    #[inline(always)]
    pub(crate) fn record_orbit_test(_succeeded: bool) {}

    #[cfg(feature = "canonical-search-instrumentation")]
    pub(crate) fn record_target_cell(cell_size: usize, representatives: usize) {
        LARGEST_TARGET_CELL.fetch_max(cell_size, Ordering::Relaxed);
        CHILDREN_PRUNED.fetch_add(cell_size.saturating_sub(representatives), Ordering::Relaxed);
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    #[inline(always)]
    pub(crate) fn record_target_cell(_cell_size: usize, _representatives: usize) {}

    #[cfg(feature = "canonical-search-instrumentation")]
    pub(crate) fn record_budget_exhaustion() {
        BUDGET_EXHAUSTIONS.fetch_add(1, Ordering::Relaxed);
    }
    #[cfg(not(feature = "canonical-search-instrumentation"))]
    #[inline(always)]
    pub(crate) fn record_budget_exhaustion() {}
}

pub use stats::CanonicalSearchStats;
pub use stats::{reset as reset_search_stats, snapshot as search_stats_snapshot};

/// Budget for the fallible, bounded canonicalization API
/// ([`canonical_smiles_with_limits`]). `None` on either field means "no cap
/// on that dimension".
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CanonicalizationLimits {
    pub max_search_nodes: Option<usize>,
    pub max_automorphism_tests: Option<usize>,
}

impl CanonicalizationLimits {
    /// No cap on either dimension. This is the policy the existing
    /// infallible `canonical_smiles`/`canonical_atom_order` use (see their
    /// doc comments): with automorphism-orbit pruning removing the
    /// combinatorial blowup that made the old fixed 10,000-branch cap
    /// necessary, an unbounded search over the (now much smaller) set of
    /// *necessary* branches is safe in practice -- termination is still
    /// mathematically guaranteed (individualize-refine recursion depth is
    /// bounded by atom count regardless of branching factor), just not
    /// bounded in wall-clock for a theoretical adversarial worst case that
    /// has not been observed on this project's fixtures or 5,000-molecule
    /// corpus.
    pub fn unbounded() -> Self {
        Self {
            max_search_nodes: None,
            max_automorphism_tests: None,
        }
    }
}

/// Errors from the bounded canonicalization API. Never silently substituted
/// for a correct answer: budget exhaustion is always `Err`, never a
/// partial-search result returned as if it were the true minimum.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CanonicalizationError {
    /// The configured `CanonicalizationLimits` were exceeded before the
    /// search could complete and prove its result correct.
    SearchBudgetExceeded {
        nodes_visited: usize,
        automorphism_tests: usize,
    },
    /// An internal invariant was violated (e.g. the search produced no leaf
    /// at all for a non-empty molecule). Should never happen; surfaced as a
    /// typed error rather than a panic or a silently-wrong string so a
    /// caller can detect and report it.
    InvalidInternalMapping { detail: String },
}

impl std::fmt::Display for CanonicalizationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SearchBudgetExceeded {
                nodes_visited,
                automorphism_tests,
            } => write!(
                f,
                "canonicalization search budget exceeded ({nodes_visited} nodes, {automorphism_tests} automorphism tests)"
            ),
            Self::InvalidInternalMapping { detail } => {
                write!(f, "internal canonicalization invariant violated: {detail}")
            }
        }
    }
}

impl std::error::Error for CanonicalizationError {}

struct Incumbent {
    ranks: Vec<u64>,
    string: String,
}

struct SearchNode {
    ranks: Vec<u64>,
    target_cell: Option<SmallVec<[usize; 8]>>,
    depth: usize,
}

#[derive(Default)]
struct SearchBudget {
    nodes_visited: usize,
    automorphism_tests: usize,
}

/// Lowest-ranked non-singleton cell, without constructing a `Vec` for every
/// rank class. Search only ever branches on this one cell.
fn first_non_singleton_cell(ranks: &[u64]) -> Option<SmallVec<[usize; 8]>> {
    let mut counts: SmallVec<[usize; 64]> = smallvec![0; ranks.len()];
    for &rank in ranks {
        counts[rank as usize] += 1;
    }
    let target = counts.iter().position(|&count| count > 1)?;
    let mut members = SmallVec::with_capacity(counts[target]);
    for (atom, &rank) in ranks.iter().enumerate() {
        if rank as usize == target {
            members.push(atom);
        }
    }
    Some(members)
}

/// Return the canonical SMILES for `mol`, respecting `limits`. `Err` on
/// budget exhaustion -- never a partial/truncated string.
pub fn canonical_smiles_with_limits(
    mol: &Molecule,
    limits: &CanonicalizationLimits,
) -> Result<String, CanonicalizationError> {
    if mol.atom_count() == 0 {
        return Ok(String::new());
    }
    let (_, s) = winning_individualized_ranks_with_limits(mol, limits)?;
    Ok(s)
}

/// Same search as [`canonical_smiles_with_limits`], but also returns the
/// winning fully-discrete rank vector (used by `canonical_atom_order`).
pub(crate) fn winning_individualized_ranks_with_limits(
    mol: &Molecule,
    limits: &CanonicalizationLimits,
) -> Result<(Vec<u64>, String), CanonicalizationError> {
    let ranks = crate::canonical::morgan_ranks(mol);

    // Most drug-like molecules have no Morgan-rank plateau after iterative
    // refinement.  In that case individualization cannot change the order,
    // so constructing the writer-visible colored graph and running the
    // automorphism machinery is redundant.  Keep the exact same writer and
    // rank vector, but skip that setup entirely.  This is a correctness
    // preserving fast path: branching is only required when two atoms share
    // a rank, and the existing exhaustive/orbit-pruned paths remain the
    // authority for tied molecules.
    let root_cell = first_non_singleton_cell(&ranks);
    if root_cell.is_none() {
        let string = crate::canonical::CanonicalWriter::new(mol, &ranks).write_all();
        return Ok((ranks, string));
    }

    let graph = CanonicalColoredGraph::new(mol);
    let mut budget = SearchBudget::default();
    let mut incumbent: Option<Incumbent> = None;
    search_canonical(
        mol,
        &graph,
        SearchNode {
            ranks,
            target_cell: root_cell,
            depth: 0,
        },
        limits,
        &mut budget,
        &mut incumbent,
    )?;
    match incumbent {
        Some(Incumbent { ranks, string }) => Ok((ranks, string)),
        None => Err(CanonicalizationError::InvalidInternalMapping {
            detail: "orbit-pruned search produced no leaf for a non-empty molecule".to_string(),
        }),
    }
}

fn search_canonical(
    mol: &Molecule,
    graph: &CanonicalColoredGraph,
    node: SearchNode,
    limits: &CanonicalizationLimits,
    budget: &mut SearchBudget,
    incumbent: &mut Option<Incumbent>,
) -> Result<(), CanonicalizationError> {
    let SearchNode {
        ranks,
        target_cell,
        depth,
    } = node;
    budget.nodes_visited += 1;
    stats::record_node(depth);
    if let Some(max) = limits.max_search_nodes
        && budget.nodes_visited > max
    {
        stats::record_budget_exhaustion();
        return Err(CanonicalizationError::SearchBudgetExceeded {
            nodes_visited: budget.nodes_visited,
            automorphism_tests: budget.automorphism_tests,
        });
    }

    // Same heuristic as the legacy enumeration (section 10: unchanged in
    // this PR): the lowest-ranked non-singleton cell.
    let Some(members) = target_cell.or_else(|| first_non_singleton_cell(&ranks)) else {
        let s = CanonicalWriter::new(mol, &ranks).write_all();
        stats::record_leaf();
        let is_better = incumbent.as_ref().is_none_or(|c| s < c.string);
        if is_better {
            *incumbent = Some(Incumbent { ranks, string: s });
        }
        return Ok(());
    };

    // Exact partition (full writer-visible chemical color, intersected with
    // the current individualization state carried by `ranks`), refined to a
    // fixpoint. Used only to decide which members of `members` are provably
    // in the same automorphism orbit -- NEVER used to alter `ranks` itself
    // (see module docs).
    let representatives = if members[1..].iter().all(|&other| {
        local_twins_without_partition(graph, AtomIdx(members[0] as u32), AtomIdx(other as u32))
    }) {
        // Every member is related to the first by an independently proven
        // fixed-point swap automorphism.  The entire target cell is
        // therefore one orbit, without needing exact partition refinement
        // or the general automorphism backtracker.
        vec![members[0]]
    } else {
        let partition = exact_refine(graph, initial_partition(graph, &ranks));
        exact_orbit_representatives(graph, &partition, &members, limits, budget)?
    };
    stats::record_target_cell(members.len(), representatives.len());

    for rep in representatives {
        let individualized = individualize(&ranks, rep);
        let re_refined = refine_ranks(mol, individualized);
        search_canonical(
            mol,
            graph,
            SearchNode {
                ranks: re_refined,
                target_cell: None,
                depth: depth + 1,
            },
            limits,
            budget,
            incumbent,
        )?;
    }
    Ok(())
}

/// Partition `members` (raw atom indices, ascending) into automorphism
/// orbits under the current node's coloring, and return one representative
/// per orbit -- the minimum-index member, ascending across orbits. This
/// exactly mirrors the legacy `enumerate_discrete_ranks`'s traversal order
/// (`for &atom_idx in members` in ascending-index order) restricted to
/// orbit representatives, so that whenever every relevant cell turns out to
/// be a genuine orbit, the very first leaf this search reaches is identical
/// to the first leaf the legacy exhaustive enumeration would have reached
/// (see the design doc for why this matters for tie-break determinism).
///
/// Union-Find; only merges pairs an actual verified automorphism test
/// proved equivalent (`has_colored_automorphism_mapping`). A false negative
/// (failing to merge a genuinely-automorphic pair) only costs performance;
/// a false positive is structurally impossible here **given `ranks` (and
/// therefore `coloring`, seeded from it in `initial_partition`) correctly
/// reflects every individualization already committed in an ancestor call**
/// -- every union really is gated on a real, independently-reverified
/// bijection over `graph`'s own color/edge data, which never merges two
/// atoms that differ in any writer-visible attribute.
///
/// That parenthetical is not free, though (independent Round-2 false-prune
/// audit, PR #193): `ranks` comes from `crate::canonical::individualize` +
/// `refine_ranks`. `individualize` itself is exact integer arithmetic --
/// the atom it distinguishes provably gets a rank no other atom in the
/// vector holds, zero collision risk. But `refine_ranks` immediately
/// re-hashes from that point (`fnv_hash_sequence`, pre-existing, unchanged
/// by this PR) and rank normalization groups by raw hash-value *equality* --
/// so a genuine 64-bit FNV-1a collision there could in principle re-merge
/// an already-individualized atom back in with a formerly-tied sibling,
/// and `coloring` would inherit that error (nothing downstream re-derives
/// "was this atom individualized" independently of `ranks` -- `VertexColor`
/// deliberately does not encode search-time individualization history, only
/// intrinsic atom attributes). This dependency is not new: `refine_ranks`
/// rank-equality is *already* the sole basis for the crate's pre-existing
/// `equivalent_atom_classes`/`are_atoms_equivalent` public APIs, with
/// identical collision exposure, unrelated to orbit pruning. What this PR
/// changes is the *consequence* of a hypothetical collision: in the legacy
/// exhaustive engine it would cause redundant-but-still-correct
/// over-exploration (every member of a wrongly-merged cell still gets
/// individualized and compared); here it could instead cause a genuinely
/// distinct branch to be silently skipped. Not observed on any fixture, the
/// n<=5 exhaustive suite, hundreds of randomized fuzz trials, or the
/// 5,000-molecule corpus; would require a correlated 64-bit hash collision
/// reconstructing an entire real symmetry's cell structure to manifest.
/// Judged not worth threading a parallel, hash-free individualization-state
/// vector through the hot path to close (a bigger change than this PR's
/// scope, defending against a risk already implicitly accepted crate-wide)
/// -- documented here and in `docs/rfcs/canonical_automorphism_pruning.md`
/// instead of silently left as an unqualified "structurally impossible"
/// claim.
fn exact_orbit_representatives(
    graph: &CanonicalColoredGraph,
    coloring: &Partition,
    members: &[usize],
    limits: &CanonicalizationLimits,
    budget: &mut SearchBudget,
) -> Result<Vec<usize>, CanonicalizationError> {
    // Union-Find is indexed by POSITION within `members` (0..members.len()),
    // not by atom index -- `members` are raw atom indices which may be
    // sparse/non-contiguous (e.g. a target cell of atoms 3 and 7), so
    // `parent` must never store a raw atom index as if it were a position.
    let mut parent: Vec<usize> = (0..members.len()).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    for i in 0..members.len() {
        for j in (i + 1)..members.len() {
            let ri = find(&mut parent, i);
            let rj = find(&mut parent, j);
            if ri == rj {
                continue;
            }
            // A local-twin swap is an exact automorphism: both vertices have
            // the same writer-visible color, and fixing every other vertex
            // while swapping this pair preserves the complete edge-colored
            // graph. This avoids invoking the much more expensive general
            // automorphism search for repeated terminal groups such as the
            // methyl arms of tBu/Boc. It is deliberately stricter than a
            // Morgan/WL rank comparison and rejects asymmetric dative edges.
            if local_twins(
                graph,
                coloring,
                AtomIdx(members[i] as u32),
                AtomIdx(members[j] as u32),
            ) {
                parent[ri] = rj;
                continue;
            }
            budget.automorphism_tests += 1;
            if let Some(max) = limits.max_automorphism_tests
                && budget.automorphism_tests > max
            {
                stats::record_budget_exhaustion();
                return Err(CanonicalizationError::SearchBudgetExceeded {
                    nodes_visited: budget.nodes_visited,
                    automorphism_tests: budget.automorphism_tests,
                });
            }
            let equivalent = has_colored_automorphism_mapping(
                graph,
                coloring,
                AtomIdx(members[i] as u32),
                AtomIdx(members[j] as u32),
            );
            stats::record_orbit_test(equivalent);
            if equivalent {
                parent[ri] = rj;
            }
        }
    }

    // One representative (minimum atom index) per orbit, orbits ordered by
    // that minimum ascending.
    let mut by_root: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    for (i, &atom_idx) in members.iter().enumerate() {
        let r = find(&mut parent, i);
        by_root
            .entry(r)
            .and_modify(|min_idx| {
                if atom_idx < *min_idx {
                    *min_idx = atom_idx;
                }
            })
            .or_insert(atom_idx);
    }
    let mut reps: Vec<usize> = by_root.into_values().collect();
    reps.sort_unstable();
    Ok(reps)
}

/// Return whether swapping `a` and `b` while fixing every other vertex is an
/// automorphism of the current colored graph. This is a sufficient, exact
/// test for local (true or false) twins, not a heuristic equivalence test.
fn local_twins(
    graph: &CanonicalColoredGraph,
    coloring: &Partition,
    a: AtomIdx,
    b: AtomIdx,
) -> bool {
    if a == b || coloring.cell_of[a.0 as usize] != coloring.cell_of[b.0 as usize] {
        return false;
    }

    local_twins_without_partition(graph, a, b)
}

/// Exact fixed-point swap automorphism, independent of a search partition.
/// If this succeeds, every other vertex is literally left in place, so no
/// partition refinement can invalidate the mapping.
fn local_twins_without_partition(graph: &CanonicalColoredGraph, a: AtomIdx, b: AtomIdx) -> bool {
    if a == b || graph.vertex_color(a) != graph.vertex_color(b) {
        return false;
    }

    let mut a_edges: SmallVec<[(AtomIdx, crate::canonical_partition::EdgeColor); 4]> = graph
        .neighbors(a)
        .filter(|&(neighbor, _)| neighbor != b)
        .map(|(neighbor, bond)| (neighbor, graph.edge_color(a, bond)))
        .collect();
    let mut b_edges: SmallVec<[(AtomIdx, crate::canonical_partition::EdgeColor); 4]> = graph
        .neighbors(b)
        .filter(|&(neighbor, _)| neighbor != a)
        .map(|(neighbor, bond)| (neighbor, graph.edge_color(b, bond)))
        .collect();
    a_edges.sort_unstable();
    b_edges.sort_unstable();
    if a_edges != b_edges {
        return false;
    }

    // If the pair is adjacent, its edge must be invariant under reversal.
    // This rejects a dative bond whose donor/acceptor direction would flip.
    let Some((_, a_to_b)) = graph.neighbors(a).find(|&(neighbor, _)| neighbor == b) else {
        return true;
    };
    let Some((_, b_to_a)) = graph.neighbors(b).find(|&(neighbor, _)| neighbor == a) else {
        return false;
    };
    graph.edge_color(a, a_to_b) == graph.edge_color(b, b_to_a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canonical::canonical_smiles_exhaustive_oracle;
    use crate::parser::parse;

    #[test]
    fn unbounded_matches_exhaustive_oracle_on_symmetric_molecules() {
        for smi in [
            "c1ccccc1",
            "C1CCCCC1",
            "CC(C)(C)C",
            "FC(F)(F)C",
            "C12CC3CC(CC(C3)C1)C2", // adamantane-shaped
            // E/Z alkene with an unrelated CF3 group elsewhere -- the
            // stereo-neighbor-pinning judgment call (see design doc) must
            // not disable pruning for symmetry unrelated to the stereo
            // system.
            "F/C=C/C(F)(F)F",
            "F/C=C\\C(F)(F)F",
            "c2c3c4c1c6c(ccc7ccc5ccc(c4c5c67)cc3)ccc1c2", // coronene (C24H12, 24 atoms, 7 aromatic rings, verified geometrically -- see canonical_orbit_perf.rs)
        ] {
            let mol = parse(smi).unwrap();
            let got =
                canonical_smiles_with_limits(&mol, &CanonicalizationLimits::unbounded()).unwrap();
            let oracle = canonical_smiles_exhaustive_oracle(&mol);
            assert_eq!(got, oracle, "mismatch for {smi}");
        }
    }

    #[test]
    fn local_twins_are_exact_but_dative_direction_is_not_symmetric() {
        let tbu = parse("CC(C)(C)C").unwrap();
        let graph = CanonicalColoredGraph::new(&tbu);
        let partition = initial_partition(&graph, &vec![0; graph.n()]);
        assert!(local_twins(&graph, &partition, AtomIdx(0), AtomIdx(2)));

        let mut b = chematic_core::MoleculeBuilder::new();
        let donor = b.add_atom(chematic_core::Atom::organic(chematic_core::Element::N));
        let acceptor = b.add_atom(chematic_core::Atom::organic(chematic_core::Element::N));
        b.add_bond(donor, acceptor, chematic_core::BondOrder::Dative)
            .unwrap();
        let dative = b.build();
        let dative_graph = CanonicalColoredGraph::new(&dative);
        let dative_partition = initial_partition(&dative_graph, &vec![0; dative_graph.n()]);
        assert!(!local_twins(
            &dative_graph,
            &dative_partition,
            donor,
            acceptor
        ));
    }

    /// Section 14 chemical fixtures: a molecule whose otherwise-symmetric
    /// skeleton is broken at exactly one position by isotope / formal
    /// charge / explicit H / tetrahedral stereo, plus an atom-mapped
    /// reaction fragment and a repeated disconnected fragment. Each must
    /// match the unbounded exhaustive oracle exactly (no false prune could
    /// survive that check) -- the single broken position must never be
    /// silently merged back into its formerly-symmetric siblings.
    #[test]
    fn single_position_symmetry_breaks_match_oracle() {
        for smi in [
            // Neopentane, one methyl isotope-labeled.
            "CC(C)(C)[13CH3]",
            // Neopentane-like center, one arm bearing a formal + charge
            // (via a quaternary ammonium analog) -- otherwise-identical
            // substituents, one distinguished by charge alone.
            "[NH4+]",
            "C[N+](C)(C)C",
            // Explicit-H bracket forcing on one of four otherwise-identical
            // substituents (still the same effective H count, but
            // `needs_bracket` differs -- see design doc's vertex-color
            // audit).
            "CC(C)(C)[CH3]",
            // Tetrahedral stereo breaking an otherwise symmetric skeleton
            // at exactly one position.
            "CC(C)(C)[C@H](N)O",
            // Atom-mapped reaction fragment (SMIRKS-style atom maps),
            // exercising the atom-map-is-writer-visible policy.
            "[CH3:1][C:2]([CH3:3])([CH3:4])[CH3:5]",
            // Repeated disconnected fragment (independent, structurally
            // identical components). Correction (independent Round-2
            // false-prune audit, PR #193): cross-component automorphism
            // pruning between different physical copies of the same
            // fragment DOES fire here, correctly and safely (measured:
            // leaves_written=1, orbit_unions=14 with
            // --features canonical-search-instrumentation) -- an earlier
            // version of this comment claimed the opposite, which was
            // factually backwards. The oracle check below still holds
            // regardless of which orbits get merged, since it only
            // constrains the final string, not the search's internal
            // merge decisions.
            "FC(F)(F)C.FC(F)(F)C.FC(F)(F)C",
        ] {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
            let got =
                canonical_smiles_with_limits(&mol, &CanonicalizationLimits::unbounded()).unwrap();
            let oracle = canonical_smiles_exhaustive_oracle(&mol);
            assert_eq!(got, oracle, "mismatch for {smi}");
        }
    }

    /// Independent Round-1 correctness audit (PR #193) found a real
    /// verification gap: every oracle cross-check above only ever compares
    /// the winning canonical *string*, never the *rank vector*
    /// `winning_individualized_ranks_with_limits` also returns and the
    /// public `canonical_atom_order` API consumes. Two branches within one
    /// automorphism orbit can legitimately share a minimal string via
    /// different rank vectors, so string equality does not imply
    /// rank-vector equality. Compares against the unbounded exhaustive
    /// oracle's rank vector (`canonical_smiles_exhaustive_oracle_with_ranks`)
    /// -- not the legacy engine's single winning leaf, which can also
    /// legitimately differ for the same reason (Round 2 flagged this too).
    #[test]
    fn rank_vector_matches_exhaustive_oracle_not_just_string() {
        for smi in [
            "c1ccccc1",
            "C1CCCCC1",
            "CC(C)(C)C",
            "FC(F)(F)C",
            "C12CC3CC(CC(C3)C1)C2",
            "F/C=C/C(F)(F)F",
            "F/C=C\\C(F)(F)F",
            "CC(C)(C)[13CH3]",
            "[NH4+]",
            "C[N+](C)(C)C",
            "CC(C)(C)[CH3]",
            "CC(C)(C)[C@H](N)O",
            "[CH3:1][C:2]([CH3:3])([CH3:4])[CH3:5]",
            "FC(F)(F)C.FC(F)(F)C.FC(F)(F)C",
            // Independent fixtures added by the Round-1 audit, not
            // previously exercised by this test suite.
            "C1CCC2(CCCC2)C1",                      // spiro[4.4]nonane
            "c1ncncn1",                             // 1,3,5-triazine
            "Nc1nc(N)nc(N)n1",                      // melamine
            "Cc1c(C)c(C)c(C)c(C)c1C",               // hexamethylbenzene
            "C1CN2CCC1CC2",                         // quinuclidine
            "C1CC2(CC1)OCCO2",                      // 1,4-dioxaspiro[4.5]decane
            "C12C3C4C1C5C4C3C25", // cubane (this repo's existing fixture spelling)
            "CC(C)(C)c1cc(C(C)(C)C)cc(C(C)(C)C)c1", // 1,3,5-tri-tert-butylbenzene
            "C1CN2CCN1CC2",       // DABCO
        ] {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
            let (got_ranks, got_string) = winning_individualized_ranks_with_limits(
                &mol,
                &CanonicalizationLimits::unbounded(),
            )
            .unwrap();
            let (oracle_ranks, oracle_string) =
                crate::canonical::canonical_smiles_exhaustive_oracle_with_ranks(&mol);
            assert_eq!(got_string, oracle_string, "string mismatch for {smi}");
            assert_eq!(got_ranks, oracle_ranks, "rank-vector mismatch for {smi}");
        }
    }

    #[test]
    fn ez_pinning_does_not_disable_unrelated_cf3_pruning() {
        // E and Z isomers of the same skeleton must still produce different
        // canonical strings (the stereo information itself must never be
        // lost or pruned away)...
        let e = parse("F/C=C/C(F)(F)F").unwrap();
        let z = parse("F/C=C\\C(F)(F)F").unwrap();
        let e_out = canonical_smiles_with_limits(&e, &CanonicalizationLimits::unbounded()).unwrap();
        let z_out = canonical_smiles_with_limits(&z, &CanonicalizationLimits::unbounded()).unwrap();
        assert_ne!(
            e_out, z_out,
            "E and Z isomers must not canonicalize to the same string"
        );

        // ...while the CF3 group's three (mutually automorphic, chemically
        // unrelated to the E/Z system) fluorines still get pruned to a
        // single explored orbit representative rather than 3 separate
        // individualize branches. Checked via the exhaustive oracle
        // cross-check above (equality holds) plus a direct node-count
        // sanity bound here: a plain CF3-only molecule needs few search
        // nodes, and adding an unrelated, already-fully-discrete E/Z system
        // must not multiply that cost by 3 (which is what NOT pruning the
        // CF3 fluorines would look like).
        let (_, budget_nodes) = {
            let mut budget = SearchBudget::default();
            let mut incumbent = None;
            let graph = CanonicalColoredGraph::new(&e);
            let ranks = crate::canonical::morgan_ranks(&e);
            search_canonical(
                &e,
                &graph,
                SearchNode {
                    ranks,
                    target_cell: None,
                    depth: 0,
                },
                &CanonicalizationLimits::unbounded(),
                &mut budget,
                &mut incumbent,
            )
            .unwrap();
            (incumbent, budget.nodes_visited)
        };
        // Without CF3 pruning, the 3 fluorines alone would force >= 3! / (orbit) branch
        // multiplicity layered on top of the E/Z system; with pruning, a
        // handful of nodes suffices (well under a naive-enumeration-scale
        // count). This is a coarse regression guard, not a precise bound.
        assert!(
            budget_nodes < 20,
            "expected CF3 pruning to keep node count small, got {budget_nodes}"
        );
    }

    /// Issue #421: on a real 94-atom ChEMBL molecule (3 near-identical
    /// Boc-protected benzylamine arms off a symmetric polyamine core)
    /// reordered into `canonical_atom_order`'s own output order,
    /// `canonical_smiles` used to hang -- observed running past 2 minutes,
    /// never confirmed to terminate. Root cause was in
    /// `canonical_automorphism::extend_mapping`, not this module: an
    /// unbounded backtracking search with no internal step cap, which the
    /// `SearchBudget` here cannot see (it only counts *calls* to
    /// `has_colored_automorphism_mapping`, not work done inside one call).
    /// Fixed there via an always-on `MAX_EXTEND_MAPPING_STEPS` ceiling that
    /// falls back to `false` (a documented-safe result -- see that module's
    /// own invariant) rather than searching unbounded. This test pins both
    /// that the fix actually bounds the search (a generous but finite node
    /// budget suffices, where before it would not terminate at all) and
    /// that atom-order invariance holds despite the internal fallback (the
    /// reordered input must still canonicalize to the exact same string as
    /// the original).
    #[test]
    fn issue421_reordered_symmetric_molecule_does_not_hang() {
        use chematic_core::{AtomIdx, MoleculeBuilder};

        let smi = "CC(C)(C)OC(=O)N(CCCCCN1CCCN(CCCCCN(Cc2ccccc2)C(=O)OC(C)(C)C)CCN(CCCCCN(Cc2ccccc2)C(=O)OC(C)(C)C)CCCN(CCCCCN(Cc2ccccc2)C(=O)OC(C)(C)C)CC1)Cc1ccccc1";
        let mol = parse(smi).unwrap();
        assert_eq!(mol.atom_count(), 94);

        // Reorder atoms into canonical_atom_order's own output order -- the
        // exact composition that triggered the hang.
        let order = crate::canonical_atom_order(&mol);
        let mut builder = MoleculeBuilder::new();
        let mut remap = std::collections::HashMap::new();
        for &old in &order {
            let old_idx = AtomIdx(old as u32);
            remap.insert(old_idx, builder.add_atom(mol.atom(old_idx).clone()));
        }
        for i in 0..mol.bond_count() {
            let b = mol.bond(chematic_core::BondIdx(i as u32));
            builder
                .add_bond(remap[&b.atom1], remap[&b.atom2], b.order)
                .unwrap();
        }
        let reordered = builder.build();

        // Generous but finite budget: before the fix, no finite budget on
        // either axis mattered because a single automorphism-test call
        // itself never returned.
        let limits = CanonicalizationLimits {
            max_search_nodes: Some(10_000),
            max_automorphism_tests: Some(1_000_000),
        };
        let reordered_result = canonical_smiles_with_limits(&reordered, &limits)
            .expect("reordered molecule must canonicalize within a bounded search");

        let original_result = canonical_smiles_with_limits(&mol, &limits)
            .expect("original-order molecule must canonicalize within a bounded search");
        assert_eq!(
            reordered_result, original_result,
            "canonical_smiles must be atom-order-invariant even when the internal \
             automorphism-search step cap falls back to `false`"
        );
    }

    #[test]
    fn tiny_node_budget_fails_closed_not_empty_string() {
        let mol = parse("c1ccccc1").unwrap();
        let limits = CanonicalizationLimits {
            max_search_nodes: Some(1),
            max_automorphism_tests: None,
        };
        let result = canonical_smiles_with_limits(&mol, &limits);
        assert!(
            matches!(
                result,
                Err(CanonicalizationError::SearchBudgetExceeded { .. })
            ),
            "expected SearchBudgetExceeded, got {result:?}"
        );
    }

    #[test]
    fn tiny_automorphism_test_budget_fails_closed() {
        let mol = parse("c1ccccc1").unwrap();
        let limits = CanonicalizationLimits {
            max_search_nodes: None,
            max_automorphism_tests: Some(0),
        };
        let result = canonical_smiles_with_limits(&mol, &limits);
        assert!(matches!(
            result,
            Err(CanonicalizationError::SearchBudgetExceeded { .. })
        ));
    }

    #[test]
    fn generous_budget_succeeds_on_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let limits = CanonicalizationLimits {
            max_search_nodes: Some(1000),
            max_automorphism_tests: Some(1000),
        };
        assert!(canonical_smiles_with_limits(&mol, &limits).is_ok());
    }

    // --- Negative control: never return Ok on budget exhaustion -----------
    #[test]
    fn budget_exceeded_is_always_err_never_ok_with_wrong_string() {
        let mol = parse("C1CC2CCC1CC2").unwrap(); // bicyclic, has real ties
        let limits = CanonicalizationLimits {
            max_search_nodes: Some(1),
            max_automorphism_tests: Some(1),
        };
        match canonical_smiles_with_limits(&mol, &limits) {
            Err(CanonicalizationError::SearchBudgetExceeded { .. }) => {}
            other => panic!("expected SearchBudgetExceeded, got {other:?}"),
        }
    }
}
