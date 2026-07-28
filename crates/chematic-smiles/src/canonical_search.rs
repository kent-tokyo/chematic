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
//! sequence -- see `docs/canonical_automorphism_pruning.md`, "why pruning
//! cannot change surviving output", for the full argument and its
//! `canonical_smiles_exhaustive_oracle` cross-check.

use chematic_core::{AtomIdx, Molecule};

use crate::canonical::{CanonicalWriter, group_by_rank, individualize, refine_ranks};
use crate::canonical_automorphism::has_colored_automorphism_mapping;
use crate::canonical_partition::{
    CanonicalColoredGraph, Partition, exact_refine, initial_partition,
};

mod stats {
    //! Feature-gated internal work counters, following the exact pattern
    //! established by `chematic-rxn`'s `perf_counters.rs`
    //! (`docs/reaction_transform_perf.md`): process-global `AtomicUsize`s,
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

#[derive(Default)]
struct SearchBudget {
    nodes_visited: usize,
    automorphism_tests: usize,
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
    let graph = CanonicalColoredGraph::new(mol);
    let ranks = crate::canonical::morgan_ranks(mol);
    let mut budget = SearchBudget::default();
    let mut incumbent: Option<Incumbent> = None;
    search_canonical(mol, &graph, ranks, 0, limits, &mut budget, &mut incumbent)?;
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
    ranks: Vec<u64>,
    depth: usize,
    limits: &CanonicalizationLimits,
    budget: &mut SearchBudget,
    incumbent: &mut Option<Incumbent>,
) -> Result<(), CanonicalizationError> {
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

    let cells = group_by_rank(&ranks);
    // Same heuristic as the legacy enumeration (section 10: unchanged in
    // this PR): the lowest-ranked non-singleton cell.
    let Some(members) = cells.iter().find(|m| m.len() > 1) else {
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
    let partition = exact_refine(graph, initial_partition(graph, &ranks));
    let representatives = exact_orbit_representatives(graph, &partition, members, limits, budget)?;
    stats::record_target_cell(members.len(), representatives.len());

    for rep in representatives {
        let individualized = individualize(&ranks, rep);
        let re_refined = refine_ranks(mol, individualized);
        search_canonical(mol, graph, re_refined, depth + 1, limits, budget, incumbent)?;
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
/// a false positive is structurally impossible here since every union is
/// gated on a real, independently-reverified bijection.
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
            // Repeated disconnected fragment (independent components,
            // no cross-component automorphism should ever be proposed
            // between DIFFERENT physical copies for pruning -- each is
            // still correctly recognized as its own separate 3-fold CF3
            // orbit within itself).
            "FC(F)(F)C.FC(F)(F)C.FC(F)(F)C",
        ] {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
            let got =
                canonical_smiles_with_limits(&mol, &CanonicalizationLimits::unbounded()).unwrap();
            let oracle = canonical_smiles_exhaustive_oracle(&mol);
            assert_eq!(got, oracle, "mismatch for {smi}");
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
                ranks,
                0,
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
