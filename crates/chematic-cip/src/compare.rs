//! Recursive, branch-by-branch CIP comparison -- Rules 1a (atomic number), 1b
//! (duplicate-node handling), and 2 (isotope). Deliberately **not** Rules 3+
//! (stereo-dependent): this module has no notion of R/S, E/Z, or pseudoasymmetry.
//!
//! This is the actual replacement for the old, approximate engine's
//! `cip_branch_spheres`/`compare_branches` (`crates/chematic-chem/src/cip.rs`), which
//! pools every atom at a given BFS depth into one sorted multiset and compares
//! shell-by-shell -- discarding exactly the branch/provenance information a correct
//! comparison needs (proven by a reverted triple-bond fix that went net negative on
//! that engine, see `docs/cip_accurate_rfc.md`). [`compare_ligands`] instead recurses
//! branch-by-branch: two nodes are compared by their own key first, then -- only if
//! that ties -- by their *ranked* children, position by position.
//!
//! # Rule 1b, scoped
//!
//! Strict IUPAC 2013 Rule 1b is a *duplicate-vs-duplicate* tiebreak (which of two
//! same-element ring-closure duplicates is closer, along the digraph, to the real atom
//! it duplicates). This module does not implement that -- it's a ring-closure-specific
//! concern, deferred to the ring/aromatic milestone. What Milestone 2 *does* give
//! duplicate nodes "for free" is atomic-number-based comparison via their effective
//! element ([`node_key`]): a duplicate is a childless leaf, so once its atomic number
//! ties with a real atom's, the real atom's (usually non-empty) child list outranks the
//! duplicate's (always empty) one during the recursive children comparison -- no
//! separate "real beats duplicate" rule needed. Worked example: an aldehyde carbon
//! (`C=O`) vs. a hydroxymethyl carbon (`C-OH`) both present a real oxygen at rank 1
//! (tied); at rank 2 the aldehyde side has an oxygen *duplicate* (atomic number 8)
//! against the alcohol side's hydrogen (atomic number 1) -- Rule 1a alone decides it.
//!
//! # Rule 2, scoped
//!
//! Rule 2 compares `Atom.isotope: Option<u16>` (the literal isotope label -- what
//! distinguishes `13C` from `12C`), never `Element::atomic_mass()` (a fixed per-element
//! default that cannot). The old engine's `atom_key`/`cmp_key` conflates isotope
//! (Rule 2) with atomic mass (Rule 4, a tiebreaker) in one tuple; this module keeps
//! them separate and implements only Rule 2, since Rule 4 is out of scope this
//! milestone.
//!
//! # Why `rank_children` never calls `sort_by`
//!
//! A raw `sort_by`/`sort_unstable_by` call fed a comparator that hasn't fully resolved
//! every pair first can silently produce a wrong order if the comparator isn't a true
//! total order across the whole slice -- reproducing the shape of bug this crate
//! exists to avoid. [`rank_children`] instead: computes the full pairwise comparison
//! matrix up front (children counts are small, O(n^2) is cheap), merges `Equal`/
//! `Unresolved` pairs into equivalence classes via union-find (so genuinely tied
//! siblings are treated as fungible by construction, never arbitrarily ordered), and
//! only then orders the *classes* using their already-computed pairwise results --
//! sorting fully-resolved, deduplicated data, not a lazily-evaluated comparator.

use std::cmp::Ordering;
use std::collections::HashMap;

use chematic_core::Molecule;

use crate::CipError;
use crate::digraph::CipDigraph;
use crate::node::{CipNodeKind, NodeId};
use crate::trace::{ComparisonTrace, DecisionStep};

/// The outcome of comparing two ligand branches under Rules 1a/1b/2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BranchComparison {
    Higher,
    Lower,
    /// Provably indistinguishable under Rules 1a/1b/2, given the budget actually
    /// available -- a confident, stable tie (not a truncation artifact).
    Equal,
    /// Reserved for "a higher rule (3+, stereo-dependent) would be needed to break
    /// this tie" -- out of scope for this module entirely. Under a Rules-1a/1b/2-only
    /// comparator over a digraph whose finiteness is structurally guaranteed
    /// (Milestone 1's ancestor-path ring-closure rule), every comparison terminates in
    /// `Higher`/`Lower`/`Equal` or errors as `CipCompareError::BudgetExceeded` --
    /// `Unresolved` is unreachable from this module's own logic. Kept in the enum (a
    /// later milestone needs it) and handled defensively wherever it could appear, but
    /// not something [`compare_ligands`] itself ever returns as of Milestone 2.
    Unresolved,
}

fn invert(cmp: BranchComparison) -> BranchComparison {
    match cmp {
        BranchComparison::Higher => BranchComparison::Lower,
        BranchComparison::Lower => BranchComparison::Higher,
        BranchComparison::Equal => BranchComparison::Equal,
        BranchComparison::Unresolved => BranchComparison::Unresolved,
    }
}

/// Errors from comparison, distinct from a successful (even if `Unresolved`) outcome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CipCompareError {
    /// The comparison's own recursion budget was exceeded (separate from, but usually
    /// triggered alongside, the underlying digraph's own node/depth/expansion budget).
    BudgetExceeded {
        expanded_nodes: usize,
        recursive_calls: usize,
    },
    InvalidDigraph(String),
    /// The underlying digraph's own budget was exceeded while expanding children.
    Digraph(CipError),
}

impl core::fmt::Display for CipCompareError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CipCompareError::BudgetExceeded {
                expanded_nodes,
                recursive_calls,
            } => write!(
                f,
                "CIP comparison budget exceeded: {expanded_nodes} nodes expanded, {recursive_calls} recursive calls"
            ),
            CipCompareError::InvalidDigraph(reason) => write!(f, "invalid digraph: {reason}"),
            CipCompareError::Digraph(e) => write!(f, "digraph error: {e}"),
        }
    }
}

impl std::error::Error for CipCompareError {}

/// Mutable state threaded through a comparison: a recursion-call budget (separate from
/// the digraph's own node budget) and an optional trace sink.
pub struct CompareContext<'t> {
    pub recursive_calls: usize,
    pub max_recursive_calls: usize,
    trace: Option<&'t mut ComparisonTrace>,
}

impl Default for CompareContext<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'t> CompareContext<'t> {
    pub fn new() -> Self {
        Self {
            recursive_calls: 0,
            max_recursive_calls: 1_000_000,
            trace: None,
        }
    }

    pub fn with_trace(trace: &'t mut ComparisonTrace) -> Self {
        Self {
            recursive_calls: 0,
            max_recursive_calls: 1_000_000,
            trace: Some(trace),
        }
    }
}

/// `(atomic_number, isotope)` for a node's *effective* atom: a real `Atom` node looks
/// itself up; a duplicate node looks up the atom it duplicates (its own atomic number,
/// per Rule 1a -- see module docs for why no separate real-vs-duplicate rule is
/// needed); the implicit-hydrogen sentinel is `(1, None)`.
fn node_key(mol: &Molecule, kind: CipNodeKind) -> (u8, Option<u16>) {
    let atom_idx = match kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        CipNodeKind::MultipleBondDuplicate {
            duplicated_atom, ..
        } => duplicated_atom,
        CipNodeKind::RingDuplicate { closure_atom, .. } => closure_atom,
        CipNodeKind::ImplicitHydrogen => return (1, None),
    };
    let atom = mol.atom(atom_idx);
    (atom.element.atomic_number(), atom.isotope)
}

/// Rule 1a (atomic number) then Rule 2 (isotope: `Some` beats `None`, higher isotope
/// number beats lower).
fn cmp_key(a: (u8, Option<u16>), b: (u8, Option<u16>)) -> Ordering {
    match a.0.cmp(&b.0) {
        Ordering::Equal => {}
        other => return other,
    }
    match (a.1, b.1) {
        (Some(x), Some(y)) => x.cmp(&y),
        (Some(_), None) => Ordering::Greater,
        (None, Some(_)) => Ordering::Less,
        (None, None) => Ordering::Equal,
    }
}

fn kind_label(kind: CipNodeKind) -> String {
    match kind {
        CipNodeKind::Atom { atom_idx } => format!("atom({})", atom_idx.0),
        CipNodeKind::MultipleBondDuplicate {
            duplicated_atom, ..
        } => {
            format!("dup-multi({})", duplicated_atom.0)
        }
        CipNodeKind::RingDuplicate { closure_atom, .. } => format!("dup-ring({})", closure_atom.0),
        CipNodeKind::ImplicitHydrogen => "H".to_string(),
    }
}

/// Compare two ligand branches recursively: own key first (Rules 1a/2), then --
/// only if that ties -- ranked children position by position. See module docs for the
/// full algorithm and its rationale.
pub fn compare_ligands(
    graph: &mut CipDigraph,
    left: NodeId,
    right: NodeId,
    ctx: &mut CompareContext,
) -> Result<BranchComparison, CipCompareError> {
    ctx.recursive_calls += 1;
    if ctx.recursive_calls > ctx.max_recursive_calls {
        return Err(CipCompareError::BudgetExceeded {
            expanded_nodes: graph.nodes().len(),
            recursive_calls: ctx.recursive_calls,
        });
    }

    let left_kind = graph.node(left).kind;
    let right_kind = graph.node(right).kind;
    let depth = graph.node(left).depth;
    let left_key = node_key(graph.molecule(), left_kind);
    let right_key = node_key(graph.molecule(), right_kind);

    let key_cmp = cmp_key(left_key, right_key);
    let (outcome, rule) = if key_cmp != Ordering::Equal {
        let cmp = if key_cmp == Ordering::Greater {
            BranchComparison::Higher
        } else {
            BranchComparison::Lower
        };
        (cmp, "1a/2")
    } else {
        let left_children = graph
            .expand_children(left)
            .map_err(CipCompareError::Digraph)?;
        let right_children = graph
            .expand_children(right)
            .map_err(CipCompareError::Digraph)?;

        if left_children.is_empty() && right_children.is_empty() {
            (BranchComparison::Equal, "leaf")
        } else {
            let left_groups = rank_children(graph, &left_children, ctx)?;
            let right_groups = rank_children(graph, &right_children, ctx)?;
            let cmp = compare_ranked(graph, &left_groups, &right_groups, ctx)?;
            (cmp, "children")
        }
    };

    if let Some(trace) = ctx.trace.as_deref_mut() {
        trace.decisions.push(DecisionStep {
            depth,
            left_kind: kind_label(left_kind),
            right_kind: kind_label(right_kind),
            outcome,
            rule,
        });
    }

    Ok(outcome)
}

/// Flatten two sides' ranked groups into priority-ordered sequences and compare
/// position by position. Safe to pick any representative from a tied group at a given
/// position: by construction, every member of a group compares `Equal` to every other
/// member, so the outcome of comparing "some position N" doesn't depend on which
/// member of its group was chosen to sit there.
fn compare_ranked(
    graph: &mut CipDigraph,
    left_groups: &[Vec<NodeId>],
    right_groups: &[Vec<NodeId>],
    ctx: &mut CompareContext,
) -> Result<BranchComparison, CipCompareError> {
    let left_flat: Vec<NodeId> = left_groups.iter().flatten().copied().collect();
    let right_flat: Vec<NodeId> = right_groups.iter().flatten().copied().collect();

    let n = left_flat.len().min(right_flat.len());
    for i in 0..n {
        match compare_ligands(graph, left_flat[i], right_flat[i], ctx)? {
            BranchComparison::Equal | BranchComparison::Unresolved => continue,
            other => return Ok(other),
        }
    }
    match left_flat.len().cmp(&right_flat.len()) {
        Ordering::Greater => Ok(BranchComparison::Higher),
        Ordering::Less => Ok(BranchComparison::Lower),
        Ordering::Equal => Ok(BranchComparison::Equal),
    }
}

/// Rank `children` into priority-ordered groups (highest first; multiple entries in
/// one group means those children are mutually tied under Rules 1a/1b/2). Never calls
/// `sort_by`/`sort_unstable_by` on a lazily-evaluated comparator -- see module docs.
// Index-based loops throughout: this fills and reads a 2D pairwise matrix by (i, j)
// index pairs while also indexing the separate `children` slice by the same `i`/`j` --
// an iterator-adapter rewrite would need to zip three collections by position and is
// less readable than the direct indices here.
#[allow(clippy::needless_range_loop)]
pub fn rank_children(
    graph: &mut CipDigraph,
    children: &[NodeId],
    ctx: &mut CompareContext,
) -> Result<Vec<Vec<NodeId>>, CipCompareError> {
    let n = children.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    if n == 1 {
        return Ok(vec![vec![children[0]]]);
    }

    // Full pairwise matrix, computed once, up front.
    let mut pairwise = vec![vec![BranchComparison::Equal; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let cmp = compare_ligands(graph, children[i], children[j], ctx)?;
            pairwise[i][j] = cmp;
            pairwise[j][i] = invert(cmp);
        }
    }

    // Union-find: Equal and Unresolved pairs merge into one equivalence class --
    // Unresolved is treated as "tied for this ranking's purposes," never used to
    // establish an order (see BranchComparison::Unresolved's doc comment).
    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            if matches!(
                pairwise[i][j],
                BranchComparison::Equal | BranchComparison::Unresolved
            ) {
                union(&mut parent, i, j);
            }
        }
    }

    let mut groups_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        groups_map.entry(root).or_default().push(i);
    }
    let mut groups: Vec<Vec<usize>> = groups_map.into_values().collect();

    // Order the (already deduplicated, already fully-resolved) classes using their
    // precomputed representative comparison -- sorting resolved data, not invoking the
    // comparator lazily mid-sort.
    groups.sort_by(|a, b| match pairwise[a[0]][b[0]] {
        BranchComparison::Higher => Ordering::Less,
        BranchComparison::Lower => Ordering::Greater,
        BranchComparison::Equal | BranchComparison::Unresolved => Ordering::Equal,
    });

    // Defensive: under Rules 1a/1b/2 alone (a genuine recursive comparator, not
    // shell-pooling) the classes should form a strict total order with no cycles. If
    // this ever fires, the comparator has a real transitivity bug -- surface it loudly
    // in debug builds rather than silently ship a possibly-wrong order.
    debug_assert!(
        groups
            .windows(2)
            .all(|w| !matches!(pairwise[w[1][0]][w[0][0]], BranchComparison::Higher)),
        "rank_children: inconsistent (non-transitive) pairwise comparison among children"
    );

    Ok(groups
        .into_iter()
        .map(|g| g.into_iter().map(|i| children[i]).collect())
        .collect())
}

fn find(parent: &mut [usize], x: usize) -> usize {
    if parent[x] != x {
        parent[x] = find(parent, parent[x]);
    }
    parent[x]
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}
