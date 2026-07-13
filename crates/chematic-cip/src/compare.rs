//! Recursive, branch-by-branch CIP comparison -- Rules 1a (atomic number), 1b
//! (duplicate-node handling), and 2 (isotope). Deliberately **not** Rules 3+
//! (stereo-dependent): this module has no notion of R/S, E/Z, or pseudoasymmetry.
//!
//! This is the actual replacement for the old, approximate engine's
//! `cip_branch_spheres`/`compare_branches` (`crates/chematic-chem/src/cip.rs`), which
//! pools every atom at a given BFS depth into one sorted multiset and compares
//! shell-by-shell -- discarding exactly the branch/provenance information a correct
//! comparison needs (proven by a reverted triple-bond fix that went net negative on
//! that engine, see `docs/cip_accurate_rfc.md`).
//!
//! # Sphere-by-sphere, not one branch to completion before the next
//!
//! [`compare_ligands`] compares two ligand roots one *level* (sphere) at a time: level 0
//! is the two roots' own keys; level 1 is their ranked children's own keys, compared
//! position by position *without recursing into any position's own children first*.
//! Only once an entire level ties completely across every position does the comparison
//! advance to the next level. An earlier version instead resolved each position
//! *fully* (recursing arbitrarily deep) before checking the next position -- depth-first
//! when CIP requires breadth-first, and wrong whenever a later, shallow difference
//! should decide the comparison before an earlier position's much deeper (and
//! irrelevant) difference is even reached.
//!
//! Concretely, this is what the corpus's lone `triple_bond_dup` case exposed: comparing
//! an ethynyl branch (`-C#CH`, whose own children are three carbons -- the real
//! terminal carbon plus two triple-bond duplicates) against a carboxymethyl branch
//! (`-CH2COOH`, whose own children are one carbon and two hydrogens) must be decided at
//! *that* level -- second child C(6) beats second child H(1) -- without ever descending
//! into the carboxyl group's oxygens several levels further down the carboxymethyl
//! branch. The depth-first version dove into those oxygens first (since the *first*
//! children on both sides tied on atomic number alone) and returned the wrong winner
//! without ever checking the second children. Fixed by making the recursion genuinely
//! level-order: each step advances *all* live position-pairs on both sides by exactly
//! one sphere, together, and only descends past a level once that whole level has tied.
//!
//! Ranking children *within* one parent (via [`rank_children`], to establish which
//! child is "position 0" vs "position 1" for the next level) is a separate,
//! locally-scoped comparison problem -- ordering priority among a handful of siblings
//! under a single node -- and legitimately recurses as deep as needed to resolve that;
//! it lacks the "a shallow sibling difference must win" hazard the top-level
//! ligand-vs-ligand comparison has, because it's ordering one node's own children
//! against each other, not advancing two whole subtrees in lockstep.
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
/// the digraph's own node budget), an optional trace sink, and the innermost
/// `rank_children` call currently in progress (see [`DecisionStep::ranking_parent`]).
pub struct CompareContext<'t> {
    pub recursive_calls: usize,
    pub max_recursive_calls: usize,
    trace: Option<&'t mut ComparisonTrace>,
    /// Which node's children `rank_children` is currently ordering, if a
    /// `rank_children` call is on the stack -- set/restored around its pairwise
    /// `compare_ligands` calls (see `rank_children`'s body) so each recorded
    /// `DecisionStep` can be tagged with the sibling group it belongs to.
    ranking_parent: Option<NodeId>,
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
            ranking_parent: None,
        }
    }

    pub fn with_trace(trace: &'t mut ComparisonTrace) -> Self {
        Self {
            recursive_calls: 0,
            max_recursive_calls: 1_000_000,
            trace: Some(trace),
            ranking_parent: None,
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

/// One slot in a level being compared: either a real digraph node, or padding for a
/// position whose counterpart (established as corresponding by the previous level's tie)
/// has more substituents than this side does at this exact spot. Padding with a
/// phantom (atomic number 0, always childless) lets a substituent-*count* mismatch at
/// any one position resolve through the same position-by-position own-key comparison as
/// everything else, localized to the position it actually occurs at -- rather than
/// flattening every parent's children into one combined list and comparing raw total
/// lengths, which shifts every position *after* the mismatch out of alignment (found via
/// the `aromatic_mancude` corpus bucket regressing hard under an earlier version of this
/// function that did exactly that).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LevelSlot {
    Node(NodeId),
    Phantom,
}

fn slot_key(graph: &CipDigraph, slot: LevelSlot) -> (u8, Option<u16>) {
    match slot {
        LevelSlot::Node(n) => node_key(graph.molecule(), graph.node(n).kind),
        LevelSlot::Phantom => (0, None),
    }
}

/// A real node's own ranked children, as plain `NodeId`s (a phantom slot's "children"
/// are always empty -- it never gets here since [`compare_ligands`] only calls this for
/// `LevelSlot::Node`).
fn ranked_child_ids(
    graph: &mut CipDigraph,
    node: NodeId,
    ctx: &mut CompareContext,
) -> Result<Vec<NodeId>, CipCompareError> {
    let children = graph
        .expand_children(node)
        .map_err(CipCompareError::Digraph)?;
    Ok(rank_children(graph, &children, ctx)?
        .into_iter()
        .flatten()
        .collect())
}

/// Compare two ligand branches sphere by sphere: own key first (Rules 1a/2), then --
/// only if that ties -- ranked children, one whole level at a time. See module docs
/// ("Sphere-by-sphere...") for why this must be level-order rather than resolving one
/// position's subtree to completion before checking its siblings.
pub fn compare_ligands(
    graph: &mut CipDigraph,
    left: NodeId,
    right: NodeId,
    ctx: &mut CompareContext,
) -> Result<BranchComparison, CipCompareError> {
    let left_kind = graph.node(left).kind;
    let right_kind = graph.node(right).kind;
    let depth = graph.node(left).depth;

    let mut left_level = vec![LevelSlot::Node(left)];
    let mut right_level = vec![LevelSlot::Node(right)];
    let mut rule = "1a/2";

    let outcome = loop {
        ctx.recursive_calls += 1;
        if ctx.recursive_calls > ctx.max_recursive_calls {
            return Err(CipCompareError::BudgetExceeded {
                expanded_nodes: graph.nodes().len(),
                recursive_calls: ctx.recursive_calls,
            });
        }

        // Invariant: left_level.len() == right_level.len() always -- any per-position
        // substituent-count mismatch is padded with a Phantom slot at the position it
        // originates (see the next_left/next_right construction below), never left to
        // shift later positions out of alignment.
        debug_assert_eq!(left_level.len(), right_level.len());
        let n = left_level.len();

        // Compare this whole level's own keys, position by position, *before*
        // descending into any position's children -- a shallow difference at position
        // 1 must win even if position 0 would only differ several spheres further down.
        // A Phantom's key (0, None) loses to any real atom, so a substituent-count
        // mismatch decides here too, via the same mechanism as an atomic-number
        // difference.
        let mut decided = None;
        for i in 0..n {
            let lk = slot_key(graph, left_level[i]);
            let rk = slot_key(graph, right_level[i]);
            match cmp_key(lk, rk) {
                Ordering::Greater => {
                    decided = Some(BranchComparison::Higher);
                    break;
                }
                Ordering::Less => {
                    decided = Some(BranchComparison::Lower);
                    break;
                }
                Ordering::Equal => {}
            }
        }
        if let Some(outcome) = decided {
            break outcome;
        }
        if n == 0 {
            rule = "leaf";
            break BranchComparison::Equal;
        }

        // Whole level tied on own keys: expand to the next sphere. Each position's
        // children are ranked *within that one parent* (a separate, locally-scoped
        // comparison -- see module docs) before being placed into the next level;
        // per-position padding (not a global flatten) keeps "position i" well-defined
        // on both sides even when a tied pair's substituent counts differ.
        rule = "children";
        let mut next_left = Vec::new();
        let mut next_right = Vec::new();
        for i in 0..n {
            let left_children = match left_level[i] {
                LevelSlot::Node(node) => ranked_child_ids(graph, node, ctx)?,
                LevelSlot::Phantom => Vec::new(),
            };
            let right_children = match right_level[i] {
                LevelSlot::Node(node) => ranked_child_ids(graph, node, ctx)?,
                LevelSlot::Phantom => Vec::new(),
            };
            let max_len = left_children.len().max(right_children.len());
            for j in 0..max_len {
                next_left.push(
                    left_children
                        .get(j)
                        .map(|&n| LevelSlot::Node(n))
                        .unwrap_or(LevelSlot::Phantom),
                );
                next_right.push(
                    right_children
                        .get(j)
                        .map(|&n| LevelSlot::Node(n))
                        .unwrap_or(LevelSlot::Phantom),
                );
            }
        }
        left_level = next_left;
        right_level = next_right;
    };

    let ranking_parent = ctx.ranking_parent;
    if let Some(trace) = ctx.trace.as_deref_mut() {
        trace.decisions.push(DecisionStep {
            depth,
            left_kind: kind_label(left_kind),
            right_kind: kind_label(right_kind),
            outcome,
            rule,
            ranking_parent,
        });
    }

    Ok(outcome)
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

    // Every child in `children` is a sibling under the same parent (by construction --
    // this is always called with one node's own `expand_children` output) -- record it
    // for the pairwise calls below so their DecisionSteps are tagged with which sibling
    // group they belong to. Restored right after the pairwise fill (on both the success
    // and the `?`-propagated error path -- an error here unwinds the entire recursive
    // call chain via `?`, so no later code in this ctx's lifetime observes a stale value
    // either way, but restoring promptly keeps that true by inspection, not by relying
    // on the unwind).
    let siblings_parent = graph.node(children[0]).parent;
    let saved_ranking_parent = ctx.ranking_parent;
    ctx.ranking_parent = siblings_parent;

    // Full pairwise matrix, computed once, up front.
    let mut pairwise = vec![vec![BranchComparison::Equal; n]; n];
    let mut fill_err = None;
    'fill: for i in 0..n {
        for j in (i + 1)..n {
            match compare_ligands(graph, children[i], children[j], ctx) {
                Ok(cmp) => {
                    pairwise[i][j] = cmp;
                    pairwise[j][i] = invert(cmp);
                }
                Err(e) => {
                    fill_err = Some(e);
                    break 'fill;
                }
            }
        }
    }
    ctx.ranking_parent = saved_ranking_parent;
    if let Some(e) = fill_err {
        return Err(e);
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
