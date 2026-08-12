//! Recursive, branch-by-branch CIP comparison -- Rules 1a (atomic number), 1b
//! (duplicate-node handling), and 2 (isotope). Deliberately **not** Rules 3+
//! (stereo-dependent): this module has no notion of R/S, E/Z, or pseudoasymmetry.
//!
//! This is the actual replacement for the old, approximate engine's
//! `cip_branch_spheres`/`compare_branches` (`crates/chematic-chem/src/cip.rs`), which
//! pools every atom at a given BFS depth into one sorted multiset and compares
//! shell-by-shell -- discarding exactly the branch/provenance information a correct
//! comparison needs (proven by a reverted triple-bond fix that went net negative on
//! that engine, see `docs/rfcs/cip_accurate_rfc.md`).
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
//! # Rule 1b: investigated in Milestone 3A, implemented, then reverted (not a defect)
//!
//! Rule 1b applies **only** to [`CipNodeKind::RingDuplicate`] nodes -- never
//! `MultipleBondDuplicate` (verified against RDKit's actual `Rule1b::compare`, which
//! checks `isSet(Node::RING_DUPLICATE)` specifically). Most real-vs-duplicate and
//! duplicate-vs-duplicate distinctions are already free under Rule 1a alone, since a
//! duplicate is a childless leaf: once its atomic number ties with a real atom's, the
//! real atom's non-empty child list outranks the duplicate's empty one one sphere deeper
//! (worked example: an aldehyde carbon `C=O` vs. a hydroxymethyl carbon `C-OH` both
//! present a real oxygen at rank 1 -- tied; at rank 2 the aldehyde side has an oxygen
//! *duplicate* against the alcohol side's hydrogen -- Rule 1a alone decides it, no Rule
//! 1b needed).
//!
//! Milestone 3A implemented Rule 1b as a second `compare_by_level` pass (per RDKit's
//! `Rule1b.{h,cpp}`/`Node.cpp`: a ring duplicate outranks a non-duplicate
//! *unconditionally*, and between two duplicates the one whose real atom is closer to
//! the root wins) and verified it against RDKit source and oracle-checked molecules --
//! the implementation itself was correct. But on the frozen corpus's 24
//! `uncharacterized`-bucket tied cases it changed **zero** outputs: instrumenting the
//! Rule 1b comparator across all 8 non-pseudoasymmetric tied cases (the other 16 are
//! pseudoasymmetric, Rule 5 territory) found 4056 real invocations, none of which ever
//! compared a duplicate against a same-element non-duplicate, or two duplicates at
//! different depths -- Rule 1a's own child-count check (the paragraph above) had already
//! resolved every such position one sphere earlier. Deliberate synthetic constructions
//! (spiro, fused, bridged, decalin-like branches) reproduced the same pattern, and
//! re-deriving RDKit's own `SequenceRule::recursiveCompare` (`Sort.cpp`/
//! `SequenceRule.cpp`) shows this isn't a chematic-specific accident: RDKit's duplicate
//! nodes are *also* unconditionally childless leaves (`EXPANDED` is set at construction
//! for any `DUPLICATE` node, so `getEdges()` never expands one), so RDKit's own Rule 1a
//! hits the identical empty-vs-nonempty child-count difference before its Rule 1b ever
//! runs. Rule 1b is architecturally shadowed by Rule 1a in both engines whenever every
//! rule runs as an unconditional full pass over the whole comparison.
//!
//! Given a correct-but-inert implementation doubles `compare_ligands`'s cost for zero
//! behavior change on this corpus, the Rule 1b pass was reverted rather than shipped --
//! see `docs/rfcs/cip_accurate_rfc.md`'s Milestone 3A entry for the finding in full. The
//! `compare_by_level<K>` generic walk stays (this paragraph's `K`-parameterization is
//! exactly what let Rule 1b be tried, measured, and cleanly removed without touching the
//! Rule 1a/2 path). **Design note for Milestone 4's Rule 3/4/5**: don't repeat the
//! "unconditional full pass per rule" scheduling that made Rule 1b's cost symmetric with
//! its (absent) benefit here -- run each rule only over the equivalence classes the
//! *previous* rule left tied (`rank_by_rule_1a` -> `refine_tied_groups(_, rule_1b)` ->
//! `refine_tied_groups(_, rule_2)` -> ...), so a comparison Rule 1a already decided never
//! pays for a later rule's pass.
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
use crate::budget::CipBudget;
use crate::digraph::CipDigraph;
use crate::node::{CipNode, CipNodeKind, NodeId};
use crate::rational::{AtomicNumberKey, cmp_atomic_number_key};
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

/// Which comparator/rule this crate's own `compare_ligands` runs -- currently the
/// only variant that exists, since this crate implements Rules 1a/1b/2 only (see
/// module docs for why Rule 3+ is out of scope). Included in [`PairwiseCacheKey`] as
/// a forward-compatible discriminator: if a future Rule 3/4/5 pass is added and ever
/// shares this same cache/context type, its own comparisons must never collide with a
/// Rule-1a/2 result cached under the same `(left, right)` NodeId pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum RuleMode {
    Rule1a2,
}

/// Cache key for issue #107's `compare_ligands` pairwise-comparison memoization,
/// scoped to exactly one [`CompareContext`]'s lifetime (i.e. one stereocenter
/// resolution -- `rank_children`/`compare_ligands` never share a `CompareContext`
/// across two different resolutions, so this cache is never global). Per that
/// issue's own "first implementation candidate" spec, includes:
///
/// - `left`/`right`: the literal compared [`NodeId`]s, normalized to a canonical
///   order (`left.0 <= right.0`) so `(a, b)` and `(b, a)` share one entry -- the
///   stored [`BranchComparison`] is always relative to *this* canonical order; a
///   query in the opposite order inverts the looked-up value before returning it
///   (`Self::canonical` / the inversion in [`compare_ligands`] itself). This is the
///   "outcome-direction normalization" the issue's spec calls for.
/// - `rule_mode`: see [`RuleMode`].
/// - `mancude_identity`/`budget`: see [`CipDigraph::mancude_identity`]/
///   [`CipDigraph::budget`] -- both are in fact constant for one `CompareContext`'s
///   whole lifetime (one digraph, one resolution), so within a single cache instance
///   these fields never actually vary; they exist so this key type is still correct
///   by construction if a future change ever widens the cache's scope, rather than
///   relying on "the cache happens to only ever see one digraph" as an unenforced
///   invariant.
///
/// Ring-closure/duplicate-node context (a `RingDuplicate`'s `closure_atom`, a
/// `MultipleBondDuplicate`'s `duplicated_atom`/`bond_order`) does **not** need its own
/// key field: a [`NodeId`] already uniquely identifies one exact, immutable digraph
/// node for the digraph's whole lifetime, so keying on the literal `NodeId` (not a
/// content-based structural signature -- issue #107's own diagnosis found that
/// distinction load-bearing, see `docs/rfcs/rank_children_heavy_tail_diagnosis_107.md`)
/// already captures it: two different nodes that happen to *look* alike (same
/// represented element, different closure/duplication provenance) get different
/// `NodeId`s and therefore different, independently-cached entries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct PairwiseCacheKey {
    left: NodeId,
    right: NodeId,
    rule_mode: RuleMode,
    mancude_identity: Option<usize>,
    budget: CipBudget,
}

impl PairwiseCacheKey {
    /// Builds the canonical (order-normalized) key for `(a, b)`, plus whether `(a, b)`
    /// itself was already in canonical order (`false` means the caller's `a`/`b` are
    /// swapped relative to the key -- the cached/stored value must be inverted before
    /// being returned to a caller who asked in that order).
    fn canonical(graph: &CipDigraph, a: NodeId, b: NodeId) -> (Self, bool) {
        let already_canonical = a.0 <= b.0;
        let (left, right) = if already_canonical { (a, b) } else { (b, a) };
        (
            Self {
                left,
                right,
                rule_mode: RuleMode::Rule1a2,
                mancude_identity: graph.mancude_identity(),
                budget: graph.budget(),
            },
            already_canonical,
        )
    }
}

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
    /// Count of Rule-1a/2 key comparisons where at least one side's `AtomicNumberKey` is
    /// `Rational` -- a pure "path reached" counter (see [`Self::fractional_decisions`] for
    /// the stricter, load-bearing one). Incremented in the instrumented `cmp_fn` closure
    /// `compare_ligands` passes to `compare_by_level`, not inside `cmp_key` itself (which
    /// stays pure and independently unit-tested) or inside `compare_by_level` (which stays
    /// fully generic over `K`, with no MANCUDE-specific knowledge).
    pub fractional_comparisons: u64,
    /// Count of Rule-1a/2 key comparisons where the fraction *itself* decided the
    /// ordering -- i.e. collapsing the `Rational` side to its integer part would have
    /// produced a different (or tied) result. Deliberately **not** "non-Equal and either
    /// side is Rational": that naive definition fires even when an element difference
    /// alone decides it (e.g. `Rational(6/1)` vs `Integral(8)`, carbon vs oxygen), which
    /// would misrepresent a fraction as load-bearing when it did nothing.
    ///
    /// At the Milestone 3B-1b closeout, the RFC asserted this would be zero on the
    /// frozen corpus, "consistent with" a byte-identical-output finding -- but that
    /// assertion was only ever checked against a single curated molecule, never
    /// measured at full-corpus scale. A later full-corpus diagnostic on the same
    /// corpus (SHA-256: `1c47371d...`) found the real value is 248, across 36
    /// stereocenters in 21 molecules (Milestone MANCUDE-Decision-A0; see
    /// `docs/rfcs/cip_accurate_rfc.md`'s tripwire closeout entry for the full breakdown).
    ///
    /// This is **not** a regression, and it does **not** contradict the byte-identical
    /// finding once fraction is properly isolated from Kekule-respelling structure (the
    /// naive "with vs without `MancudeContext`" contrast conflates both): for all 36
    /// affected centers, including the 3 where the naive contrast's *final* Pass-1
    /// ranking also differs, holding structure fixed (`CipDigraph::new` on the same
    /// Kekule-respelled molecule, with vs without the attached `MancudeContext`) shows
    /// **identical** root-child partitions -- the fraction touches individual
    /// sub-branch comparisons (hence nonzero `fractional_decisions`) without ever
    /// changing the resolved label. All 36 are classification D ("fraction locally
    /// load-bearing, final-label-inert"); **zero** are classification E ("fraction
    /// changes the final label"). The 3 that flip vs the un-Kekulized baseline flip
    /// because of Kekule-respelling alone, matching modern and legacy RDKit both --
    /// frozen as regression fixtures in `tests/mancude_decision_regression.rs`.
    ///
    /// A nonzero value is therefore a tripwire requiring classification, not evidence
    /// of a problem by itself:
    /// - `fractional_decisions > 0` alone -> diagnosis required (is it D or E?).
    /// - isolate fraction from structure before concluding either way -- a bundled
    ///   with/without-`MancudeContext` contrast is not sufficient, as this milestone's
    ///   own first attempt at that contrast found the hard way.
    /// - a genuine (isolated) fractional-driven final-label change -> oracle
    ///   classification required.
    /// - one that *reduces* oracle agreement -> correctness blocker.
    pub fractional_decisions: u64,
    /// Issue #107's pairwise-comparison memoization cache -- see
    /// [`PairwiseCacheKey`]'s doc for the key design. Bounded by construction: every
    /// entry is only ever inserted as the direct result of one real
    /// `compare_ligands` call, and the number of such calls in one resolution is
    /// itself already bounded by `max_recursive_calls` (default 1,000,000) -- this
    /// cache can therefore never hold more entries than that ceiling, and a fresh,
    /// empty cache is created (and dropped) with every fresh `CompareContext`, i.e.
    /// once per stereocenter resolution, never shared or reused across resolutions
    /// or molecules. Measured in practice (5,000-molecule corpus, `docs/
    /// cip_perf_memoization_107.md`) to hold a few hundred entries per resolution
    /// even on the corpus's own worst-case fixtures -- nowhere near that ceiling.
    cache: HashMap<PairwiseCacheKey, BranchComparison>,
    /// Diagnostic only (never read to decide a comparison): how many
    /// `compare_ligands` calls were satisfied from `cache` instead of recomputed.
    pub cache_hits: u64,
    /// Diagnostic only: how many calls computed and inserted a fresh cache entry.
    pub cache_misses: u64,
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
            fractional_comparisons: 0,
            fractional_decisions: 0,
            cache: HashMap::new(),
            cache_hits: 0,
            cache_misses: 0,
        }
    }

    pub fn with_trace(trace: &'t mut ComparisonTrace) -> Self {
        Self {
            recursive_calls: 0,
            max_recursive_calls: 1_000_000,
            trace: Some(trace),
            ranking_parent: None,
            cache: HashMap::new(),
            cache_hits: 0,
            cache_misses: 0,
            fractional_comparisons: 0,
            fractional_decisions: 0,
        }
    }
}

/// `(atomic_number, isotope)` for a node's *effective* atom, Rule 1a's comparison key.
/// The atomic-number half is simply `node.atomic_number` -- already computed once at
/// digraph-construction time (see `crate::digraph`'s `atomic_number_for`), plain integer
/// for every node except a MANCUDE-affected `MultipleBondDuplicate`, which carries a
/// `Rational` (see `crate::mancude`'s module docs for why the owner's, not the
/// represented atom's, fraction is correct here).
///
/// The isotope half needs its own atom lookup and, for a MANCUDE-affected duplicate,
/// its own owner-vs-represented decision: `duplicated_atom` is a *specific* Kekulé-form-
/// dependent partner (exactly the value `atomic_number` deliberately stops depending on
/// for this same node), so reading isotope from it *would* reintroduce Kekulé-form
/// dependence into Rule 2 while Rule 1a stays invariant -- by the same argument that
/// motivates the atomic-number owner decision. Whenever `atomic_number` is `Rational`,
/// isotope is read from `source_atom` (the owner) instead. **Unverified against RDKit,
/// corpus-inert**, not "correct by the same reasoning" -- unlike the atomic-number
/// decision (checked against a concrete RDKit-reimplemented value), this branch has never
/// been exercised: 0 isotope-labeled atoms exist anywhere in the frozen 155-row corpus,
/// and the one isotope label in the full ~5,000-molecule verification corpus
/// (`[3H]` on a plain aliphatic ring, unrelated to any MANCUDE system) never reaches this
/// code path. Structurally consistent with the atomic-number decision, not empirically
/// confirmed by it. An ordinary (non-MANCUDE) duplicate keeps today's existing behavior
/// unchanged: `duplicated_atom`'s real isotope.
fn node_key(mol: &Molecule, node: &CipNode) -> (AtomicNumberKey, Option<u16>) {
    let isotope_atom = match node.kind {
        CipNodeKind::Atom { atom_idx } => Some(atom_idx),
        CipNodeKind::MultipleBondDuplicate {
            source_atom,
            duplicated_atom,
            ..
        } => Some(match node.atomic_number {
            AtomicNumberKey::Rational(_) => source_atom,
            AtomicNumberKey::Integral(_) => duplicated_atom,
        }),
        CipNodeKind::RingDuplicate { closure_atom, .. } => Some(closure_atom),
        CipNodeKind::ImplicitHydrogen => None,
    };
    let isotope = isotope_atom.and_then(|idx| mol.atom(idx).isotope);
    (node.atomic_number, isotope)
}

/// Rule 1a (atomic number, via `cmp_atomic_number_key` -- promotes a plain integer to a
/// `RationalAtomicNumber` when compared against a MANCUDE fraction) then Rule 2 (isotope:
/// `Some` beats `None`, higher isotope number beats lower).
fn cmp_key(a: &(AtomicNumberKey, Option<u16>), b: &(AtomicNumberKey, Option<u16>)) -> Ordering {
    match cmp_atomic_number_key(a.0, b.0) {
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

/// `key`, with a `Rational` atomic-number half collapsed to its integer part -- used only
/// to detect whether a fraction *decided* a comparison (see
/// `CompareContext::fractional_decisions`), never for a real ranking decision.
fn collapse_to_integer(key: (AtomicNumberKey, Option<u16>)) -> (AtomicNumberKey, Option<u16>) {
    let atomic_number = match key.0 {
        AtomicNumberKey::Integral(n) => AtomicNumberKey::Integral(n),
        AtomicNumberKey::Rational(r) => {
            AtomicNumberKey::Integral((r.numerator() / r.denominator()) as u8)
        }
    };
    (atomic_number, key.1)
}

/// `cmp_key`, instrumented to update [`CompareContext::fractional_comparisons`]/
/// [`CompareContext::fractional_decisions`] -- kept as a thin wrapper at this one call
/// site (not inside `cmp_key` itself, which stays pure and independently unit-tested, and
/// not inside `compare_by_level`, which stays fully generic over `K` with no
/// MANCUDE-specific knowledge).
fn cmp_key_instrumented(
    a: &(AtomicNumberKey, Option<u16>),
    b: &(AtomicNumberKey, Option<u16>),
    ctx: &mut CompareContext,
) -> Ordering {
    let ordering = cmp_key(a, b);
    let fractional_involved =
        matches!(a.0, AtomicNumberKey::Rational(_)) || matches!(b.0, AtomicNumberKey::Rational(_));
    if fractional_involved {
        ctx.fractional_comparisons += 1;
        if cmp_key(&collapse_to_integer(*a), &collapse_to_integer(*b)) != ordering {
            ctx.fractional_decisions += 1;
        }
    }
    ordering
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

fn rule1a2_slot_key(graph: &CipDigraph, slot: LevelSlot) -> (AtomicNumberKey, Option<u16>) {
    match slot {
        LevelSlot::Node(n) => node_key(graph.molecule(), graph.node(n)),
        LevelSlot::Phantom => (AtomicNumberKey::Integral(0), None),
    }
}

/// A real node's own ranked children, as plain `NodeId`s (a phantom slot's "children"
/// are always empty -- it never gets here since [`compare_by_level`] only calls this for
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

/// The result of one key-parameterized sphere-by-sphere walk (see [`compare_by_level`]):
/// either a definite winner (with the two keys that decided it, for tracing), or a full
/// tie -- every live position, at every sphere, matched exactly under this pass's key,
/// all the way to simultaneous exhaustion.
enum LevelOutcome<K> {
    Decided(BranchComparison, K, K),
    FullyTied,
}

/// Walk two ligand branches sphere by sphere under a single comparison rule, given as
/// `key_fn` (what to compare at each position) and `cmp_fn` (how to compare two keys --
/// kept separate from `K: Ord` so a rule whose priority order isn't a plain
/// ascending/descending sort on its own key type doesn't need a newtype wrapper).
/// Compares a whole level's keys, position by position, *before* descending into any
/// position's children -- a shallow difference at position 1 must win even if position 0
/// would only differ several spheres further down (see module docs, "Sphere-by-sphere").
/// Generic so a future sequence rule (Milestone 4's Rule 3/4/5) can reuse this walk
/// without duplicating the phantom-padding/level-advance machinery -- see module docs'
/// "Rule 1b" section for why a naive "run every rule as an unconditional second pass"
/// scheduling was tried and reverted (correct but ~2x cost for zero behavior change on
/// the current corpus).
fn compare_by_level<K>(
    graph: &mut CipDigraph,
    left: NodeId,
    right: NodeId,
    ctx: &mut CompareContext,
    key_fn: impl Fn(&CipDigraph, LevelSlot) -> K + Copy,
    cmp_fn: impl Fn(&K, &K, &mut CompareContext) -> Ordering + Copy,
) -> Result<LevelOutcome<K>, CipCompareError> {
    let mut left_level = vec![LevelSlot::Node(left)];
    let mut right_level = vec![LevelSlot::Node(right)];

    loop {
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

        let mut decided = None;
        for i in 0..n {
            let lk = key_fn(graph, left_level[i]);
            let rk = key_fn(graph, right_level[i]);
            match cmp_fn(&lk, &rk, ctx) {
                Ordering::Greater => {
                    decided = Some((BranchComparison::Higher, lk, rk));
                    break;
                }
                Ordering::Less => {
                    decided = Some((BranchComparison::Lower, lk, rk));
                    break;
                }
                Ordering::Equal => {}
            }
        }
        if let Some((outcome, lk, rk)) = decided {
            return Ok(LevelOutcome::Decided(outcome, lk, rk));
        }
        if n == 0 {
            return Ok(LevelOutcome::FullyTied);
        }

        // Whole level tied on this pass's key: expand to the next sphere. Each
        // position's children are ranked *within that one parent* (a separate,
        // locally-scoped comparison using the *full* two-pass comparator -- see
        // `rank_children` and module docs) before being placed into the next level;
        // per-position padding (not a global flatten) keeps "position i" well-defined on
        // both sides even when a tied pair's substituent counts differ.
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
    }
}

/// Compare two ligand branches under Rule 1a/2 (see module docs for why Rule 1b is not
/// wired in as a second pass here despite being correctly implementable -- it was tried
/// in Milestone 3A and found to change zero outputs at ~2x cost, since Rule 1a's own
/// child-count check always resolves a duplicate-vs-real-atom position first).
pub fn compare_ligands(
    graph: &mut CipDigraph,
    left: NodeId,
    right: NodeId,
    ctx: &mut CompareContext,
) -> Result<BranchComparison, CipCompareError> {
    let left_kind = graph.node(left).kind;
    let right_kind = graph.node(right).kind;
    let depth = graph.node(left).depth;

    // Issue #107: same-(left,right) pairwise memoization, scoped to this
    // `CompareContext` alone (see `PairwiseCacheKey`'s doc for the key design and
    // why it's safe). Checked/populated here, immediately around the actual
    // recursive comparison -- the trace still records one `DecisionStep` per
    // logical `compare_ligands` call regardless of hit/miss, so anything reading
    // the trace (this crate's own diagnostics, e.g. `examples/
    // rank_children_heavy_tail_diagnosis.rs`) sees a complete, faithful log of
    // every comparison a caller actually asked for, not just the ones that did
    // fresh recursive work.
    let (cache_key, already_canonical) = PairwiseCacheKey::canonical(graph, left, right);
    let (outcome, rule) = if let Some(&cached) = ctx.cache.get(&cache_key) {
        ctx.cache_hits += 1;
        let outcome = if already_canonical {
            cached
        } else {
            invert(cached)
        };
        (outcome, "1a/2 (cached)".to_string())
    } else {
        let (outcome, rule) = match compare_by_level(
            graph,
            left,
            right,
            ctx,
            rule1a2_slot_key,
            cmp_key_instrumented,
        )? {
            LevelOutcome::Decided(c, lk, rk) => {
                let rule = if ctx.trace.is_some() {
                    format!("1a/2 ({lk:?} vs {rk:?})")
                } else {
                    String::new()
                };
                (c, rule)
            }
            LevelOutcome::FullyTied => (BranchComparison::Equal, "leaf".to_string()),
        };
        ctx.cache_misses += 1;
        let canonical_outcome = if already_canonical {
            outcome
        } else {
            invert(outcome)
        };
        ctx.cache.insert(cache_key, canonical_outcome);
        (outcome, rule)
    };

    let ranking_parent = ctx.ranking_parent;
    if let Some(trace) = ctx.trace.as_deref_mut() {
        trace.decisions.push(DecisionStep {
            depth,
            left_node: left,
            right_node: right,
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
