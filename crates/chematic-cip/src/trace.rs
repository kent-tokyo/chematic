//! Human-readable per-comparison decision traces.
//!
//! Accuracy percentages alone don't explain *why* a comparator picked one branch over
//! another -- this project's own history bears that out (the `d0e726b` root cause took
//! hand-tracing exact substituent orders to find). A [`ComparisonTrace`], recorded by
//! [`crate::compare::compare_ligands`] when a caller opts in via
//! [`crate::compare::CompareContext::with_trace`], is meant to make future mismatch
//! triage (Milestone 3's aromatic-ring-adjacent bucket, in particular) tractable
//! without re-deriving the decision path by hand each time.

use crate::compare::BranchComparison;
use crate::node::NodeId;

/// One step in a comparison's decision path.
///
/// Every [`crate::compare::compare_ligands`] call -- whether it's the comparison a
/// caller asked for directly, or a sibling-pairwise sub-comparison
/// [`crate::compare::rank_children`] makes internally to order *some other* node's own
/// children -- pushes exactly one `DecisionStep`. Without `ranking_parent`, those two
/// kinds of step interleave indistinguishably in trace order (found the hard way while
/// root-causing a wrong-but-confident corpus case: a step reading "atom(0) vs
/// dup-multi(2) -> Lower" gives no way to tell whether that's part of the comparison you
/// asked about, or a sub-step ranking some unrelated node's children three levels
/// removed). `ranking_parent` names which node's children this particular step's `left`
/// and `right` are siblings under, making the interleaving legible again.
#[derive(Debug, Clone)]
pub struct DecisionStep {
    pub depth: u32,
    /// The two compared nodes' own `NodeId`s -- not needed by any display/reporting
    /// logic (which reads `left_kind`/`right_kind` instead), but required by
    /// diagnostics that need real node identity, e.g. detecting whether the same
    /// `NodeId` pair (or, via [`crate::digraph::CipDigraph::branch_signature`], the
    /// same *isomorphic* subtree pair) gets re-compared redundantly across a single
    /// resolution. See issue #107 (CIP-Perf-A1).
    pub left_node: NodeId,
    pub right_node: NodeId,
    pub left_kind: String,
    pub right_kind: String,
    pub outcome: BranchComparison,
    /// Which rule decided this step, and the specific keys compared: `"1a/2 (<left key>
    /// vs <right key>)"`, or `"leaf"` (both sides fully tied under Rule 1a/2 -- would
    /// need Rule 1b/3+ to go further; see `compare.rs`'s module docs for why Rule 1b
    /// isn't wired in as an automatic second pass).
    pub rule: String,
    /// The node whose children `left`/`right` are siblings under, i.e. whose
    /// `rank_children` call this step belongs to. `None` only if the digraph root itself
    /// somehow appears as a comparison operand (not expected in practice -- every
    /// comparison happens between two siblings under a real parent).
    pub ranking_parent: Option<NodeId>,
}

impl core::fmt::Display for DecisionStep {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "depth {}: {} vs {} -> {:?} (rule {}) [ranking children of {}]",
            self.depth,
            self.left_kind,
            self.right_kind,
            self.outcome,
            self.rule,
            self.ranking_parent
                .map(|n| n.0.to_string())
                .unwrap_or_else(|| "?".to_string()),
        )
    }
}

/// The full decision path recorded for one top-level [`compare_ligands`] call
/// (including every nested recursive comparison it made along the way).
///
/// [`compare_ligands`]: crate::compare::compare_ligands
#[derive(Debug, Clone)]
pub struct ComparisonTrace {
    pub left_root: NodeId,
    pub right_root: NodeId,
    pub decisions: Vec<DecisionStep>,
}

impl ComparisonTrace {
    pub fn new(left_root: NodeId, right_root: NodeId) -> Self {
        Self {
            left_root,
            right_root,
            decisions: Vec::new(),
        }
    }
}

impl core::fmt::Display for ComparisonTrace {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        writeln!(
            f,
            "comparing root {} vs root {}:",
            self.left_root.0, self.right_root.0
        )?;
        for step in &self.decisions {
            writeln!(f, "  {step}")?;
        }
        Ok(())
    }
}
