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
#[derive(Debug, Clone)]
pub struct DecisionStep {
    pub depth: u32,
    pub left_kind: String,
    pub right_kind: String,
    pub outcome: BranchComparison,
    /// Which rule decided this step: `"1a/2"` (own-key comparison), `"leaf"` (both
    /// sides childless), or `"children"` (recursed into ranked children).
    pub rule: &'static str,
}

impl core::fmt::Display for DecisionStep {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "depth {}: {} vs {} -> {:?} (rule {})",
            self.depth, self.left_kind, self.right_kind, self.outcome, self.rule
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
