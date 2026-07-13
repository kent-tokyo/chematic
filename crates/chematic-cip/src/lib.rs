//! `chematic-cip` — hierarchical digraph foundation for an accurate CIP
//! (Cahn-Ingold-Prelog) stereochemistry engine.
//!
//! **Milestone 1 only** (see `docs/cip_accurate_rfc.md` at the workspace root): this
//! crate builds a structural digraph rooted at a candidate stereocenter. It does
//! **not** rank substituents, assign R/S, or otherwise replace
//! `chematic_chem::assign_cip()`, which is unaffected by this crate's existence.
//!
//! The motivation is a real, proven limitation in `chematic-chem`'s existing engine:
//! `cip_branch_spheres`/`compare_branches` pool every atom at a given BFS depth into
//! one sorted multiset and compare shell-by-shell. That's an approximation of CIP, not
//! the real algorithm (which recursively compares branch-by-branch, following the
//! highest-priority sub-branch first) -- and the approximation is unsafe to patch
//! piecemeal: extending a correct double-bond duplication rule to triple bonds went net
//! negative (16 newly-wrong stereocenters vs. 1 newly-fixed) because the pooled
//! representation had already discarded the branch/provenance information a correct
//! comparison needs. This crate's [`digraph::CipDigraph`] keeps that information as
//! explicit, provenance-carrying nodes instead.
#![forbid(unsafe_code)]

pub mod assign;
pub mod budget;
pub mod compare;
pub mod debug;
pub mod digraph;
pub mod edge;
pub mod node;
pub mod trace;

#[cfg(test)]
mod tests;

pub use assign::{AccurateCipAssignment, SkipReason, assign_cip_accurate_experimental};
pub use budget::CipBudget;
pub use compare::{
    BranchComparison, CipCompareError, CompareContext, compare_ligands, rank_children,
};
pub use digraph::{CipDigraph, DigraphExpander};
pub use edge::{CipEdge, EdgeId};
pub use node::{CipNode, CipNodeKind, NodeId};
pub use trace::{ComparisonTrace, DecisionStep};

/// Errors from digraph construction/expansion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CipError {
    /// Expansion exceeded a [`CipBudget`] limit (node count, depth, or expansion
    /// count). Never returned as a substitute for a wrong-but-plausible answer --
    /// exceeding budget always surfaces as this error, not a silently truncated or
    /// guessed structure.
    BudgetExceeded { reason: String },
}

impl core::fmt::Display for CipError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CipError::BudgetExceeded { reason } => {
                write!(f, "CIP digraph budget exceeded: {reason}")
            }
        }
    }
}

impl std::error::Error for CipError {}
