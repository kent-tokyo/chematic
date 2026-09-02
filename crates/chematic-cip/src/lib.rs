//! `chematic-cip` — low-level, hierarchical-digraph CIP (Cahn-Ingold-Prelog)
//! stereochemistry engine for chematic.
//!
//! **Most applications should use [`chematic_chem::assign_cip_with_mode`] rather
//! than depending on this crate directly.** That's the stable, supported entry point
//! (opt-in `CipMode::Accurate`, merged with legacy E/Z/allene handling; see
//! `docs/rfcs/cip_accurate_rfc.md`'s Milestone 5A). This crate is the low-level engine
//! behind it, published separately so `chematic-chem` can depend on it normally.
//! The accurate engine remains experimental and may receive breaking API revisions
//! in a future 0.x minor release.
//!
//! Full milestone history in `docs/rfcs/cip_accurate_rfc.md` at the workspace root: this
//! crate now builds a provenance-carrying digraph, ranks substituents (Rules 1a/1b/2/4b/5),
//! and assigns R/S — the current residual-corpus report estimates 99.74% agreement
//! against modern RDKit `rdCIPLabeler` on the project's full validation corpus,
//! with the remaining phosphorus cases kept unresolved when independent respellings
//! disagree. It does not replace
//! `chematic_chem::assign_cip()`'s default (legacy) path, which is unaffected.
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
mod auxiliary;
pub mod budget;
pub mod compare;
pub mod debug;
pub mod digraph;
pub mod digraph_diff;
pub mod edge;
pub mod mancude;
pub mod node;
pub mod rational;
mod resolver;
mod rule4b;
pub mod trace;

#[cfg(test)]
mod tests;

pub use assign::{
    AccurateCipAssignment, SkipReason, assign_cip_accurate_experimental,
    assign_cip_accurate_experimental_without_mancude,
};
pub use budget::CipBudget;
pub use compare::{
    BranchComparison, CipCompareError, CompareContext, compare_ligands, rank_children,
};
pub use digraph::{CipDigraph, DigraphExpander};
pub use edge::{CipEdge, EdgeId};
pub use mancude::{
    MancudeBudget, MancudeComponentId, MancudeContext, MancudeError, effective_atomic_number,
    enumerate_kekule_matchings, prepare_kekule_form,
};
pub use node::{CipNode, CipNodeKind, NodeId};
pub use rational::{AtomicNumberKey, RationalAtomicNumber};
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
