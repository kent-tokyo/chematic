//! Expansion budget: an explicit backstop against pathological node-count blow-up.
//!
//! The digraph's *finiteness* (every root-to-leaf path terminates) is guaranteed
//! structurally by the ancestor-path ring-closure rule in [`crate::digraph`], not by
//! this budget -- a molecule with many fused/bridged rings can still have a very large
//! *total node count* even though every individual path is finite (many branches, each
//! short). The budget exists to catch that case and fail loudly
//! (`CipError::BudgetExceeded`) rather than let expansion run away, and to give an
//! honest "I don't know" instead of ever silently truncating or guessing a structure.

/// Limits on digraph expansion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CipBudget {
    pub max_nodes: usize,
    pub max_depth: usize,
    pub max_expansions: usize,
}

impl CipBudget {
    /// A generous default, sized for real molecules rather than pathological inputs.
    pub fn default_budget() -> Self {
        Self {
            max_nodes: 100_000,
            max_depth: 64,
            max_expansions: 100_000,
        }
    }
}
