//! Digraph edge types.

use chematic_core::BondOrder;

use crate::node::NodeId;

/// Index of a [`CipEdge`] in a [`crate::digraph::CipDigraph`]'s arena.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EdgeId(pub u32);

/// The edge connecting a node to the child that was reached to create it.
///
/// Kept as its own type (rather than folding `bond_order` onto [`crate::node::CipNode`])
/// because Milestone 2's ranking comparator is expected to need bond-order-sensitive
/// comparison along an edge, not just at a node -- if this stays a pure
/// `(parent, child)` pass-through through Milestone 2, it's a fold-into-node candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipEdge {
    pub id: EdgeId,
    pub parent: NodeId,
    pub child: NodeId,
    pub bond_order: BondOrder,
}
