//! Pairwise substituent ranking at an arbitrary centre (double-bond ends).
//!
//! [`assign_cip_accurate_experimental`](crate::assign_cip_accurate_experimental)
//! only labels tetrahedral centres. E/Z labels need the same hierarchical
//! digraph ranking (Rules 1a/2 with MANCUDE duplicate atomic numbers) applied
//! to the two substituents of each double-bond end; [`SubstituentRanker`]
//! exposes exactly that, computing the molecule's Kekulé form and MANCUDE
//! context once for all ends.

use std::cmp::Ordering;

use chematic_core::{AtomIdx, Molecule};

use crate::budget::CipBudget;
use crate::compare::{CipCompareError, CompareContext, rank_children};
use crate::digraph::CipDigraph;
use crate::mancude::{MancudeContext, prepare_kekule_form};
use crate::node::CipNodeKind;

/// Ranks substituents of one centre against each other with the accurate
/// engine's comparator. Build once per molecule.
pub struct SubstituentRanker {
    kekule: Option<(Molecule, MancudeContext)>,
    budget: CipBudget,
}

impl SubstituentRanker {
    /// Prepare ranking for `mol` (Kekulé form and MANCUDE context are computed
    /// once; when `mol` cannot be kekulized the plain digraph is used, as the
    /// tetrahedral engine does).
    pub fn new(mol: &Molecule, budget: CipBudget) -> Self {
        Self {
            kekule: prepare_kekule_form(mol).ok(),
            budget,
        }
    }

    /// CIP priority of substituent `a` relative to substituent `b`, both
    /// neighbours of `center`: `Greater` when `a` ranks higher. `Ok(None)`
    /// when the two are tied under the implemented rules or either atom is not
    /// a neighbour of `center` — never a guess.
    pub fn compare(
        &self,
        mol: &Molecule,
        center: AtomIdx,
        a: AtomIdx,
        b: AtomIdx,
    ) -> Result<Option<Ordering>, CipCompareError> {
        let mut graph = match &self.kekule {
            Some((kekule_mol, ctx)) => {
                CipDigraph::new_with_mancude(kekule_mol, center, self.budget, ctx)
            }
            None => CipDigraph::new(mol, center, self.budget),
        }
        .map_err(CipCompareError::Digraph)?;
        let root = graph.root();
        let children = graph
            .expand_children(root)
            .map_err(CipCompareError::Digraph)?;
        let find = |graph: &CipDigraph, atom: AtomIdx| {
            children.iter().copied().find(|&id| {
                matches!(graph.node(id).kind, CipNodeKind::Atom { atom_idx } if atom_idx == atom)
            })
        };
        let (Some(node_a), Some(node_b)) = (find(&graph, a), find(&graph, b)) else {
            return Ok(None);
        };
        let mut ctx = CompareContext::new();
        let groups = rank_children(&mut graph, &children, &mut ctx)?;
        let group_of = |node| groups.iter().position(|g| g.contains(&node));
        match (group_of(node_a), group_of(node_b)) {
            (Some(ga), Some(gb)) if ga != gb => Ok(Some(gb.cmp(&ga))),
            _ => Ok(None),
        }
    }
}
