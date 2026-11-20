//! Pairwise substituent ranking at an arbitrary centre (double-bond ends).
//!
//! [`assign_cip_accurate_experimental`](crate::assign_cip_accurate_experimental)
//! only labels tetrahedral centres. E/Z labels need the same hierarchical
//! digraph ranking (Rules 1a/2 with MANCUDE and Rule 3 for configured alkene ties) applied
//! to the two substituents of each double-bond end; [`SubstituentRanker`]
//! exposes exactly that, computing the molecule's Kekulé form and MANCUDE
//! context once for all ends.

use std::cmp::Ordering;

use chematic_core::{AtomIdx, Molecule};

use crate::budget::CipBudget;
use crate::compare::{BranchComparison, CipCompareError, CompareContext, compare_ligands};
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
        self.compare_nodes(
            mol,
            center,
            CipNodeKind::Atom { atom_idx: a },
            CipNodeKind::Atom { atom_idx: b },
        )
    }

    /// CIP priority of an explicit ligand relative to the implicit hydrogen
    /// on `center`. `Ok(None)` if a ligand is absent or the priorities tie.
    pub fn compare_to_implicit_hydrogen(
        &self,
        mol: &Molecule,
        center: AtomIdx,
        atom: AtomIdx,
    ) -> Result<Option<Ordering>, CipCompareError> {
        self.compare_nodes(
            mol,
            center,
            CipNodeKind::Atom { atom_idx: atom },
            CipNodeKind::ImplicitHydrogen,
        )
    }

    fn compare_nodes(
        &self,
        mol: &Molecule,
        center: AtomIdx,
        a: CipNodeKind,
        b: CipNodeKind,
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
            .expand_child_ids(root)
            .map_err(CipCompareError::Digraph)?;
        let find = |kind| {
            children
                .iter()
                .copied()
                .find(|&id| graph.node(id).kind == kind)
        };
        let (Some(node_a), Some(node_b)) = (find(a), find(b)) else {
            return Ok(None);
        };
        let mut ctx = CompareContext::new();
        match compare_ligands(&mut graph, node_a, node_b, &mut ctx)? {
            BranchComparison::Higher => Ok(Some(Ordering::Greater)),
            BranchComparison::Lower => Ok(Some(Ordering::Less)),
            BranchComparison::Equal => {
                let order = crate::rule3::compare(&mut graph, node_a, node_b, &mut ctx)?;
                Ok((order != Ordering::Equal).then_some(order))
            }
            BranchComparison::Unresolved => Ok(None),
        }
    }
}
