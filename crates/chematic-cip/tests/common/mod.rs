//! Shared helpers for `tests/mancude_*.rs` -- a `tests/common/mod.rs` file is the
//! standard way to share code between separate integration-test binaries in Rust (unlike
//! a plain `tests/foo.rs`, this path is not itself compiled as a test crate).

use chematic_cip::{CipBudget, CipDigraph, CipNodeKind, NodeId, assign_cip_accurate_experimental};
use chematic_core::{AtomIdx, CipCode, Molecule};

pub fn code_str(code: CipCode) -> &'static str {
    match code {
        CipCode::R => "R",
        CipCode::S => "S",
        CipCode::E => "E",
        CipCode::Z => "Z",
    }
}

/// True if `node`'s subtree (up to `budget` nodes) contains an aromatic ring atom that
/// is fully substituted (3 real molecular neighbors, no hydrogen) but has only 2 digraph
/// children -- the "missing mancude duplicate" signature `uncharacterized_diagnosis.rs`
/// (Milestone 2.5) uses to diagnose its 2 `BucketMisclassified` cases.
pub fn subtree_has_undercounted_aromatic(
    mol: &Molecule,
    graph: &mut CipDigraph,
    node: NodeId,
    budget: &mut usize,
) -> bool {
    *budget = budget.saturating_sub(1);
    if *budget == 0 {
        return false;
    }
    let atom_idx = match graph.node(node).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        CipNodeKind::MultipleBondDuplicate { .. }
        | CipNodeKind::RingDuplicate { .. }
        | CipNodeKind::ImplicitHydrogen => {
            return false;
        }
    };
    let Ok(children) = graph.expand_children(node) else {
        return false;
    };
    let no_hydrogen_child = !children
        .iter()
        .any(|&c| matches!(graph.node(c).kind, CipNodeKind::ImplicitHydrogen));
    if mol.atom(atom_idx).aromatic
        && mol.neighbors(atom_idx).count() == 3
        && children.len() == 2
        && no_hydrogen_child
    {
        return true;
    }
    children
        .into_iter()
        .any(|c| subtree_has_undercounted_aromatic(mol, graph, c, budget))
}

/// Mirrors `uncharacterized_diagnosis.rs`'s exact gating: the aromatic-undercounting
/// check only applies to cases `AccurateExperimental` got *wrong* (produced an incorrect
/// assignment), never to *tied* cases (a different, unrelated cause,
/// `NeedsLaterSequenceRule`). Applying it unconditionally over-matches: many tied cases
/// also happen to contain a 3-neighbor aromatic atom somewhere, incidentally, without
/// that being why they're tied.
pub fn is_bucket_misclassified(mol: &Molecule, atom_idx: u32, modern: &str) -> bool {
    let Ok(accurate) = assign_cip_accurate_experimental(mol, CipBudget::default_budget()) else {
        return false;
    };
    let assigned = accurate
        .assignments
        .iter()
        .find(|(i, _)| i.0 == atom_idx)
        .map(|(_, code)| *code);
    let Some(code) = assigned else {
        return false; // tied, budget-exceeded, or no_assign -- not a "wrong" case
    };
    if code_str(code) == modern {
        return false; // correctly matched, not wrong
    }

    let Ok(mut graph) = CipDigraph::new(mol, AtomIdx(atom_idx), CipBudget::default_budget()) else {
        return false;
    };
    let Ok(root_children) = graph.expand_children(graph.root()) else {
        return false;
    };
    root_children.iter().any(|&n| {
        let mut budget = 5_000;
        subtree_has_undercounted_aromatic(mol, &mut graph, n, &mut budget)
    })
}

/// True iff this corpus row is in Milestone 3B-0's 98-case scope: the 96
/// `aromatic_mancude`-bucket cases plus the 2 `uncharacterized` cases mis-tagged by the
/// corpus's own `classify_bucket()` heuristic (Milestone 2.5's `BucketMisclassified`).
pub fn in_mancude_scope(mol: &Molecule, atom_idx: u32, bucket: Option<&str>, modern: &str) -> bool {
    bucket == Some("aromatic_mancude")
        || (bucket.is_none() && is_bucket_misclassified(mol, atom_idx, modern))
}
