//! Human-readable and JSON dumps of a [`CipDigraph`].
//!
//! CIP bugs are hard to reason about without seeing the actual path taken to reach a
//! node -- this project's own history bears that out (the `d0e726b` root-cause hunt
//! required hand-tracing exact substituent orders through a debugger). A debug
//! representation is part of Milestone 1, not a later add-on, so every future milestone
//! can dump and diff digraphs from day one.
//!
//! The JSON dump is hand-rolled rather than built on `serde_json` (a dev-only
//! dependency of this crate, used solely by the residual-corpus test to parse its
//! *input* fixture) -- the schema here is small and fixed, so a real JSON library adds
//! a runtime dependency for no real benefit. Both dumps assume the subtree has already
//! been fully expanded (e.g. via [`CipDigraph::expand_all`]); they don't expand
//! anything themselves.

use crate::digraph::CipDigraph;
use crate::node::{CipNodeKind, NodeId};

/// Indented tree dump, depth-first, children in arena order.
pub fn to_tree_string(digraph: &CipDigraph) -> String {
    let mut out = String::new();
    let root = digraph.root();
    out.push_str(&kind_label(digraph.node(root).kind));
    out.push('\n');
    write_children(digraph, root, "", &mut out);
    out
}

fn write_children(digraph: &CipDigraph, id: NodeId, prefix: &str, out: &mut String) {
    let children = children_of(digraph, id);
    for (i, &child) in children.iter().enumerate() {
        let last = i + 1 == children.len();
        let branch = if last { "└─ " } else { "├─ " };
        out.push_str(prefix);
        out.push_str(branch);
        out.push_str(&kind_label(digraph.node(child).kind));
        out.push('\n');
        let child_prefix = if last {
            format!("{prefix}   ")
        } else {
            format!("{prefix}│  ")
        };
        write_children(digraph, child, &child_prefix, out);
    }
}

fn children_of(digraph: &CipDigraph, id: NodeId) -> Vec<NodeId> {
    digraph
        .nodes()
        .iter()
        .filter(|n| n.parent == Some(id))
        .map(|n| n.id)
        .collect()
}

fn kind_label(kind: CipNodeKind) -> String {
    match kind {
        CipNodeKind::Atom { atom_idx } => format!("atom={}", atom_idx.0),
        CipNodeKind::MultipleBondDuplicate {
            source_atom,
            duplicated_atom,
            bond_order,
        } => format!(
            "duplicate(multiple-bond, source={}, atom={}, order={})",
            source_atom.0, duplicated_atom.0, bond_order
        ),
        CipNodeKind::RingDuplicate {
            source_atom,
            closure_atom,
        } => format!(
            "ring-duplicate(source={}, atom={})",
            source_atom.0, closure_atom.0
        ),
        CipNodeKind::ImplicitHydrogen => "implicit-H".to_string(),
    }
}

/// JSON dump: `{"root_atom": N, "nodes": [...]}`, one object per node, in
/// arena/[`NodeId`] order -- never a `HashMap`/`HashSet` iteration order, so this is
/// deterministic by construction, not by luck of an unspecified hasher.
pub fn to_json(digraph: &CipDigraph) -> String {
    let root_atom = match digraph.node(digraph.root()).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx.0,
        _ => unreachable!("digraph root is always an Atom node"),
    };

    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("  \"root_atom\": {root_atom},\n"));
    out.push_str("  \"nodes\": [\n");
    let nodes = digraph.nodes();
    for (i, node) in nodes.iter().enumerate() {
        out.push_str("    ");
        out.push_str(&node_json(node.id, node.kind, node.parent, node.depth));
        if i + 1 < nodes.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ]\n}");
    out
}

fn node_json(id: NodeId, kind: CipNodeKind, parent: Option<NodeId>, depth: u32) -> String {
    let parent_json = match parent {
        Some(p) => p.0.to_string(),
        None => "null".to_string(),
    };
    match kind {
        CipNodeKind::Atom { atom_idx } => format!(
            "{{\"id\": {}, \"kind\": \"atom\", \"atom_idx\": {}, \"parent\": {}, \"depth\": {}}}",
            id.0, atom_idx.0, parent_json, depth
        ),
        CipNodeKind::MultipleBondDuplicate {
            source_atom,
            duplicated_atom,
            bond_order,
        } => format!(
            "{{\"id\": {}, \"kind\": \"multiple_bond_duplicate\", \"source_atom\": {}, \"duplicated_atom\": {}, \"bond_order\": {}, \"parent\": {}, \"depth\": {}}}",
            id.0, source_atom.0, duplicated_atom.0, bond_order, parent_json, depth
        ),
        CipNodeKind::RingDuplicate {
            source_atom,
            closure_atom,
        } => format!(
            "{{\"id\": {}, \"kind\": \"ring_duplicate\", \"source_atom\": {}, \"closure_atom\": {}, \"parent\": {}, \"depth\": {}}}",
            id.0, source_atom.0, closure_atom.0, parent_json, depth
        ),
        CipNodeKind::ImplicitHydrogen => format!(
            "{{\"id\": {}, \"kind\": \"implicit_hydrogen\", \"parent\": {}, \"depth\": {}}}",
            id.0, parent_json, depth
        ),
    }
}
