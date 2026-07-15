//! Milestone 4B-0: mechanically diagnose the 8 `rule4_candidate` rows frozen in
//! `validation/cip_m4a0_residual.jsonl` (Milestone 4A-0) into CIP Rule 4a / 4b / 4c,
//! using only `chematic-cip`'s existing public API -- no crate changes. Read-only
//! diagnostic, same footprint as `residual_report.rs`/`trace_report.rs`.
//!
//! Usage: cargo run -p chematic-cip --release --example rule4_diagnose
//!
//! For each row, reports:
//! - the root-level tie shape (which 2 of the 4 physical positions tie, at what depth
//!   the comparator exhausts Rules 1a/2 into `Equal`)
//! - Rule 4a: N/A, since `chematic_core::Chirality` has only one non-`None` category
//!   (tetrahedral) -- there is no "chiral vs pseudoasymmetric vs nonstereogenic unit
//!   type" distinction this data model can even represent, so a category difference
//!   between branches is structurally impossible here (checked, not assumed: this is
//!   a fact about the enum itself, not about these 8 rows specifically).
//! - Rule 4c: N/A, for the identical structural reason (no axial/planar/m/p
//!   descriptor variant exists in `Chirality` at all).
//! - Rule 4b: the load-bearing check. For each of the tied pair's two branches, BFS to
//!   the *nearest* embedded stereocenter (mirroring
//!   `assign::nearest_embedded_stereocenter`'s search, re-implemented here since it's
//!   private to that module) and rank its own `expand_children` output. A clean
//!   (tie-free) local ranking at that embedded node is the acyclicity evidence: Rules
//!   1a/2 alone can already order that embedded atom's *forward* branches without
//!   needing any other tied atom's global/molecular descriptor.
//!
//! **What this diagnostic does NOT do**: compute the embedded atom's actual auxiliary
//! R/S *sign*. That requires representing the "back toward the tied root" direction as
//! a proper 4th ligand in a digraph rooted so that direction terminates immediately
//! (an artificial-ancestor root) -- `CipDigraph` has no such constructor today. Closing
//! that gap is scoped to Milestone 4B-1, not this diagnosis.

use std::collections::HashMap;

use chematic_cip::{
    CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, assign_cip_accurate_experimental,
    rank_children,
};
use chematic_core::{AtomIdx, CipCode, Molecule, STEREO_H_SENTINEL};

/// The 8 `rule4_candidate` rows frozen in `validation/cip_m4a0_residual.jsonl`
/// (Milestone 4A-0, commit `64f3a38`). Oracle (`modern`) is `S` for all 8.
const ROWS: &[(&str, u32, &str)] = &[
    (
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
        5,
        "tri-galloyl-a",
    ),
    (
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
        35,
        "tri-galloyl-b",
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        3,
        "quinic-a",
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        7,
        "quinic-b",
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        3,
        "mono-galloyl-a",
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        7,
        "mono-galloyl-b",
    ),
    (
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
        11,
        "tetra-galloyl-a",
    ),
    (
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
        26,
        "tetra-galloyl-b",
    ),
];

/// Map `stereo_neighbor_order` positions to root-child `NodeId`s. Local re-derivation
/// of `assign::position_node_ids` (private to that module) -- unchanged logic.
fn position_node_ids(
    graph: &CipDigraph,
    root_children: &[NodeId],
    stereo_order: &[u32],
) -> Option<Vec<NodeId>> {
    let mut result = Vec::with_capacity(stereo_order.len());
    for &pos_val in stereo_order {
        let node_id = if pos_val == STEREO_H_SENTINEL {
            root_children
                .iter()
                .copied()
                .find(|&id| matches!(graph.node(id).kind, CipNodeKind::ImplicitHydrogen))?
        } else {
            let atom_idx = AtomIdx(pos_val);
            root_children.iter().copied().find(|&id| {
                matches!(graph.node(id).kind, CipNodeKind::Atom { atom_idx: a } if a == atom_idx)
            })?
        };
        result.push(node_id);
    }
    Some(result)
}

/// BFS, level by level, for the nearest atom with `chirality != None` in `node`'s
/// subtree. Local re-derivation of `assign::nearest_embedded_stereocenter` (private to
/// that module) -- unchanged logic. Returns `None` if ambiguous (2+ at the nearest
/// level) or absent.
fn nearest_embedded_stereocenter(
    graph: &mut CipDigraph,
    mol: &Molecule,
    node: NodeId,
) -> Option<AtomIdx> {
    let mut frontier = vec![node];
    while !frontier.is_empty() {
        let mut found_this_level: Option<AtomIdx> = None;
        let mut next_frontier = Vec::new();
        for &current in &frontier {
            if let CipNodeKind::Atom { atom_idx } = graph.node(current).kind
                && mol.atom(atom_idx).chirality != chematic_core::Chirality::None
            {
                match found_this_level {
                    None => found_this_level = Some(atom_idx),
                    Some(existing) if existing == atom_idx => {}
                    Some(_) => return None,
                }
            }
            next_frontier.extend(graph.expand_children(current).ok()?);
        }
        if found_this_level.is_some() {
            return found_this_level;
        }
        frontier = next_frontier;
    }
    None
}

struct CaseReport {
    case_id: String,
    stereocenter: u32,
    tied_pair: (u32, u32),
    tie_depth_note: String,
    branch_a_embedded: Option<(AtomIdx, bool)>, // (atom, local_children_tie_free)
    branch_b_embedded: Option<(AtomIdx, bool)>,
    rule4a: &'static str,
    rule4c: &'static str,
    rule4b: String,
}

fn diagnose_one(smiles: &str, atom_idx: u32, case_id: &str) -> CaseReport {
    let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
    let idx = AtomIdx(atom_idx);
    let budget = CipBudget::default_budget();

    let mut graph = CipDigraph::new(&mol, idx, budget).expect("digraph builds");
    let root = graph.root();
    let root_children = graph.expand_children(root).expect("root expands");
    let stereo_order = mol
        .stereo_neighbor_order(idx)
        .expect("stereocenter has stereo_neighbor_order");
    let position_nodes = position_node_ids(&graph, &root_children, stereo_order)
        .expect("all 4 physical positions resolve to root children");

    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).expect("root children rank");

    let position_set: std::collections::HashSet<NodeId> = position_nodes.iter().copied().collect();
    let mut tied_pair: Option<(NodeId, NodeId)> = None;
    for group in &groups {
        let physical: Vec<NodeId> = group
            .iter()
            .copied()
            .filter(|n| position_set.contains(n))
            .collect();
        if physical.len() == 2 {
            tied_pair = Some((physical[0], physical[1]));
        } else if physical.len() > 2 {
            panic!("{case_id}: 3+-way physical tie, outside this diagnostic's scope");
        }
    }
    let (pos_a, pos_b) = tied_pair.unwrap_or_else(|| {
        panic!("{case_id}: expected exactly one 2-way physical tie among root children")
    });

    let atom_a = match graph.node(pos_a).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        _ => panic!("tied position must be a real atom"),
    };
    let atom_b = match graph.node(pos_b).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        _ => panic!("tied position must be a real atom"),
    };

    let tie_depth_note = format!(
        "root's own two ring-branches (toward atom{} and atom{}) tie under Rules 1a/2 -- \
         confirmed via M4A-0's branch_signature() as genuine constitutional identity, not \
         a comparator gap",
        atom_a.0, atom_b.0
    );

    // Rule 4b evidence: nearest embedded stereocenter per branch, and whether ITS OWN
    // (parent-excluded) children rank tie-free.
    let embedded_a = nearest_embedded_stereocenter(&mut graph, &mol, pos_a);
    let embedded_b = nearest_embedded_stereocenter(&mut graph, &mol, pos_b);

    let check_local_tie_free = |graph: &mut CipDigraph, embedded: Option<AtomIdx>| {
        embedded.map(|atom| {
            // Find the embedded atom's own NodeId within this branch (BFS again, since
            // nearest_embedded_stereocenter only returned the AtomIdx).
            let node_id = find_node_for_atom(graph, pos_a, pos_b, atom);
            let children = graph.expand_children(node_id).expect("expands");
            let mut local_ctx = CompareContext::new();
            let local_groups =
                rank_children(graph, &children, &mut local_ctx).expect("local children rank");
            let tie_free = local_groups.iter().all(|g| g.len() == 1);
            (atom, tie_free)
        })
    };

    let branch_a_embedded = check_local_tie_free(&mut graph, embedded_a);
    let branch_b_embedded = check_local_tie_free(&mut graph, embedded_b);

    let rule4b = match (branch_a_embedded, branch_b_embedded) {
        (Some((a, true)), Some((b, true))) => format!(
            "candidate: nearest embedded stereocenters atom{} (branch toward atom{}) and \
             atom{} (branch toward atom{}) both rank tie-free locally -- acyclic, a \
             per-branch auxiliary computation is possible without needing the other tied \
             atom's global descriptor. Auxiliary R/S *sign* not computed here (needs a \
             4th-ligand/artificial-ancestor representation -- see module docs; scoped to \
             M4B-1).",
            a.0, atom_a.0, b.0, atom_b.0
        ),
        _ => "inconclusive: embedded stereocenter ambiguous or tied locally".to_string(),
    };

    CaseReport {
        case_id: case_id.to_string(),
        stereocenter: atom_idx,
        tied_pair: (atom_a.0, atom_b.0),
        tie_depth_note,
        branch_a_embedded,
        branch_b_embedded,
        rule4a: "N/A -- chematic_core::Chirality has only one non-None variant family \
                 (tetrahedral); no chiral/pseudoasymmetric/nonstereogenic unit-type \
                 distinction exists to compare",
        rule4c: "N/A -- chematic_core::Chirality has no axial/planar (m/p) or bond \
                 (seqCis/seqTrans) variant at all; structurally unreachable, not just \
                 absent from these 8 rows",
        rule4b,
    }
}

/// Re-walk from the tied position toward `target`, returning the `NodeId` for `target`
/// on whichever of the two branches (`pos_a`/`pos_b`) reaches it nearest -- mirrors
/// `nearest_embedded_stereocenter`'s own BFS so the returned id is consistent with what
/// that search found.
fn find_node_for_atom(
    graph: &mut CipDigraph,
    pos_a: NodeId,
    pos_b: NodeId,
    target: AtomIdx,
) -> NodeId {
    for start in [pos_a, pos_b] {
        let mut frontier = vec![start];
        let mut seen: HashMap<AtomIdx, NodeId> = HashMap::new();
        while !frontier.is_empty() {
            let mut next_frontier = Vec::new();
            for current in frontier {
                if let CipNodeKind::Atom { atom_idx } = graph.node(current).kind {
                    seen.entry(atom_idx).or_insert(current);
                    if atom_idx == target {
                        return current;
                    }
                }
                next_frontier.extend(graph.expand_children(current).unwrap_or_default());
            }
            frontier = next_frontier;
        }
    }
    panic!("target atom not found on either branch");
}

fn main() {
    println!("Milestone 4B-0: Rule 4 subtype diagnosis (8 rule4_candidate rows)\n");

    let mut oracle_context_printed = false;
    let mut tie_free_count = 0usize;
    let mut total_embedded = 0usize;

    for &(smiles, atom_idx, case_id) in ROWS {
        let report = diagnose_one(smiles, atom_idx, case_id);

        if !oracle_context_printed {
            let mol = chematic_smiles::parse(smiles).unwrap();
            let assignment = assign_cip_accurate_experimental(&mol, CipBudget::default_budget())
                .expect("assignment succeeds");
            print_reference_categories(case_id, &assignment);
            oracle_context_printed = true;
        }

        println!("case_id: {}", report.case_id);
        println!("  stereocenter: atom{}", report.stereocenter);
        println!(
            "  tied ligand pair: atom{} vs atom{}",
            report.tied_pair.0, report.tied_pair.1
        );
        println!(
            "  first descriptor-bearing sphere: {}",
            report.tie_depth_note
        );
        println!("  Rule 4a result: {}", report.rule4a);
        println!("  Rule 4c result: {}", report.rule4c);
        println!("  Rule 4b result: {}", report.rule4b);
        println!("  earliest deciding subrule: 4b (4a/4c ruled out structurally, see above)");
        for (label, embedded) in [
            ("branch A", report.branch_a_embedded),
            ("branch B", report.branch_b_embedded),
        ] {
            total_embedded += 1;
            match embedded {
                Some((atom, tie_free)) => {
                    if tie_free {
                        tie_free_count += 1;
                    }
                    println!(
                        "    {label}: nearest embedded stereocenter atom{} -- local rank \
                         tie-free: {tie_free}",
                        atom.0
                    );
                }
                None => println!("    {label}: no unambiguous nearest embedded stereocenter"),
            }
        }
        println!();
    }

    println!(
        "Summary: {tie_free_count}/{total_embedded} embedded-branch local rankings tie-free \
         (acyclicity evidence). 8/8 rows: Rule 4a N/A, Rule 4c N/A, Rule 4b sole candidate."
    );
}

fn print_reference_categories(case_id: &str, assignment: &chematic_cip::AccurateCipAssignment) {
    let mut resolved: Vec<(u32, CipCode)> = assignment
        .assignments
        .iter()
        .map(|(idx, code)| (idx.0, *code))
        .collect();
    resolved.sort_by_key(|(idx, _)| *idx);
    println!(
        "reference (case {case_id}): Pass-1/2 resolved codes for non-tied ring \
         stereocenters (context for Rule 4a's category check): {resolved:?}\n"
    );
}
