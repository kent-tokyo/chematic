//! Milestone 4B-1, gating deliverable: compute each embedded reference stereocenter's
//! true AUXILIARY R/S sign via [`CipDigraph::new_with_artificial_ancestor`] (new this
//! milestone -- see `crates/chematic-cip/src/digraph.rs`), then check whether the
//! simplest Rule-4b-shaped combination (R-precedes-S tie-break on the two branches'
//! auxiliary references, mirroring Milestone 4A's own Rule 5 convention in
//! `assign::assign_one_with_rule5`) reproduces the oracle `S` label for all 8
//! `rule4_candidate` rows (Milestone 4A-0/4B-0).
//!
//! This is the go/no-go gate for the rest of Milestone 4B: Milestone 4B-0's diagnosis
//! confirmed the embedded reference's *local* (parent-excluded) children rank tie-free,
//! but explicitly could not confirm its auxiliary *sign* without this construct --
//! which is exactly where a hidden circularity (needing the other tied atom's own
//! still-unresolved global value) would have to surface, if it exists at all.
//!
//! Usage: cargo run -p chematic-cip --release --example rule4b_aux_sign

use std::collections::{HashMap, HashSet};

use chematic_cip::{
    CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, assign_cip_accurate_experimental,
    rank_children,
};
use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

/// Same 8 `rule4_candidate` rows as `examples/rule4_diagnose.rs` (Milestone 4B-0).
/// Oracle (`modern`) is `S` for all 8.
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

/// Map `stereo_neighbor_order` positions to root-child `NodeId`s, matching either a
/// real `Atom` or -- for the artificial-ancestor slot -- a `RingDuplicate` whose
/// `closure_atom` is the target. Generalizes `assign::position_node_ids` (private to
/// that module, and written before this milestone's artificial-ancestor construct
/// existed) to also recognize that duplicate as a valid physical position.
fn position_node_ids_generalized(
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
            root_children
                .iter()
                .copied()
                .find(|&id| match graph.node(id).kind {
                    CipNodeKind::Atom { atom_idx: a } => a == atom_idx,
                    CipNodeKind::RingDuplicate { closure_atom, .. } => closure_atom == atom_idx,
                    _ => false,
                })?
        };
        result.push(node_id);
    }
    Some(result)
}

/// Identical to `assign::resolve_is_r_from_groups` (private to that module) --
/// duplicated here rather than exposed, since this is diagnostic/validation code, not
/// production wiring (Milestone 4B-1's own scope limit).
fn resolve_is_r_from_groups(
    groups: &[Vec<NodeId>],
    position_nodes: &[NodeId],
    chirality: Chirality,
) -> Option<bool> {
    let n = groups.len() as u8;
    let mut rank_of: HashMap<NodeId, u8> = HashMap::new();
    for (group_idx, group) in groups.iter().enumerate() {
        for &node in group {
            rank_of.insert(node, n - group_idx as u8);
        }
    }

    let raw_ranks: Vec<u8> = position_nodes.iter().map(|node| rank_of[node]).collect();
    let mut distinct_ranks = raw_ranks.clone();
    distinct_ranks.sort_unstable();
    distinct_ranks.dedup();
    let ranks: Vec<u8> = raw_ranks
        .iter()
        .map(|&r| distinct_ranks.iter().position(|&x| x == r).unwrap() as u8 + 1)
        .collect();

    let lowest_pos = ranks.iter().position(|&r| r == 1)?;
    let parity_odd = lowest_pos % 2 == 1;
    let smiles_cw = chirality == Chirality::Clockwise;
    let cw_from_lowest = smiles_cw ^ parity_odd;

    let remaining_ranks: Vec<u8> = (0..4usize)
        .filter(|&i| i != lowest_pos)
        .map(|i| ranks[i])
        .collect();
    let remaining_swaps_odd = swap_parity(&remaining_ranks)?;

    Some(cw_from_lowest ^ remaining_swaps_odd)
}

/// Identical to `assign::swap_parity` (private to that module).
fn swap_parity(remaining_ranks: &[u8]) -> Option<bool> {
    let mut r = remaining_ranks.to_vec();
    let target = [2u8, 3, 4];
    let mut swaps = 0usize;
    for i in 0..3 {
        if r[i] != target[i] {
            let j_rel = r[i + 1..].iter().position(|&x| x == target[i])?;
            r.swap(i, j_rel + i + 1);
            swaps += 1;
        }
    }
    Some(swaps % 2 == 1)
}

/// BFS, level by level, for the nearest atom with `chirality != None` in `start`'s
/// subtree, returning `(embedded_atom, arrival_atom)` -- `arrival_atom` is the real
/// atom immediately preceding `embedded_atom` on this specific root-to-node path, i.e.
/// exactly the atom [`CipDigraph::new_with_artificial_ancestor`] needs as its seed to
/// reproduce this same branch-local view when rooted fresh at `embedded_atom`. `None`
/// if ambiguous (2+ distinct atoms at the nearest level) or absent -- mirrors
/// `assign::nearest_embedded_stereocenter`'s own ambiguity handling.
fn nearest_embedded_with_arrival(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
) -> Option<(AtomIdx, AtomIdx)> {
    let mut frontier = vec![start];
    while !frontier.is_empty() {
        let mut found_this_level: Option<(AtomIdx, NodeId)> = None;
        let mut next_frontier = Vec::new();
        for &current in &frontier {
            if let CipNodeKind::Atom { atom_idx } = graph.node(current).kind
                && mol.atom(atom_idx).chirality != Chirality::None
            {
                match found_this_level {
                    None => found_this_level = Some((atom_idx, current)),
                    Some((existing, _)) if existing == atom_idx => {}
                    Some(_) => return None,
                }
            }
            next_frontier.extend(graph.expand_children(current).ok()?);
        }
        if let Some((embedded_atom, node_id)) = found_this_level {
            let parent_id = graph.node(node_id).parent?;
            let arrival_atom = match graph.node(parent_id).kind {
                CipNodeKind::Atom { atom_idx } => atom_idx,
                _ => return None,
            };
            return Some((embedded_atom, arrival_atom));
        }
        frontier = next_frontier;
    }
    None
}

fn atom_of(graph: &CipDigraph, node: NodeId) -> AtomIdx {
    match graph.node(node).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        other => panic!("expected a physical Atom node, got {other:?}"),
    }
}

/// The embedded reference atom's true auxiliary R/S sign, computed via a fresh
/// [`CipDigraph::new_with_artificial_ancestor`] root -- `Some(true)` = R, `Some(false)`
/// = S, `None` if the auxiliary digraph itself can't resolve it (still tied, or budget
/// exceeded).
fn auxiliary_sign(
    mol: &Molecule,
    embedded_atom: AtomIdx,
    arrival_atom: AtomIdx,
    budget: CipBudget,
) -> Option<bool> {
    let mut graph =
        CipDigraph::new_with_artificial_ancestor(mol, embedded_atom, arrival_atom, budget).ok()?;
    let root = graph.root();
    let root_children = graph.expand_children(root).ok()?;
    let stereo_order = mol.stereo_neighbor_order(embedded_atom)?;
    let position_nodes = position_node_ids_generalized(&graph, &root_children, stereo_order)?;
    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).ok()?;
    resolve_is_r_from_groups(&groups, &position_nodes, mol.atom(embedded_atom).chirality)
}

struct RowResult {
    case_id: String,
    stereocenter: u32,
    tied_pair: (u32, u32),
    embedded_a: Option<(AtomIdx, AtomIdx, Option<bool>)>,
    embedded_b: Option<(AtomIdx, AtomIdx, Option<bool>)>,
    global_a: Option<CipCode>,
    global_b: Option<CipCode>,
    predicted: Option<CipCode>,
}

fn diagnose_and_resolve(smiles: &str, atom_idx: u32, case_id: &str) -> RowResult {
    let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
    let idx = AtomIdx(atom_idx);
    let budget = CipBudget::default_budget();

    let assignment = assign_cip_accurate_experimental(&mol, budget)
        .expect("Pass-1/2 assignment succeeds (for the reference global codes only)");
    let global: HashMap<AtomIdx, CipCode> = assignment.assignments.into_iter().collect();

    let mut graph = CipDigraph::new(&mol, idx, budget).expect("digraph builds");
    let root = graph.root();
    let root_children = graph.expand_children(root).expect("root expands");
    let stereo_order = mol
        .stereo_neighbor_order(idx)
        .expect("stereocenter has stereo_neighbor_order");
    let position_nodes = position_node_ids_generalized(&graph, &root_children, stereo_order)
        .expect("all 4 physical positions resolve to root children");

    let mut ctx = CompareContext::new();
    let groups = rank_children(&mut graph, &root_children, &mut ctx).expect("root children rank");

    let position_set: HashSet<NodeId> = position_nodes.iter().copied().collect();
    let mut tied_group_idx = None;
    let mut tied_pair_nodes = None;
    for (gi, group) in groups.iter().enumerate() {
        let physical: Vec<NodeId> = group
            .iter()
            .copied()
            .filter(|n| position_set.contains(n))
            .collect();
        if physical.len() == 2 {
            tied_group_idx = Some(gi);
            tied_pair_nodes = Some((physical[0], physical[1]));
        } else if physical.len() > 2 {
            panic!("{case_id}: 3+-way physical tie, outside this validation's scope");
        }
    }
    let tied_group_idx =
        tied_group_idx.unwrap_or_else(|| panic!("{case_id}: expected exactly one physical tie"));
    let (pos_a, pos_b) = tied_pair_nodes.unwrap();
    let atom_a = atom_of(&graph, pos_a);
    let atom_b = atom_of(&graph, pos_b);

    let embedded_a_ids = nearest_embedded_with_arrival(&mut graph, &mol, pos_a);
    let embedded_b_ids = nearest_embedded_with_arrival(&mut graph, &mol, pos_b);

    let embedded_a = embedded_a_ids.map(|(embedded, arrival)| {
        (
            embedded,
            arrival,
            auxiliary_sign(&mol, embedded, arrival, budget),
        )
    });
    let embedded_b = embedded_b_ids.map(|(embedded, arrival)| {
        (
            embedded,
            arrival,
            auxiliary_sign(&mol, embedded, arrival, budget),
        )
    });

    let global_a = embedded_a_ids.and_then(|(e, _)| global.get(&e).copied());
    let global_b = embedded_b_ids.and_then(|(e, _)| global.get(&e).copied());

    // Candidate Rule-4b-shaped combination: R precedes S on the two branches'
    // auxiliary references (mirrors Rule 5's own tie-break convention in
    // `assign::assign_one_with_rule5`). Only distinguishing if the two auxiliary
    // signs actually differ -- if they match, this convention has no discriminating
    // power and the row is left unpredicted rather than guessed.
    let aux_a = embedded_a.and_then(|(_, _, s)| s);
    let aux_b = embedded_b.and_then(|(_, _, s)| s);
    let predicted = match (aux_a, aux_b) {
        (Some(ra), Some(rb)) if ra != rb => {
            let (higher, lower) = if ra { (pos_a, pos_b) } else { (pos_b, pos_a) };
            let mut resolved_groups: Vec<Vec<NodeId>> = Vec::with_capacity(groups.len() + 1);
            for (gi, group) in groups.iter().enumerate() {
                if gi == tied_group_idx {
                    resolved_groups.push(vec![higher]);
                    let mut lower_group = vec![lower];
                    lower_group
                        .extend(group.iter().copied().filter(|&n| n != higher && n != lower));
                    resolved_groups.push(lower_group);
                } else {
                    resolved_groups.push(group.clone());
                }
            }
            resolve_is_r_from_groups(&resolved_groups, &position_nodes, mol.atom(idx).chirality)
                .map(|is_r| if is_r { CipCode::R } else { CipCode::S })
        }
        _ => None,
    };

    let _ = (atom_a, atom_b); // used only via tied_pair below

    RowResult {
        case_id: case_id.to_string(),
        stereocenter: atom_idx,
        tied_pair: (atom_a.0, atom_b.0),
        embedded_a,
        embedded_b,
        global_a,
        global_b,
        predicted,
    }
}

fn code_str(code: Option<CipCode>) -> String {
    match code {
        Some(CipCode::R) => "R".to_string(),
        Some(CipCode::S) => "S".to_string(),
        Some(CipCode::LowerR) => "r".to_string(),
        Some(CipCode::LowerS) => "s".to_string(),
        Some(CipCode::E) => "E".to_string(),
        Some(CipCode::Z) => "Z".to_string(),
        None => "?".to_string(),
    }
}

fn aux_str(sign: Option<bool>) -> &'static str {
    match sign {
        Some(true) => "R",
        Some(false) => "S",
        None => "?",
    }
}

fn main() {
    println!("Milestone 4B-1 gate: auxiliary-sign construct vs oracle (8 rule4_candidate rows)\n");

    let mut matched = 0usize;
    let mut mismatched = 0usize;
    let mut inconclusive = 0usize;

    for &(smiles, atom_idx, case_id) in ROWS {
        let r = diagnose_and_resolve(smiles, atom_idx, case_id);

        println!("case_id: {}", r.case_id);
        println!("  stereocenter: atom{}", r.stereocenter);
        println!(
            "  tied ligand pair: atom{} vs atom{}",
            r.tied_pair.0, r.tied_pair.1
        );
        if let Some((embedded, arrival, sign)) = r.embedded_a {
            println!(
                "  branch A: embedded atom{} (arrival atom{}) -- auxiliary sign: {}, \
                 global (Pass-1) code: {}",
                embedded.0,
                arrival.0,
                aux_str(sign),
                code_str(r.global_a)
            );
        }
        if let Some((embedded, arrival, sign)) = r.embedded_b {
            println!(
                "  branch B: embedded atom{} (arrival atom{}) -- auxiliary sign: {}, \
                 global (Pass-1) code: {}",
                embedded.0,
                arrival.0,
                aux_str(sign),
                code_str(r.global_b)
            );
        }
        match r.predicted {
            Some(CipCode::S) => {
                matched += 1;
                println!("  predicted: S -- MATCHES oracle (S)");
            }
            Some(other) => {
                mismatched += 1;
                println!(
                    "  predicted: {} -- MISMATCHES oracle (S)",
                    code_str(Some(other))
                );
            }
            None => {
                inconclusive += 1;
                println!(
                    "  predicted: inconclusive (auxiliary signs did not distinguish the \
                     two branches, or one side unresolved)"
                );
            }
        }
        println!();
    }

    println!(
        "Summary: {matched}/8 matched oracle, {mismatched}/8 mismatched, {inconclusive}/8 \
         inconclusive."
    );
}
