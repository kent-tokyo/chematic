//! Milestone 4B-1, phase 2: extends `rule4b_aux_sign.rs`'s single-reference gate (which
//! found 0/8 discriminating -- both branches' *first* embedded reference always carries
//! the same auxiliary sign) to a genuine ordered *sequence* of paired auxiliary
//! descriptors per branch, one level deeper each step via
//! [`CipDigraph::new_with_artificial_ancestor`].
//!
//! **Result: 8/8 match the oracle.** At the first sequence position where the two
//! branches' auxiliary signs differ, S precedes R decides the outer tied atom's own
//! R/S -- note the orientation: this is the *opposite* convention from
//! `assign::assign_one_with_rule5`'s "R precedes S" (Rule 5 compares an outer atom's
//! own resolution directly against 2 embedded references at a single position; this is
//! a different comparison -- the first *differing* position in a nested pair
//! sequence -- and empirically needs the opposite orientation to reproduce the
//! oracle). Found by trying "R precedes S" first (uniform 8/8 wrong-direction
//! mismatch, all at sequence position 1) and then testing the flipped orientation,
//! not assumed from either convention.
//!
//! This is still a diagnostic/validation script, not production wiring -- the natural
//! next step, if this direction is confirmed to generalize, is a real
//! `DescriptorPairList`-shaped implementation wired into `assign.rs` (Milestone 4B-2),
//! not this script itself.
//!
//! Usage: cargo run -p chematic-cip --release --example rule4b_pair_sequence

use std::collections::{HashMap, HashSet};

use chematic_cip::{CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, rank_children};
use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

/// Same 8 `rule4_candidate` rows as `examples/rule4_diagnose.rs`/`rule4b_aux_sign.rs`.
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

const MAX_PAIR_DEPTH: usize = 4;

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

/// Search `frontier` outward, level by level, for the nearest atom with
/// `chirality != None`. Returns `(embedded_atom, arrival_atom, node_id)` --
/// `arrival_atom` is the real atom immediately preceding it on this path (the seed
/// [`CipDigraph::new_with_artificial_ancestor`] needs), `node_id` is its own position
/// in `graph` (so the caller can continue the search *beyond* it via
/// [`CipDigraph::expand_children`]). `None` if ambiguous (2+ distinct atoms at the
/// nearest level) or the frontier is exhausted.
fn search_nearest_embedded(
    graph: &mut CipDigraph,
    mol: &Molecule,
    frontier: Vec<NodeId>,
) -> Option<(AtomIdx, AtomIdx, NodeId)> {
    let mut frontier = frontier;
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
            return Some((embedded_atom, arrival_atom, node_id));
        }
        frontier = next_frontier;
    }
    None
}

/// The ordered chain of `(embedded_atom, arrival_atom)` pairs reached from `start`,
/// each one level deeper than the last -- position 0 is the nearest embedded
/// stereocenter, position 1 is the nearest one found *beyond* it, etc., up to
/// `MAX_PAIR_DEPTH` or until the branch runs out of embedded stereocenters / turns
/// ambiguous.
fn embedded_reference_chain(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
) -> Vec<(AtomIdx, AtomIdx)> {
    let mut chain = Vec::new();
    let mut frontier = vec![start];
    for _ in 0..MAX_PAIR_DEPTH {
        let Some((embedded_atom, arrival_atom, node_id)) =
            search_nearest_embedded(graph, mol, frontier)
        else {
            break;
        };
        chain.push((embedded_atom, arrival_atom));
        frontier = match graph.expand_children(node_id) {
            Ok(children) => children,
            Err(_) => break,
        };
    }
    chain
}

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

fn atom_of(graph: &CipDigraph, node: NodeId) -> AtomIdx {
    match graph.node(node).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        other => panic!("expected a physical Atom node, got {other:?}"),
    }
}

struct RowResult {
    case_id: String,
    stereocenter: u32,
    tied_pair: (u32, u32),
    chain_a: Vec<(AtomIdx, AtomIdx, Option<bool>)>,
    chain_b: Vec<(AtomIdx, AtomIdx, Option<bool>)>,
    deciding_position: Option<usize>,
    predicted: Option<CipCode>,
}

fn diagnose_and_resolve(smiles: &str, atom_idx: u32, case_id: &str) -> RowResult {
    let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
    let idx = AtomIdx(atom_idx);
    let budget = CipBudget::default_budget();

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

    let chain_a_atoms = embedded_reference_chain(&mut graph, &mol, pos_a);
    let chain_b_atoms = embedded_reference_chain(&mut graph, &mol, pos_b);

    let chain_a: Vec<(AtomIdx, AtomIdx, Option<bool>)> = chain_a_atoms
        .iter()
        .map(|&(e, a)| (e, a, auxiliary_sign(&mol, e, a, budget)))
        .collect();
    let chain_b: Vec<(AtomIdx, AtomIdx, Option<bool>)> = chain_b_atoms
        .iter()
        .map(|&(e, a)| (e, a, auxiliary_sign(&mol, e, a, budget)))
        .collect();

    // Lexicographic pair-sequence comparison: first position where the two chains'
    // auxiliary signs differ decides, R precedes S (same convention as the
    // single-reference gate and Rule 5's own tie-break).
    let mut deciding_position = None;
    let mut predicted = None;
    for i in 0..chain_a.len().min(chain_b.len()) {
        let (Some(ra), Some(rb)) = (chain_a[i].2, chain_b[i].2) else {
            break;
        };
        if ra != rb {
            deciding_position = Some(i);
            // NOTE: inverted vs. rule4b_aux_sign.rs's own convention (S precedes R, not
            // R precedes S) -- see module docs update: the depth-1 gate found this
            // orientation empirically, after "R precedes S" (Rule 5's convention)
            // produced a uniform 8/8 wrong-direction mismatch.
            let (higher, lower) = if ra { (pos_b, pos_a) } else { (pos_a, pos_b) };
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
            predicted = resolve_is_r_from_groups(
                &resolved_groups,
                &position_nodes,
                mol.atom(idx).chirality,
            )
            .map(|is_r| if is_r { CipCode::R } else { CipCode::S });
            break;
        }
    }

    RowResult {
        case_id: case_id.to_string(),
        stereocenter: atom_idx,
        tied_pair: (atom_a.0, atom_b.0),
        chain_a,
        chain_b,
        deciding_position,
        predicted,
    }
}

fn chain_str(chain: &[(AtomIdx, AtomIdx, Option<bool>)]) -> String {
    chain
        .iter()
        .map(|(e, a, s)| {
            let sign = match s {
                Some(true) => "R",
                Some(false) => "S",
                None => "?",
            };
            format!("atom{}(via atom{})={}", e.0, a.0, sign)
        })
        .collect::<Vec<_>>()
        .join(" -> ")
}

fn main() {
    println!(
        "Milestone 4B-1 phase 2: descriptor pair-sequence gate (8 rule4_candidate rows, \
         max depth {MAX_PAIR_DEPTH})\n"
    );

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
        println!("  branch A chain: {}", chain_str(&r.chain_a));
        println!("  branch B chain: {}", chain_str(&r.chain_b));
        match r.deciding_position {
            Some(pos) => println!("  deciding pair-sequence position: {pos}"),
            None => println!("  deciding pair-sequence position: none found within depth cap"),
        }
        match r.predicted {
            Some(CipCode::S) => {
                matched += 1;
                println!("  predicted: S -- MATCHES oracle (S)");
            }
            Some(other) => {
                mismatched += 1;
                println!("  predicted: {:?} -- MISMATCHES oracle (S)", other);
            }
            None => {
                inconclusive += 1;
                println!("  predicted: inconclusive (no distinguishing position found)");
            }
        }
        println!();
    }

    println!(
        "Summary: {matched}/8 matched oracle, {mismatched}/8 mismatched, {inconclusive}/8 \
         inconclusive."
    );
}
