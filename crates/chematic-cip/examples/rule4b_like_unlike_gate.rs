//! Milestone 4B-2, pre-implementation instrumentation: resolve a real contradiction
//! before writing any production code. `rule4b_pair_sequence.rs`'s empirical "S
//! precedes R" tie-break (8/8 vs oracle) is flagged by the user as likely a
//! coincidental overfit to this 8-row family, not the true Rule 4b mechanism -- the
//! real rule is reference-relative: choose a reference `DescriptorFamily` (majority at
//! the first descriptor-bearing level), convert each branch's sequence to Like/Unlike
//! relative to that reference, and Like precedes Unlike (never an absolute S>R or
//! R>S).
//!
//! Hand-tracing this predicts a *disagreement* on the 4 rows whose shared position-0
//! value is `R` (Like there = R, so a faithful Like/Unlike ranks the R-branch higher,
//! opposite of the empirical "S precedes R" rule) -- meaning if "S precedes R"
//! genuinely reproduces oracle 8/8, a faithful Like/Unlike implementation, fed through
//! the exact same `resolve_is_r_from_groups` parity mapping, should only get 4/8.
//! Rather than trust that hand-trace, this script prints both mappings' predicted
//! final label side by side against the oracle for all 8 rows, settling which mapping
//! (and which orientation) is actually correct by direct instrumentation.
//!
//! Usage: cargo run -p chematic-cip --release --example rule4b_like_unlike_gate

use std::collections::{HashMap, HashSet};

use chematic_cip::{CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, rank_children};
use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

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

/// Mirror a SMILES string: swap every `@` <-> `@@` tetrahedral tag, textually. A
/// well-defined, purely mechanical transformation (no RDKit needed) -- the resulting
/// molecule is the enantiomer of the input, so CIP theory *guarantees* every
/// stereocenter's oracle label flips (S<->R) uniformly, without needing to look the
/// mirrored oracle up. This is what makes the inversion test decisive without a fresh
/// RDKit run: the mirrored oracle is derived, not measured.
fn mirror_smiles(smiles: &str) -> String {
    const PLACEHOLDER: char = '\u{0}';
    smiles
        .replace("@@", &PLACEHOLDER.to_string())
        .replace('@', "@@")
        .replace(PLACEHOLDER, "@")
}

/// One of 4 candidate conventions for turning a deciding-position pair
/// `(branch_a_is_r, branch_b_is_r, reference_is_r)` into which branch (A or B) should
/// be treated as higher priority -- brute-forced rather than hand-derived, since two
/// prior hand-derivations in this exact investigation were wrong.
#[derive(Debug, Clone, Copy)]
#[allow(clippy::enum_variant_names)]
enum Convention {
    /// Branch matching the (position-0-shared) reference wins.
    LikeWinsSharedRef,
    /// Branch NOT matching the (position-0-shared) reference wins.
    UnlikeWinsSharedRef,
    /// Branch matching the OPPOSITE of the position-0-shared value wins.
    LikeWinsOppositeRef,
    /// Branch NOT matching the opposite-of-shared value wins (= matches shared value).
    UnlikeWinsOppositeRef,
}

const CONVENTIONS: [Convention; 4] = [
    Convention::LikeWinsSharedRef,
    Convention::UnlikeWinsSharedRef,
    Convention::LikeWinsOppositeRef,
    Convention::UnlikeWinsOppositeRef,
];

fn higher_branch(
    convention: Convention,
    ra: bool,
    reference_is_r: bool,
    pos_a: NodeId,
    pos_b: NodeId,
) -> (NodeId, NodeId) {
    let effective_ref = match convention {
        Convention::LikeWinsSharedRef | Convention::UnlikeWinsSharedRef => reference_is_r,
        Convention::LikeWinsOppositeRef | Convention::UnlikeWinsOppositeRef => !reference_is_r,
    };
    let a_matches = ra == effective_ref;
    let a_wins = match convention {
        Convention::LikeWinsSharedRef | Convention::LikeWinsOppositeRef => a_matches,
        Convention::UnlikeWinsSharedRef | Convention::UnlikeWinsOppositeRef => !a_matches,
    };
    if a_wins {
        (pos_a, pos_b)
    } else {
        (pos_b, pos_a)
    }
}

fn run_gate(
    smiles: &str,
    atom_idx: u32,
    case_id: &str,
    oracle: CipCode,
) -> HashMap<&'static str, Option<CipCode>> {
    let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
    let idx = AtomIdx(atom_idx);
    let budget = CipBudget::default_budget();

    let mut graph = CipDigraph::new(&mol, idx, budget).expect("digraph builds");
    let root = graph.root();
    let root_children = graph.expand_children(root).expect("root expands");
    let stereo_order = mol
        .stereo_neighbor_order(idx)
        .expect("has stereo_neighbor_order");
    let position_nodes = position_node_ids_generalized(&graph, &root_children, stereo_order)
        .expect("all 4 physical positions resolve");

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
        }
    }
    let tied_group_idx = tied_group_idx.expect("expected exactly one physical tie");
    let (pos_a, pos_b) = tied_pair_nodes.unwrap();

    let chain_a_atoms = embedded_reference_chain(&mut graph, &mol, pos_a);
    let chain_b_atoms = embedded_reference_chain(&mut graph, &mol, pos_b);
    let chain_a: Vec<Option<bool>> = chain_a_atoms
        .iter()
        .map(|&(e, a)| auxiliary_sign(&mol, e, a, budget))
        .collect();
    let chain_b: Vec<Option<bool>> = chain_b_atoms
        .iter()
        .map(|&(e, a)| auxiliary_sign(&mol, e, a, budget))
        .collect();

    let reference_is_r = match (
        chain_a.first().copied().flatten(),
        chain_b.first().copied().flatten(),
    ) {
        (Some(ra), Some(rb)) if ra == rb => ra,
        _ => {
            println!("{case_id}: reference ambiguous/unavailable, skipping");
            return HashMap::new();
        }
    };

    let mut deciding: Option<(bool, bool)> = None; // (ra, rb) at the first differing position
    for i in 0..chain_a.len().min(chain_b.len()) {
        let (Some(ra), Some(rb)) = (chain_a[i], chain_b[i]) else {
            break;
        };
        if ra != rb {
            deciding = Some((ra, rb));
            break;
        }
    }

    let predict = |higher: NodeId, lower: NodeId| -> Option<CipCode> {
        let mut resolved_groups: Vec<Vec<NodeId>> = Vec::with_capacity(groups.len() + 1);
        for (gi, group) in groups.iter().enumerate() {
            if gi == tied_group_idx {
                resolved_groups.push(vec![higher]);
                let mut lower_group = vec![lower];
                lower_group.extend(group.iter().copied().filter(|&n| n != higher && n != lower));
                resolved_groups.push(lower_group);
            } else {
                resolved_groups.push(group.clone());
            }
        }
        resolve_is_r_from_groups(&resolved_groups, &position_nodes, mol.atom(idx).chirality)
            .map(|is_r| if is_r { CipCode::R } else { CipCode::S })
    };

    let mut results = HashMap::new();
    let Some((ra, _rb)) = deciding else {
        println!("{case_id}: no deciding position found");
        return results;
    };
    let names = [
        "LikeWinsSharedRef",
        "UnlikeWinsSharedRef",
        "LikeWinsOppositeRef",
        "UnlikeWinsOppositeRef",
    ];
    print!("  {case_id}: oracle={oracle:?} ");
    for (conv, name) in CONVENTIONS.iter().zip(names.iter()) {
        let (higher, lower) = higher_branch(*conv, ra, reference_is_r, pos_a, pos_b);
        let pred = predict(higher, lower);
        results.insert(*name, pred);
        print!(
            " {name}={}",
            pred.map(|c| format!("{c:?}")).unwrap_or_else(|| "?".into())
        );
    }
    println!();
    results
}

fn main() {
    println!(
        "Milestone 4B-2 pre-implementation gate: Like/Unlike (reference-relative) vs \
         S>R (absolute, empirical) -- side by side against oracle, then the decisive \
         check: same two mechanisms on each molecule's ENANTIOMER (every @/@@ flipped),\n\
         whose oracle label is *derived*, not measured: mirroring a molecule inverts \
         every stereocenter's CIP label, always -- so the mirrored oracle is R for all \
         8 rows.\n"
    );

    let names = [
        "LikeWinsSharedRef",
        "UnlikeWinsSharedRef",
        "LikeWinsOppositeRef",
        "UnlikeWinsOppositeRef",
    ];
    let mut orig_scores: HashMap<&str, usize> = HashMap::new();
    let mut mirror_scores: HashMap<&str, usize> = HashMap::new();

    println!("-- original molecules (oracle = S for all 8) --");
    for &(smiles, atom_idx, case_id) in ROWS {
        let results = run_gate(smiles, atom_idx, case_id, CipCode::S);
        for name in names {
            if results.get(name).copied().flatten() == Some(CipCode::S) {
                *orig_scores.entry(name).or_insert(0) += 1;
            }
        }
    }
    println!("-- mirrored molecules (derived oracle = R for all 8) --");
    for &(smiles, atom_idx, case_id) in ROWS {
        let mirrored = mirror_smiles(smiles);
        let results = run_gate(&mirrored, atom_idx, case_id, CipCode::R);
        for name in names {
            if results.get(name).copied().flatten() == Some(CipCode::R) {
                *mirror_scores.entry(name).or_insert(0) += 1;
            }
        }
    }

    println!("\nSummary (original/8, mirrored/8):");
    for name in names {
        println!(
            "  {name}: {}/8, {}/8",
            orig_scores.get(name).copied().unwrap_or(0),
            mirror_scores.get(name).copied().unwrap_or(0)
        );
    }
    println!(
        "\nA correct, invariant convention scores 8/8 on BOTH lines. 4/4 (mid) is \
         invariant but still wrong; 8/0 or 0/8 is the absolute-rule failure mode."
    );
}
