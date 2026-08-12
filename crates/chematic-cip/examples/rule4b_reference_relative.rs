//! Milestone 4B-1.5: the reference-relative Rule 4b diagnostic, run faithfully to spec
//! (no S>R, no reference inversion) against two corpora at once, with the discarded
//! absolute rule kept alongside as a negative control -- per advisor guidance after the
//! prior round's "8/8 both sets" turned out to be a hollow fit (see
//! `docs/rfcs/cip_accurate_rfc.md`'s reduction-proof section).
//!
//! Reference selection here is exactly "the branch's own nearest embedded stereocenter's
//! auxiliary descriptor" (`chain[0]`) -- no majority-across-multiple-tied-siblings step
//! is exercised, because every row in both corpora (mechanically confirmed below, not
//! assumed) has exactly one chirality-bearing atom at every BFS level along each branch;
//! `search_nearest_embedded` (reused from `rule4b_pair_sequence.rs`) already returns
//! `None`/terminates the chain early if it ever finds 2+ distinct atoms tied at one
//! level, so a full-length chain on both branches *is* the mechanical proof of "no
//! sibling tie", not an assumption.
//!
//! Two corpora:
//! - `ROWS_ORIGINAL`: the 8 frozen `rule4_candidate` rows (Milestone 4A-0/4B-0). Kept for
//!   continuity, but proven this round to be reference-selection-degenerate (both
//!   branches' own reference is uniformly R in all 16 branch positions) -- it cannot by
//!   itself validate reference selection, only the operator.
//! - `ROWS_DISCRIMINATING`: 16 new synthetic rows (`validation/cip_rule4b_discriminating_corpus.jsonl`,
//!   RDKit `rdCIPLabeler`-labeled), specifically constructed so the two branches'
//!   references differ (`independent_reference`) or agree at position 0 but diverge only
//!   at position 1/2 (`same_reference_deep_divergence`) -- exactly the two properties the
//!   original 8 rows cannot exercise.
//!
//! Usage: cargo run -p chematic-cip --release --example rule4b_reference_relative

use std::collections::{HashMap, HashSet};

use chematic_cip::{CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, rank_children};
use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

const MAX_PAIR_DEPTH: usize = 4;

/// The 8 frozen `rule4_candidate` rows (Milestone 4A-0/4B-0). Oracle is `S` for all 8.
const ROWS_ORIGINAL: &[(&str, u32, &str, CipCode)] = &[
    (
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
        5,
        "tri-galloyl-a",
        CipCode::S,
    ),
    (
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
        35,
        "tri-galloyl-b",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        3,
        "quinic-a",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        7,
        "quinic-b",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        3,
        "mono-galloyl-a",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        7,
        "mono-galloyl-b",
        CipCode::S,
    ),
    (
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
        11,
        "tetra-galloyl-a",
        CipCode::S,
    ),
    (
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
        26,
        "tetra-galloyl-b",
        CipCode::S,
    ),
];

/// 16 synthetic rows from `validation/cip_rule4b_discriminating_corpus.jsonl` (frozen
/// there via RDKit `rdCIPLabeler`; hardcoded here to keep this example self-contained
/// like its siblings). Scaffold: `CH3-C1-C2-C3-C4(root)-C5-C6-C7-CH3`, all-OH
/// substituted -- root's two branches are constitutionally identical, so Rules 1a/2 tie
/// and resolution is Rule-4b-only. Case tag is informational (see module docs).
const ROWS_DISCRIMINATING: &[(&str, u32, &str, CipCode)] = &[
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C",
        7,
        "indep-ref-1",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)C",
        7,
        "indep-ref-2",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@@H](O)C",
        7,
        "indep-ref-3",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
        7,
        "indep-ref-4",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C",
        7,
        "indep-ref-5",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)C",
        7,
        "indep-ref-6",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)C",
        7,
        "indep-ref-7",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)C",
        7,
        "indep-ref-8",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-1",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C",
        7,
        "deep-div-2",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)C",
        7,
        "deep-div-3",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-4",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-5",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C",
        7,
        "deep-div-6",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)C",
        7,
        "deep-div-7",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-8",
        CipCode::R,
    ),
];

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
                    // 2+ distinct chirality-bearing atoms tied at the same BFS level:
                    // this IS the sibling-tie case Step B's re-sorter would need to
                    // resolve. Returning None here (chain terminates early) is the
                    // mechanical proof, not an assumption, that no such tie occurred
                    // whenever a chain reaches its full expected length below.
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

fn atom_of(graph: &CipDigraph, node: NodeId) -> AtomIdx {
    match graph.node(node).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        other => panic!("expected a physical Atom node, got {other:?}"),
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Winner {
    A,
    B,
}

/// Faithful reference-relative Rule 4b: each branch's reference is its OWN nearest
/// (position-0) auxiliary descriptor -- never a shared/cross-branch value, never
/// inverted. Like precedes Unlike at the first differing relation.
fn faithful_winner(chain_a: &[Option<bool>], chain_b: &[Option<bool>]) -> Option<Winner> {
    let ref_a = (*chain_a.first()?)?;
    let ref_b = (*chain_b.first()?)?;
    for i in 0..chain_a.len().min(chain_b.len()) {
        let (Some(a), Some(b)) = (chain_a[i], chain_b[i]) else {
            break;
        };
        let like_a = a == ref_a;
        let like_b = b == ref_b;
        if like_a != like_b {
            return Some(if like_a { Winner::A } else { Winner::B });
        }
    }
    None
}

/// The discarded absolute rule, kept only as a negative control: raw S precedes R at
/// the first differing raw auxiliary sign, no reference, no branch-relativity at all.
fn s_beats_r_winner(chain_a: &[Option<bool>], chain_b: &[Option<bool>]) -> Option<Winner> {
    for i in 0..chain_a.len().min(chain_b.len()) {
        let (Some(a), Some(b)) = (chain_a[i], chain_b[i]) else {
            break;
        };
        if a != b {
            // a==false means S; S beats R means the S-valued branch wins.
            return Some(if !a { Winner::A } else { Winner::B });
        }
    }
    None
}

fn predict_root(
    winner: Winner,
    groups: &[Vec<NodeId>],
    tied_group_idx: usize,
    pos_a: NodeId,
    pos_b: NodeId,
    position_nodes: &[NodeId],
    chirality: Chirality,
) -> Option<CipCode> {
    let (higher, lower) = match winner {
        Winner::A => (pos_a, pos_b),
        Winner::B => (pos_b, pos_a),
    };
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
    resolve_is_r_from_groups(&resolved_groups, position_nodes, chirality)
        .map(|is_r| if is_r { CipCode::R } else { CipCode::S })
}

struct RowOutcome {
    case_id: String,
    oracle: CipCode,
    chain_a: Vec<Option<bool>>,
    chain_b: Vec<Option<bool>>,
    faithful: Option<CipCode>,
    s_beats_r: Option<CipCode>,
}

fn diagnose(smiles: &str, atom_idx: u32, case_id: &str, oracle: CipCode) -> RowOutcome {
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
            panic!("{case_id}: 3+-way physical tie, outside this diagnostic's scope");
        }
    }
    let tied_group_idx =
        tied_group_idx.unwrap_or_else(|| panic!("{case_id}: expected exactly one physical tie"));
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

    let chirality = mol.atom(idx).chirality;
    let faithful = faithful_winner(&chain_a, &chain_b).and_then(|w| {
        predict_root(
            w,
            &groups,
            tied_group_idx,
            pos_a,
            pos_b,
            &position_nodes,
            chirality,
        )
    });
    let s_beats_r = s_beats_r_winner(&chain_a, &chain_b).and_then(|w| {
        predict_root(
            w,
            &groups,
            tied_group_idx,
            pos_a,
            pos_b,
            &position_nodes,
            chirality,
        )
    });

    let _ = atom_of(&graph, pos_a); // keep helper used across both branches, symmetry with pos_b below
    let _ = atom_of(&graph, pos_b);

    RowOutcome {
        case_id: case_id.to_string(),
        oracle,
        chain_a,
        chain_b,
        faithful,
        s_beats_r,
    }
}

fn chain_fmt(chain: &[Option<bool>]) -> String {
    chain
        .iter()
        .map(|s| match s {
            Some(true) => "R",
            Some(false) => "S",
            None => "?",
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn run_corpus(name: &str, rows: &[(&str, u32, &str, CipCode)]) -> (usize, usize, usize) {
    println!("== {name} ==");
    let (mut faithful_match, mut sr_match, mut diverge) = (0usize, 0usize, 0usize);
    for &(smiles, atom_idx, case_id, oracle) in rows {
        let r = diagnose(smiles, atom_idx, case_id, oracle);
        let f_ok = r.faithful == Some(r.oracle);
        let s_ok = r.s_beats_r == Some(r.oracle);
        if r.faithful != r.s_beats_r {
            diverge += 1;
        }
        if f_ok {
            faithful_match += 1;
        }
        if s_ok {
            sr_match += 1;
        }
        println!(
            "  {:16} A=[{}] B=[{}] oracle={:?} faithful={:?}[{}] s_beats_r={:?}[{}]{}",
            r.case_id,
            chain_fmt(&r.chain_a),
            chain_fmt(&r.chain_b),
            r.oracle,
            r.faithful,
            if f_ok { "OK" } else { "MISS" },
            r.s_beats_r,
            if s_ok { "OK" } else { "MISS" },
            if r.faithful != r.s_beats_r {
                "  <- diverges from S>R"
            } else {
                ""
            },
        );
    }
    println!(
        "  -> faithful {faithful_match}/{len} , s_beats_r(negative control) {sr_match}/{len} , \
         diverge {diverge}/{len}\n",
        len = rows.len()
    );
    (faithful_match, sr_match, diverge)
}

fn main() {
    println!(
        "Milestone 4B-1.5: reference-relative Rule 4b, faithful vs. discarded-absolute \
         negative control\n"
    );

    let (orig_f, orig_sr, orig_div) =
        run_corpus("ROWS_ORIGINAL (8 rows, known degenerate)", ROWS_ORIGINAL);
    let (disc_f, disc_sr, disc_div) = run_corpus(
        "ROWS_DISCRIMINATING (16 rows, non-degenerate)",
        ROWS_DISCRIMINATING,
    );

    println!("=== Summary ===");
    println!(
        "ROWS_ORIGINAL:        faithful {orig_f}/8  s_beats_r {orig_sr}/8  diverge {orig_div}/8 \
         (expect diverge=0 -- this corpus cannot distinguish the two rules)"
    );
    println!(
        "ROWS_DISCRIMINATING:  faithful {disc_f}/16  s_beats_r {disc_sr}/16  diverge {disc_div}/16 \
         (expect diverge>0 -- proves this corpus has teeth; faithful should beat s_beats_r here)"
    );
}
