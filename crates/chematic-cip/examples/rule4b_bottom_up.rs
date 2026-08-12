//! Milestone 4B-2 (scope grown from "M4B-1.5 diagnostic" per the primary source, see
//! `docs/rfcs/cip_accurate_rfc.md`): a from-scratch, in-place auxiliary-descriptor engine
//! that replaces the per-atom `new_with_artificial_ancestor` re-rooting architecture
//! with the one Hanson, Musacchio, Mayfield et al. 2018 (*J. Chem. Inf. Model.* 58(9),
//! 1755-1765) actually specify:
//!
//! - **One digraph**, rooted at the true outer stereocenter (`CipDigraph::new`,
//!   unmodified). Every embedded stereocenter is resolved *in place*, using the node
//!   already present in this digraph -- never by building a fresh digraph rooted at
//!   that embedded atom. Per the paper, re-rooting per atom gives ring-closure
//!   semantics relative to the *wrong* root, which is exactly what produced the
//!   pair-antisymmetry bug this milestone found in `rule4b_reference_relative.rs`.
//! - **Bottom-up**: `resolve_chirality` is a plain memoized recursion (memoized by
//!   `NodeId`, which is already a unique identity within one digraph -- ring atoms that
//!   recur via closure get *duplicate* nodes, never a second real node, so no
//!   `(atom, arrival)` composite key is needed). Recursion only ever descends into
//!   strictly-deeper nodes, matching the paper's postulate that an auxiliary center's
//!   descriptor never depends on anything between it and the root -- so this always
//!   terminates and never cycles.
//! - **Back-to-root ligand, Rule 1a only**: per the paper, "the priority of a ligand
//!   leading back to the digraph root will always be ranked by Rule 1a, with no need to
//!   consider auxiliary centers... the path back to the root is always unique in
//!   connectivity and atomic numbers." Implemented as [`BackItem`]/[`expand_back_item`]:
//!   a synthetic frontier that walks *up* through the existing parent chain (reusing
//!   already-built off-path subtrees as-is -- no rebuilding, no fresh digraph), compared
//!   against a real forward ligand via plain sorted-atomic-number shell comparison
//!   ([`compare_rule1a_only`]). Shell-pooling is unsound for the *full* rule cascade
//!   (see `compare.rs` module docs) but is exactly Rule 1a's own classical definition,
//!   so it is sound restricted to this one rule.
//!
//! Validation order per advisor guidance: smallest, simplest coupled case first
//! (`VS196`, hexachlorocyclohexane, tag `4b` only, all-uppercase -- no pseudoasymmetry,
//! so no need for this diagnostic's `bool` (`is_r`) representation to express `r`/`s`
//! yet), then the original 8-row `rule4_candidate` residual (still *not* proof by
//! itself -- it is all-`S`, see the reduction-proof section of the RFC -- but a
//! regression floor), then the 16-row discriminating corpus (must hold as regression).
//!
//! Usage: cargo run -p chematic-cip --release --example rule4b_bottom_up

use std::collections::HashMap;

const DEBUG: bool = false;

fn atom_of_node(graph: &CipDigraph, n: NodeId) -> String {
    match graph.node(n).kind {
        CipNodeKind::Atom { atom_idx } => format!("atom{}", atom_idx.0),
        CipNodeKind::RingDuplicate { closure_atom, .. } => format!("dup(atom{})", closure_atom.0),
        CipNodeKind::ImplicitHydrogen => "H".to_string(),
        CipNodeKind::MultipleBondDuplicate {
            duplicated_atom, ..
        } => {
            format!("mdup(atom{})", duplicated_atom.0)
        }
    }
}

use chematic_cip::{
    CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, compare::CipCompareError,
    rank_children, rational::cmp_atomic_number_key,
};
use chematic_core::{AtomIdx, Chirality, CipCode, Molecule, STEREO_H_SENTINEL};

/// `VS196`: hexachlorocyclohexane, one ring, every ring carbon a stereocenter bearing
/// one Cl + one H, tag `4b` only (no Rule 5 needed -- all reference labels are
/// uppercase). Suite reference: `2R 3S 4S 5R 6S 7R` (1-indexed); RDKit `rdCIPLabeler`
/// matches exactly (spot-checked, see RFC). 0-indexed here.
const VS196_SMILES: &str = "Cl[C@H]1[C@H]([C@H]([C@@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl";
const VS196_ROWS: &[(u32, CipCode)] = &[
    (1, CipCode::R),
    (2, CipCode::S),
    (3, CipCode::S),
    (4, CipCode::R),
    (5, CipCode::S),
    (6, CipCode::R),
];

/// `VS197`: same scaffold, different substitution pattern. Suite: `2R 3R 4R 5R 6S 7S`.
const VS197_SMILES: &str = "Cl[C@H]1[C@@H]([C@H]([C@@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl";
const VS197_ROWS: &[(u32, CipCode)] = &[
    (1, CipCode::R),
    (2, CipCode::R),
    (3, CipCode::R),
    (4, CipCode::R),
    (5, CipCode::S),
    (6, CipCode::S),
];

/// The frozen 8-row `rule4_candidate` residual (Milestone 4A-0/4B-0). All-`S` --
/// regression floor only, never proof by itself (see RFC reduction-proof section).
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

/// `validation/cip_rule4b_discriminating_corpus.jsonl` (16 rows, RDKit-oracled,
/// acyclic, independent/differing references + same-reference-deep-divergence). The
/// old per-atom-re-rooted engine got 16/16 here (`rule4b_reference_relative.rs`) --
/// mandatory regression floor for the new in-place engine, since this corpus still
/// exercises a real back-to-root ligand (toward the outer root) with no ring/duplicate
/// complications, and localizes whether a ring failure is plumbing or operator.
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

/// A position in the back-to-root ligand's growing frontier.
#[derive(Clone, Copy, PartialEq, Eq)]
enum BackItem {
    /// A node already present in the outer digraph -- normal forward expansion from
    /// here on via `expand_children`, no special handling.
    Reused(NodeId),
    /// The walk continues upward: `ancestor` is this sphere's atom; `came_from` is the
    /// specific child of `ancestor` to exclude when expanding (the one leading back
    /// down toward the atom under resolution -- never re-enter its subtree).
    Ascending { ancestor: NodeId, came_from: NodeId },
}

fn back_item_atomic_number(graph: &CipDigraph, item: BackItem) -> chematic_cip::AtomicNumberKey {
    let node = match item {
        BackItem::Reused(n) => n,
        BackItem::Ascending { ancestor, .. } => ancestor,
    };
    graph.node(node).atomic_number
}

fn expand_back_item(
    graph: &mut CipDigraph,
    item: BackItem,
) -> Result<Vec<BackItem>, CipCompareError> {
    match item {
        BackItem::Reused(n) => Ok(graph
            .expand_children(n)
            .map_err(CipCompareError::Digraph)?
            .into_iter()
            .map(BackItem::Reused)
            .collect()),
        BackItem::Ascending {
            ancestor,
            came_from,
        } => {
            let mut result: Vec<BackItem> = graph
                .expand_children(ancestor)
                .map_err(CipCompareError::Digraph)?
                .into_iter()
                .filter(|&c| c != came_from)
                .map(BackItem::Reused)
                .collect();
            if let Some(grandparent) = graph.node(ancestor).parent {
                result.push(BackItem::Ascending {
                    ancestor: grandparent,
                    came_from: ancestor,
                });
            }
            Ok(result)
        }
    }
}

/// Rule 1a only, shell-pooled (sound for Rule 1a alone -- it is the classical
/// sum/sorted-multiset-of-atomic-numbers rule; shell-pooling is only unsound once
/// branch-provenance-aware rules like 1b/4b/5 are mixed in, which this function never
/// invokes). `Greater` means the back-to-root ligand (`lhs`) outranks `rhs`.
fn compare_rule1a_only(
    graph: &mut CipDigraph,
    mut lhs: Vec<BackItem>,
    mut rhs: Vec<NodeId>,
) -> Result<std::cmp::Ordering, CipCompareError> {
    use std::cmp::Ordering;
    loop {
        let mut lhs_keys: Vec<_> = lhs
            .iter()
            .map(|&i| back_item_atomic_number(graph, i))
            .collect();
        let mut rhs_keys: Vec<_> = rhs.iter().map(|&n| graph.node(n).atomic_number).collect();
        lhs_keys.sort_unstable_by(|a, b| cmp_atomic_number_key(*b, *a));
        rhs_keys.sort_unstable_by(|a, b| cmp_atomic_number_key(*b, *a));
        if DEBUG {
            eprintln!("    compare_rule1a_only sphere: lhs={lhs_keys:?} rhs={rhs_keys:?}");
        }
        let n = lhs_keys.len().max(rhs_keys.len());
        for i in 0..n {
            let l = lhs_keys.get(i).copied();
            let r = rhs_keys.get(i).copied();
            let ord = match (l, r) {
                (Some(a), Some(b)) => cmp_atomic_number_key(a, b),
                (Some(_), None) => Ordering::Greater,
                (None, Some(_)) => Ordering::Less,
                (None, None) => Ordering::Equal,
            };
            if ord != Ordering::Equal {
                if DEBUG {
                    eprintln!("    compare_rule1a_only DECIDED: {ord:?}");
                }
                return Ok(ord);
            }
        }
        if lhs.is_empty() && rhs.is_empty() {
            return Ok(Ordering::Equal);
        }
        let mut next_lhs = Vec::new();
        for item in lhs {
            next_lhs.extend(expand_back_item(graph, item)?);
        }
        let mut next_rhs = Vec::new();
        for node in rhs {
            next_rhs.extend(
                graph
                    .expand_children(node)
                    .map_err(CipCompareError::Digraph)?,
            );
        }
        lhs = next_lhs;
        rhs = next_rhs;
    }
}

/// BFS, in place within the outer digraph (no re-rooting), for the nearest
/// chirality-bearing atom reachable from `start`'s children onward. `None` if 2+
/// distinct atoms tie at the nearest level (ambiguous -- mechanical proof of "no
/// sibling tie" when this never fires, same discipline as
/// `rule4b_reference_relative.rs`'s `search_nearest_embedded`).
fn nearest_embedded(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
) -> Result<Option<NodeId>, CipCompareError> {
    let mut frontier = graph
        .expand_children(start)
        .map_err(CipCompareError::Digraph)?;
    loop {
        if frontier.is_empty() {
            return Ok(None);
        }
        let mut found: Option<(AtomIdx, NodeId)> = None;
        for &n in &frontier {
            if let CipNodeKind::Atom { atom_idx } = graph.node(n).kind
                && mol.atom(atom_idx).chirality != Chirality::None
            {
                match found {
                    None => found = Some((atom_idx, n)),
                    Some((existing, _)) if existing == atom_idx => {}
                    Some(_) => return Ok(None),
                }
            }
        }
        if let Some((_, node_id)) = found {
            return Ok(Some(node_id));
        }
        let mut next = Vec::new();
        for &n in &frontier {
            next.extend(graph.expand_children(n).map_err(CipCompareError::Digraph)?);
        }
        frontier = next;
    }
}

const MAX_CHAIN_DEPTH: usize = 4;

/// The ordered chain of nearest-embedded-stereocenter `NodeId`s reached from `start`,
/// one level deeper each step -- in place, no re-rooting (contrast with
/// `rule4b_reference_relative.rs`'s `embedded_reference_chain`, which built a fresh
/// digraph per element via `auxiliary_sign`).
fn embedded_chain(
    graph: &mut CipDigraph,
    mol: &Molecule,
    start: NodeId,
) -> Result<Vec<NodeId>, CipCompareError> {
    // `start` itself is chain position 0 if it's chirality-bearing (it always is here
    // -- `start` is one of a tied group's own members, i.e. the outer atom's own
    // forward child, and every case in this residual family has the tied children
    // themselves be embedded stereocenters). `nearest_embedded` starts searching from
    // `start`'s *children*, so it is only correct for continuing the chain past an
    // already-found element, never for the first one.
    let mut chain = Vec::new();
    if let CipNodeKind::Atom { atom_idx } = graph.node(start).kind
        && mol.atom(atom_idx).chirality != Chirality::None
    {
        chain.push(start);
    }
    let mut current = start;
    for _ in 0..MAX_CHAIN_DEPTH {
        match nearest_embedded(graph, mol, current)? {
            Some(node_id) => {
                chain.push(node_id);
                current = node_id;
            }
            None => break,
        }
    }
    Ok(chain)
}

/// Rule 4b tiebreak for an exactly-2-way tie (every case observed in this corpus family
/// so far; 3+-way ties are out of scope and reported rather than guessed). Faithful
/// reference-relative Like/Unlike: each branch's own reference is its own nearest
/// (position-0) resolved chirality, never shared, never inverted -- the operator
/// confirmed correct by `rule4b_operator_tests.rs` (7/7). The chain elements'
/// chirality is resolved via *recursive, in-place* `resolve_chirality`, not the old
/// per-atom re-rooted `auxiliary_sign`.
fn break_tie_rule4b(
    graph: &mut CipDigraph,
    mol: &Molecule,
    a: NodeId,
    b: NodeId,
    budget: CipBudget,
    cache: &mut HashMap<NodeId, Option<bool>>,
) -> Result<Option<std::cmp::Ordering>, CipCompareError> {
    use std::cmp::Ordering;
    let chain_a = embedded_chain(graph, mol, a)?;
    let chain_b = embedded_chain(graph, mol, b)?;
    let signs_a: Vec<Option<bool>> = chain_a
        .iter()
        .map(|&n| resolve_chirality(graph, mol, n, budget, cache))
        .collect::<Result<_, _>>()?;
    let signs_b: Vec<Option<bool>> = chain_b
        .iter()
        .map(|&n| resolve_chirality(graph, mol, n, budget, cache))
        .collect::<Result<_, _>>()?;

    let (Some(Some(ref_a)), Some(Some(ref_b))) = (signs_a.first(), signs_b.first()) else {
        return Ok(None);
    };
    for i in 0..signs_a.len().min(signs_b.len()) {
        let (Some(sa), Some(sb)) = (signs_a[i], signs_b[i]) else {
            break;
        };
        let like_a = sa == *ref_a;
        let like_b = sb == *ref_b;
        if like_a != like_b {
            return Ok(Some(if like_a {
                Ordering::Greater
            } else {
                Ordering::Less
            }));
        }
    }
    Ok(None)
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
    let mut r = remaining_ranks.clone();
    let target = [2u8, 3, 4];
    let mut swaps = 0usize;
    for i in 0..3 {
        if r[i] != target[i] {
            let j_rel = r[i + 1..].iter().position(|&x| x == target[i])?;
            r.swap(i, j_rel + i + 1);
            swaps += 1;
        }
    }
    Some(cw_from_lowest ^ (swaps % 2 == 1))
}

/// `stereo_order`'s 4 physical positions matched against `candidates` (forward children
/// plus the back-to-root parent node) by underlying `AtomIdx`.
fn position_node_ids(
    graph: &CipDigraph,
    stereo_order: &[u32],
    candidates: &[NodeId],
) -> Option<Vec<NodeId>> {
    let mut result = Vec::with_capacity(stereo_order.len());
    for &pos_val in stereo_order {
        let node_id = if pos_val == STEREO_H_SENTINEL {
            candidates
                .iter()
                .copied()
                .find(|&id| matches!(graph.node(id).kind, CipNodeKind::ImplicitHydrogen))?
        } else {
            let atom_idx = AtomIdx(pos_val);
            candidates
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

/// The core recursive, memoized, in-place resolver. `node_id` must be a real `Atom`
/// node (not root -- the root has no "back to root" ligand and is handled by the
/// caller, `resolve_outer_root`, using the plain 4-forward-children path).
fn resolve_chirality(
    graph: &mut CipDigraph,
    mol: &Molecule,
    node_id: NodeId,
    budget: CipBudget,
    cache: &mut HashMap<NodeId, Option<bool>>,
) -> Result<Option<bool>, CipCompareError> {
    if let Some(&cached) = cache.get(&node_id) {
        return Ok(cached);
    }

    let atom_idx = match graph.node(node_id).kind {
        CipNodeKind::Atom { atom_idx } => atom_idx,
        _ => {
            cache.insert(node_id, None);
            return Ok(None);
        }
    };
    if mol.atom(atom_idx).chirality == Chirality::None {
        cache.insert(node_id, None);
        return Ok(None);
    }
    let Some(parent_id) = graph.node(node_id).parent else {
        // True digraph root -- use resolve_outer_root instead.
        cache.insert(node_id, None);
        return Ok(None);
    };

    let forward_children = graph
        .expand_children(node_id)
        .map_err(CipCompareError::Digraph)?;
    let mut ctx = CompareContext::new();
    let mut groups = rank_children(graph, &forward_children, &mut ctx)?;

    // Rule 4b tiebreak: only a clean 2-way tie is handled (scope match with every case
    // seen in this residual family). 3+-way ties are reported, not guessed.
    for gi in 0..groups.len() {
        if groups[gi].len() == 2 {
            let (a, b) = (groups[gi][0], groups[gi][1]);
            if let Some(ord) = break_tie_rule4b(graph, mol, a, b, budget, cache)? {
                let (higher, lower) = match ord {
                    std::cmp::Ordering::Greater => (a, b),
                    std::cmp::Ordering::Less => (b, a),
                    std::cmp::Ordering::Equal => (a, b),
                };
                if ord != std::cmp::Ordering::Equal {
                    groups[gi] = vec![lower];
                    groups.insert(gi, vec![higher]);
                }
            }
        } else if groups[gi].len() > 2 {
            cache.insert(node_id, None);
            return Ok(None); // out of scope this round
        }
    }

    // Insert the back-to-root ligand at its Rule-1a-only rank among `groups`. Compare
    // sphere-1-to-sphere-1: the back ligand's sphere 1 is `ancestor` itself; the
    // forward ligand's sphere 1 is `rep` itself -- NOT `rep`'s children (that would be
    // off by one sphere).
    let mut insert_at = groups.len();
    for (gi, group) in groups.iter().enumerate() {
        let rep = group[0];
        let back_frontier = vec![BackItem::Ascending {
            ancestor: parent_id,
            came_from: node_id,
        }];
        let ord = compare_rule1a_only(graph, back_frontier, vec![rep])?;
        if ord == std::cmp::Ordering::Greater {
            insert_at = gi;
            break;
        }
    }
    let mut final_groups = groups.clone();
    final_groups.insert(insert_at, vec![parent_id]);
    if DEBUG {
        eprintln!(
            "resolve_chirality({}): parent={} forward_groups={:?} insert_at={} final_groups={:?}",
            atom_of_node(graph, node_id),
            atom_of_node(graph, parent_id),
            groups
                .iter()
                .map(|g| g
                    .iter()
                    .map(|&n| atom_of_node(graph, n))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            insert_at,
            final_groups
                .iter()
                .map(|g| g
                    .iter()
                    .map(|&n| atom_of_node(graph, n))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        );
    }

    let stereo_order = match mol.stereo_neighbor_order(atom_idx) {
        Some(o) => o,
        None => {
            cache.insert(node_id, None);
            return Ok(None);
        }
    };
    let mut candidates = forward_children.clone();
    candidates.push(parent_id);
    let position_nodes = match position_node_ids(graph, stereo_order, &candidates) {
        Some(p) => p,
        None => {
            cache.insert(node_id, None);
            return Ok(None);
        }
    };
    let result =
        resolve_is_r_from_groups(&final_groups, &position_nodes, mol.atom(atom_idx).chirality);
    if DEBUG {
        eprintln!(
            "  -> resolve_chirality(atom{}) = {:?}",
            atom_idx.0,
            result.map(|r| if r { "R" } else { "S" })
        );
    }
    cache.insert(node_id, result);
    Ok(result)
}

/// The true outer root has 4 real forward children and no "back to root" ligand at
/// all -- the ordinary case, just with the same Rule-4b-aware `rank_children`.
fn resolve_outer_root(
    graph: &mut CipDigraph,
    mol: &Molecule,
    root: NodeId,
    root_atom: AtomIdx,
    budget: CipBudget,
    cache: &mut HashMap<NodeId, Option<bool>>,
) -> Result<Option<CipCode>, CipCompareError> {
    let children = graph
        .expand_children(root)
        .map_err(CipCompareError::Digraph)?;
    let mut ctx = CompareContext::new();
    let mut groups = rank_children(graph, &children, &mut ctx)?;

    for gi in 0..groups.len() {
        if groups[gi].len() == 2 {
            let (a, b) = (groups[gi][0], groups[gi][1]);
            if let Some(ord) = break_tie_rule4b(graph, mol, a, b, budget, cache)? {
                let (higher, lower) = match ord {
                    std::cmp::Ordering::Greater => (a, b),
                    _ => (b, a),
                };
                if ord != std::cmp::Ordering::Equal {
                    groups[gi] = vec![lower];
                    groups.insert(gi, vec![higher]);
                }
            }
        } else if groups[gi].len() > 2 {
            return Ok(None);
        }
    }

    if DEBUG {
        eprintln!(
            "resolve_outer_root(atom{}): final_groups={:?}",
            root_atom.0,
            groups
                .iter()
                .map(|g| g
                    .iter()
                    .map(|&n| atom_of_node(graph, n))
                    .collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        );
    }
    let Some(stereo_order) = mol.stereo_neighbor_order(root_atom) else {
        return Ok(None);
    };
    let Some(position_nodes) = position_node_ids(graph, stereo_order, &children) else {
        return Ok(None);
    };
    Ok(
        resolve_is_r_from_groups(&groups, &position_nodes, mol.atom(root_atom).chirality)
            .map(|is_r| if is_r { CipCode::R } else { CipCode::S }),
    )
}

fn run_case(
    smiles: &str,
    atom_idx: u32,
    case_id: &str,
    oracle: CipCode,
) -> (Option<CipCode>, bool) {
    let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
    let idx = AtomIdx(atom_idx);
    let budget = CipBudget::default_budget();
    let mut graph = CipDigraph::new(&mol, idx, budget).expect("digraph builds");
    let root = graph.root();
    let mut cache = HashMap::new();
    let predicted =
        resolve_outer_root(&mut graph, &mol, root, idx, budget, &mut cache).unwrap_or(None);
    let ok = predicted == Some(oracle);
    println!(
        "  {case_id:16} oracle={oracle:?} predicted={predicted:?} [{}]",
        if ok { "OK" } else { "MISS" }
    );
    (predicted, ok)
}

/// Textually swaps every `@`/`@@` -- the true enantiomer, whose oracle labels are all
/// flipped R<->S/r<->s (never re-derived from RDKit; this is the same decisive
/// falsification tool that caught the "S precedes R" and "opposite-of-global-code"
/// hollow fits earlier this milestone).
fn mirror_smiles(smiles: &str) -> String {
    let mut out = String::with_capacity(smiles.len());
    let mut chars = smiles.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '@' {
            if chars.peek() == Some(&'@') {
                chars.next();
                out.push('@');
            } else {
                out.push_str("@@");
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn mirror_code(c: CipCode) -> CipCode {
    match c {
        CipCode::R => CipCode::S,
        CipCode::S => CipCode::R,
        other => other,
    }
}

fn run_mirror_check(name: &str, rows: &[(&str, u32, &str, CipCode)]) -> (usize, usize) {
    println!("== {name} mirrored (enantiomer, oracle labels flipped) ==");
    let mut ok = 0;
    for &(smiles, atom_idx, case_id, oracle) in rows {
        let mirrored_smiles = mirror_smiles(smiles);
        let (_, row_ok) = run_case(&mirrored_smiles, atom_idx, case_id, mirror_code(oracle));
        ok += row_ok as usize;
    }
    println!("  -> {ok}/{}\n", rows.len());
    (ok, rows.len())
}

fn main() {
    println!("Milestone 4B-2: bottom-up, in-place, single-digraph Rule 4b engine\n");

    println!("== VS196 (Hanson 2018 suite, tag `4b` only) ==");
    let mut vs196_ok = 0;
    for &(atom_idx, oracle) in VS196_ROWS {
        let (_, ok) = run_case(
            VS196_SMILES,
            atom_idx,
            &format!("VS196-atom{atom_idx}"),
            oracle,
        );
        vs196_ok += ok as usize;
    }
    println!("  -> {vs196_ok}/{}\n", VS196_ROWS.len());

    println!("== VS196 mirrored ==");
    let vs196_mirrored_smiles = mirror_smiles(VS196_SMILES);
    let mut vs196_mirror_ok = 0;
    for &(atom_idx, oracle) in VS196_ROWS {
        let (_, ok) = run_case(
            &vs196_mirrored_smiles,
            atom_idx,
            &format!("VS196m-atom{atom_idx}"),
            mirror_code(oracle),
        );
        vs196_mirror_ok += ok as usize;
    }
    println!("  -> {vs196_mirror_ok}/{}\n", VS196_ROWS.len());

    println!("== VS197 (Hanson 2018 suite, tag `4b` only) ==");
    let mut vs197_ok = 0;
    for &(atom_idx, oracle) in VS197_ROWS {
        let (_, ok) = run_case(
            VS197_SMILES,
            atom_idx,
            &format!("VS197-atom{atom_idx}"),
            oracle,
        );
        vs197_ok += ok as usize;
    }
    println!("  -> {vs197_ok}/{}\n", VS197_ROWS.len());

    println!("== VS197 mirrored ==");
    let vs197_mirrored_smiles = mirror_smiles(VS197_SMILES);
    let mut vs197_mirror_ok = 0;
    for &(atom_idx, oracle) in VS197_ROWS {
        let (_, ok) = run_case(
            &vs197_mirrored_smiles,
            atom_idx,
            &format!("VS197m-atom{atom_idx}"),
            mirror_code(oracle),
        );
        vs197_mirror_ok += ok as usize;
    }
    println!("  -> {vs197_mirror_ok}/{}\n", VS197_ROWS.len());

    println!("== ROWS_ORIGINAL (8 rows, all-S -- regression floor, never proof alone) ==");
    let mut orig_ok = 0;
    for &(smiles, atom_idx, case_id, oracle) in ROWS_ORIGINAL {
        let (_, ok) = run_case(smiles, atom_idx, case_id, oracle);
        orig_ok += ok as usize;
    }
    println!("  -> {orig_ok}/{}\n", ROWS_ORIGINAL.len());
    run_mirror_check("ROWS_ORIGINAL", ROWS_ORIGINAL);

    println!(
        "== ROWS_DISCRIMINATING (16 rows, acyclic, non-degenerate -- mandatory regression) =="
    );
    let mut disc_ok = 0;
    for &(smiles, atom_idx, case_id, oracle) in ROWS_DISCRIMINATING {
        let (_, ok) = run_case(smiles, atom_idx, case_id, oracle);
        disc_ok += ok as usize;
    }
    println!("  -> {disc_ok}/{}\n", ROWS_DISCRIMINATING.len());
    run_mirror_check("ROWS_DISCRIMINATING", ROWS_DISCRIMINATING);
}
