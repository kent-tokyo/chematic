//! Milestone 2.5's core gate: every one of the frozen corpus's `uncharacterized`-bucket
//! cases (the one bucket with no known out-of-scope excuse, unlike `aromatic_mancude`/
//! `phosphorus`) that `AccurateExperimental` doesn't resolve correctly must carry an
//! automated, evidence-backed diagnosis -- never left as an unexplained wrong answer or
//! an unexplained tie. This test *is* that gate: each diagnosis is a structural check run
//! against the actual digraph, not a hand-typed guess table, so it fails loudly if a case
//! stops matching its claimed explanation (e.g. after a future milestone changes the
//! digraph's aromatic-bond representation).
//!
//! Two mechanisms account for all 26 cases (2 wrong + 24 tied) as of Milestone 2.5:
//!
//! - **`NeedsLaterSequenceRule`** (all 24 tied cases): every member of the tied group,
//!   when fully expanded, contains at least one `RingDuplicate` node -- confirming the
//!   tie genuinely stems from two branches that are structurally indistinguishable under
//!   Rules 1a/2 alone and only diverge via ring-closure duplicates. Milestone 2.5 guessed
//!   this meant Rule 1b's root-distance duplicate tiebreak would resolve most of these;
//!   Milestone 3A implemented Rule 1b and found it resolves **0/24** -- confirmed both
//!   empirically (zero decisive Rule 1b comparisons across all 8 non-pseudoasymmetric
//!   cases) and structurally (Rule 1a's own child-count check shadows Rule 1b in both
//!   chematic and RDKit; see `compare.rs`'s module docs). Of the 24: 16 have a lowercase
//!   `modern` expected value (pseudoasymmetric centers -- Rule 5, structurally
//!   unrepresentable by this crate's uppercase-only `CipCode` today) and 8 have an
//!   uppercase expected value but remain tied after Rule 1b -- their actual deciding
//!   mechanism (likely Rule 3 and/or Rule 4's auxiliary-descriptor comparison, both
//!   requiring cross-stereocenter information this crate doesn't compute yet) is
//!   unconfirmed and out of scope for this diagnosis.
//! - **`BucketMisclassified`** (both wrong cases): both stereocenters have, among their
//!   substituent branches, an aromatic ring atom that is fully substituted (3 real
//!   neighbors: 2 ring bonds + 1 exocyclic bond, no hydrogen) yet has only 2 digraph
//!   children -- i.e. no duplicate node stands in for the "extra" aromatic bond order a
//!   full mancude-ring CIP treatment would give it. That gap is precisely Milestone 3's
//!   scope (aromatic/mancude ring duplicate representation), not a Rules 1a/1b/2 defect;
//!   these two cases were mis-tagged `uncharacterized` by the corpus's own
//!   `classify_bucket()` heuristic (which evidently keys off whether the *stereocenter
//!   itself* sits on an aromatic ring, not whether one appears one level into a
//!   substituent branch) and structurally belong in `aromatic_mancude`.

use chematic_cip::{
    CipBudget, CipDigraph, CipNodeKind, CompareContext, NodeId, SkipReason,
    assign_cip_accurate_experimental, rank_children,
};
use chematic_core::{AtomIdx, CipCode, Molecule};

const CORPUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cip_label_corpus.jsonl"
));

fn code_str(code: CipCode) -> &'static str {
    match code {
        CipCode::R => "R",
        CipCode::S => "S",
        CipCode::E => "E",
        CipCode::Z => "Z",
        CipCode::LowerR => "r",
        CipCode::LowerS => "s",
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Diagnosis {
    NeedsLaterSequenceRule,
    BucketMisclassified,
}

/// True if `node`'s subtree (expanded up to `budget` nodes visited) contains a
/// `RingDuplicate` anywhere. Stops descending at duplicate/hydrogen leaves.
fn subtree_has_ring_duplicate(graph: &mut CipDigraph, node: NodeId, budget: &mut usize) -> bool {
    *budget = budget.saturating_sub(1);
    if *budget == 0 {
        return false;
    }
    match graph.node(node).kind {
        CipNodeKind::RingDuplicate { .. } => return true,
        CipNodeKind::MultipleBondDuplicate { .. } | CipNodeKind::ImplicitHydrogen => return false,
        CipNodeKind::Atom { .. } => {}
    }
    let Ok(children) = graph.expand_children(node) else {
        return false;
    };
    children
        .into_iter()
        .any(|c| subtree_has_ring_duplicate(graph, c, budget))
}

/// True if `node`'s subtree (up to `budget` nodes) contains an aromatic ring atom that
/// is fully substituted (3 real molecular neighbors: no hydrogen) but has only 2 digraph
/// children -- the "missing mancude duplicate" signature.
fn subtree_has_undercounted_aromatic(
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

fn diagnose_tied(graph: &mut CipDigraph, groups: &[Vec<NodeId>]) -> Option<Diagnosis> {
    for group in groups {
        if group.len() < 2 {
            continue;
        }
        let all_ring_duplicated = group.iter().all(|&n| {
            let mut budget = 5_000;
            subtree_has_ring_duplicate(graph, n, &mut budget)
        });
        if all_ring_duplicated {
            return Some(Diagnosis::NeedsLaterSequenceRule);
        }
    }
    None
}

fn diagnose_wrong(
    mol: &Molecule,
    graph: &mut CipDigraph,
    root_children: &[NodeId],
) -> Option<Diagnosis> {
    let found = root_children.iter().any(|&n| {
        let mut budget = 5_000;
        subtree_has_undercounted_aromatic(mol, graph, n, &mut budget)
    });
    found.then_some(Diagnosis::BucketMisclassified)
}

#[test]
fn all_uncharacterized_cases_are_diagnosed() {
    let mut total = 0usize;
    let mut undiagnosed: Vec<String> = Vec::new();
    let mut tally: std::collections::BTreeMap<&'static str, usize> = Default::default();

    for line in CORPUS.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let value: serde_json::Value = serde_json::from_str(line).unwrap();
        let Some(smiles) = value.get("smiles").and_then(|v| v.as_str()) else {
            continue; // manifest line
        };
        if value.get("bucket").and_then(|v| v.as_str()).is_some() {
            continue; // only the "uncharacterized" (bucket=null) set is this gate's scope
        }
        let atom_idx = value.get("atom_idx").and_then(|v| v.as_u64()).unwrap() as u32;
        let modern = value.get("modern").and_then(|v| v.as_str()).unwrap();
        let idx = AtomIdx(atom_idx);

        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse: {smiles}: {e:?}"));
        let accurate = assign_cip_accurate_experimental(&mol, CipBudget::default_budget())
            .unwrap_or_else(|e| panic!("accurate assignment errored: {smiles}: {e}"));

        total += 1;

        let matched = accurate
            .assignments
            .iter()
            .find(|(i, _)| i.0 == atom_idx)
            .map(|(_, code)| code_str(*code) == modern)
            .unwrap_or(false);
        if matched {
            *tally.entry("OK").or_default() += 1;
            continue;
        }

        let is_wrong = accurate.assignments.iter().any(|(i, _)| i.0 == atom_idx);
        let is_tied = accurate
            .skipped
            .iter()
            .any(|(i, r)| i.0 == atom_idx && *r == SkipReason::Tied);
        assert!(
            is_wrong || is_tied,
            "{smiles} atom {atom_idx}: expected wrong-or-tied, got neither -- outcome \
             bucketing drifted from what this test assumes (budget-exceeded or \
             no-assignment cases in this bucket would need their own diagnosis)"
        );

        let mut graph = CipDigraph::new(&mol, idx, CipBudget::default_budget()).unwrap();
        let root_children = graph.expand_children(graph.root()).unwrap();

        let diagnosis = if is_tied {
            let mut ctx = CompareContext::new();
            let groups = rank_children(&mut graph, &root_children, &mut ctx).unwrap();
            diagnose_tied(&mut graph, &groups)
        } else {
            diagnose_wrong(&mol, &mut graph, &root_children)
        };

        match diagnosis {
            Some(tag) => {
                let label = match tag {
                    Diagnosis::NeedsLaterSequenceRule => "NeedsLaterSequenceRule",
                    Diagnosis::BucketMisclassified => "BucketMisclassified",
                };
                *tally.entry(label).or_default() += 1;
            }
            None => {
                undiagnosed.push(format!(
                    "{smiles} atom {atom_idx}: {} -- no automated diagnosis matched",
                    if is_wrong { "WRONG" } else { "TIED" }
                ));
            }
        }
    }

    println!("\n=== uncharacterized-bucket diagnosis tally ===");
    for (label, count) in &tally {
        println!("  {label:24} {count:3}");
    }
    println!("  {:24} {:3}", "checked", total);

    assert!(
        undiagnosed.is_empty(),
        "Milestone 2.5 gate: every uncharacterized wrong/tied case must carry an \
         automated diagnosis. {} case(s) did not:\n{}",
        undiagnosed.len(),
        undiagnosed.join("\n")
    );
}
