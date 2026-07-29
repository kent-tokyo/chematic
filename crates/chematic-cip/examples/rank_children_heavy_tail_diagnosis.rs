//! CIP-Perf-A1 (issue #107): diagnosis-only instrumentation of the heavy tail
//! `cip_perf_diagnosis.rs` (CIP-Perf-A0) found -- a small minority of stereocenters
//! that resolve "trivially" at Rules 1a/1b/2 do 40-3700x the median comparator work.
//! Answers the issue's own suggested investigation order, in order, using two named
//! worst-case fixtures plus a full-corpus aggregate. Ships **no production behavior
//! change** -- the only non-test-only production edit this required was adding
//! `left_node`/`right_node: NodeId` to `DecisionStep` (trace.rs), a pure data-carrying
//! addition to an already debug-only struct; every counting/aggregation policy lives
//! here.
//!
//! Usage:
//!   cargo run -p chematic-cip --release --example rank_children_heavy_tail_diagnosis
//!   cargo run -p chematic-cip --release --example rank_children_heavy_tail_diagnosis -- <SMILES.csv>

use std::collections::HashMap;
use std::env;
use std::fs;
use std::time::Instant;

use chematic_cip::{
    BranchComparison, CipBudget, CipDigraph, CipNodeKind, CompareContext, ComparisonTrace, NodeId,
    prepare_kekule_form, rank_children,
};
use chematic_core::{AtomIdx, Chirality, Molecule};

/// The two worst-offender fixtures named directly in issue #107.
const WORST_PASS1: &str = "O=C1OCC2OC(=O)c3cc(O)c(O)c(O)c3-c3c(O)c(O)c(O)c4c3C(=O)OC(C2OC(=O)c2cc(O)c(O)c(O)c2-c2c1cc(O)c(O)c2O)C1OC(=O)c2c-4c(O)c(O)c(O)c2[C@@H]1c1c(O)cc(O)c2c1O[C@H](c1ccc(O)c(O)c1)[C@@H](O)C2";
const WORST_NEEDS_MORE: &str =
    "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3";

fn is_leaf_kind(kind: CipNodeKind) -> bool {
    matches!(
        kind,
        CipNodeKind::MultipleBondDuplicate { .. }
            | CipNodeKind::RingDuplicate { .. }
            | CipNodeKind::ImplicitHydrogen
    )
}

/// Builds the digraph for one stereocenter and returns the full trace of every
/// pairwise `compare_ligands` call made while ranking its root's children (mirrors
/// `cip_perf_diagnosis.rs`'s Q2 setup exactly).
fn trace_for_center<'a>(
    mol: &'a Molecule,
    idx: AtomIdx,
    budget: CipBudget,
    kekule: &'a Option<(Molecule, chematic_cip::MancudeContext)>,
) -> Option<(CipDigraph<'a>, ComparisonTrace, u128, u128)> {
    let build_start = Instant::now();
    let mut graph = match kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).ok()?
        }
        None => CipDigraph::new(mol, idx, budget).ok()?,
    };
    let root = graph.root();
    let root_children = graph.expand_children(root).ok()?;
    let build_nanos = build_start.elapsed().as_nanos();

    // Finding 5: cost split. Fully materialize the whole reachable digraph *before*
    // timing the comparison pass, so `expand_children`'s lazy node-materialization
    // cost (molecule traversal, Kekule-context lookups, duplicate-node synthesis) is
    // charged to construction, not to comparison -- every node touched by
    // `rank_children`'s recursive descent is already resident by the time comparison
    // starts.
    let expand_start = Instant::now();
    let _ = graph.expand_all(root);
    let expand_nanos = expand_start.elapsed().as_nanos();

    let mut trace = ComparisonTrace::new(root, root);
    let mut ctx = CompareContext::with_trace(&mut trace);
    let compare_start = Instant::now();
    rank_children(&mut graph, &root_children, &mut ctx).ok()?;
    let compare_nanos = compare_start.elapsed().as_nanos();

    Some((graph, trace, build_nanos + expand_nanos, compare_nanos))
}

/// Unordered key so (a, b) and (b, a) collapse to one bucket.
fn unordered<T: Ord + Copy>(a: T, b: T) -> (T, T) {
    if a <= b { (a, b) } else { (b, a) }
}

/// `NodeId` doesn't derive `Ord` (it's an opaque arena index, not meant to be
/// compared for anything but identity elsewhere in the crate) -- order on the raw
/// `u32` here, purely for this diagnostic's own unordered-pair bucketing.
fn unordered_node(a: NodeId, b: NodeId) -> (NodeId, NodeId) {
    if a.0 <= b.0 { (a, b) } else { (b, a) }
}

fn invert_outcome(cmp: BranchComparison) -> BranchComparison {
    match cmp {
        BranchComparison::Higher => BranchComparison::Lower,
        BranchComparison::Lower => BranchComparison::Higher,
        BranchComparison::Equal => BranchComparison::Equal,
        BranchComparison::Unresolved => BranchComparison::Unresolved,
    }
}

struct HeavyTailReport {
    total_comparisons: usize,
    distinct_nodeid_pairs: usize,
    max_nodeid_pair_repeats: usize,
    distinct_signature_pairs: usize,
    max_signature_pair_repeats: usize,
    comparisons_saved_if_signature_memoized: usize,
    // `branch_signature` hashes atomic number + isotope only -- it does not encode
    // MANCUDE fractional atomic numbers or which ancestor a `RingDuplicate` closes
    // back to, both of which are real inputs to `compare_ligands`'s actual outcome.
    // Equal signatures are therefore a *candidate* cache key, not a proven-safe one.
    // These two counts are the actual discriminating check: among repeated
    // (signature, signature) buckets, how many always produced the same canonical
    // outcome (a signature-keyed cache would have been correct) vs. produced more
    // than one distinct outcome across their repeats (a signature-keyed cache would
    // have been wrong at least once).
    signature_buckets_with_repeats: usize,
    signature_buckets_outcome_homogeneous: usize,
    signature_buckets_outcome_mixed: usize,
    comparisons_in_mixed_buckets: usize,
    leaf_involved_comparisons: usize,
    both_leaf_comparisons: usize,
    distinct_ranking_parents: usize,
    parents_sharing_a_signature_set_with_another_parent: usize,
}

/// Runs findings 1, 2, 3, and 4 against one already-collected trace.
fn analyze(graph: &mut CipDigraph, trace: &ComparisonTrace) -> HeavyTailReport {
    let total_comparisons = trace.decisions.len();

    // Finding 1: same-NodeId-pair re-comparison.
    let mut nodeid_pair_counts: HashMap<(NodeId, NodeId), usize> = HashMap::new();
    for step in &trace.decisions {
        let key = unordered_node(step.left_node, step.right_node);
        *nodeid_pair_counts.entry(key).or_default() += 1;
    }
    let distinct_nodeid_pairs = nodeid_pair_counts.len();
    let max_nodeid_pair_repeats = nodeid_pair_counts.values().copied().max().unwrap_or(0);

    // Finding 2: isomorphic-subtree-pair re-comparison, via the already-existing,
    // already-correct `branch_signature` (order/numbering-invariant structural hash).
    // Memoize per NodeId -- each node's signature is queried at most once regardless
    // of how many comparisons involve it.
    let mut sig_cache: HashMap<NodeId, u64> = HashMap::new();
    let mut signature_of = |graph: &mut CipDigraph, n: NodeId| -> u64 {
        *sig_cache
            .entry(n)
            .or_insert_with(|| graph.branch_signature(n).unwrap_or(0))
    };
    let mut signature_pair_counts: HashMap<(u64, u64), usize> = HashMap::new();
    // Canonical (order-normalized) outcomes observed per signature-pair bucket, to
    // check whether `branch_signature` alone is actually a safe cache key (see the
    // report struct's doc comment) rather than assuming it from the raw repeat count.
    let mut signature_pair_outcomes: HashMap<(u64, u64), Vec<BranchComparison>> = HashMap::new();
    let mut leaf_involved_comparisons = 0usize;
    let mut both_leaf_comparisons = 0usize;
    let mut parent_signature_sets: HashMap<NodeId, Vec<u64>> = HashMap::new();
    for step in &trace.decisions {
        let ls = signature_of(graph, step.left_node);
        let rs = signature_of(graph, step.right_node);
        *signature_pair_counts.entry(unordered(ls, rs)).or_default() += 1;
        let canonical_outcome = if ls <= rs {
            step.outcome
        } else {
            invert_outcome(step.outcome)
        };
        signature_pair_outcomes
            .entry(unordered(ls, rs))
            .or_default()
            .push(canonical_outcome);

        let left_leaf = is_leaf_kind(graph.node(step.left_node).kind);
        let right_leaf = is_leaf_kind(graph.node(step.right_node).kind);
        if left_leaf || right_leaf {
            leaf_involved_comparisons += 1;
        }
        if left_leaf && right_leaf {
            both_leaf_comparisons += 1;
        }

        if let Some(parent) = step.ranking_parent {
            let set = parent_signature_sets.entry(parent).or_default();
            if !set.contains(&ls) {
                set.push(ls);
            }
            if !set.contains(&rs) {
                set.push(rs);
            }
        }
    }
    let distinct_signature_pairs = signature_pair_counts.len();
    let max_signature_pair_repeats = signature_pair_counts.values().copied().max().unwrap_or(0);
    // Upper bound only -- see the report struct's doc comment. Real savings depend on
    // `signature_buckets_outcome_mixed` being zero; a mixed bucket means this key is
    // not by itself a correct cache key for those repeats.
    let comparisons_saved_if_signature_memoized = signature_pair_counts
        .values()
        .map(|&c| c.saturating_sub(1))
        .sum();

    let mut signature_buckets_with_repeats = 0usize;
    let mut signature_buckets_outcome_homogeneous = 0usize;
    let mut signature_buckets_outcome_mixed = 0usize;
    let mut comparisons_in_mixed_buckets = 0usize;
    for (key, count) in &signature_pair_counts {
        if *count <= 1 {
            continue;
        }
        signature_buckets_with_repeats += 1;
        let outcomes = &signature_pair_outcomes[key];
        let distinct: std::collections::HashSet<_> =
            outcomes.iter().map(|o| format!("{o:?}")).collect();
        if distinct.len() <= 1 {
            signature_buckets_outcome_homogeneous += 1;
        } else {
            signature_buckets_outcome_mixed += 1;
            comparisons_in_mixed_buckets += count;
        }
    }

    // Finding 3: does the same *set* of sibling signatures recur across different
    // `rank_children` calls (different parents) within this one center's resolution?
    // If so, the whole pairwise matrix for that sibling group is redundant work the
    // second (and later) time it's built, not just individual pair comparisons.
    for sigs in parent_signature_sets.values_mut() {
        sigs.sort_unstable();
    }
    let mut set_counts: HashMap<Vec<u64>, usize> = HashMap::new();
    for sigs in parent_signature_sets.values() {
        if sigs.len() > 1 {
            *set_counts.entry(sigs.clone()).or_default() += 1;
        }
    }
    let parents_sharing_a_signature_set_with_another_parent =
        set_counts.values().filter(|&&c| c > 1).copied().sum();

    HeavyTailReport {
        total_comparisons,
        distinct_nodeid_pairs,
        max_nodeid_pair_repeats,
        distinct_signature_pairs,
        max_signature_pair_repeats,
        comparisons_saved_if_signature_memoized,
        signature_buckets_with_repeats,
        signature_buckets_outcome_homogeneous,
        signature_buckets_outcome_mixed,
        comparisons_in_mixed_buckets,
        leaf_involved_comparisons,
        both_leaf_comparisons,
        distinct_ranking_parents: parent_signature_sets.len(),
        parents_sharing_a_signature_set_with_another_parent,
    }
}

fn print_report(label: &str, r: &HeavyTailReport, build_nanos: u128, compare_nanos: u128) {
    println!("--- {label} ---");
    println!(
        "total_comparisons={}  build+expand_all={:.2}ms  rank_children(compare only)={:.2}ms",
        r.total_comparisons,
        build_nanos as f64 / 1_000_000.0,
        compare_nanos as f64 / 1_000_000.0,
    );
    println!(
        "[1] distinct (NodeId,NodeId) pairs compared={}  max repeats of one pair={}  \
         -- expected/structural: `rank_children`'s own pairwise fill visits each (i,j) \
         exactly once per call, so a repeat here means the SAME descendant pair was \
         reached again from a DIFFERENT ancestor comparison elsewhere in the digraph, \
         not a bug in the fill loop itself. This is the same phenomenon [2] measures \
         with a coarser (isomorphism, not identity) key.",
        r.distinct_nodeid_pairs, r.max_nodeid_pair_repeats,
    );
    println!(
        "[2] distinct isomorphic (signature,signature) pairs={}  max repeats={}  \
         UPPER BOUND comparisons_saved_if_memoized={} ({:.1}% of total) -- \
         branch_signature hashes atomic number + isotope only, NOT MANCUDE fractional \
         atomic numbers or ring-closure ancestor identity, so equal signatures are a \
         *candidate* cache key, not a proven-safe one. Discriminating check: of {} \
         signature-pair buckets that repeat, {} always produced the SAME canonical \
         outcome (safe to memoize) and {} produced MORE THAN ONE distinct outcome \
         across their repeats ({} comparisons live in those unsafe buckets -- a \
         signature-only cache would have returned a wrong answer for at least one of \
         them).",
        r.distinct_signature_pairs,
        r.max_signature_pair_repeats,
        r.comparisons_saved_if_signature_memoized,
        100.0 * r.comparisons_saved_if_signature_memoized as f64
            / r.total_comparisons.max(1) as f64,
        r.signature_buckets_with_repeats,
        r.signature_buckets_outcome_homogeneous,
        r.signature_buckets_outcome_mixed,
        r.comparisons_in_mixed_buckets,
    );
    println!(
        "[3] distinct rank_children calls (ranking_parents)={}  parents whose own \
         sibling-signature SET recurs elsewhere in this trace={}",
        r.distinct_ranking_parents, r.parents_sharing_a_signature_set_with_another_parent,
    );
    println!(
        "[4] comparisons with >=1 leaf-kind (dup/H) operand={} ({:.1}%)  both-leaf={} \
         ({:.1}%) -- both-leaf comparisons are decidable by atomic number alone, no \
         recursion needed",
        r.leaf_involved_comparisons,
        100.0 * r.leaf_involved_comparisons as f64 / r.total_comparisons.max(1) as f64,
        r.both_leaf_comparisons,
        100.0 * r.both_leaf_comparisons as f64 / r.total_comparisons.max(1) as f64,
    );
    println!();
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let budget = CipBudget::default_budget();

    println!("=== Named worst-case fixtures (from issue #107 directly) ===\n");
    for (label, smi) in [
        (
            "worst pass1_resolved (89,250 comparisons reported)",
            WORST_PASS1,
        ),
        (
            "worst needs_pass2_or_3 (36,198 comparisons reported)",
            WORST_NEEDS_MORE,
        ),
    ] {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(e) => {
                println!("{label}: FAILED TO PARSE ({e:?}) -- skipping\n");
                continue;
            }
        };
        let kekule = prepare_kekule_form(&mol).ok();
        let mut found_any = false;
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if mol.atom(idx).chirality == Chirality::None {
                continue;
            }
            let Some(stereo_order) = mol.stereo_neighbor_order(idx) else {
                continue;
            };
            if stereo_order.len() != 4 {
                continue;
            }
            let Some((mut graph, trace, build_nanos, compare_nanos)) =
                trace_for_center(&mol, idx, budget, &kekule)
            else {
                continue;
            };
            // Only the worst center per molecule is interesting for a targeted trace.
            if trace.decisions.len() < 1000 {
                continue;
            }
            found_any = true;
            let report = analyze(&mut graph, &trace);
            print_report(
                &format!("{label} [atom {}]", idx.0),
                &report,
                build_nanos,
                compare_nanos,
            );
        }
        if !found_any {
            println!(
                "{label}: no stereocenter in this molecule crossed the 1000-comparison threshold this run (budget/heuristics may differ from the original A0 run) -- skipping\n"
            );
        }
    }

    // Full-corpus aggregate, same denominator as cip_perf_diagnosis.rs, if a corpus
    // path was given -- otherwise the two named fixtures above already answer the
    // issue's investigation order on its own worst-case evidence.
    let Some(csv_path) = args.get(1) else {
        println!(
            "(no corpus path given -- pass `~/Downloads/SMILES.csv` or similar as arg 1 \
             for a full-corpus aggregate; named-fixture findings above stand on their own)"
        );
        return;
    };
    let content = match fs::read_to_string(csv_path) {
        Ok(c) => c,
        Err(e) => {
            println!("could not read corpus at {csv_path}: {e}");
            return;
        }
    };
    let smis: Vec<&str> = content
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .collect();

    let mut agg_total_comparisons = 0u64;
    let mut agg_saved_if_memoized = 0u64;
    let mut agg_leaf_involved = 0u64;
    let mut agg_both_leaf = 0u64;
    let mut agg_build_nanos: u128 = 0;
    let mut agg_compare_nanos: u128 = 0;
    let mut centers_examined = 0u64;
    // Homogeneity check, aggregated per-center (matching the realistic cache scope:
    // one CompareContext/graph/MancudeContext per stereocenter resolution -- issue
    // #107's own "First implementation candidate" already scopes the cache key to
    // include MancudeContext identity, i.e. no cross-molecule/cross-center cache is
    // being proposed, so pooling buckets *within* one center's own trace, then summing
    // across centers, is the correct check -- not pooling signatures globally).
    let mut agg_buckets_with_repeats = 0u64;
    let mut agg_buckets_homogeneous = 0u64;
    let mut agg_buckets_mixed = 0u64;
    let mut agg_comparisons_in_mixed_buckets = 0u64;
    let mut centers_with_any_mixed_bucket = 0u64;

    for smi in &smis {
        let Ok(mol) = chematic_smiles::parse(smi) else {
            continue;
        };
        let has_chirality =
            (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None);
        if !has_chirality {
            continue;
        }
        let kekule = prepare_kekule_form(&mol).ok();
        for i in 0..mol.atom_count() {
            let idx = AtomIdx(i as u32);
            if mol.atom(idx).chirality == Chirality::None {
                continue;
            }
            let Some(stereo_order) = mol.stereo_neighbor_order(idx) else {
                continue;
            };
            if stereo_order.len() != 4 {
                continue;
            }
            let Some((mut graph, trace, build_nanos, compare_nanos)) =
                trace_for_center(&mol, idx, budget, &kekule)
            else {
                continue;
            };
            centers_examined += 1;
            agg_build_nanos += build_nanos;
            agg_compare_nanos += compare_nanos;
            let report = analyze(&mut graph, &trace);
            agg_total_comparisons += report.total_comparisons as u64;
            agg_saved_if_memoized += report.comparisons_saved_if_signature_memoized as u64;
            agg_leaf_involved += report.leaf_involved_comparisons as u64;
            agg_both_leaf += report.both_leaf_comparisons as u64;
            agg_buckets_with_repeats += report.signature_buckets_with_repeats as u64;
            agg_buckets_homogeneous += report.signature_buckets_outcome_homogeneous as u64;
            agg_buckets_mixed += report.signature_buckets_outcome_mixed as u64;
            agg_comparisons_in_mixed_buckets += report.comparisons_in_mixed_buckets as u64;
            if report.signature_buckets_outcome_mixed > 0 {
                centers_with_any_mixed_bucket += 1;
            }
        }
    }

    println!("=== Full-corpus aggregate (centers_examined={centers_examined}) ===");
    println!(
        "total_comparisons={agg_total_comparisons}  \
         UPPER BOUND comparisons_saved_if_signature_memoized={agg_saved_if_memoized} ({:.1}%)  \
         leaf_involved={agg_leaf_involved} ({:.1}%)  both_leaf={agg_both_leaf} ({:.1}%)",
        100.0 * agg_saved_if_memoized as f64 / agg_total_comparisons.max(1) as f64,
        100.0 * agg_leaf_involved as f64 / agg_total_comparisons.max(1) as f64,
        100.0 * agg_both_leaf as f64 / agg_total_comparisons.max(1) as f64,
    );
    println!(
        "[2, corpus-wide] discriminating check, summed per-center (matching the \
         realistic per-resolution cache scope, not a pooled-across-molecules cache): \
         {agg_buckets_with_repeats} repeating signature-pair buckets total, \
         {agg_buckets_homogeneous} homogeneous (safe) / {agg_buckets_mixed} mixed \
         (unsafe) -- {agg_comparisons_in_mixed_buckets} comparisons \
         ({:.2}% of total_comparisons) live in a mixed bucket somewhere. \
         {centers_with_any_mixed_bucket}/{centers_examined} centers have >=1 mixed \
         bucket.",
        100.0 * agg_comparisons_in_mixed_buckets as f64 / agg_total_comparisons.max(1) as f64,
    );
    println!(
        "[5] construction (digraph build + expand_all) = {:.1}ms total   \
         comparison (rank_children pairwise fill on an already-fully-expanded graph) = {:.1}ms total   \
         ratio compare/construct = {:.2}x   \
         CAVEAT: 'construction' here forces expand_all() on the WHOLE reachable digraph \
         before timing comparison, specifically so lazy materialization isn't charged \
         to the comparison timer. Production's actual expand_children() is lazy and \
         only ever materializes nodes the comparator truly visits -- this split is 'cost \
         if eagerly built' vs 'compare cost given an already-fully-built graph', NOT two \
         disjoint slices of one real production run. Do not read the ratio as \
         'construction dominates production cost' without measuring lazy production \
         cost separately.",
        agg_build_nanos as f64 / 1_000_000.0,
        agg_compare_nanos as f64 / 1_000_000.0,
        agg_compare_nanos as f64 / agg_build_nanos.max(1) as f64,
    );
}
