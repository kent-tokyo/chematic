//! CIP-Perf-A0: diagnosis-only structural profiling of the AccurateExperimental CIP
//! engine's Rules-1a/1b/2 comparator. Motivated by RDKit PR #9171 ("CIP labeller
//! performance: Don't calculate auxiliary descriptors unnecessarily", merged
//! 2026-05-06): RDKit was computing Rule-4b-equivalent auxiliary stereodescriptors
//! for *every* center, when only centers still tied after the constitutional rules
//! need them -- PR #9171's own headline case (PDB 4AXM) went 14s -> 0.036s. Verified
//! against the PR body directly (`gh pr view 9171 --repo rdkit/rdkit`), not paraphrased
//! secondhand. This diagnostic checks whether chematic has the same bug.
//!
//! Ships **no production changes**: every call below is existing `pub` API, the same
//! combination `trace_report.rs` already exercises (`CipDigraph` + `CompareContext` +
//! `rank_children`) plus the two existing corpus-facing entry points. Answers two
//! questions, per stereocenter:
//!
//! 1. Resolution level -- `_without_mancude` is used here as a proxy for whether a
//!    center is resolved during the constitutional-rules (1a/1b/2) pass, or needs the
//!    live engine's later Rule 4b/5 passes -- both already gated on `SkipReason::Tied`
//!    from the prior pass (see `assign.rs`/`resolver.rs`), i.e. chematic already gates
//!    auxiliary-descriptor computation at the center level, the same granularity
//!    #9171 fixed in RDKit. The proxy does **not** require `fractional_decisions` to
//!    be zero: a MANCUDE fraction may change a *local* recursive comparison somewhere
//!    in Pass 1's tree without changing whether the center resolves or what its final
//!    label is -- Q3 (below) verifies this in-run and classifies nonzero counts by
//!    final-label impact, rather than treating any nonzero count as invalidating this
//!    split. See `examples/mancude_decision_diagnosis.rs` and
//!    `docs/rfcs/cip_accurate_rfc.md`'s MANCUDE-Decision-A0 entry for the full
//!    classification methodology and result (D=36/E=0 on this corpus). Rule 4b and
//!    Rule 5 are **not** split apart here (both `pub(crate)`, no Pass-1+4b-only public
//!    entry point) -- reported together as `resolved_pass2_or_3`; a future run could
//!    split them with a small `pub(crate)` exposure if that distinction becomes
//!    decision-relevant.
//! 2. Comparator size -- for the Rules-1a/1b/2 ranking itself, how many digraph nodes
//!    get materialized (and what fraction are `MultipleBondDuplicate`/`RingDuplicate`
//!    phantom nodes, not real atoms), how many pairwise comparisons `rank_children`
//!    makes, and how deep recursion goes? `rank_children` computes a full pairwise
//!    matrix up front for every sibling group at every visited depth (by design, see
//!    its own module docs) -- this measures whether that matrix stays small in
//!    practice or is a real cost center, split by whether the center turned out to be
//!    Pass-1-trivial.
//! 3. Fractional-decision accounting -- accumulates `CompareContext::fractional_decisions`
//!    (a fraction actually deciding a ranking, not just being present) across every
//!    Q2 comparison. A nonzero value triggers D/E classification (see
//!    `mancude_decision_diagnosis.rs`); it does not by itself invalidate Q1's
//!    resolution-level accounting -- printed here for corpus-scale awareness, not as
//!    a pass/fail check on this tool's own output.
//!
//! elapsed_us is corroborating only, not authoritative -- single-threaded local run,
//! read alongside the structural counts, not as a regression gate (see this project's
//! own criterion-gate pseudo-replication finding, issue #70).
//!
//! Usage:
//!   cargo run -p chematic-cip --release --example cip_perf_diagnosis -- <SMILES.csv>

use std::env;
use std::fs;
use std::time::Instant;

use chematic_cip::{
    CipBudget, CipDigraph, CipNodeKind, CompareContext, ComparisonTrace,
    assign_cip_accurate_experimental, assign_cip_accurate_experimental_without_mancude,
    prepare_kekule_form, rank_children,
};
use chematic_core::{AtomIdx, Chirality};

#[derive(Default)]
struct SizeStats {
    nodes: Vec<usize>,
    duplicate_nodes: Vec<usize>,
    comparisons: Vec<usize>,
    max_depth: Vec<u32>,
}

impl SizeStats {
    fn push(&mut self, nodes: usize, duplicate_nodes: usize, comparisons: usize, max_depth: u32) {
        self.nodes.push(nodes);
        self.duplicate_nodes.push(duplicate_nodes);
        self.comparisons.push(comparisons);
        self.max_depth.push(max_depth);
    }

    fn summarize(&self, label: &str) {
        let n = self.nodes.len();
        if n == 0 {
            println!("{label}: n=0");
            return;
        }
        let total_nodes: usize = self.nodes.iter().sum();
        let total_dupes: usize = self.duplicate_nodes.iter().sum();
        println!(
            "{label}: n={n}  nodes[median/p95/max]={}/{}/{}  dup_nodes%={:.1}  comparisons[median/p95/max]={}/{}/{}  depth[median/p95/max]={}/{}/{}",
            pct(&self.nodes, 50),
            pct(&self.nodes, 95),
            self.nodes.iter().max().unwrap(),
            100.0 * total_dupes as f64 / total_nodes.max(1) as f64,
            pct(&self.comparisons, 50),
            pct(&self.comparisons, 95),
            self.comparisons.iter().max().unwrap(),
            pct(&self.max_depth, 50),
            pct(&self.max_depth, 95),
            self.max_depth.iter().max().unwrap(),
        );
    }
}

fn pct<T: Ord + Copy>(values: &[T], p: usize) -> T {
    let mut sorted = values.to_vec();
    sorted.sort_unstable();
    let idx = (sorted.len() * p / 100).min(sorted.len() - 1);
    sorted[idx]
}

/// The `SMILES.csv` this tool's own aggregate counts (99.3%/0.7% resolution split,
/// comparator-size percentiles) were measured against and quoted in
/// `docs/rfcs/cip_accurate_rfc.md`'s MANCUDE-Decision-A0 entry -- see [`corpus_sha256`]'s
/// call site in `main` for the runtime check against whatever corpus is actually
/// passed in.
const EXPECTED_CORPUS_SHA256: &str =
    "1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e";

/// SHA-256 of `path`'s contents, via the platform `shasum` binary -- `sha2` is not a
/// workspace dependency and this is a one-shot diagnostic tool, not production code.
fn corpus_sha256(path: &str) -> String {
    use std::process::Command;
    let output = Command::new("shasum")
        .arg("-a")
        .arg("256")
        .arg(path)
        .output()
        .expect("run `shasum -a 256` (required to verify corpus identity)");
    String::from_utf8(output.stdout)
        .expect("utf8 shasum output")
        .split_whitespace()
        .next()
        .expect("shasum hash field")
        .to_string()
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| format!("{}/Downloads/SMILES.csv", env::var("HOME").unwrap()));

    let content = fs::read_to_string(&csv_path).expect("read SMILES.csv");
    let actual_sha256 = corpus_sha256(&csv_path);
    println!("corpus: {csv_path}");
    println!("corpus_sha256: {actual_sha256}");
    if actual_sha256 == EXPECTED_CORPUS_SHA256 {
        println!(
            "classification A (corpus drift): RULED OUT -- corpus_sha256 matches this \
             tool's recorded MANCUDE-Decision-A0 closeout value exactly."
        );
    } else {
        println!(
            "*** WARNING: corpus_sha256 does NOT match the expected \
             {EXPECTED_CORPUS_SHA256} -- the aggregate counts below are NOT directly \
             comparable to docs/rfcs/cip_accurate_rfc.md's recorded closeout numbers. ***"
        );
    }
    println!();
    let smis: Vec<&str> = content
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .collect();

    let budget = CipBudget::default_budget();

    let mut resolved_pass1 = 0usize;
    let mut resolved_pass2_or_3 = 0usize;
    let mut still_tied = 0usize;
    let mut skip_other = 0usize;

    let mut stats_pass1 = SizeStats::default();
    let mut stats_needs_more = SizeStats::default();

    let mut total_size_pass_nanos: u128 = 0;
    let mut total_resolution_pass_nanos: u128 = 0;

    // Worst offender by comparison count, tracked separately for the two groups so a
    // single pathological molecule doesn't hide what the *other* group's tail looks
    // like.
    let mut worst_pass1: Option<(usize, String, u32)> = None;
    let mut worst_needs_more: Option<(usize, String, u32)> = None;

    // Q3: does a MANCUDE fraction ever actually decide a ranking on this corpus? A
    // nonzero value triggers D/E classification; it does not by itself invalidate the
    // Q1 resolution-level accounting above (see mancude_decision_diagnosis.rs for the
    // structure-isolated classification that answers whether any final label changes).
    let mut fractional_decisions_total: u64 = 0;

    for smi in &smis {
        let Ok(mol) = chematic_smiles::parse(smi) else {
            continue;
        };
        let has_chirality =
            (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None);
        if !has_chirality {
            continue;
        }

        // Question 1: resolution level (whole-molecule calls, existing pub entry points).
        let start = Instant::now();
        let pass1 = assign_cip_accurate_experimental_without_mancude(&mol, budget);
        let final_result = assign_cip_accurate_experimental(&mol, budget);
        total_resolution_pass_nanos += start.elapsed().as_nanos();

        let (Ok(pass1), Ok(final_result)) = (pass1, final_result) else {
            continue;
        };
        let pass1_resolved: std::collections::HashSet<AtomIdx> =
            pass1.assignments.iter().map(|(idx, _)| *idx).collect();

        for (idx, _) in &final_result.assignments {
            if pass1_resolved.contains(idx) {
                resolved_pass1 += 1;
            } else {
                resolved_pass2_or_3 += 1;
            }
        }
        for (idx, reason) in &final_result.skipped {
            match reason {
                chematic_cip::SkipReason::Tied => still_tied += 1,
                _ => {
                    let _ = idx;
                    skip_other += 1;
                }
            }
        }

        // Question 2: comparator size, per candidate stereocenter (4-substituent atoms
        // only -- mirrors assign.rs's own Pass-1 filter).
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

            let start = Instant::now();
            let mut graph = match &kekule {
                Some((kekule_mol, ctx)) => {
                    match CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx) {
                        Ok(g) => g,
                        Err(_) => continue,
                    }
                }
                None => match CipDigraph::new(&mol, idx, budget) {
                    Ok(g) => g,
                    Err(_) => continue,
                },
            };
            let root = graph.root();
            let Ok(root_children) = graph.expand_children(root) else {
                continue;
            };
            let mut trace = ComparisonTrace::new(root, root);
            let mut ctx = CompareContext::with_trace(&mut trace);
            if rank_children(&mut graph, &root_children, &mut ctx).is_err() {
                continue;
            }
            total_size_pass_nanos += start.elapsed().as_nanos();
            fractional_decisions_total += ctx.fractional_decisions;

            let node_count = graph.nodes().len();
            let duplicate_node_count = graph
                .nodes()
                .iter()
                .filter(|n| {
                    matches!(
                        n.kind,
                        CipNodeKind::MultipleBondDuplicate { .. }
                            | CipNodeKind::RingDuplicate { .. }
                    )
                })
                .count();
            let comparison_count = trace.decisions.len();
            let max_depth = trace.decisions.iter().map(|d| d.depth).max().unwrap_or(0);

            if pass1_resolved.contains(&idx) {
                stats_pass1.push(
                    node_count,
                    duplicate_node_count,
                    comparison_count,
                    max_depth,
                );
                if worst_pass1
                    .as_ref()
                    .is_none_or(|(c, _, _)| comparison_count > *c)
                {
                    worst_pass1 = Some((comparison_count, smi.to_string(), idx.0));
                }
            } else {
                stats_needs_more.push(
                    node_count,
                    duplicate_node_count,
                    comparison_count,
                    max_depth,
                );
                if worst_needs_more
                    .as_ref()
                    .is_none_or(|(c, _, _)| comparison_count > *c)
                {
                    worst_needs_more = Some((comparison_count, smi.to_string(), idx.0));
                }
            }
        }
    }

    let total_candidates = resolved_pass1 + resolved_pass2_or_3 + still_tied;
    println!(
        "=== Q1: resolution level (Pass 1 = Rules 1a/1b/2 alone, via _without_mancude proxy) ==="
    );
    println!(
        "resolved_pass1={resolved_pass1} ({:.1}%)  resolved_pass2_or_3={resolved_pass2_or_3} ({:.1}%)  \
         still_tied_after_all={still_tied} ({:.1}%)  [of {total_candidates} R/S-eligible centers; \
         skip_other(not4/budget)={skip_other} excluded from this denominator]",
        100.0 * resolved_pass1 as f64 / total_candidates.max(1) as f64,
        100.0 * resolved_pass2_or_3 as f64 / total_candidates.max(1) as f64,
        100.0 * still_tied as f64 / total_candidates.max(1) as f64,
    );
    println!(
        "resolution-pass elapsed: {:.1}ms total (corroborating only, not a gate)",
        total_resolution_pass_nanos as f64 / 1_000_000.0
    );
    println!();
    println!(
        "=== Q2: Rules-1a/1b/2 comparator size, split by whether Pass 1 resolved the center ==="
    );
    stats_pass1.summarize("pass1_resolved  ");
    stats_needs_more.summarize("needs_pass2_or_3");
    println!(
        "size-pass elapsed: {:.1}ms total (corroborating only, not a gate)",
        total_size_pass_nanos as f64 / 1_000_000.0
    );
    if let Some((c, smi, atom)) = worst_pass1 {
        println!("worst pass1_resolved   ({c} comparisons): {smi}  atom {atom}");
    }
    if let Some((c, smi, atom)) = worst_needs_more {
        println!("worst needs_pass2_or_3 ({c} comparisons): {smi}  atom {atom}");
    }
    println!();
    println!("=== Q3: fractional-decision accounting (does not by itself invalidate Q1) ===");
    println!(
        "fractional_decisions_total={fractional_decisions_total} -- {}",
        if fractional_decisions_total == 0 {
            "0; Q1's proxy is sound for this run"
        } else {
            "NONZERO -- a fraction was decision-involved somewhere in Pass 1's recursive \
             comparison. This does NOT by itself mean Q1's resolved_pass1 vs \
             resolved_pass2_or_3 split is wrong -- 'decision-involved' and 'changes the \
             resolved label' are different claims (see \
             examples/mancude_decision_diagnosis.rs and \
             docs/rfcs/cip_accurate_rfc.md's MANCUDE-Decision-A0 entry for the full, \
             structure-isolated classification: on this corpus, 0 stereocenters have \
             their final label changed by the fraction)."
        }
    );
}
