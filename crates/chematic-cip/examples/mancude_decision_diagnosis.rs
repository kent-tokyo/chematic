//! MANCUDE-Decision-A0: classify the nonzero `fractional_decisions` count found while
//! running CIP-Perf-A0 (`cip_perf_diagnosis`) against the frozen `SMILES.csv` corpus
//! (`corpus_sha256=1c47371d...`, confirmed identical to the corpus
//! `docs/rfcs/cip_accurate_rfc.md`'s Milestone 3B-1b closeout used for its "byte-identical,
//! fractional_decisions expected 0" claim -- see that file, and
//! `compare.rs`'s own `CompareContext::fractional_decisions` doc comment, which names a
//! future nonzero value as Milestone 3B-2's own resumption condition #1).
//!
//! Ships **no production changes** -- reuses the same existing `pub` API combination
//! `cip_perf_diagnosis.rs` and `trace_report.rs` already exercise. Breaks the nonzero
//! count down into the four units requested before any classification (A: corpus
//! drift -- already ruled out by the SHA match above; B: instrumentation/counter bug;
//! C: engine evolution since M3B-1b, e.g. Rule 4b/5 reaching new comparisons; D: real
//! decision, but final-label-inert; E: real decision that changes a final label):
//!
//! - comparison events: `fractional_comparisons` / `fractional_decisions` totals
//!   (the same two counters `cip_perf_diagnosis` already reads, reprinted here for a
//!   single source of truth on the raw numbers).
//! - unique stereocenters / unique molecules: how many distinct (molecule, atom_idx)
//!   pairs and distinct molecules have `fractional_decisions > 0` on their own
//!   Rules-1a/1b/2 root-children ranking.
//! - unique `ranking_parent` nodes touched: a **proxy**, not the strict
//!   `fractional_decisions` count -- counts `DecisionStep`s whose recorded rule string
//!   contains a `Rational(...)` key (i.e. a MANCUDE fraction was *involved* in the
//!   decisive comparison, `compare.rs`'s own weaker `fractional_comparisons` notion),
//!   grouped by `DecisionStep::ranking_parent`. This is a proxy because
//!   `fractional_decisions` itself has no per-event node identity recorded anywhere
//!   (a private `u64` counter, incremented deep inside `compare_by_level`) -- adding
//!   that would be a production change to `compare.rs`, out of scope for a
//!   diagnosis-only tool. Reported as an upper bound / concentration signal, not an
//!   exact count.
//! - final-assignment impact, **structure-isolated**: compares
//!   `assign_cip_accurate_experimental_without_mancude` run on the *Kekule-respelled*
//!   molecule (no `MancudeContext` attached -- integer atomic numbers, same digraph
//!   structure the live engine uses) against `assign_cip_accurate_experimental`
//!   (Kekule-respelled + `MancudeContext` + Rule 4b/5). This is the real D-vs-E
//!   classification: holding structure fixed and toggling only the fraction. A
//!   **naive** contrast -- `_without_mancude` on the *original, un-Kekulized* molecule
//!   vs the live engine -- is also reported, labeled clearly as naive: it bundles
//!   Kekule-respelling structure with the fraction, and an earlier version of this
//!   tool used *only* the naive contrast, wrongly classifying 3 centers as E when
//!   they're D (see `docs/rfcs/cip_accurate_rfc.md`'s MANCUDE-Decision-A0 entry for the
//!   full correction). Do not conclude E from the naive numbers alone.
//!
//! elapsed_us is corroborating only, not a gate (see CIP-Perf-A0's own note, and this
//! project's criterion pseudo-replication finding, issue #70).
//!
//! Usage:
//!   cargo run -p chematic-cip --release --example mancude_decision_diagnosis -- <SMILES.csv>

use std::collections::HashSet;
use std::env;
use std::fs;

use chematic_cip::{
    CipBudget, CipDigraph, CompareContext, ComparisonTrace, assign_cip_accurate_experimental,
    assign_cip_accurate_experimental_without_mancude, prepare_kekule_form, rank_children,
};
use chematic_core::{AtomIdx, Chirality, CipCode};

/// The `SMILES.csv` this tool's own aggregate counts (and
/// `docs/rfcs/cip_accurate_rfc.md`'s MANCUDE-Decision-A0 entry, which quotes them) were
/// measured against -- see [`corpus_sha256`]'s call site in `main` for the runtime
/// check against whatever corpus is actually passed in.
const EXPECTED_CORPUS_SHA256: &str =
    "1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e";

fn main() {
    let args: Vec<String> = env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| format!("{}/Downloads/SMILES.csv", env::var("HOME").unwrap()));

    let content = fs::read_to_string(&csv_path).expect("read SMILES.csv");
    let actual_sha256 = corpus_sha256(&csv_path);
    let smis: Vec<&str> = content
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .collect();

    let budget = CipBudget::default_budget();

    let mut total_fractional_comparisons: u64 = 0;
    let mut total_fractional_decisions: u64 = 0;
    let mut centers_with_decisions: usize = 0;
    let mut molecules_with_decisions: HashSet<usize> = HashSet::new();
    let mut proxy_ranking_parents_touched: usize = 0;

    let mut pass1_resolved_final_label_changed: usize = 0;
    let mut pass1_resolved_final_label_same: usize = 0;
    let mut naive_final_label_changed: usize = 0;
    let mut naive_final_label_same: usize = 0;
    let mut changed_rows: Vec<(String, u32, CipCode, CipCode, bool)> = Vec::new();

    let mut worst: Option<(u64, String, u32)> = None;

    for (mol_idx, smi) in smis.iter().enumerate() {
        let Ok(mol) = chematic_smiles::parse(smi) else {
            continue;
        };
        let has_chirality =
            (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None);
        if !has_chirality {
            continue;
        }

        let kekule = prepare_kekule_form(&mol).ok();

        // Naive baseline (bundles Kekule-respelling structure with the fraction --
        // reported for narrative reproducibility only, NOT the classification; see
        // module docs).
        let naive_pass1 = assign_cip_accurate_experimental_without_mancude(&mol, budget);
        // Structure-isolated baseline: same Kekule-respelled molecule the live engine
        // uses, no MancudeContext attached. This is the real "integer-collapsed
        // control" -- structure held fixed, only the fraction toggled.
        let structure_pass1 = match &kekule {
            Some((kekule_mol, _)) => {
                assign_cip_accurate_experimental_without_mancude(kekule_mol, budget)
            }
            None => assign_cip_accurate_experimental_without_mancude(&mol, budget),
        };
        let final_result = assign_cip_accurate_experimental(&mol, budget);
        let (Ok(naive_pass1), Ok(structure_pass1), Ok(final_result)) =
            (naive_pass1, structure_pass1, final_result)
        else {
            continue;
        };
        let naive_pass1_codes: std::collections::HashMap<AtomIdx, CipCode> =
            naive_pass1.assignments.iter().copied().collect();
        let structure_pass1_codes: std::collections::HashMap<AtomIdx, CipCode> =
            structure_pass1.assignments.iter().copied().collect();
        let final_codes: std::collections::HashMap<AtomIdx, CipCode> =
            final_result.assignments.iter().copied().collect();

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
            let Ok(manc_groups) = rank_children(&mut graph, &root_children, &mut ctx) else {
                continue;
            };

            total_fractional_comparisons += ctx.fractional_comparisons;
            total_fractional_decisions += ctx.fractional_decisions;

            if ctx.fractional_decisions > 0 {
                centers_with_decisions += 1;
                molecules_with_decisions.insert(mol_idx);

                if ctx.fractional_decisions > worst.as_ref().map(|(c, _, _)| *c).unwrap_or(0) {
                    worst = Some((ctx.fractional_decisions, smi.to_string(), idx.0));
                }

                // Proxy: DecisionSteps whose rule string names a Rational key, grouped
                // by ranking_parent -- see module docs for why this is an upper bound,
                // not the strict fractional_decisions count.
                let touched_parents: HashSet<Option<chematic_cip::NodeId>> = trace
                    .decisions
                    .iter()
                    .filter(|d| d.rule.contains("Rational("))
                    .map(|d| d.ranking_parent)
                    .collect();
                proxy_ranking_parents_touched += touched_parents.len();

                // Naive final-assignment impact (narrative only -- bundles structure
                // with fraction, see module docs).
                if let (Some(&naive_code), Some(&final_code)) =
                    (naive_pass1_codes.get(&idx), final_codes.get(&idx))
                {
                    if naive_code == final_code {
                        naive_final_label_same += 1;
                    } else {
                        naive_final_label_changed += 1;
                    }
                }

                // Structure-isolated final-assignment impact -- the real D-vs-E
                // classification. `structure_pass1_codes` and `final_codes` share the
                // same Kekule-respelled digraph structure; the only thing toggled is
                // whether a `MancudeContext` is attached.
                if let (Some(&structure_code), Some(&final_code)) =
                    (structure_pass1_codes.get(&idx), final_codes.get(&idx))
                {
                    if structure_code == final_code {
                        pass1_resolved_final_label_same += 1;
                    } else {
                        pass1_resolved_final_label_changed += 1;
                        // Corroborating structural check: same comparison at the
                        // digraph-partition level, not just the resolved label --
                        // rebuild a plain (no MancudeContext) root digraph on the SAME
                        // Kekule-respelled molecule `manc_groups` used, and compare
                        // root-child partitions by index (NodeIds aren't comparable
                        // across separately-built graphs, but root-child expansion
                        // order is identical on both since the structural divergence,
                        // if any, is deeper in the tree).
                        let partitions_match = (|| {
                            let (kekule_mol, _) = kekule.as_ref()?;
                            let mut plain_graph = CipDigraph::new(kekule_mol, idx, budget).ok()?;
                            let plain_root = plain_graph.root();
                            let plain_children = plain_graph.expand_children(plain_root).ok()?;
                            let mut plain_ctx = CompareContext::new();
                            let plain_groups =
                                rank_children(&mut plain_graph, &plain_children, &mut plain_ctx)
                                    .ok()?;

                            let to_index_partition =
                                |groups: &[Vec<chematic_cip::NodeId>],
                                 children: &[chematic_cip::NodeId]|
                                 -> Vec<Vec<usize>> {
                                    groups
                                        .iter()
                                        .map(|g| {
                                            g.iter()
                                                .filter_map(|n| {
                                                    children.iter().position(|c| c == n)
                                                })
                                                .collect()
                                        })
                                        .collect()
                                };
                            let plain_partition =
                                to_index_partition(&plain_groups, &plain_children);
                            let manc_partition = to_index_partition(&manc_groups, &root_children);
                            Some(plain_partition == manc_partition)
                        })()
                        .unwrap_or(false);
                        changed_rows.push((
                            smi.to_string(),
                            idx.0,
                            structure_code,
                            final_code,
                            partitions_match,
                        ));
                    }
                }
            }
        }
    }

    println!("=== MANCUDE-Decision-A0 ===");
    println!("corpus: {csv_path}");
    println!("corpus_sha256: {actual_sha256}");
    if actual_sha256 == EXPECTED_CORPUS_SHA256 {
        println!(
            "classification A (corpus drift): RULED OUT -- corpus_sha256 matches \
             docs/rfcs/cip_accurate_rfc.md's recorded MANCUDE-Decision-A0/M3B-1b closeout \
             value exactly."
        );
    } else {
        println!(
            "*** WARNING: corpus_sha256 does NOT match the expected \
             {EXPECTED_CORPUS_SHA256} -- classification A (corpus drift) is LIVE for \
             this run. The aggregate counts below are NOT directly comparable to \
             docs/rfcs/cip_accurate_rfc.md's recorded MANCUDE-Decision-A0 closeout numbers. \
             ***"
        );
    }
    println!();
    println!("--- comparison events (strict counters, same as CIP-Perf-A0's Q3) ---");
    println!("fractional_comparisons_total: {total_fractional_comparisons}");
    println!("fractional_decisions_total:   {total_fractional_decisions}");
    println!();
    println!("--- unique stereocenters / molecules ---");
    println!("stereocenters_with_decisions: {centers_with_decisions}");
    println!(
        "molecules_with_decisions:     {}",
        molecules_with_decisions.len()
    );
    println!();
    println!("--- unique ranking_parent nodes touched (PROXY, upper bound -- see module docs) ---");
    println!(
        "proxy_ranking_parents_touched (summed per affected center): {proxy_ranking_parents_touched}"
    );
    println!();
    println!(
        "--- NAIVE final-assignment impact (narrative only -- un-Kekulized baseline vs live, \
         bundles structure+fraction, do NOT classify D/E from this) ---"
    );
    println!("naive_final_label_same:    {naive_final_label_same}");
    println!(
        "naive_final_label_changed: {naive_final_label_changed}  (these differ from the un-Kekulized \
         baseline -- see the structure-isolated numbers below for why: Kekule-respelling \
         structure, not the fraction, per docs/rfcs/cip_accurate_rfc.md's MANCUDE-Decision-A0 entry)"
    );
    println!();
    println!("--- STRUCTURE-ISOLATED final-assignment impact (real classification D vs E) ---");
    println!(
        "pass1_resolved_final_label_same:    {pass1_resolved_final_label_same}  (fraction was decision-involved but final R/S unchanged once structure is held fixed -> D)"
    );
    println!(
        "pass1_resolved_final_label_changed: {pass1_resolved_final_label_changed}  (fraction changed the resolved R/S with structure held fixed -> E)"
    );
    if let Some((c, smi, atom)) = worst {
        println!();
        println!("worst center by fractional_decisions ({c}): {smi}  atom {atom}");
    }
    if !changed_rows.is_empty() {
        println!();
        println!(
            "--- structure-isolated rows where structure_pass1 != final (should be empty; if not, re-classify) ---"
        );
        println!(
            "partitions_match=true  -> plain (no MancudeContext) and mancude root-child \
             partitions on the SAME Kekule-respelled structure agree; a label difference here \
             would need a different explanation than the partition, investigate further."
        );
        println!(
            "partitions_match=false -> the fraction's own root-level partition differs from \
             the structure-only partition -- classification E, on the SAME Kekule-respelled \
             structure (not confounded with Kekule-respelling itself)."
        );
        for (smi, atom, structure_code, fin, partitions_match) in &changed_rows {
            println!(
                "{smi}\tatom {atom}\tstructure_only(no_mancude)={}\twith_mancude(final)={}\tpartitions_match={partitions_match}",
                code_str(*structure_code),
                code_str(*fin)
            );
        }
    }
}

fn code_str(c: CipCode) -> &'static str {
    match c {
        CipCode::R => "R",
        CipCode::S => "S",
        CipCode::E => "E",
        CipCode::Z => "Z",
        CipCode::LowerR => "r",
        CipCode::LowerS => "s",
    }
}

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
