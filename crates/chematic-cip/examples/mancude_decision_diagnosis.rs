//! MANCUDE-Decision-A0: classify the nonzero `fractional_decisions` count found while
//! running CIP-Perf-A0 (`cip_perf_diagnosis`) against the frozen `SMILES.csv` corpus
//! (`corpus_sha256=1c47371d...`, confirmed identical to the corpus
//! `docs/cip_accurate_rfc.md`'s Milestone 3B-1b closeout used for its "byte-identical,
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
//! - final-assignment impact: for stereocenters where Pass 1 (Rules 1a/1b/2, both with
//!   and without `MancudeContext`) alone resolves the center, does the resolved
//!   `CipCode` actually differ between `assign_cip_accurate_experimental_without_mancude`
//!   and `assign_cip_accurate_experimental`? This isolates the fraction's effect from
//!   Rule 4b/5 ever running at all (which `_without_mancude` never reaches), directly
//!   answering classification D vs E for the Pass-1-resolved subset.
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

fn main() {
    let args: Vec<String> = env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| format!("{}/Downloads/SMILES.csv", env::var("HOME").unwrap()));

    let content = fs::read_to_string(&csv_path).expect("read SMILES.csv");
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

        let pass1 = assign_cip_accurate_experimental_without_mancude(&mol, budget);
        let final_result = assign_cip_accurate_experimental(&mol, budget);
        let (Ok(pass1), Ok(final_result)) = (pass1, final_result) else {
            continue;
        };
        let pass1_codes: std::collections::HashMap<AtomIdx, CipCode> =
            pass1.assignments.iter().copied().collect();
        let final_codes: std::collections::HashMap<AtomIdx, CipCode> =
            final_result.assignments.iter().copied().collect();

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

                // Final-assignment impact, isolated to the Pass-1-resolved subset.
                if let (Some(&p1_code), Some(&final_code)) =
                    (pass1_codes.get(&idx), final_codes.get(&idx))
                {
                    if p1_code == final_code {
                        pass1_resolved_final_label_same += 1;
                    } else {
                        pass1_resolved_final_label_changed += 1;
                        // Discriminating test (isolates Pass 1's own mancude effect from
                        // Rule 4b/5, which `pass1_codes` vs `final_codes` alone conflates
                        // -- `final_codes` includes atoms Rule 4b/5 resolved, not just
                        // Pass-1-with-mancude): rebuild a *plain* root digraph (no
                        // MancudeContext, on the original `mol` -- exactly what
                        // `_without_mancude` uses internally) and compare its root-child
                        // group partition against `manc_groups` above, by root-child
                        // *index* (NodeIds aren't comparable across separate graphs, but
                        // root-child expansion order is identical on both since the
                        // structural divergence is deeper in the tree).
                        let partitions_match = (|| {
                            let mut plain_graph = CipDigraph::new(&mol, idx, budget).ok()?;
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
                            p1_code,
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
    println!(
        "classification A (corpus drift): RULED OUT -- `shasum -a 256` on this file was \
         independently verified (outside this tool) to match docs/cip_accurate_rfc.md's \
         recorded M3B-1b closeout corpus_sha256=1c47371d... exactly."
    );
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
        "--- final-assignment impact, Pass-1-resolved subset only (classification D vs E) ---"
    );
    println!(
        "pass1_resolved_final_label_same:    {pass1_resolved_final_label_same}  (fraction was decision-involved but final R/S unchanged -> D)"
    );
    println!(
        "pass1_resolved_final_label_changed: {pass1_resolved_final_label_changed}  (fraction actually changed the resolved R/S -> E, needs RDKit-agreement check)"
    );
    if let Some((c, smi, atom)) = worst {
        println!();
        println!("worst center by fractional_decisions ({c}): {smi}  atom {atom}");
    }
    if !changed_rows.is_empty() {
        println!();
        println!(
            "--- rows where without_mancude != final, with the Pass-1-only discriminating test ---"
        );
        println!(
            "partitions_match=true  -> plain-Pass-1 and mancude-Pass-1 rank identically; the \
             flip is Rule 4b/5's doing, not mancude's -- classification D, not E."
        );
        println!(
            "partitions_match=false -> mancude-Pass-1's OWN ranking differs from plain-Pass-1's \
             -- classification E stands for this row."
        );
        for (smi, atom, p1, fin, partitions_match) in &changed_rows {
            println!(
                "{smi}\tatom {atom}\twithout_mancude={}\twith_mancude(final)={}\tpartitions_match={partitions_match}",
                code_str(*p1),
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
