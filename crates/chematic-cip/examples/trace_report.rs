//! Milestone 4A-0: dump the full Rules-1a/1b/2 decision path (a [`ComparisonTrace`]) for
//! one stereocenter's root-children ranking. Existing tooling (`residual_report`,
//! `scripts/cip_accurate_full_corpus_report.py`) says *whether* a row is right, wrong, or
//! tied; this says *why*, in the same shape `d0e726b`'s root-cause work needed by hand --
//! see `trace.rs`'s module docs.
//!
//! Usage:
//!   cargo run -p chematic-cip --release --example trace_report -- '<smiles>' <atom_idx>
//!
//! Builds the exact same digraph `assign_cip_accurate_experimental` uses (Kekule form +
//! MANCUDE fractional numbers when available), then ranks the stereocenter's root
//! children with a trace-enabled [`CompareContext`], printing every decision step.

use std::env;

use chematic_cip::{CipBudget, CipDigraph, ComparisonTrace};
use chematic_core::AtomIdx;

fn main() {
    let args: Vec<String> = env::args().collect();
    let (Some(smi), Some(atom_idx)) = (args.get(1), args.get(2)) else {
        eprintln!("usage: trace_report '<smiles>' <atom_idx>");
        std::process::exit(64);
    };
    let atom_idx: u32 = atom_idx.parse().expect("atom_idx must be a number");

    let mol = chematic_smiles::parse(smi).expect("valid SMILES");
    let idx = AtomIdx(atom_idx);
    let budget = CipBudget::default_budget();

    let kekule = chematic_cip::prepare_kekule_form(&mol).ok();
    let mut graph = match &kekule {
        Some((kekule_mol, ctx)) => {
            CipDigraph::new_with_mancude(kekule_mol, idx, budget, ctx).expect("digraph builds")
        }
        None => CipDigraph::new(&mol, idx, budget).expect("digraph builds"),
    };
    let root = graph.root();
    let root_children = graph.expand_children(root).expect("root expands");

    let mut trace = ComparisonTrace::new(root, root);
    let mut cmp_ctx = chematic_cip::CompareContext::with_trace(&mut trace);
    let groups = chematic_cip::rank_children(&mut graph, &root_children, &mut cmp_ctx)
        .expect("ranking succeeds");

    println!("{smi}  atom {atom_idx}");
    println!(
        "groups (highest priority first): {:?}",
        groups
            .iter()
            .map(|g| g.iter().map(|n| n.0).collect::<Vec<_>>())
            .collect::<Vec<_>>()
    );
    println!();
    print!("{trace}");
}
