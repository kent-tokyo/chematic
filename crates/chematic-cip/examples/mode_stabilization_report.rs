//! Milestone 5B: Accurate CIP opt-in stabilization measurements. Measurement only --
//! no behavior changes. Runs `chematic_chem::assign_cip_with_mode` in both
//! `CipMode::LegacyFast` and `CipMode::Accurate` across a full SMILES corpus and
//! reports: (1) per-molecule wall-clock time for each mode, (2) the Accurate-mode
//! unresolved rate broken down by cause (`Tied` vs `BudgetExceeded`), (3) how many
//! molecules have both an R/S and an E/Z assignment (the merge path Milestone 5A
//! added), and a byte-for-byte parity check between `CipMode::LegacyFast` and the
//! pre-existing `chematic_chem::assign_cip()` (must be identical by construction --
//! this is the regression backstop, not a new measurement).
//!
//! `SMILES.csv` is not checked into this repo (same convention as
//! `scripts/bench5k.py`/`corpus_snapshot.rs`) -- supply your own corpus file.
//!
//! Usage: cargo run -p chematic-cip --release --example mode_stabilization_report -- <SMILES.csv>

use std::env;
use std::fs;
use std::time::Instant;

use chematic_chem::{CipMode, CipUnresolvedReason, assign_cip, assign_cip_with_mode};

fn main() {
    let args: Vec<String> = env::args().collect();
    let csv_path = args
        .get(1)
        .unwrap_or_else(|| panic!("usage: mode_stabilization_report <SMILES.csv>"));

    let content = fs::read_to_string(csv_path).expect("read SMILES.csv");
    let smis: Vec<&str> = content
        .lines()
        .skip(1)
        .filter_map(|line| line.split(',').next())
        .filter(|s| !s.is_empty())
        .collect();

    let mut legacy_fast_mismatch = 0usize;
    let mut legacy_total_ns = 0u128;
    let mut accurate_total_ns = 0u128;
    let mut molecules_with_stereo = 0usize;
    let mut molecules_with_ez_and_rs = 0usize;
    let mut atoms_resolved_accurate = 0usize;
    let mut atoms_tied = 0usize;
    let mut atoms_budget_exceeded = 0usize;
    let mut atoms_oracle_unstable = 0usize;
    let mut engine_errors = 0usize;
    let mut parse_errors = 0usize;
    let mut tied_examples: Vec<(String, u32, String)> = Vec::new();

    for smi in &smis {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(_) => {
                parse_errors += 1;
                continue;
            }
        };

        let t0 = Instant::now();
        let legacy = assign_cip(&mol);
        legacy_total_ns += t0.elapsed().as_nanos();

        let t1 = Instant::now();
        let via_mode_legacy = assign_cip_with_mode(&mol, CipMode::LegacyFast);
        let mode_legacy_ns = t1.elapsed().as_nanos();
        let _ = mode_legacy_ns; // included in legacy_total_ns via the direct call above

        match &via_mode_legacy {
            Ok(r) if r.assignments == legacy.assignments && r.unresolved.is_empty() => {}
            _ => legacy_fast_mismatch += 1,
        }

        let t2 = Instant::now();
        let accurate = assign_cip_with_mode(&mol, CipMode::Accurate);
        accurate_total_ns += t2.elapsed().as_nanos();

        if !legacy.assignments.is_empty() {
            molecules_with_stereo += 1;
        }

        match accurate {
            Ok(result) => {
                let has_rs = result.assignments.iter().any(|(_, c)| {
                    matches!(c, chematic_core::CipCode::R | chematic_core::CipCode::S)
                });
                let has_ez = result.assignments.iter().any(|(_, c)| {
                    matches!(c, chematic_core::CipCode::E | chematic_core::CipCode::Z)
                });
                if has_rs && has_ez {
                    molecules_with_ez_and_rs += 1;
                }
                atoms_resolved_accurate += result.assignments.len();
                for (idx, reason) in &result.unresolved {
                    match reason {
                        CipUnresolvedReason::Tied => {
                            atoms_tied += 1;
                            let elem = mol.atom(*idx).element.symbol();
                            tied_examples.push((smi.to_string(), idx.0, elem.to_string()));
                        }
                        CipUnresolvedReason::BudgetExceeded => atoms_budget_exceeded += 1,
                        CipUnresolvedReason::OracleUnstable => atoms_oracle_unstable += 1,
                    }
                }
            }
            Err(_) => engine_errors += 1,
        }
    }

    let unresolved_total = atoms_tied + atoms_budget_exceeded + atoms_oracle_unstable;
    let unresolved_rate = if atoms_resolved_accurate + unresolved_total > 0 {
        100.0 * unresolved_total as f64 / (atoms_resolved_accurate + unresolved_total) as f64
    } else {
        0.0
    };

    println!("=== Milestone 5B stabilization report ===");
    println!("input molecules:                    {}", smis.len());
    println!("parse errors (excluded):             {parse_errors}");
    println!("molecules with any legacy stereo:    {molecules_with_stereo}");
    println!(
        "molecules with R/S AND E/Z (Accurate, the Milestone 5A merge path): {molecules_with_ez_and_rs}"
    );
    println!();
    println!("--- CipMode::LegacyFast byte-for-byte parity with assign_cip() ---");
    println!("mismatches (should be 0):            {legacy_fast_mismatch}");
    println!();
    println!("--- Accurate-mode unresolved rate ---");
    println!("atoms resolved:                      {atoms_resolved_accurate}");
    println!("atoms unresolved (Tied):             {atoms_tied}");
    println!("atoms unresolved (BudgetExceeded):   {atoms_budget_exceeded}");
    println!("atoms unresolved (OracleUnstable):    {atoms_oracle_unstable}");
    println!("unresolved rate:                     {unresolved_rate:.3}%");
    println!("engine errors (should be 0):          {engine_errors}");
    println!();
    println!("--- Tied atoms, breakdown by element ---");
    let mut by_elem: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for (_, _, elem) in &tied_examples {
        *by_elem.entry(elem.clone()).or_insert(0) += 1;
    }
    let mut elems: Vec<_> = by_elem.into_iter().collect();
    elems.sort_by_key(|(_, n)| std::cmp::Reverse(*n));
    for (elem, n) in &elems {
        println!("  {elem}: {n}");
    }
    println!("all tied (smiles, atom_idx, element):");
    for (smi, idx, elem) in &tied_examples {
        println!("  {smi}  atom={idx}  element={elem}");
    }
    println!();
    println!("--- Perf: total wall-clock across corpus ---");
    println!(
        "assign_cip() (legacy, direct call):  {:.2} ms  ({:.2} us/molecule)",
        legacy_total_ns as f64 / 1e6,
        legacy_total_ns as f64 / 1e3 / smis.len().max(1) as f64
    );
    println!(
        "assign_cip_with_mode(Accurate):      {:.2} ms  ({:.2} us/molecule)",
        accurate_total_ns as f64 / 1e6,
        accurate_total_ns as f64 / 1e3 / smis.len().max(1) as f64
    );
    println!(
        "Accurate / Legacy ratio:             {:.2}x",
        accurate_total_ns as f64 / legacy_total_ns.max(1) as f64
    );
}
