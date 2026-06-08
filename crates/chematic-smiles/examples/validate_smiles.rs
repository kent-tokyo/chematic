//! Large-scale SMILES validation: reads one SMILES per line from stdin (or a file),
//! parses each, writes canonical SMILES, re-parses, and confirms atom/bond counts match.
//!
//! Usage:
//!   cargo run -p chematic-smiles --example validate_smiles --release < chembl.smi
//!   cargo run -p chematic-smiles --example validate_smiles --release -- chembl.smi
//!
//! Output (stderr): progress every 10000 molecules
//! Output (stdout): final summary + any failures printed as  FAIL\t<smiles>\t<reason>

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use chematic_smiles::{canonical_smiles, parse};

fn main() {
    let reader: Box<dyn BufRead> = match env::args().nth(1) {
        Some(path) => Box::new(BufReader::new(
            File::open(&path).unwrap_or_else(|e| panic!("cannot open {path}: {e}")),
        )),
        None => Box::new(BufReader::new(io::stdin())),
    };

    let mut total: u64 = 0;
    let mut ok: u64 = 0;
    let mut fail_parse: u64 = 0;
    let mut fail_roundtrip: u64 = 0;
    let mut error_counts: HashMap<String, u64> = HashMap::new();
    let mut first_failures: Vec<(String, String)> = Vec::new(); // (smiles, reason)

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        let smiles = line.trim();
        if smiles.is_empty() || smiles.starts_with('#') {
            continue;
        }
        // Some SMILES files have <smiles>\t<name> format — take the first field
        let smiles = smiles.split('\t').next().unwrap_or(smiles);

        total += 1;
        #[allow(clippy::manual_is_multiple_of)]
        if total % 10_000 == 0 {
            eprintln!(
                "  [{total}] ok={ok} fail_parse={fail_parse} fail_roundtrip={fail_roundtrip}"
            );
        }

        // Step 1: parse
        let mol = match parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                fail_parse += 1;
                let key = error_key(&e.to_string());
                *error_counts.entry(key).or_default() += 1;
                if first_failures.len() < 50 {
                    first_failures.push((smiles.to_string(), format!("parse: {e}")));
                }
                continue;
            }
        };

        let atoms1 = mol.atom_count();
        let bonds1 = mol.bond_count();

        // Step 2: write canonical SMILES
        let canon = canonical_smiles(&mol);

        // Step 3: re-parse
        let mol2 = match parse(&canon) {
            Ok(m) => m,
            Err(e) => {
                fail_roundtrip += 1;
                let key = format!("roundtrip_parse: {}", error_key(&e.to_string()));
                *error_counts.entry(key).or_default() += 1;
                if first_failures.len() < 50 {
                    first_failures.push((
                        smiles.to_string(),
                        format!("roundtrip parse of '{canon}': {e}"),
                    ));
                }
                continue;
            }
        };

        // Step 4: atom/bond count consistency
        let atoms2 = mol2.atom_count();
        let bonds2 = mol2.bond_count();

        if atoms1 != atoms2 || bonds1 != bonds2 {
            fail_roundtrip += 1;
            let reason = format!(
                "count mismatch: atoms {atoms1}→{atoms2}, bonds {bonds1}→{bonds2}  canon='{canon}'"
            );
            *error_counts
                .entry("roundtrip_count_mismatch".to_string())
                .or_default() += 1;
            if first_failures.len() < 50 {
                first_failures.push((smiles.to_string(), reason));
            }
            continue;
        }

        ok += 1;
    }

    // ── Summary ─────────────────────────────────────────────────────────────
    let fail_total = fail_parse + fail_roundtrip;
    let success_rate = if total > 0 {
        ok as f64 / total as f64 * 100.0
    } else {
        0.0
    };

    println!("=== chematic SMILES validation ===");
    println!("Total       : {total}");
    println!("OK          : {ok}  ({success_rate:.3}%)");
    println!("Fail parse  : {fail_parse}");
    println!("Fail rt     : {fail_roundtrip}");
    println!("Fail total  : {fail_total}");
    println!();
    if !error_counts.is_empty() {
        println!("--- Error breakdown ---");
        let mut sorted: Vec<_> = error_counts.iter().collect();
        sorted.sort_by_key(|&(_, v)| std::cmp::Reverse(v));
        for (k, v) in sorted.iter().take(20) {
            println!("  {v:6}  {k}");
        }
        println!();
    }
    if !first_failures.is_empty() {
        println!("--- First failures (up to 50) ---");
        for (smi, reason) in &first_failures {
            println!("  FAIL  {smi}");
            println!("        {reason}");
        }
    }
}

/// Collapse an error message to a short category key.
fn error_key(msg: &str) -> String {
    if msg.contains("unknown element") {
        // Extract just the symbol
        if let Some(start) = msg.find('\'')
            && let Some(end) = msg[start + 1..].find('\'')
        {
            return format!("unknown_element_{}", &msg[start + 1..start + 1 + end]);
        }
        "unknown_element".to_string()
    } else if msg.contains("unmatched ring") {
        "unmatched_ring_closure".to_string()
    } else if msg.contains("mismatched paren") {
        "mismatched_parentheses".to_string()
    } else if msg.contains("invalid bracket") {
        "invalid_bracket_atom".to_string()
    } else if msg.contains("conflicting bond") {
        "conflicting_ring_bond".to_string()
    } else if msg.contains("unexpected end") {
        "unexpected_end".to_string()
    } else {
        msg.chars().take(40).collect()
    }
}
