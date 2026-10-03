use std::collections::BTreeMap;
use std::fs;

use chematic_chem::add_hydrogens;
use chematic_ff::assign_mmff94_numeric_types;
use chematic_smiles::parse;

// Explicit integration gate: generate the TSV using scripts/mmff94_atom_type_oracle_tsv.py
// with RDKit 2026.03.6, then set SCHEMATIC_MMFF94_ORACLE_TSV before running.
#[test]
#[ignore = "requires pinned RDKit 2026.03.6 oracle TSV"]
fn compare_full_corpus_to_pinned_rdkit_oracle() {
    let corpus = fs::read_to_string("validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")
        .or_else(|_| {
            fs::read_to_string("../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")
        })
        .unwrap();
    let oracle = fs::read_to_string(std::env::var("SCHEMATIC_MMFF94_ORACLE_TSV").unwrap()).unwrap();
    let mut buckets = BTreeMap::<(u8, u8), usize>::new();
    let mut status = BTreeMap::<&str, usize>::new();
    let mut atom_rows = String::new();
    let mut oracle_rows = 0usize;
    for ((i, line), expected) in corpus.lines().enumerate().zip(oracle.lines()) {
        oracle_rows += 1;
        let (oracle_i, value) = expected.split_once('\t').unwrap();
        assert_eq!(oracle_i.parse::<usize>().unwrap(), i);
        if matches!(value, "parse" | "unsupported") {
            continue;
        }
        let reference: Vec<u8> = value.split(',').map(|x| x.parse().unwrap()).collect();
        let smiles = line.split_whitespace().next().unwrap();
        let Ok(mol) = parse(smiles) else {
            *status.entry("parse_error").or_default() += 1;
            continue;
        };
        let with_h = add_hydrogens(&mol);
        let Ok(actual) = assign_mmff94_numeric_types(&with_h) else {
            *status.entry("typing_error").or_default() += 1;
            continue;
        };
        if actual.len() != reference.len() {
            *status.entry("length_mismatch").or_default() += 1;
            continue;
        }
        *status.entry("compared").or_default() += 1;
        for (atom_i, (a, b)) in actual.into_iter().zip(reference).enumerate() {
            atom_rows.push_str(&format!("{i}\t{atom_i}\t{b}\t{a}\n"));
            if a != b {
                *buckets.entry((b, a)).or_default() += 1;
            }
        }
    }
    println!("STATUS {status:?}");
    println!("DIFFS {buckets:?}");
    assert_eq!(oracle_rows, 10_000, "the benchmark corpus must stay pinned");
    assert_eq!(corpus.lines().count(), oracle_rows);
    assert_eq!(oracle.lines().count(), oracle_rows);
    assert_eq!(status.get("compared"), Some(&9_774));
    assert_eq!(status.get("typing_error"), Some(&22));
    assert!(
        buckets.values().sum::<usize>() <= 459,
        "MMFF94 type parity regressed"
    );
    if let Ok(path) = std::env::var("SCHEMATIC_MMFF94_ATOM_ROWS_OUT") {
        fs::write(path, atom_rows).unwrap();
    }
}
