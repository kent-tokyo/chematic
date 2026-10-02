//! Published-crate chemistry lane: every input and SMARTS cell is accounted for.
//! Run from a Cargo.lock resolving chematic =1.0.30 from crates.io, not a path.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use chematic::{chem, core, fp, smarts, smiles};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn read_jsonl(path: &Path) -> Result<Vec<Value>, Box<dyn std::error::Error>> {
    BufReader::new(File::open(path)?)
        .lines()
        .map(|line| Ok(serde_json::from_str(&line?)?))
        .collect()
}

fn code_string(code: core::CipCode) -> &'static str {
    match code {
        core::CipCode::R => "R",
        core::CipCode::S => "S",
        core::CipCode::E => "E",
        core::CipCode::Z => "Z",
        core::CipCode::LowerR => "r",
        core::CipCode::LowerS => "s",
    }
}

fn cip_maps(
    mol: &core::Molecule,
) -> Result<
    (
        BTreeMap<String, String>,
        BTreeMap<String, String>,
        BTreeMap<String, String>,
    ),
    Box<dyn std::error::Error>,
> {
    let accurate = chem::assign_cip_with_mode(mol, chem::CipMode::Accurate)?;
    let mut ez = chem::assign_ez_bonds_with_mode(mol, chem::CipMode::Accurate);
    let mut atoms = BTreeMap::new();
    let mut bonds = BTreeMap::new();
    for (index, code) in accurate.assignments {
        match code {
            core::CipCode::E | core::CipCode::Z => {
                let key = if let Some(pos) = ez
                    .iter()
                    .position(|(bond_idx, c)| *c == code && mol.bond(*bond_idx).atom1 == index)
                {
                    let (bond_idx, _) = ez.remove(pos);
                    let bond = mol.bond(bond_idx);
                    let (a, b) = (
                        bond.atom1.0.min(bond.atom2.0),
                        bond.atom1.0.max(bond.atom2.0),
                    );
                    format!("{a}-{b}")
                } else {
                    format!("ambiguous-ez-atom-{}", index.0)
                };
                bonds.insert(key, code_string(code).to_string());
            }
            _ => {
                atoms.insert(index.0.to_string(), code_string(code).to_string());
            }
        }
    }
    let unresolved = accurate
        .unresolved
        .into_iter()
        .map(|(index, reason)| {
            let label = match reason {
                chem::CipUnresolvedReason::Tied => "tied",
                chem::CipUnresolvedReason::BudgetExceeded => "budget_exceeded",
                chem::CipUnresolvedReason::OracleUnstable => "oracle_unstable",
                chem::CipUnresolvedReason::LonePairCenter => "lone_pair_center",
            };
            (index.0.to_string(), label.to_string())
        })
        .collect();
    Ok((atoms, bonds, unresolved))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 7 {
        return Err("usage: published-chematic-gate CORPUS QUERIES PYTHON_ROWS SMARTS_ORACLE ROWS_OUT SUMMARY_OUT".into());
    }
    let corpus_path = Path::new(&args[1]);
    let queries_path = Path::new(&args[2]);
    let reference_path = Path::new(&args[3]);
    let oracle_path = Path::new(&args[4]);
    let rows_path = Path::new(&args[5]);
    let summary_path = Path::new(&args[6]);
    let corpus_bytes = fs::read(corpus_path)?;
    let queries_bytes = fs::read(queries_path)?;
    let corpus: Vec<String> = String::from_utf8(corpus_bytes.clone())?
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.trim().starts_with('#'))
        .map(|line| line.split_whitespace().next().unwrap_or("").to_string())
        .collect();
    let queries_doc: Value = serde_json::from_slice(&queries_bytes)?;
    let queries: Vec<String> = queries_doc["queries"]
        .as_array()
        .ok_or("missing queries")?
        .iter()
        .map(|item| item.as_str().unwrap_or("").to_string())
        .collect();
    let prepared: Vec<_> = queries
        .iter()
        .map(|query| smarts::parse_smarts(query))
        .collect::<Result<_, _>>()?;
    let reference = read_jsonl(reference_path)?;
    let mut oracle = read_jsonl(oracle_path)?;
    let manifest = oracle.remove(0);
    if manifest["rdkit_version"] != "2026.03.6"
        || manifest["corpus_sha256"] != digest(&corpus_bytes)
        || manifest["queries_sha256"] != digest(&queries_bytes)
        || corpus.len() != reference.len()
        || corpus.len() != oracle.len()
    {
        return Err("oracle/corpus/query correspondence failed".into());
    }
    let mut counts = BTreeMap::<&str, usize>::new();
    let mut add = |key: &'static str, n: usize| {
        *counts.entry(key).or_default() += n;
    };
    let mut output = BufWriter::new(File::create(rows_path)?);
    for (index, smile) in corpus.iter().enumerate() {
        let baseline = &reference[index];
        let truth = &oracle[index];
        if baseline["input_index"] != index
            || truth["input_index"] != index
            || baseline["smiles"] != *smile
            || truth["smiles"] != *smile
        {
            return Err(format!("input correspondence failed at {index}").into());
        }
        let mut row = json!({"input_index": index, "smiles": smile, "status": "completed", "differences": []});
        let mol = match smiles::parse(smile) {
            Ok(value) => value,
            Err(err) => {
                add("parse_failure", 1);
                row["status"] = json!("parse_failure");
                row["error"] = json!(err.to_string());
                writeln!(output, "{row}")?;
                continue;
            }
        };
        row["canonical"] = json!(smiles::canonical_smiles(&mol));
        if row["canonical"] != baseline["smiles_parse_write"]["chematic_canonical"] {
            add("canonical_cross_binding_differences", 1);
        }

        let (atoms, bonds, unresolved) = cip_maps(&mol)?;
        let expected_atoms = &baseline["cip"]["rdkit_atoms"];
        let expected_bonds = &baseline["cip"]["rdkit_bonds"];
        let wrong_label = atoms.iter().any(|(idx, code)| expected_atoms[idx] != *code);
        let missing: Vec<_> = expected_atoms
            .as_object()
            .ok_or("missing RDKit CIP atoms")?
            .keys()
            .filter(|idx| {
                atoms.get(*idx).map(|label| label.as_str()) != expected_atoms[*idx].as_str()
            })
            .collect();
        row["cip"] = json!({"atoms": atoms, "bonds": bonds, "unresolved": unresolved});
        if wrong_label
            || row["cip"]["bonds"] != *expected_bonds
            || missing.iter().any(|idx| !unresolved.contains_key(*idx))
        {
            add("cip_wrong_confident", 1);
        } else if missing.is_empty() {
            add("cip_exact", 1);
        } else {
            add("cip_typed_refusal", 1);
        }

        match fp::rdkit_morgan_ecfp4_bitvec(&mol) {
            Ok(bits) => {
                let fp_sha = digest(&bits.to_le_bytes());
                row["morgan_sha256"] = json!(fp_sha);
                if row["morgan_sha256"] == baseline["morgan"]["rdkit_sha256"] {
                    add("morgan_exact", 1);
                } else {
                    add("morgan_wrong_confident", 1);
                }
            }
            Err(err) => {
                row["morgan_error"] = json!(err.to_string());
                add("morgan_typed_refusal", 1);
            }
        }

        let mut matches = Vec::with_capacity(prepared.len());
        for (query_index, query) in prepared.iter().enumerate() {
            let mut actual = smarts::find_match_atom_sets_perceived(
                query,
                &mol,
                &smarts::MatchConfig::default(),
            );
            actual.sort();
            actual.dedup();
            if json!(actual) != truth["matches"][query_index] {
                add("smarts_differences", 1);
                row["differences"]
                    .as_array_mut()
                    .unwrap()
                    .push(json!({"operation": "smarts", "query_index": query_index}));
            }
            add("smarts_cells", 1);
            matches.push(actual);
        }
        row["smarts_matches"] = json!(matches);
        add("completed", 1);
        writeln!(output, "{row}")?;
    }
    output.flush()?;
    let summary = json!({
        "schema": "published-rust-chemistry-v1",
        "crate": {"name": "chematic", "version": "1.0.30", "origin": "crates.io exact dependency; verify Cargo.lock"},
        "oracle": {"rdkit_version": "2026.03.6", "corpus_sha256": digest(&corpus_bytes),
            "queries_sha256": digest(&queries_bytes), "smarts_all_cells_sha256": digest(&fs::read(oracle_path)?),
            "python_reference_sha256": digest(&fs::read(reference_path)?)},
        "counts": {"input": corpus.len(), "completed": counts.get("completed").copied().unwrap_or(0),
            "parse_failure": counts.get("parse_failure").copied().unwrap_or(0),
            "canonical_cross_binding_differences": counts.get("canonical_cross_binding_differences").copied().unwrap_or(0),
            "cip_exact": counts.get("cip_exact").copied().unwrap_or(0),
            "cip_typed_refusal": counts.get("cip_typed_refusal").copied().unwrap_or(0),
            "cip_wrong_confident": counts.get("cip_wrong_confident").copied().unwrap_or(0),
            "morgan_exact": counts.get("morgan_exact").copied().unwrap_or(0),
            "morgan_typed_refusal": counts.get("morgan_typed_refusal").copied().unwrap_or(0),
            "morgan_wrong_confident": counts.get("morgan_wrong_confident").copied().unwrap_or(0),
            "smarts_cells": counts.get("smarts_cells").copied().unwrap_or(0),
            "smarts_differences": counts.get("smarts_differences").copied().unwrap_or(0)},
        "rows": {"path": rows_path, "sha256": digest(&fs::read(rows_path)?)},
        "reaction_transform": {"status": "not_yet_measured"},
        "operation_matrix": {"status": "not_yet_measured"},
    });
    fs::write(summary_path, serde_json::to_string_pretty(&summary)? + "\n")?;
    println!("{}", summary["counts"]);
    Ok(())
}
