//! Dump RDKit-compatible top-k search results for JSONL search cases.

use chematic_fp::FpType;
use serde_json::json;
use std::io::{self, BufRead};

fn main() {
    let mut cached: Option<(Vec<String>, chematic_fp::PreparedFingerprintIndex)> = None;
    for line in io::stdin().lock().lines() {
        let line = line.expect("read search case");
        if line.trim().is_empty() {
            continue;
        }
        let case: serde_json::Value = serde_json::from_str(&line).expect("valid search case JSON");
        let query_smiles = case["query"].as_str().expect("query string");
        let db_smiles = case["db"].as_array().expect("database array");
        let db_strings = db_smiles
            .iter()
            .map(|smiles| smiles.as_str().expect("database SMILES string").to_string())
            .collect::<Vec<_>>();
        let k = case["k"].as_u64().expect("k integer") as usize;
        let threshold = case.get("threshold").and_then(|value| value.as_f64());
        let query = chematic_smiles::parse(query_smiles).expect("valid query SMILES");
        if cached
            .as_ref()
            .is_none_or(|(previous, _)| previous != &db_strings)
        {
            let db = db_strings
                .iter()
                .map(|smiles| chematic_smiles::parse(smiles).expect("valid database SMILES"))
                .collect::<Vec<_>>();
            let index = chematic_fp::PreparedFingerprintIndex::try_new(&db, FpType::RdkitEcfp4)
                .expect("RDKit-compatible search index must build");
            cached = Some((db_strings, index));
        }
        let results = if let Some(threshold) = threshold {
            cached
                .as_ref()
                .unwrap()
                .1
                .try_search_threshold(&query, threshold, k)
        } else {
            cached.as_ref().unwrap().1.try_search(&query, k)
        }
        .expect("RDKit-compatible search must not fail for the fixture");
        println!(
            "{}",
            json!({"results": results.iter().map(|(i, score)| json!({"index": i, "tanimoto": score})).collect::<Vec<_>>() })
        );
    }
}
