//! Published Rust crate reaction output, retaining typed transform failures.

use std::fs;

use chematic::{rxn, smiles};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 {
        return Err("usage: reactions BASE_CASES_JSON STRATA_JSON ROWS_JSON SUMMARY_JSON".into());
    }
    let base_bytes = fs::read(&args[1])?;
    let strata_bytes = fs::read(&args[2])?;
    let mut cases: Vec<Value> = serde_json::from_slice(&base_bytes)?;
    if cases.len() != 57 {
        return Err("historical base fixture must have 57 cases".into());
    }
    let strata: Value = serde_json::from_slice(&strata_bytes)?;
    if strata["base_fixture_sha256"] != digest(&base_bytes) {
        return Err("strata manifest does not pin base fixture hash".into());
    }
    cases.extend(
        strata["cases"]
            .as_array()
            .ok_or("missing strata cases")?
            .iter()
            .cloned(),
    );
    let mut rows = Vec::with_capacity(cases.len());
    for case in &cases {
        let mut reactants = Vec::new();
        let mut parse_error = None;
        for text in case["reactants"].as_array().ok_or("missing reactants")? {
            match smiles::parse(text.as_str().ok_or("non-string reactant")?) {
                Ok(mol) => reactants.push(mol),
                Err(err) => {
                    parse_error = Some(err.to_string());
                    break;
                }
            }
        }
        if let Some(detail) = parse_error {
            rows.push(json!({"id": case["id"], "status": "typed_refusal", "stage": "reactant_parse", "reason": "smiles_parse", "detail": detail}));
            continue;
        }
        let refs: Vec<_> = reactants.iter().collect();
        let smirks = case["smirks"].as_str().ok_or("non-string SMIRKS")?;
        match rxn::run_reactants_with_diagnostics(
            smirks,
            &refs,
            &rxn::ReactionTransformLimits::default(),
        ) {
            Ok(report) => {
                let raw_count = report.products.len();
                let mut sets: Vec<Vec<String>> = report
                    .products
                    .into_iter()
                    .map(|set| {
                        let mut values: Vec<String> =
                            set.iter().map(smiles::canonical_smiles).collect();
                        values.sort();
                        values
                    })
                    .collect();
                sets.sort();
                sets.dedup();
                let diagnostics = report.diagnostics;
                let status = if sets.is_empty() && diagnostics.valence_rejected_matches > 0 {
                    "diagnosed_valence_refusal"
                } else {
                    "products"
                };
                rows.push(json!({"id": case["id"], "status": status, "sets": sets, "raw_product_sets": raw_count,
                    "diagnostics": {"accepted_matches": diagnostics.accepted_matches,
                        "applied_products": diagnostics.applied_products,
                        "valence_rejected_matches": diagnostics.valence_rejected_matches,
                        "truncated_matches": diagnostics.truncated_matches}}));
            }
            Err(err) => {
                let reason = match &err {
                    rxn::TransformError::ResourceLimit { .. } => "resource_limit",
                    rxn::TransformError::ReactantCountMismatch { .. } => "reactant_count_mismatch",
                    rxn::TransformError::SmirksParse(_) => "smirks_parse",
                };
                rows.push(json!({"id": case["id"], "status": "typed_refusal", "stage": "transform", "reason": reason, "detail": err.to_string()}));
            }
        }
    }
    let raw = serde_json::to_vec_pretty(&rows)?;
    fs::write(&args[3], &raw)?;
    let summary = json!({"schema": "published-rust-reactions/v1", "crate": "chematic 1.0.30 from crates.io",
        "base_cases_sha256": digest(&base_bytes), "strata_sha256": digest(&strata_bytes), "rows_sha256": digest(&raw),
        "input_count": cases.len(), "typed_refusals": rows.iter().filter(|r| r["status"] == "typed_refusal").count(),
        "diagnosed_valence_refusals": rows.iter().filter(|r| r["status"] == "diagnosed_valence_refusal").count(),
        "product_rows": rows.iter().filter(|r| r["status"] == "products").count()});
    fs::write(&args[4], serde_json::to_vec_pretty(&summary)?)?;
    println!("{summary}");
    Ok(())
}
