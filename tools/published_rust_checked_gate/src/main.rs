//! Checked RDKit-profile reaction rows from the published `chematic` crate.
//!
//! Mirrors Python `run_smirks_checked(..., rdkit_compat=True)` on the 83
//! pinned fixtures and writes the candidates schema read by
//! `scripts/reaction_python_checked_provenance_gate.py --candidates`.
//!
//! Usage: published-chematic-checked-gate BASE_JSON STRATA_JSON OUT_JSON

use std::fs;

use chematic::{rxn, smiles};
use serde_json::{Map, Value, json};
use sha2::{Digest, Sha256};

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn diagnostics(accepted: usize, applied: usize, rejected: usize, truncated: bool) -> Value {
    json!({"accepted_matches": accepted, "applied_products": applied,
           "valence_rejected_matches": rejected, "truncated_matches": truncated})
}

fn refusal(status: &str, reason: &str, detail: String) -> Value {
    json!({"status": status, "reason": reason, "detail": detail,
           "diagnostics": diagnostics(0, 0, 0, false)})
}

fn row(case: &Value) -> Result<Value, Box<dyn std::error::Error>> {
    let mut reactants = Vec::new();
    for text in case["reactants"].as_array().ok_or("missing reactants")? {
        match smiles::parse(text.as_str().ok_or("non-string reactant")?) {
            Ok(mol) => reactants.push(mol),
            Err(err) => return Ok(refusal("typed_refusal", "reactant_parse", err.to_string())),
        }
    }
    if reactants.iter().any(|mol| mol.atom_count() > 300) {
        return Ok(refusal(
            "typed_refusal",
            "reactant_too_large",
            "max 300 atoms per reactant".into(),
        ));
    }
    let refs: Vec<_> = reactants.iter().collect();
    let smirks = case["smirks"].as_str().ok_or("non-string SMIRKS")?;
    let limits = rxn::ReactionTransformLimits::default();
    let report = match rxn::run_reactants_traced_rdkit_2026_03_6(smirks, &refs, &limits) {
        Ok(rxn::RdkitProfileOutcome::Report(report)) => report,
        Ok(rxn::RdkitProfileOutcome::Unsupported(unsupported)) => {
            return Ok(refusal(
                "typed_unsupported",
                unsupported.reason_code(),
                "RDKit 2026.03.6 reaction semantics cannot be reproduced for this input".into(),
            ));
        }
        Err(error) => {
            let reason = match &error {
                rxn::TransformError::ResourceLimit { .. } => "resource_limit",
                rxn::TransformError::ReactantCountMismatch { .. } => "reactant_count_mismatch",
                rxn::TransformError::SmirksParse(_) => "smirks_parse",
            };
            return Ok(refusal("typed_refusal", reason, error.to_string()));
        }
    };
    let d = &report.diagnostics;
    let status = if d.valence_rejected_matches > 0 || d.truncated_matches {
        if report.products.is_empty() {
            "typed_refusal"
        } else {
            "partial_products"
        }
    } else if report.products.is_empty() {
        "no_match"
    } else {
        "products"
    };
    let diag = diagnostics(
        d.accepted_matches,
        d.applied_products,
        d.valence_rejected_matches,
        d.truncated_matches,
    );
    if status != "products" && status != "no_match" {
        let reason = if d.truncated_matches {
            "truncated_matches"
        } else {
            "product_valence"
        };
        return Ok(
            json!({"status": status, "reason": reason, "detail": Value::Null,
                         "diagnostics": diag}),
        );
    }
    let mut sets = Vec::new();
    for set in report.products {
        let mut items = Vec::new();
        for traced in set {
            let (spelling, order) = smiles::canonical_smiles_with_atom_order(&traced.molecule);
            let sources: Vec<Value> = order
                .iter()
                .map(|&atom| match traced.atom_sources[atom.0 as usize] {
                    Some(source) => json!([source.reactant, source.atom.0]),
                    None => Value::Null,
                })
                .collect();
            let maps: Vec<Value> = order
                .iter()
                .map(|&atom| json!(traced.template_maps[atom.0 as usize]))
                .collect();
            items.push(json!({"smiles": spelling, "atom_sources": sources,
                              "template_map_numbers": maps}));
        }
        sets.push(Value::Array(items));
    }
    Ok(json!({"status": status, "sets": sets, "diagnostics": diag}))
}

/// The published crate version pinned in Cargo.toml.
const CHEMATIC_VERSION: &str = "1.0.35";

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 4 {
        return Err("usage: published-chematic-checked-gate BASE_JSON STRATA_JSON OUT_JSON".into());
    }
    let base_bytes = fs::read(&args[1])?;
    let strata_bytes = fs::read(&args[2])?;
    let mut cases: Vec<Value> = serde_json::from_slice(&base_bytes)?;
    if cases.len() != 57 {
        return Err("historical base fixture must have 57 cases".into());
    }
    let strata: Value = serde_json::from_slice(&strata_bytes)?;
    if strata["schema_version"] != 2 || strata["base_fixture_sha256"] != digest(&base_bytes) {
        return Err("strata manifest does not pin base fixture hash".into());
    }
    cases.extend(
        strata["cases"]
            .as_array()
            .ok_or("missing strata cases")?
            .iter()
            .cloned(),
    );
    if cases.len() != 83 {
        return Err("expected 83 cases".into());
    }
    let mut rows = Map::new();
    for case in &cases {
        let id = case["id"].as_str().ok_or("non-string id")?.to_owned();
        rows.insert(id, row(case)?);
    }
    let lock = fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.lock"))?;
    let checksum = lock
        .split("[[package]]")
        .find(|block| block.contains("name = \"chematic\"\n"))
        .and_then(|block| block.lines().find(|l| l.starts_with("checksum = ")))
        .map(|l| {
            l.trim_start_matches("checksum = ")
                .trim_matches('"')
                .to_owned()
        });
    let report = json!({
        "schema": "python-checked-reaction-candidates/v1",
        "fixtures": {"base_sha256": digest(&base_bytes), "strata_sha256": digest(&strata_bytes)},
        "artifact": {"kind": "published_crate", "published": true,
                     "published_from": format!("https://crates.io/crates/chematic/{CHEMATIC_VERSION}"),
                     "chematic_version": CHEMATIC_VERSION, "crate_checksum": checksum},
        "rows": rows,
    });
    fs::write(&args[3], serde_json::to_string_pretty(&report)? + "\n")?;
    println!("{}: {} rows", args[3], cases.len());
    Ok(())
}
