//! Exposed reaction map-label diagnostic from the published v1.0.30 crate.
//! Map numbers are recovered from the product template's atom order. This
//! includes newly created mapped atoms, which have no traced reactant origin.
//! The ordinary product molecule intentionally clears map annotations.

use std::collections::HashMap;
use std::fs;

use chematic::{core::AtomIdx, rxn, smiles};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5 && args.len() != 6 {
        return Err(
            "usage: reaction_template_maps BASE_JSON STRATA_JSON [SUPPLEMENT_JSON] ROWS_JSON SUMMARY_JSON".into(),
        );
    }
    let output_offset = usize::from(args.len() == 6);
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
    let supplement_sha256 = if output_offset == 1 {
        let bytes = fs::read(&args[3])?;
        let supplement: Value = serde_json::from_slice(&bytes)?;
        if supplement["base_fixture_sha256"] != digest(&base_bytes)
            || supplement["strata_fixture_sha256"] != digest(&strata_bytes)
        {
            return Err("supplement does not pin base and strata fixture hashes".into());
        }
        cases.extend(
            supplement["cases"]
                .as_array()
                .ok_or("missing supplement cases")?
                .iter()
                .cloned(),
        );
        Some(digest(&bytes))
    } else {
        None
    };

    let mut rows = Vec::with_capacity(cases.len());
    'case_loop: for case in &cases {
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
            rows.push(json!({"id": case["id"], "status": "typed_refusal",
                "stage": "reactant_parse", "reason": "smiles_parse", "detail": detail}));
            continue;
        }
        let smirks = case["smirks"].as_str().ok_or("non-string SMIRKS")?;
        let variants = match rxn::expand_atomic_number_primitives(smirks) {
            Ok(value) => value,
            Err(err) => {
                rows.push(json!({"id": case["id"], "status": "typed_refusal",
                    "stage": "smirks_parse", "reason": "smirks_parse", "detail": err.to_string()}));
                continue;
            }
        };
        let prepared_variants = match variants
            .iter()
            .map(|variant| rxn::PreparedReaction::new(variant))
            .collect::<Result<Vec<_>, _>>()
        {
            Ok(value) => value,
            Err(err) => {
                rows.push(json!({"id": case["id"], "status": "typed_refusal",
                    "stage": "smirks_parse", "reason": "smirks_parse", "detail": err.to_string()}));
                continue;
            }
        };
        let refs: Vec<_> = reactants.iter().collect();
        let mut accepted_matches = 0;
        let mut valence_rejected_matches = 0;
        let mut product_sets = Vec::new();
        for (variant, prepared) in variants.iter().zip(prepared_variants.iter()) {
            let template_products = rxn::parse_reaction(variant)?.products;
            let matches = match prepared.find_matches(&refs) {
                Ok(value) => value,
                Err(err) => {
                    rows.push(json!({"id": case["id"], "status": "typed_refusal",
                        "stage": "match", "reason": "match_error", "detail": err.to_string()}));
                    continue 'case_loop;
                }
            };
            accepted_matches += matches.len();
            for matched in &matches {
                let map_positions = matched.atom_map_positions(variant)?;
                let mut map_by_source: HashMap<(usize, AtomIdx), u16> = HashMap::new();
                for (map_number, source) in map_positions {
                    if map_by_source.insert(source, map_number).is_some() {
                        return Err(format!(
                            "{}: two template maps resolve to one source atom",
                            case["id"]
                        )
                        .into());
                    }
                }
                let Some(products) = prepared.apply_match_traced(&refs, matched, true)? else {
                    valence_rejected_matches += 1;
                    continue;
                };
                let mut product_set = Vec::new();
                if products.len() != template_products.len() {
                    return Err(format!("{}: product/template count differs", case["id"]).into());
                }
                for (product, template) in products.into_iter().zip(&template_products) {
                    let (spelling, order) =
                        smiles::canonical_smiles_with_atom_order(&product.molecule);
                    if order.len() != product.atom_sources.len() {
                        return Err(
                            format!("{}: atom output order is incomplete", case["id"]).into()
                        );
                    }
                    let mut atom_sources = Vec::with_capacity(order.len());
                    let mut template_map_numbers = Vec::with_capacity(order.len());
                    for atom in order {
                        let source = product.atom_sources[atom.0 as usize];
                        let map_number = if (atom.0 as usize) < template.atom_count() {
                            template.atom(atom).atom_map
                        } else {
                            None
                        };
                        if let (Some(label), Some(value)) = (map_number, source)
                            && map_by_source.get(&(value.reactant, value.atom)) != Some(&label)
                        {
                            return Err(
                                format!("{}: template/source map differs", case["id"]).into()
                            );
                        }
                        atom_sources.push(
                            source
                                .map(|value| json!([value.reactant, value.atom.0]))
                                .unwrap_or(Value::Null),
                        );
                        template_map_numbers.push(map_number);
                    }
                    product_set.push(json!({"smiles": spelling, "atom_sources": atom_sources,
                        "template_map_numbers": template_map_numbers}));
                }
                product_set.sort_by_key(|value| value["smiles"].as_str().unwrap_or("").to_owned());
                product_sets.push(product_set);
            }
        }
        let status = if product_sets.is_empty() && valence_rejected_matches > 0 {
            "diagnosed_valence_refusal"
        } else {
            "products"
        };
        rows.push(
            json!({"id": case["id"], "status": status, "sets": product_sets,
            "accepted_matches": accepted_matches,
            "valence_rejected_matches": valence_rejected_matches}),
        );
    }
    let raw = serde_json::to_vec_pretty(&rows)?;
    fs::write(&args[3 + output_offset], &raw)?;
    let mut summary = json!({"schema": "published-rust-reaction-template-maps/v1",
        "crate": "chematic 1.0.30 from crates.io", "input_count": cases.len(),
        "base_cases_sha256": digest(&base_bytes), "strata_sha256": digest(&strata_bytes),
        "rows_sha256": digest(&raw)});
    if let Some(hash) = supplement_sha256 {
        summary["supplement_sha256"] = json!(hash);
    }
    fs::write(
        &args[4 + output_offset],
        serde_json::to_vec_pretty(&summary)?,
    )?;
    println!("{summary}");
    Ok(())
}
