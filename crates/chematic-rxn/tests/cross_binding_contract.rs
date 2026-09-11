use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn reaction_document_contract_round_trips_through_typed_model() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["rxn_document_contract"];
    assert_eq!(contract["schema_version"], 1);
    let authored = serde_json::to_string(&contract["document"]).unwrap();
    let parsed = chematic_rxn::ReactionDocument::from_json_str(&authored).unwrap();
    let observed: Vec<Value> = parsed.steps[0]
        .components
        .iter()
        .map(|component| serde_json::json!({"role": component.role, "smiles": component.smiles}))
        .collect();
    let expected: Vec<Value> =
        serde_json::from_value(contract["expected_components"].clone()).unwrap();
    assert_eq!(observed, expected);
}

#[test]
fn reaction_application_contract_expands_atomic_number_primitives() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["reaction_application_contract"];
    let reactant = chematic_smiles::parse(contract["reactants"][0].as_str().unwrap()).unwrap();
    let products =
        chematic_rxn::run_reactants(contract["smirks"].as_str().unwrap(), &[&reactant]).unwrap();
    assert!(
        products.len()
            >= contract["expected"]["minimum_product_sets"]
                .as_u64()
                .unwrap() as usize
    );
    assert_eq!(
        products[0][0].atom_count(),
        contract["expected"]["product_atom_count"]
    );
    for case in contract["negative_cases"].as_array().unwrap() {
        let reactants: Vec<_> = case["reactants"]
            .as_array()
            .unwrap()
            .iter()
            .map(|smiles| chematic_smiles::parse(smiles.as_str().unwrap()).unwrap())
            .collect();
        let refs: Vec<_> = reactants.iter().collect();
        assert!(chematic_rxn::run_reactants(case["smirks"].as_str().unwrap(), &refs).is_err());
    }
    for case in contract["additional_cases"].as_array().unwrap() {
        let reactants: Vec<_> = case["reactants"]
            .as_array()
            .unwrap()
            .iter()
            .map(|smiles| chematic_smiles::parse(smiles.as_str().unwrap()).unwrap())
            .collect();
        let refs: Vec<_> = reactants.iter().collect();
        let products =
            chematic_rxn::run_reactants(case["smirks"].as_str().unwrap(), &refs).unwrap();
        if let Some(expected) = case["expected_product_sets"].as_u64() {
            assert_eq!(products.len(), expected as usize);
            continue;
        }
        assert!(products.len() >= case["minimum_product_sets"].as_u64().unwrap() as usize);
        assert_eq!(products[0][0].atom_count(), case["product_atom_count"]);
    }
}

#[test]
fn reaction_smarts_contract_matches_shared_fixture() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    for case in document["reaction_smarts_contract"]["cases"]
        .as_array()
        .expect("reaction SMARTS cases must be an array")
    {
        let reaction = chematic_rxn::parse_reaction(case["reaction"].as_str().unwrap()).unwrap();
        let query = chematic_rxn::parse_reaction_query(case["smarts"].as_str().unwrap()).unwrap();
        assert_eq!(
            chematic_rxn::has_reaction_substructure_match(&reaction, &query),
            case["matches"],
            "reaction SMARTS case {} against {}",
            case["smarts"],
            case["reaction"]
        );
    }
}

#[test]
fn reaction_balance_contract_matches_shared_fixture() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    for case in document["reaction_balance_contract"]["cases"]
        .as_array()
        .expect("balance cases must be an array")
    {
        let reaction = chematic_rxn::parse_reaction(case["reaction"].as_str().unwrap()).unwrap();
        let result = chematic_rxn::balance_check(&reaction);
        assert_eq!(result.balanced, case["balanced"]);
        let expected: Vec<String> = serde_json::from_value(case["diff"].clone()).unwrap();
        assert_eq!(result.diff(), expected);
    }
}

#[test]
fn reaction_center_contract_matches_shared_fixture() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    for case in document["reaction_center_contract"]["cases"]
        .as_array()
        .expect("center cases must be an array")
    {
        let reaction = chematic_rxn::parse_reaction(case["reaction"].as_str().unwrap()).unwrap();
        let center = chematic_rxn::find_reaction_center(&reaction);
        let changed: Vec<usize> = center
            .changed_atoms
            .iter()
            .map(|atom| atom.0 as usize)
            .collect();
        assert_eq!(
            center.broken_bonds.len(),
            case["broken_bonds"].as_array().unwrap().len()
        );
        assert_eq!(
            center.formed_bonds.len(),
            case["formed_bonds"].as_array().unwrap().len()
        );
        let expected: Vec<usize> = serde_json::from_value(case["changed_atoms"].clone()).unwrap();
        assert_eq!(changed, expected);
    }
}
