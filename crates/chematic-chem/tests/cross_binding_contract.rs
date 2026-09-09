use chematic_chem::TransformationRecord;
use serde_json::Value;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cross_binding_contract.json"
));

#[test]
fn shared_fragment_parent_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["fragment_parent_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let (parent, _) = chematic_chem::fragment_parent(&molecule);
        assert_eq!(
            chematic_smiles::canonical_smiles(&parent),
            case["canonical_smiles"],
            "fragment parent mismatch for {smiles}"
        );
    }
}

fn assert_parent_contract(
    document: &Value,
    section: &str,
    transform: fn(&chematic_core::Molecule) -> (chematic_core::Molecule, TransformationRecord),
) {
    let contract = &document[section];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let (parent, _) = transform(&molecule);
        assert_eq!(
            chematic_smiles::canonical_smiles(&parent),
            case["canonical_smiles"],
            "parent mismatch for {section}/{smiles}"
        );
    }
}

#[test]
fn shared_parent_transform_contracts_match() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    assert_parent_contract(
        &document,
        "charge_parent_contract",
        chematic_chem::charge_parent,
    );
    assert_parent_contract(
        &document,
        "isotope_parent_contract",
        chematic_chem::isotope_parent,
    );
    assert_parent_contract(
        &document,
        "stereo_parent_contract",
        chematic_chem::stereo_parent,
    );
}

#[test]
fn shared_canonical_tautomer_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["canonical_tautomer_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let tautomer = chematic_chem::canonical_tautomer(&molecule);
        assert_eq!(
            chematic_smiles::canonical_smiles(&tautomer),
            case["canonical_smiles"],
            "canonical tautomer mismatch for {smiles}"
        );
    }
}

#[test]
fn shared_neutralize_charges_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["neutralize_charges_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let neutral = chematic_chem::neutralize_charges(&molecule);
        assert_eq!(
            chematic_smiles::canonical_smiles(&neutral),
            case["canonical_smiles"],
            "neutralization mismatch for {smiles}"
        );
    }
}

#[test]
fn shared_largest_fragment_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["largest_fragment_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let fragment = chematic_chem::largest_fragment(&molecule);
        assert_eq!(
            chematic_smiles::canonical_smiles(&fragment),
            case["canonical_smiles"],
            "largest fragment mismatch for {smiles}"
        );
    }
}

fn assert_budgeted_parent_contract(
    document: &Value,
    section: &str,
    transform: fn(
        &chematic_core::Molecule,
        &chematic_chem::TautomerLimits,
    ) -> chematic_chem::ParentResult,
) {
    let contract = &document[section];
    assert_eq!(contract["schema_version"], 1);
    let limits_json = &contract["limits"];
    let mut limits = chematic_chem::TautomerLimits::default();
    limits.max_transforms = limits_json["max_transforms"].as_u64().unwrap() as usize;
    limits.max_tautomers = limits_json["max_tautomers"].as_u64().unwrap() as usize;
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let result = transform(&molecule, &limits);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result.molecule),
            case["canonical_smiles"],
            "parent mismatch for {section}/{smiles}"
        );
        assert_eq!(format!("{:?}", result.status), case["status"]);
    }
}

#[test]
fn shared_budgeted_parent_contracts_match() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    assert_budgeted_parent_contract(
        &document,
        "tautomer_parent_contract",
        chematic_chem::tautomer_parent,
    );
    assert_budgeted_parent_contract(
        &document,
        "super_parent_contract",
        chematic_chem::super_parent,
    );
}

#[test]
fn shared_hydrogen_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["hydrogen_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let added = chematic_chem::add_hydrogens(&molecule);
        assert_eq!(
            chematic_smiles::canonical_smiles(&added),
            case["added_canonical_smiles"],
            "add hydrogens mismatch for {smiles}"
        );
        let removed = chematic_chem::remove_hydrogens(&added);
        assert_eq!(
            chematic_smiles::canonical_smiles(&removed),
            case["removed_canonical_smiles"],
            "remove hydrogens mismatch for {smiles}"
        );
    }
}

#[test]
fn shared_super_parent_report_contract_matches() {
    let document: Value = serde_json::from_str(FIXTURE).expect("fixture JSON must parse");
    let contract = &document["super_parent_report_contract"];
    assert_eq!(contract["schema_version"], 1);
    for case in contract["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let molecule = chematic_smiles::parse(smiles).expect("fixture SMILES must parse");
        let mut limits = chematic_chem::TautomerLimits::default();
        limits.max_transforms = contract["limits"]["max_transforms"].as_u64().unwrap() as usize;
        limits.max_tautomers = contract["limits"]["max_tautomers"].as_u64().unwrap() as usize;
        let result = chematic_chem::super_parent(&molecule, &limits);
        assert_eq!(
            chematic_smiles::canonical_smiles(&result.molecule),
            case["canonical_smiles"],
            "super-parent report result mismatch for {smiles}"
        );
        assert_eq!(format!("{:?}", result.status), case["status"]);
        match result.audit {
            chematic_chem::ParentAudit::Composed(stages) => {
                assert_eq!(
                    stages.len(),
                    contract["stage_names"].as_array().unwrap().len()
                );
            }
            other => panic!("expected composed parent audit, got {other:?}"),
        }
    }
}
