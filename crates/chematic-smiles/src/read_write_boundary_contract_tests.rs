use crate::rdkit_mol_from_pdb_block;

#[test]
fn pdb_fixed_width_records_and_conect_deduplication_match_native_contracts() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-pdb-records-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let result =
            rdkit_mol_from_pdb_block(row["block"].as_str().unwrap(), true, false, 0, false);
        if let Some(reason) = row["unsupported"].as_str() {
            // The public API explicitly declines zero-order bonds and radicals.
            match result {
                Err(error) if error.to_string().contains(reason) => {}
                other => failures.push(format!(
                    "{}: expected {reason} refusal, got {}",
                    row["label"],
                    if other.is_ok() {
                        "success"
                    } else {
                        "another error"
                    }
                )),
            }
        } else if row["valid"].as_bool().unwrap() {
            match result {
                Ok(Some(actual))
                    if actual.smiles == row["smiles"].as_str().unwrap()
                        && actual.molecule.atom_count()
                            == row["atoms"].as_u64().unwrap() as usize
                        && actual.molecule.bond_count()
                            == row["bonds"].as_u64().unwrap() as usize => {}
                Ok(Some(actual)) => failures.push(format!(
                    "{}: {} != {}",
                    row["label"], actual.smiles, row["smiles"]
                )),
                Err(error) => failures.push(format!("{}: {error}", row["label"])),
                Ok(None) => failures.push(format!("{}: unexpectedly absent", row["label"])),
            }
        } else if matches!(result, Ok(Some(_))) {
            failures.push(format!("{}: native reader rejects input", row["label"]));
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 106);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn query_boolean_precedence_ranges_and_reaction_serialization_match_native() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-query-serialization-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for (key, function) in [
        (
            "queries",
            crate::rdkit_smarts_to_smarts as fn(&str) -> Result<String, crate::RdkitSmilesError>,
        ),
        ("reactions", crate::rdkit_reaction_to_smarts),
    ] {
        for row in fixture[key].as_array().unwrap() {
            let input = row["input"].as_str().unwrap();
            let result = function(input);
            match (result, row["expected"].as_str()) {
                (Ok(actual), Some(expected)) if actual == expected => {}
                (Err(_), None) => {}
                (actual, expected) => {
                    failures.push(format!("{key}: {input}: {actual:?} != {expected:?}"))
                }
            }
        }
    }
    assert_eq!(fixture["queries"].as_array().unwrap().len(), 171);
    assert_eq!(fixture["reactions"].as_array().unwrap().len(), 11);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn cleanup_reports_native_charge_bond_and_dative_donor_edits_without_mutation() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-cleanup-edits-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut changed = 0;
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = crate::parse(source).unwrap();
        let before_charges: Vec<_> = mol.atoms().map(|(_, a)| a.charge).collect();
        let before_bonds: Vec<_> = mol
            .bonds()
            .map(|(_, b)| (b.atom1, b.atom2, b.order))
            .collect();
        let actual = crate::rdkit_cleanup_edits(&mol).unwrap();
        let expected_charges: Vec<_> = row["charges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| {
                (
                    r[0].as_u64().unwrap() as usize,
                    r[1].as_i64().unwrap() as i32,
                )
            })
            .collect();
        let expected_bonds: Vec<_> = row["bonds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|r| {
                let order = match r[1].as_str().unwrap() {
                    "SINGLE" => chematic_core::BondOrder::Single,
                    "DATIVE" => chematic_core::BondOrder::Dative,
                    other => panic!("unhandled native edit {other}"),
                };
                (
                    r[0].as_u64().unwrap() as usize,
                    order,
                    r[2].as_u64().map(|i| i as usize),
                )
            })
            .collect();
        if expected_charges.is_empty() && expected_bonds.is_empty() {
            assert_eq!(actual, None, "{source}");
        } else {
            changed += 1;
            let actual = actual.unwrap();
            assert_eq!(actual.charges, expected_charges, "{source}");
            assert_eq!(actual.bonds, expected_bonds, "{source}");
        }
        assert_eq!(
            mol.atoms().map(|(_, a)| a.charge).collect::<Vec<_>>(),
            before_charges
        );
        assert_eq!(
            mol.bonds()
                .map(|(_, b)| (b.atom1, b.atom2, b.order))
                .collect::<Vec<_>>(),
            before_bonds
        );
    }
    assert_eq!(changed, 6);
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 13);
}

#[test]
fn molecular_hashes_of_ions_radicals_and_unusual_orders_match_native_defaults() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-hash-ions-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let mol = if source.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            crate::parse(source).unwrap()
        };
        let function =
            crate::RdkitHashFunction::from_name(row["function"].as_str().unwrap()).unwrap();
        let result = crate::rdkit_mol_hash(&mol, function, false);
        match result {
            Ok(actual) if actual == row["expected"].as_str().unwrap() => {}
            other => failures.push(format!(
                "{source},{function:?}: {other:?} != {}",
                row["expected"]
            )),
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 874);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
