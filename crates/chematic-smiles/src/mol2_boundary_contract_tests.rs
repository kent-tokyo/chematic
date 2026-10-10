use crate::rdkit_mol_from_mol2_block;
use serde_json::Value;
#[test]
fn mol2_corina_cleanup_and_double_bond_directions_match_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-mol2-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let block = row["block"].as_str().unwrap();
        let expected = row["smiles"].as_str().unwrap();
        let result = rdkit_mol_from_mol2_block(block, true, true, true);
        match result {
            Ok(result) if result.smiles == expected => {}
            Ok(result) => failures.push(format!(
                "{}: {} != {expected}",
                row["source_smiles"], result.smiles
            )),
            Err(error) => failures.push(format!("{}: {error}", row["source_smiles"])),
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn mol2_sanitization_hydrogen_and_cleanup_flags_match_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-mol2-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut checked = 0;
    for row in fixture["rows"].as_array().unwrap() {
        let variants = row["variants"].as_array().unwrap();
        assert_eq!(variants.len(), 8);
        for expected in variants {
            let result = rdkit_mol_from_mol2_block(
                row["block"].as_str().unwrap(),
                expected["sanitize"].as_bool().unwrap(),
                expected["remove_hs"].as_bool().unwrap(),
                expected["cleanup"].as_bool().unwrap(),
            );
            let context = format!("{} {expected}", row["source_smiles"]);
            if !expected["accepted"].as_bool().unwrap() {
                assert!(result.is_err(), "{context}");
            } else {
                let actual = result.unwrap_or_else(|e| panic!("{context}: {e}"));
                assert_eq!(
                    actual.smiles,
                    expected["smiles"].as_str().unwrap(),
                    "{context}"
                );
                assert_eq!(
                    actual.molecule.atom_count() as u64,
                    expected["atom_count"].as_u64().unwrap(),
                    "{context}"
                );
                assert_eq!(
                    actual.molecule.bond_count() as u64,
                    expected["bond_count"].as_u64().unwrap(),
                    "{context}"
                );
                assert_eq!(
                    actual.coords.len(),
                    actual.molecule.atom_count(),
                    "{context}"
                );
                assert!(
                    actual.coords.iter().flatten().all(|v| v.is_finite()),
                    "{context}"
                );
                let numbers: Vec<_> = actual
                    .molecule
                    .atoms()
                    .map(|(_, a)| a.element.atomic_number())
                    .collect();
                let charges: Vec<_> = actual.molecule.atoms().map(|(_, a)| a.charge).collect();
                assert_eq!(
                    serde_json::json!(numbers),
                    expected["atomic_numbers"],
                    "{context}"
                );
                assert_eq!(serde_json::json!(charges), expected["charges"], "{context}");
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 672);
}

#[test]
fn mol2_phosphate_quaternary_n_and_invalid_sybyl_types_match_native_cleanup() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-mol2-ionic-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let result = rdkit_mol_from_mol2_block(
            row["block"].as_str().unwrap(),
            true,
            true,
            row["cleanup"].as_bool().unwrap(),
        );
        match (result, row["expected"].as_str()) {
            (Ok(actual), Some(expected)) if actual.smiles == expected => {}
            (Err(_), None) => {}
            (Ok(actual), expected) => failures.push(format!(
                "{} {} cleanup={}: {} != {expected:?}",
                row["smiles"], row["mode"], row["cleanup"], actual.smiles
            )),
            (Err(error), expected) => failures.push(format!(
                "{} {} cleanup={}: {error} != {expected:?}",
                row["smiles"], row["mode"], row["cleanup"]
            )),
        }
    }
    assert_eq!(fixture["rows"].as_array().unwrap().len(), 36);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
