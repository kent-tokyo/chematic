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
