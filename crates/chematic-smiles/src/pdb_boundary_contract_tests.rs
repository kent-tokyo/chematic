use crate::rdkit_mol_from_pdb_block;
#[test]
fn pdb_coordination_geometries_and_malformed_fields_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-pdb-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let result = rdkit_mol_from_pdb_block(
            row["block"].as_str().unwrap(),
            row["sanitize"].as_bool().unwrap(),
            row["remove_hs"].as_bool().unwrap(),
            row["flavor"].as_u64().unwrap() as u32,
            row["proximity"].as_bool().unwrap(),
        );
        let valid = row["valid"].as_bool().unwrap();
        match result {
            Ok(Some(p)) if valid => {
                if p.smiles != row["smiles"].as_str().unwrap()
                    || p.molecule.atom_count() != row["atoms"].as_u64().unwrap() as usize
                    || p.molecule.bond_count() != row["bonds"].as_u64().unwrap() as usize
                {
                    failures.push(format!(
                        "{}: {} atoms, {} bonds, {} != {row}",
                        row["label"],
                        p.molecule.atom_count(),
                        p.molecule.bond_count(),
                        p.smiles
                    ));
                }
            }
            Ok(None) | Err(_) if !valid => {}
            Ok(_) => failures.push(format!(
                "{}: unexpected outcome (valid={valid})",
                row["label"]
            )),
            Err(e) => failures.push(format!("{}: {e}", row["label"])),
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
