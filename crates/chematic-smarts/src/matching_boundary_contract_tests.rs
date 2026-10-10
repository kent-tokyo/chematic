use crate::{RdkitParityConfig, find_match_atom_sets_rdkit_parity, parse_smarts};
#[test]
fn parity_matcher_recursive_and_atom_primitives_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-smarts-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let source = row["smiles"].as_str().unwrap();
        let smarts = row["smarts"].as_str().unwrap();
        let mol = if source.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            chematic_smiles::parse(source).unwrap()
        };
        let query = parse_smarts(smarts).unwrap();
        let expected: Vec<Vec<u32>> = serde_json::from_value(row["sets"].clone()).unwrap();
        for shared in [false, true] {
            let config = RdkitParityConfig {
                use_shared_symmetrized_sssr: shared,
                ..Default::default()
            };
            let (actual, exhausted) =
                find_match_atom_sets_rdkit_parity(&query, &mol, &config).unwrap();
            assert!(!exhausted);
            if actual != expected {
                failures.push(format!(
                    "{source}: {smarts} shared={shared}: {actual:?} != {expected:?}"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
