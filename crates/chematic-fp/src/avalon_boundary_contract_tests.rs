use crate::rdkit_avalon::{on_bits, rdkit_avalon_fp};
use serde_json::Value;
#[test]
fn avalon_boundary_fingerprints_match_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-avalon-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let text = row["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(text).unwrap();
        let expected: Vec<usize> = row["on_bits"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_u64().unwrap() as usize)
            .collect();
        let actual = on_bits(&rdkit_avalon_fp(&mol, 512).unwrap());
        if actual != expected {
            failures.push(format!("{text}: actual {actual:?}, expected {expected:?}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
