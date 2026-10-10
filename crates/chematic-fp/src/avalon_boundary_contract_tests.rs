use crate::rdkit_avalon::{avalon_fp_bytes, on_bits, read_molblock};
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
        let block = chematic_smiles::rdkit_mol_block_2d(&mol).unwrap();
        let av = read_molblock(&block).unwrap();
        let flags = row["bit_flags"].as_u64().unwrap() as u32;
        let actual = on_bits(&avalon_fp_bytes(&av, 512, flags));
        if actual != expected {
            failures.push(format!(
                "{text},flags={flags}: actual {actual:?}, expected {expected:?}"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
