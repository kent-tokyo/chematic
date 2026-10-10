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
        let bits = row["n_bits"].as_u64().unwrap_or(512) as usize;
        let actual = on_bits(&avalon_fp_bytes(&av, bits, flags));
        if actual != expected {
            failures.push(format!(
                "{text},bits={bits},flags={flags}: actual {actual:?}, expected {expected:?}"
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn avalon_v3000_charge_radical_and_isotope_records_match_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-avalon-molblock-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    for row in fixture["rows"].as_array().unwrap() {
        let block = row["block"].as_str().unwrap();
        let parsed = read_molblock(block).unwrap();
        let actual = on_bits(&avalon_fp_bytes(&parsed, 512, 15761407));
        let expected: Vec<usize> = serde_json::from_value(row["on_bits"].clone()).unwrap();
        assert_eq!(actual, expected, "{}", row["smiles"]);
        let lines: Vec<_> = block.lines().collect();
        let atom_line = lines
            .iter()
            .position(|s| s.starts_with("M  V30 1 "))
            .unwrap();
        let mut corrupt = lines.clone();
        corrupt[atom_line] = "M  V30 1 C";
        assert!(read_molblock(&corrupt.join("\n")).is_none());
    }
}
