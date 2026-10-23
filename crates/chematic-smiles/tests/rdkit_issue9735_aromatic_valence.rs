//! Minimum explicit valence must not disappear in aromatic correction (#9735).
use chematic_smiles::{
    RdkitSmilesError, parse, rdkit_canonical_smiles, rdkit_detect_chemistry_problems,
    rdkit_sanitized_model,
};

#[test]
fn aromatic_overvalence_is_refused_before_kekulization() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-9735-2026-10-24/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        // Parsing retains the graph; strict chemistry validation is a separate
        // compatibility surface and must also handle manually built graphs.
        let mol = parse(smiles).unwrap();
        if case["expected"] == "accept" {
            let first = rdkit_canonical_smiles(&mol).unwrap_or_else(|e| panic!("{smiles}: {e}"));
            assert!(rdkit_sanitized_model(&mol).is_ok(), "{smiles}");
            assert_eq!(
                rdkit_canonical_smiles(&parse(&first).unwrap()).unwrap(),
                first,
                "{smiles}"
            );
            assert!(
                rdkit_detect_chemistry_problems(&mol).unwrap().is_empty(),
                "{smiles}"
            );
            continue;
        }
        let error = rdkit_canonical_smiles(&mol).unwrap_err();
        assert!(
            matches!(error, RdkitSmilesError::Sanitization(_)),
            "{smiles}: {error}"
        );
        assert!(
            error.to_string().contains("greater than permitted"),
            "{smiles}: {error}"
        );
        assert!(rdkit_sanitized_model(&mol).is_err(), "{smiles}");
        let expected = case["atom"].as_u64().unwrap() as usize;
        let problems = rdkit_detect_chemistry_problems(&mol).unwrap();
        assert!(
            problems
                .iter()
                .any(|p| p.kind == "AtomValenceException" && p.atoms == [expected]),
            "{smiles}: {problems:?}"
        );
    }
}
