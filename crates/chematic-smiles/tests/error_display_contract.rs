use chematic_smiles::RdkitSmilesError;

#[test]
fn rdkit_smiles_errors_distinguish_scope_from_invalid_chemistry() {
    assert_eq!(
        RdkitSmilesError::Unsupported("query bond".into()).to_string(),
        "RDKit-compatible SMILES: unsupported input: query bond"
    );
    assert_eq!(
        RdkitSmilesError::Sanitization("valence".into()).to_string(),
        "RDKit-compatible SMILES: sanitization failed: valence"
    );
}
