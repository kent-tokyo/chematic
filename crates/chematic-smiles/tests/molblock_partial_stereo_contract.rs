use chematic_smiles::{RdkitSmilesError, parse, rdkit_mol_block_2d};

#[test]
fn mol_block_refuses_partial_conjugated_stereo_that_would_be_invented() {
    for input in [r"F/C=C/C=CF", r"C/C=C/C=CC"] {
        let error = rdkit_mol_block_2d(&parse(input).unwrap()).unwrap_err();
        assert!(matches!(&error, RdkitSmilesError::Unsupported(_)));
        assert!(
            error.to_string().contains(
                "an unspecified double bond next to directed bonds would read back with stereo"
            ),
            "{input}: {error}"
        );
    }
    assert!(rdkit_mol_block_2d(&parse(r"F/C=C/C=C/F").unwrap()).is_ok());
}
