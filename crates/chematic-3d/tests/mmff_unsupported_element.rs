use chematic_3d::rdkit_embed::{
    RdkitEmbedMmffError, RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_best_gradient,
    etkdg_embed_multiple_confs_mmff_corrected, rdkit_embed_multiple_confs_mmff,
    rdkit_embed_multiple_confs_mmff_diagnostics,
};

#[test]
fn unsupported_mmff_elements_report_force_field_failure() {
    let molecule = chematic_smiles::parse("[He]").unwrap();
    let options = RdkitEmbedOptions::default();
    let errors = [
        rdkit_embed_multiple_confs_mmff(&molecule, 1, &options, 1).unwrap_err(),
        etkdg_embed_multiple_confs_mmff_corrected(&molecule, 1, &options, 1).unwrap_err(),
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 1, &options, 1).unwrap_err(),
        rdkit_embed_multiple_confs_mmff_diagnostics(&molecule, 1, &options, 1).unwrap_err(),
    ];
    for error in errors {
        assert!(matches!(&error, RdkitEmbedMmffError::ForceField(_)));
        assert!(error.to_string().contains("unsupported element"), "{error}");
    }
}

#[test]
fn requesting_no_conformers_does_not_require_mmff_parameters() {
    let molecule = chematic_smiles::parse("[He]").unwrap();
    let options = RdkitEmbedOptions::default();
    assert!(
        rdkit_embed_multiple_confs_mmff(&molecule, 0, &options, 1)
            .unwrap()
            .is_empty()
    );
    assert!(
        etkdg_embed_multiple_confs_mmff_corrected(&molecule, 0, &options, 1)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_diagnostics(&molecule, 0, &options, 1)
            .unwrap()
            .is_empty()
    );
}
