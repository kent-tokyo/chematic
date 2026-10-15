//! Empty ensembles, invalid controls, and disabled rescue must have stable results.
use chematic_3d::rdkit_embed::{
    RdkitEmbedError, RdkitEmbedMmffError, RdkitEmbedOptions, rdkit_embed_multiple_confs,
    rdkit_embed_multiple_confs_mmff, rdkit_embed_multiple_confs_mmff_adaptive,
    rdkit_embed_multiple_confs_mmff_diagnostics, rdkit_embed_multiple_confs_mmff_gradient_rescue,
    rdkit_embed_multiple_confs_mmff_rescue,
};

#[test]
fn zero_conformers_stay_empty_through_every_optimizer() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("CO").unwrap());
    let options = RdkitEmbedOptions::default();
    assert!(
        rdkit_embed_multiple_confs_mmff(&mol, 0, &options, 10)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_diagnostics(&mol, 0, &options, 10)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_adaptive(&mol, 0, &options, 1, 10, 1)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_rescue(&mol, 0, &options, 1, 10, 1)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_gradient_rescue(&mol, 0, &options, 1, 10, 1, 0.0)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn invalid_controls_remain_typed_errors_for_empty_requests() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("CO").unwrap());
    let options = RdkitEmbedOptions::default();
    for threshold in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let error =
            rdkit_embed_multiple_confs_mmff_gradient_rescue(&mol, 0, &options, 1, 10, 1, threshold)
                .unwrap_err();
        assert!(matches!(error, RdkitEmbedMmffError::ForceField(_)));
        assert_eq!(
            error.to_string(),
            "RDKit-compatible MMFF94: minimum gradient test must be finite and non-negative"
        );
    }
    let invalid = RdkitEmbedOptions {
        random_seed: -1,
        ..options
    };
    let error = rdkit_embed_multiple_confs_mmff(&mol, 0, &invalid, 10).unwrap_err();
    assert!(matches!(
        error,
        RdkitEmbedMmffError::Embed(RdkitEmbedError::Unsupported(_))
    ));
    assert!(error.to_string().contains("randomSeed must be >= 0"));
}

#[test]
fn conformer_seed_overflow_is_reported_instead_of_wrapping() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("CO").unwrap());
    let options = RdkitEmbedOptions {
        random_seed: i32::MAX,
        ..RdkitEmbedOptions::default()
    };
    let error = rdkit_embed_multiple_confs(&mol, 2, &options).unwrap_err();
    assert!(matches!(error, RdkitEmbedError::Unsupported(_)));
    assert!(error.to_string().contains("conformer seed overflows"));
}

#[test]
fn disabling_gradient_rescue_preserves_initial_optimization() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let initial = rdkit_embed_multiple_confs_mmff(&mol, 2, &options, 1).unwrap();
    assert_eq!(initial.len(), 2);
    for (iterations, candidates) in [(10, 0), (0, 1)] {
        let result = rdkit_embed_multiple_confs_mmff_gradient_rescue(
            &mol, 2, &options, 1, iterations, candidates, 0.0,
        )
        .unwrap();
        assert_eq!(result.len(), initial.len());
        for (actual, expected) in result.iter().zip(&initial) {
            assert_eq!(actual.coords, expected.coords);
            assert_eq!(actual.status, expected.status);
            assert_eq!(actual.energy, expected.energy);
        }
    }
}
