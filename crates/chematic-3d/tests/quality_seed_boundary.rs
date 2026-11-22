use chematic_3d::rdkit_embed::{
    RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_best_gradient,
    rdkit_embed_multiple_confs_mmff,
};

#[test]
fn quality_ensemble_preserves_first_batch_when_next_seed_would_overflow() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let options = RdkitEmbedOptions {
        random_seed: i32::MAX,
        ..Default::default()
    };
    let compatibility = rdkit_embed_multiple_confs_mmff(&mol, 1, &options, 0).unwrap();
    let quality = etkdg_embed_multiple_confs_mmff_best_gradient(&mol, 1, &options, 0).unwrap();
    assert_eq!(compatibility.len(), 1);
    assert_eq!(quality.len(), 1);
    assert_eq!(quality[0].coords, compatibility[0].coords);
    assert_eq!(quality[0].energy, compatibility[0].energy);
}

#[test]
fn quality_ensemble_with_zero_requested_conformers_is_empty() {
    let mol = chematic_smiles::parse("CCO").unwrap();
    assert!(
        etkdg_embed_multiple_confs_mmff_best_gradient(&mol, 0, &RdkitEmbedOptions::default(), 10)
            .unwrap()
            .is_empty()
    );
}
