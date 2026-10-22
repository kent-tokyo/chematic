use super::*;
use crate::stereo_constraints::verify_stereo;
use chematic_core::{Coords3D, Point3};

#[test]
fn quality_profile_uses_new_seeds_when_gradients_match() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("CC(C)Cc1ccc(cc1)C(C)C(=O)O").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let prepared = prepare_fragments(&molecule, &options).unwrap();
    assert!(!corrected_torsion_gradient_can_differ(&prepared));

    let first = rdkit_embed_multiple_confs_mmff(&molecule, 4, &options, 200).unwrap();
    let next_coords = embed_prepared_multiple_confs(&prepared, 4, 4, &options, false).unwrap();
    let next = optimize_mmff_molecule(&molecule, next_coords, 200).unwrap();
    let third_coords = embed_prepared_multiple_confs(&prepared, 4, 8, &options, false).unwrap();
    let third = optimize_mmff_molecule(&molecule, third_coords, 200).unwrap();
    assert_eq!(first.len(), 4);
    assert_eq!(next.len(), 4);
    assert_eq!(third.len(), 4);
    assert!(first.iter().zip(&next).any(|(a, b)| a.coords != b.coords));
    assert!(next.iter().zip(&third).any(|(a, b)| a.coords != b.coords));

    let selected =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 4, &options, 200).unwrap();
    let best = |ensemble: &[RdkitMmffConformer]| {
        ensemble
            .iter()
            .map(|candidate| candidate.energy)
            .min_by(f64::total_cmp)
            .unwrap()
    };
    assert_eq!(selected.len(), 4);
    assert_eq!(
        best(&selected),
        best(&first).min(best(&next)).min(best(&third))
    );
    assert!(selected.iter().all(|chosen| {
        first
            .iter()
            .chain(&next)
            .chain(&third)
            .any(|original| chosen.coords == original.coords)
    }));
}

#[test]
fn quality_profile_retries_best_unconverged_candidate() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("COc1cc(N2CCN(c3cc(OC)c(OC)c(OC)c3)CC2)cc(OC)c1OC").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let initial = rdkit_embed_multiple_confs_mmff(&molecule, 40, &options, 200).unwrap();
    assert_eq!(initial.len(), 40);
    assert!(initial.iter().all(|candidate| candidate.status != 0));
    let best_initial = initial
        .iter()
        .map(|candidate| candidate.energy)
        .min_by(f64::total_cmp)
        .unwrap();

    let selected =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 20, &options, 200).unwrap();
    assert_eq!(selected.len(), 20);
    assert_eq!(selected[0].status, 0);
    assert!(selected[0].energy <= best_initial);
    let coordinates = Coords3D {
        points: selected[0]
            .coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let stereo = verify_stereo(&molecule, &coordinates);
    assert_eq!(stereo.n_violations(), 0);
    assert_eq!(stereo.n_unevaluable(), 0);
}

#[test]
fn quality_profile_retries_unfinished_minimum_even_with_other_converged_candidates() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("O=C(C1CCCCN1)N1CCN(Cc2cccnc2)CC1").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let initial = rdkit_embed_multiple_confs_mmff(&molecule, 40, &options, 200).unwrap();
    let best = initial
        .iter()
        .min_by(|left, right| left.energy.total_cmp(&right.energy))
        .unwrap();
    let previous_converged_energy = initial
        .iter()
        .filter(|candidate| candidate.status == 0)
        .map(|candidate| candidate.energy)
        .min_by(f64::total_cmp)
        .unwrap();
    assert_ne!(best.status, 0);
    assert!(best.energy < previous_converged_energy);

    let selected =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 20, &options, 200).unwrap();
    assert_eq!(selected[0].status, 0);
    assert!(selected[0].energy < previous_converged_energy - 0.1);
}

#[test]
fn quality_profile_retries_new_seed_batch_minimum() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("COc1cc2nc(N3CCN(/C(S)=N/c4ccccc4)CC3)nc(N)c2cc1OC").unwrap(),
    );
    let selected = etkdg_embed_multiple_confs_mmff_best_gradient(
        &molecule,
        20,
        &RdkitEmbedOptions::default(),
        200,
    )
    .unwrap();
    assert_eq!(selected.len(), 20);
    let best = selected
        .iter()
        .min_by(|left, right| left.energy.total_cmp(&right.energy))
        .unwrap();
    assert_eq!(best.status, 0);
    assert!(best.energy < -46.667, "energy: {}", best.energy);
    let coords = Coords3D {
        points: best
            .coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let stereo = verify_stereo(&molecule, &coords);
    assert_eq!(stereo.n_violations(), 0);
    assert_eq!(stereo.n_unevaluable(), 0);
}

#[test]
fn quality_pool_retains_best_energy_and_converged_candidates_from_both_batches() {
    let candidate = |energy, status, marker| RdkitMmffConformer {
        coords: vec![[marker, 0.0, 0.0]],
        status,
        energy,
    };
    let compatibility = vec![
        candidate(-8.0, 0, 1.0),
        candidate(-7.0, 0, 2.0),
        candidate(-6.0, 1, 3.0),
    ];
    let alternative = vec![
        candidate(-10.0, 1, 4.0),
        candidate(-9.0, 1, 5.0),
        candidate(-5.0, 1, 6.0),
    ];
    let selected = select_quality_candidates(compatibility, alternative, 3);
    assert_eq!(selected.len(), 3);
    assert_eq!(
        selected
            .iter()
            .map(|candidate| candidate.energy)
            .collect::<Vec<_>>(),
        vec![-10.0, -8.0, -7.0]
    );
    assert_eq!(
        selected
            .iter()
            .filter(|candidate| candidate.status == 0)
            .count(),
        2
    );
}

#[test]
fn corrected_ensemble_preserves_single_seed_and_opt_in_boundary() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CS(=O)(=O)NCC").unwrap());
    let options = RdkitEmbedOptions::default();
    let single = etkdg_embed_molecule_corrected(&molecule, &options).unwrap();
    let corrected = etkdg_embed_multiple_confs_corrected(&molecule, 2, &options).unwrap();
    let compatibility = rdkit_embed_multiple_confs(&molecule, 2, &options).unwrap();
    assert_eq!(corrected.len(), 2);
    assert_eq!(compatibility.len(), 2);
    assert_eq!(corrected[0], single);
    assert!(
        corrected[0]
            .iter()
            .zip(&compatibility[0])
            .any(|(a, b)| a.iter().zip(b).any(|(x, y)| (x - y).abs() > 1e-8))
    );
    let optimized = etkdg_embed_multiple_confs_mmff_corrected(&molecule, 2, &options, 200).unwrap();
    assert_eq!(optimized.len(), 2);
    assert!(
        optimized
            .iter()
            .all(|candidate| candidate.energy.is_finite())
    );
    let selected =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 2, &options, 200).unwrap();
    let compatibility = rdkit_embed_multiple_confs_mmff(&molecule, 2, &options, 200).unwrap();
    let best = |ensemble: &[RdkitMmffConformer]| {
        ensemble
            .iter()
            .map(|candidate| candidate.energy)
            .fold(f64::INFINITY, f64::min)
    };
    assert_eq!(selected.len(), 2);
    assert!(best(&selected) <= best(&compatibility));
    assert!(best(&selected) <= best(&optimized));
}
