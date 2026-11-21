use super::*;
use crate::stereo_constraints::verify_stereo;
use chematic_core::{Coords3D, Point3};

#[test]
fn zero_seed_quality_profile_searches_new_batches() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCCCCO").unwrap());
    let options = RdkitEmbedOptions {
        random_seed: 0,
        ..Default::default()
    };
    let prepared = prepare_fragments(&molecule, &options).unwrap();
    assert!(!corrected_torsion_gradient_can_differ(&prepared));
    let mut batches = Vec::new();
    for offset in [0, 3, 6] {
        let coords = embed_prepared_multiple_confs(&prepared, 3, offset, &options, false).unwrap();
        batches.push(optimize_mmff_molecule(&molecule, coords, 0).unwrap());
    }
    assert_ne!(batches[0][0].coords, batches[1][0].coords);
    assert_ne!(batches[1][0].coords, batches[2][0].coords);
    let best_energy = batches
        .iter()
        .flatten()
        .map(|c| c.energy)
        .min_by(f64::total_cmp)
        .unwrap();
    let selected =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 3, &options, 0).unwrap();
    assert_eq!(selected.len(), 3);
    assert_eq!(selected[0].energy, best_energy);
    assert!(
        selected
            .iter()
            .all(|chosen| batches.iter().flatten().any(|c| chosen.coords == c.coords))
    );
}

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
    let fourth_coords = embed_prepared_multiple_confs(&prepared, 4, 12, &options, false).unwrap();
    let fourth = optimize_mmff_molecule(&molecule, fourth_coords, 200).unwrap();
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
    assert!(best(&selected) <= best(&first).min(best(&next)).min(best(&third)));
    assert!(selected.iter().all(|chosen| {
        first
            .iter()
            .chain(&next)
            .chain(&third)
            .chain(&fourth)
            .any(|original| chosen.coords == original.coords)
    }));
}

#[test]
fn torsion_search_recovers_lower_converged_minimum() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("COc1cc2nc(N3CCN(S(=O)(=O)c4cccs4)CC3)nc(N)c2cc1OC").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let baseline =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 20, &options, 200).unwrap();
    let improved =
        etkdg_embed_multiple_confs_mmff_torsion_search(&molecule, 20, &options, 200).unwrap();
    let best_converged = |rows: &[RdkitMmffConformer]| {
        rows.iter()
            .filter(|row| row.status == 0)
            .map(|row| row.energy)
            .min_by(f64::total_cmp)
            .unwrap()
    };
    assert_eq!(baseline.len(), improved.len());
    let before = best_converged(&baseline);
    let after = best_converged(&improved);
    assert!(
        after <= before + 1e-6,
        "torsion search: {before} -> {after}"
    );
    // ETKDG and BFGS can enter different basins across targets. The measured
    // macOS arm64 regression has a large margin; the portable contract is
    // finite, converged, stereo-checked non-regression.
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    assert!(after < before - 1.2, "torsion search: {before} -> {after}");
}

#[test]
fn torsion_search_pool_preserves_output_count_and_legacy_result() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCO").unwrap());
    let options = RdkitEmbedOptions::default();
    let legacy =
        etkdg_embed_multiple_confs_mmff_torsion_search(&molecule, 4, &options, 200).unwrap();
    let same =
        etkdg_embed_multiple_confs_mmff_torsion_search_with_pool(&molecule, 4, 4, &options, 200)
            .unwrap();
    assert_eq!(
        legacy
            .iter()
            .map(|row| (row.status, row.energy))
            .collect::<Vec<_>>(),
        same.iter()
            .map(|row| (row.status, row.energy))
            .collect::<Vec<_>>()
    );
    let pooled =
        etkdg_embed_multiple_confs_mmff_torsion_search_with_pool(&molecule, 4, 8, &options, 200)
            .unwrap();
    assert_eq!(pooled.len(), 4);
    assert!(
        pooled
            .iter()
            .any(|row| row.status == 0 && row.energy.is_finite())
    );
    assert!(
        etkdg_embed_multiple_confs_mmff_torsion_search_with_pool(&molecule, 4, 3, &options, 200,)
            .is_err()
    );
}

#[test]
fn quality_profile_extends_search_after_third_batch_improves() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse(
            "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
        )
        .unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let prepared = prepare_fragments(&molecule, &options).unwrap();
    assert!(!corrected_torsion_gradient_can_differ(&prepared));

    let first_three = rdkit_embed_multiple_confs_mmff(&molecule, 24, &options, 1000).unwrap();
    let before = lowest_finite_energy(&first_three).unwrap();
    let selected =
        etkdg_embed_multiple_confs_mmff_best_gradient(&molecule, 8, &options, 1000).unwrap();
    let after = lowest_finite_energy(&selected).unwrap();
    assert!(selected.iter().all(|candidate| candidate.status == 0));
    assert!(after < before - 0.5, "fourth batch: {before} -> {after}");
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

#[test]
fn corrected_checked_ensemble_verifies_the_returned_gradient() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("CC(C)([NH3+])CNc1nc(-c2ccncc2)cc2cnccc12").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let original = etkdg_embed_multiple_confs_mmff_corrected(&molecule, 20, &options, 200).unwrap();
    let checked =
        etkdg_embed_multiple_confs_mmff_corrected_checked(&molecule, 20, &options, 200).unwrap();
    assert_eq!(original.len(), checked.len());
    for (before, after) in original.iter().zip(&checked) {
        assert_eq!(before.coords, after.coords);
        assert_eq!(before.energy, after.energy);
        if after.status == 0 {
            let field =
                chematic_ff::rdkit_mmff::RdkitMmffField::new(&molecule, &after.coords, 100.0, true)
                    .unwrap();
            let flat = after.coords.as_flattened();
            let mut gradient = vec![0.0; flat.len()];
            field.gradient(flat, &mut gradient);
            assert!(gradient.iter().all(|v| v.is_finite() && v.abs() <= 0.01));
        }
    }
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    assert!(
        original
            .iter()
            .zip(&checked)
            .any(|(before, after)| before.status == 0 && after.status == 1)
    );
}

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[test]
fn corrected_checked_ensemble_preserves_macrocycle_double_bond_stereo() {
    // Exposed Prime development case: candidate 43 crosses an E double bond
    // during unconstrained MMFF optimization on the Linux reference runner.
    let smiles = "CO[C@H]1C[C@H](CC[C@@H](C)[C@H](O)[C@H](C)[C@H]2OC(=O)/C=C/C(C)=C/C[C@H](O)C[C@@H]3C=CC[C@@H](C[C@H](OC)[C@@H](C)[C@@H](O)C[C@@H](O)[C@H](C)[C@@H]([C@@H](C)[C@@H](O)[C@@H](C)CC[C@H]4C[C@H](OC)C[C@H](C)O4)OC(=O)/C=C/C(C)=C\\C[C@@H](O)C[C@@H]4C=CC[C@@H](C[C@H](OC)[C@@H](C)[C@@H](O)C[C@@H](O)[C@@H]2C)O4)O3)O[C@@H](C)C1";
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
    let options = RdkitEmbedOptions::default();
    let embedded = etkdg_embed_multiple_confs_corrected(&molecule, 44, &options).unwrap();
    let unchecked =
        etkdg_embed_multiple_confs_mmff_corrected(&molecule, 44, &options, 200).unwrap();
    let checked =
        etkdg_embed_multiple_confs_mmff_corrected_checked(&molecule, 44, &options, 200).unwrap();
    assert_eq!(embedded.len(), 44);
    assert_eq!(unchecked.len(), 44);
    assert_eq!(checked.len(), 44);
    assert!(stereo_safe_retry(&molecule, &embedded[43]));
    assert!(!stereo_safe_retry(&molecule, &unchecked[43].coords));
    assert_eq!(checked[43].coords, embedded[43]);
    assert_eq!(checked[43].status, 1);
    assert!(
        checked
            .iter()
            .all(|candidate| stereo_safe_retry(&molecule, &candidate.coords))
    );
    let field =
        chematic_ff::rdkit_mmff::RdkitMmffField::new(&molecule, &embedded[0], 100.0, true).unwrap();
    assert_eq!(
        checked[43].energy,
        field.energy(embedded[43].as_flattened())
    );
}

#[test]
fn corrected_repair_is_bounded_and_stereo_safe() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("CC(C)([NH3+])CNc1nc(-c2ccncc2)cc2cnccc12").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let checked =
        etkdg_embed_multiple_confs_mmff_corrected_checked(&molecule, 20, &options, 200).unwrap();
    let disabled =
        etkdg_embed_multiple_confs_mmff_corrected_repair(&molecule, 20, &options, 200, 0).unwrap();
    let repaired =
        etkdg_embed_multiple_confs_mmff_corrected_repair(&molecule, 20, &options, 200, 400)
            .unwrap();
    assert_eq!(checked.len(), repaired.len());
    for ((initial, no_retry), after) in checked.iter().zip(&disabled).zip(&repaired) {
        assert_eq!(initial.coords, no_retry.coords);
        assert_eq!(initial.energy, no_retry.energy);
        assert_eq!(initial.status, no_retry.status);
        if initial.coords != after.coords {
            assert_eq!(after.status, 0);
            assert!(after.energy <= initial.energy);
            assert!(stereo_safe_retry(&molecule, &after.coords));
        }
    }
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    assert!(
        checked
            .iter()
            .zip(&repaired)
            .any(|(a, b)| a.coords != b.coords)
    );
}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
#[test]
fn near_linear_mmff_consistent_gradient_matches_energy_difference() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("CC(C)([NH3+])CNc1nc(-c2ccncc2)cc2cnccc12").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let rows =
        etkdg_embed_multiple_confs_mmff_corrected_checked(&molecule, 20, &options, 200).unwrap();
    for index in [9, 16, 17] {
        let coords = &rows[index].coords;
        let field = chematic_ff::rdkit_mmff::RdkitMmffField::with_terms(
            &molecule,
            coords,
            chematic_ff::rdkit_mmff::MmffVariant::Mmff94,
            100.0,
            true,
            [false, false, false, false, true, false, false],
        )
        .unwrap();
        let mut flat = coords.as_flattened().to_vec();
        let mut gradient = vec![0.0; flat.len()];
        field.gradient_with_consistent_torsions(&flat, &mut gradient);
        for axis in 0..flat.len() {
            let step = 1e-10;
            flat[axis] += step;
            let plus = field.energy(&flat);
            flat[axis] -= 2.0 * step;
            let minus = field.energy(&flat);
            flat[axis] += step;
            let finite = (plus - minus) / (2.0 * step);
            assert!(gradient[axis].is_finite() && finite.is_finite());
            assert!(
                (gradient[axis] - finite).abs() < 0.02 * (1.0 + finite.abs()),
                "conformer {index} axis {axis}: analytic {} finite {}",
                gradient[axis],
                finite
            );
        }
    }
}

#[cfg(all(target_os = "macos", target_arch = "aarch64"))]
#[test]
fn near_linear_mmff_plane_gradient_matches_energy_difference() {
    let molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse("CC(C)([NH3+])CNc1nc(-c2ccncc2)cc2cnccc12").unwrap(),
    );
    let options = RdkitEmbedOptions::default();
    let rows =
        etkdg_embed_multiple_confs_mmff_corrected_checked(&molecule, 20, &options, 200).unwrap();
    for index in [9, 16, 17] {
        let coords = &rows[index].coords;
        let field = chematic_ff::rdkit_mmff::RdkitMmffField::with_terms(
            &molecule,
            coords,
            chematic_ff::rdkit_mmff::MmffVariant::Mmff94,
            100.0,
            true,
            [false, false, false, true, false, false, false],
        )
        .unwrap();
        let mut flat = coords.as_flattened().to_vec();
        let mut gradient = vec![0.0; flat.len()];
        field.gradient_with_consistent_planes(&flat, &mut gradient);
        for axis in 0..flat.len() {
            let step = 1e-8;
            flat[axis] += step;
            let plus = field.energy(&flat);
            flat[axis] -= 2.0 * step;
            let minus = field.energy(&flat);
            flat[axis] += step;
            let finite = (plus - minus) / (2.0 * step);
            assert!(gradient[axis].is_finite() && finite.is_finite());
            assert!(
                (gradient[axis] - finite).abs() < 0.02 * (1.0 + finite.abs()),
                "conformer {index} axis {axis}: analytic {} finite {}",
                gradient[axis],
                finite
            );
        }
    }
}
