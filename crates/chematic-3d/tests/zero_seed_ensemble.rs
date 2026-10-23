//! RDKit #9738: seed zero must not collapse an ensemble into repeated attempts.
use chematic_3d::rdkit_embed::{
    RdkitEmbedOptions, etkdg_embed_multiple_confs_corrected, rdkit_embed_molecule,
    rdkit_embed_multiple_confs, rdkit_embed_multiple_confs_mmff,
};

fn options(seed: i32) -> RdkitEmbedOptions {
    RdkitEmbedOptions {
        random_seed: seed,
        ..Default::default()
    }
}

fn distance_change(a: &[[f64; 3]], b: &[[f64; 3]], heavy_atoms: usize) -> f64 {
    let distance = |p: [f64; 3], q: [f64; 3]| {
        p.iter()
            .zip(q)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    (0..heavy_atoms)
        .flat_map(|i| (0..i).map(move |j| (i, j)))
        .map(|(i, j)| (distance(a[i], a[j]) - distance(b[i], b[j])).abs())
        .fold(0.0, f64::max)
}

#[test]
fn zero_seed_ensemble_is_reproducible_nonrigid_and_prefix_stable() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCCCCO").unwrap());
    let conformers = rdkit_embed_multiple_confs(&molecule, 5, &options(0)).unwrap();
    assert_eq!(conformers.len(), 5);
    assert!(distance_change(&conformers[0], &conformers[1], 7) > 0.1);
    assert_eq!(
        conformers,
        rdkit_embed_multiple_confs(&molecule, 5, &options(0)).unwrap()
    );
    assert_eq!(
        conformers[..3],
        rdkit_embed_multiple_confs(&molecule, 3, &options(0)).unwrap()
    );
    for (coords, seed) in conformers.iter().zip([0, 2, 3, 4, 5]) {
        assert_eq!(
            *coords,
            rdkit_embed_molecule(&molecule, &options(seed)).unwrap()
        );
        assert_eq!(coords.len(), molecule.atom_count());
        assert!(coords.iter().flatten().all(|value| value.is_finite()));
    }
}

#[test]
fn positive_seeds_keep_the_existing_multiplicative_schedule() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCCCCO").unwrap());
    for seed in [1, 42] {
        let conformers = rdkit_embed_multiple_confs(&molecule, 3, &options(seed)).unwrap();
        assert_eq!(conformers.len(), 3);
        for (index, coords) in conformers.iter().enumerate() {
            assert_eq!(
                *coords,
                rdkit_embed_molecule(&molecule, &options(seed * (index as i32 + 1))).unwrap()
            );
        }
    }
}

#[test]
fn corrected_and_mmff_ensembles_inherit_zero_seed_variation() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCCCCO").unwrap());
    let corrected = etkdg_embed_multiple_confs_corrected(&molecule, 3, &options(0)).unwrap();
    assert_eq!(corrected.len(), 3);
    assert!(distance_change(&corrected[0], &corrected[1], 7) > 0.1);
    let optimized = rdkit_embed_multiple_confs_mmff(&molecule, 3, &options(0), 0).unwrap();
    assert_eq!(optimized.len(), 3);
    assert!(distance_change(&optimized[0].coords, &optimized[1].coords, 7) > 0.1);
}

#[test]
fn thread_count_does_not_change_zero_seed_candidate_order() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCCCCO").unwrap());
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| rdkit_embed_multiple_confs(&molecule, 3, &options(0)).unwrap())
    };
    assert_eq!(run(1), run(2));
}

#[test]
fn impossible_seed_ranges_fail_before_constructing_ensembles() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    assert!(rdkit_embed_multiple_confs(&molecule, usize::MAX, &options(0)).is_err());
    assert!(rdkit_embed_multiple_confs(&molecule, 2, &options(i32::MAX)).is_err());
    assert!(rdkit_embed_multiple_confs(&molecule, 1, &options(-1)).is_err());
    assert!(
        rdkit_embed_multiple_confs(&molecule, 0, &options(0))
            .unwrap()
            .is_empty()
    );
}
