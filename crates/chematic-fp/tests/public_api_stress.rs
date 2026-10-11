//! Broad fingerprint coverage over the checked-in special-chemistry corpus.

use chematic_fp::*;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn fingerprint_surface_handles_special_chemistry() {
    let mut exercised = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = chematic_smiles::parse(smiles) else {
            continue;
        };
        if mol.atom_count() == 0 || mol.atom_count() > 120 {
            continue;
        }
        exercised += 1;

        let _ = ecfp4(&mol);
        let _ = ecfp6(&mol);
        let _ = fcfp4(&mol);
        let _ = fcfp6(&mol);
        let _ = maccs(&mol);
        let _ = pattern_fp(&mol);
        let _ = layered_fp(&mol);
        let _ = layered_fp_by_layer(&mol);
        let _ = rdkit_pattern_fp(&mol);
        let _ = rdkit_layered_fp(&mol);
        let _ = rdkit_rdk_fp(&mol);
        let _ = rdkit_torsion_fp(&mol);
        let _ = rdkit_torsion_counts(&mol, 2_048);
        let _ = rdkit_torsion_sparse_counts(&mol);
        let _ = rdkit_atom_pair_fp(&mol);
        let _ = rdkit_atom_pair_counts(&mol, 2_048);
        let _ = rdkit_atom_pair_sparse_counts(&mol);
        let _ = rdkit_morgan_bitvec(&mol, &RdkitMorganConfig::default());
        let _ = rdkit_morgan_count_simulation(&mol, &RdkitMorganConfig::default());
        let _ = pharmacophore_fp_2d(&mol);
        let _ = pharmacophore_feature_counts(&mol);
        let _ = atom_pair_fp(&mol);
        let _ = torsion_fp(&mol);
        let _ = topo_path(&mol, &TopoPathConfig::default());
        let _ = rdkit_path_fp(&mol);

        if exercised <= 1_024 {
            let _ = avalon_fp(&mol);
            let _ = rdkit_avalon_fp(&mol, 2_048);
            let _ = erg(&mol);
            let _ = erg_extended(&mol);
            let _ = erg_vec(&mol);
            let _ = hdf_default(&mol);
            let _ = map4_default(&mol);
            let _ = mhfp(&mol);
            let _ = mhfp_128(&mol);
        }
    }
    assert!(
        exercised >= 4_000,
        "stress corpus unexpectedly shrank: {exercised}"
    );
}
