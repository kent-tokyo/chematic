//! Public-surface robustness over the checked-in special-chemistry corpus.
//!
//! The parity tests pin individual outputs. This test has a different job: it
//! sends the same difficult structures through every bounded RDKit-compatible
//! writer/inspection surface and asserts that none of those public entry
//! points panic. Unsupported chemistry remains an explicit `Err`.

use chematic_smiles::{
    RdkitHashFunction, RdkitSmilesParams, parse, rdkit_2d_coords, rdkit_atom_stereocenter_counts,
    rdkit_canonical_smiles, rdkit_canonical_tautomer, rdkit_chiral_centers, rdkit_cleanup_edits,
    rdkit_crippen_contribs, rdkit_crippen_contribs_no_hs, rdkit_crippen_logp_mr, rdkit_cx_smarts,
    rdkit_cx_smiles, rdkit_detect_chemistry_problems, rdkit_detect_chemistry_problems_sanitized,
    rdkit_distance_matrix, rdkit_distance_matrix_3d, rdkit_enumerate_tautomers,
    rdkit_fragment_smiles, rdkit_hybridizations, rdkit_hydrogen_suppressed,
    rdkit_hydrogen_suppressed_with_map, rdkit_legacy_stereo, rdkit_mol_block_2d, rdkit_mol_hash,
    rdkit_mol_view, rdkit_murcko_scaffold, rdkit_num_atoms, rdkit_num_pi_electrons,
    rdkit_parsed_molecule, rdkit_pdb_block, rdkit_random_smiles, rdkit_sanitized_model,
    rdkit_smarts, rdkit_smiles, rdkit_stereoisomer_count, rdkit_stereoisomer_smiles,
    rdkit_tautomer_score,
};

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");
const BROWSER_10K: &str =
    include_str!("../../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi");

#[test]
fn rdkit_compatible_public_surface_is_panic_free_on_stress_corpus() {
    let mut parsed = 0usize;
    for (row, line) in STRESS.lines().enumerate() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        parsed += 1;

        let _ = rdkit_canonical_smiles(&mol);
        let _ = rdkit_atom_stereocenter_counts(&mol);
        let _ = rdkit_cleanup_edits(&mol);
        let _ = rdkit_legacy_stereo(&mol);
        let _ = rdkit_smarts(&mol, true, None);
        let _ = rdkit_cx_smarts(&mol);
        let _ = rdkit_cx_smiles(&mol, &RdkitSmilesParams::default());
        let _ = rdkit_num_atoms(&mol);
        let _ = rdkit_distance_matrix(&mol, false, false);
        let _ = rdkit_distance_matrix(&mol, true, true);
        let _ = rdkit_detect_chemistry_problems(&mol);
        let _ = rdkit_detect_chemistry_problems_sanitized(&mol);
        let _ = rdkit_murcko_scaffold(&mol);
        let _ = rdkit_stereoisomer_count(&mol);
        let _ = rdkit_chiral_centers(&mol, true);
        let _ = rdkit_chiral_centers(&mol, false);
        let _ = rdkit_hydrogen_suppressed(&mol);
        let _ = rdkit_hydrogen_suppressed_with_map(&mol);
        let _ = rdkit_hybridizations(&mol);
        let _ = rdkit_num_pi_electrons(&mol);
        let _ = rdkit_crippen_contribs(&mol);
        let _ = rdkit_crippen_contribs_no_hs(&mol);
        let _ = rdkit_crippen_logp_mr(&mol);
        let _ = rdkit_parsed_molecule(&mol);
        let _ = rdkit_sanitized_model(&mol);
        let _ = rdkit_mol_view(&mol);

        if row < 256 {
            let params = RdkitSmilesParams {
                isomeric: row % 2 == 0,
                kekule: row % 3 == 0,
                canonical: row % 5 != 0,
                all_bonds_explicit: row % 7 == 0,
                all_hs_explicit: row % 11 == 0,
                rooted_at_atom: (mol.atom_count() > 0).then_some(row % mol.atom_count().max(1)),
            };
            let _ = rdkit_smiles(&mol, &params);
            let _ = rdkit_random_smiles(&mol, 2, 42 + row as u32, &params);
            let _ = rdkit_pdb_block(&mol, None);
            let _ = rdkit_2d_coords(&mol);
            let _ = rdkit_mol_block_2d(&mol);

            if let Ok(n) = rdkit_num_atoms(&mol) {
                let coords = vec![[0.0, 0.0, 0.0]; n];
                let _ = rdkit_distance_matrix_3d(&mol, &coords, row % 2 == 0);
            }
            if mol.atom_count() > 0 {
                let atoms: Vec<usize> = (0..mol.atom_count().min(4)).collect();
                let _ = rdkit_fragment_smiles(&mol, &atoms, None, &params);
            }
        }

        if row < 128 {
            let _ = rdkit_stereoisomer_smiles(&mol, 64);
            let _ = rdkit_canonical_tautomer(&mol);
            let _ = rdkit_enumerate_tautomers(&mol);
            let _ = rdkit_tautomer_score(&mol);
            for hash in RdkitHashFunction::ALL {
                let _ = rdkit_mol_hash(&mol, hash, row % 2 == 0);
            }
        }
    }
    assert!(parsed > 3_500, "stress corpus unexpectedly stopped parsing");
}

#[test]
fn canonical_writer_handles_the_full_browser_corpus_in_multiple_profiles() {
    let mut parsed = 0usize;
    for (row, line) in BROWSER_10K.lines().enumerate() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        parsed += 1;

        let _ = rdkit_canonical_smiles(&mol);
        let params = RdkitSmilesParams {
            isomeric: row % 2 == 0,
            kekule: row % 3 == 0,
            canonical: row % 5 != 0,
            all_bonds_explicit: row % 7 == 0,
            all_hs_explicit: row % 11 == 0,
            rooted_at_atom: (mol.atom_count() > 0).then_some(row % mol.atom_count().max(1)),
        };
        let _ = rdkit_smiles(&mol, &params);
    }
    assert!(
        parsed >= 9_900,
        "browser corpus unexpectedly stopped parsing"
    );
}
