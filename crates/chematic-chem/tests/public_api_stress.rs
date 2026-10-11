//! Broad public-API exercise over the checked-in special-chemistry corpus.
//!
//! These calls deliberately use integration tests: coverage reflects behavior
//! available to downstream crates instead of test-only implementation hooks.

use chematic_chem::*;
use chematic_smiles::parse;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn descriptor_and_standardization_surface_handles_special_chemistry() {
    let mut exercised = 0usize;
    for line in STRESS.lines() {
        if exercised == 1_024 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 120 {
            continue;
        }
        exercised += 1;

        let _ = molecular_weight(&mol);
        let _ = rdkit_molecular_weight(&mol);
        let _ = exact_mass(&mol);
        let _ = calc_mol_formula(&mol);
        let _ = heavy_atom_count(&mol);
        let _ = hba_count(&mol);
        let _ = hba_count_lipinski(&mol);
        let _ = rdkit_hba_count(&mol);
        let _ = hbd_count(&mol);
        let _ = tpsa(&mol);
        let _ = rdkit_tpsa(&mol);
        let _ = tpsa_per_atom(&mol);
        let _ = fsp3(&mol);
        let _ = fraction_rotatable_bonds(&mol);
        let _ = rotatable_bond_count(&mol);
        let _ = rotatable_bond_atom_pairs(&mol);
        let _ = ring_count(&mol);
        let _ = rdkit_num_rings(&mol);
        let _ = aromatic_ring_count(&mol);
        let _ = rdkit_aromatic_ring_count(&mol);
        let _ = ring_system_count(&mol);
        let _ = ring_bundle(&mol);
        let _ = num_aliphatic_rings(&mol);
        let _ = num_aromatic_heterocycles(&mol);
        let _ = num_aliphatic_heterocycles(&mol);
        let _ = num_saturated_rings(&mol);
        let _ = num_saturated_heterocycles(&mol);
        let _ = num_bridgehead_atoms(&mol);
        let _ = num_spiro_atoms(&mol);
        let _ = num_amide_bonds(&mol);
        let _ = num_ester_bonds(&mol);
        let _ = num_stereocenters(&mol);
        let _ = num_unspecified_stereocenters(&mol);
        let _ = potential_stereocenter_indices(&mol);
        let _ = implicit_hcount_per_atom(&mol);
        let _ = formal_charge_per_atom(&mol);
        let _ = formal_charge_sum(&mol);
        let _ = hybridization_per_atom(&mol);
        let _ = carbon_types(&mol);
        let _ = mqn(&mol);
        let _ = information_content(&mol);
        let _ = distance_descriptor_bundle(&mol);
        let _ = bcut2d(&mol);
        let _ = autocorr_2d(&mol);
        let _ = moran_autocorr(&mol);
        let _ = geary_autocorr(&mol);
        let _ = mde_carbon(&mol);
        let _ = usrcat(&mol);

        let _ = logp_and_mr(&mol);
        let _ = logp_crippen(&mol);
        let _ = logp_crippen_per_atom(&mol);
        let _ = molar_refractivity(&mol);
        let _ = mr_per_atom(&mol);
        let _ = qed(&mol);
        let _ = sa_score(&mol);
        let _ = esol_solubility(&mol);
        let _ = drug_score(&mol);
        let _ = xlogp3(&mol);
        let _ = xlogp3_per_atom(&mol);
        let _ = gasteiger_charges(&mol);
        let _ = mmff94_charges(&mol);
        let _ = mmff94_charges_bci(&mol);
        let _ = try_mmff94_charges_bci(&mol);
        let _ = estate_indices(&mol);
        let _ = estate_all(&mol);
        let _ = min_estate(&mol);
        let _ = max_estate(&mol);
        let _ = sum_estate(&mol);
        let _ = peoe_vsa(&mol);
        let _ = slogp_vsa(&mol);
        let _ = smr_vsa(&mol);
        let _ = estate_vsa(&mol);

        let _ = topology_bundle(&mol);
        let _ = topological_distance_matrix(&mol);
        let _ = graph_eccentricities(&mol);
        let _ = graph_diameter(&mol);
        let _ = graph_radius(&mol);
        let _ = wiener_index(&mol);
        let _ = balaban_j(&mol);
        let _ = try_balaban_j(&mol);
        let _ = ipc(&mol);
        let _ = try_ipc(&mol);
        let _ = bertz_ct(&mol);
        let _ = hall_kier_alpha(&mol);
        let _ = chi_all(&mol);
        let _ = kappa_all(&mol);
        let _ = labute_asa(&mol);
        let _ = labute_asa_per_atom(&mol);
        let _ = num_valence_electrons(&mol);
        let _ = randic_index(&mol);
        let _ = zagreb_index_m1(&mol);
        let _ = zagreb_index_m2(&mol);
        let _ = gravitational_index(&mol);
        let _ = eccentric_connectivity_index(&mol);
        let _ = petitjean_index(&mol);
        let _ = vabc(&mol);

        let _ = lipinski_passes(&mol);
        let _ = veber_passes(&mol);
        let _ = ro3_passes(&mol);
        let _ = ghose_passes(&mol);
        let _ = egan_passes(&mol);
        let _ = reos_passes(&mol);
        let _ = mcf_passes(&mol);
        let _ = pfizer_3_75_passes(&mol);
        let _ = lead_like_passes(&mol);
        let _ = pains_matches(&mol);
        let _ = brenk_matches(&mol);
        let _ = identify_functional_groups(&mol);
        let _ = detect_named_functional_groups(&mol);
        let _ = genotox_reactivity(&mol);
        let _ = detect_atropisomers(&mol);
        let _ = assign_atropisomer_chirality(&mol);
        let _ = assign_cip(&mol);
        let _ = assign_ez_bonds(&mol);
        let _ = ez_completeness(&mol);

        if exercised <= 512 {
            let _ = standardize(&mol, &StandardizeOptions::default());
            let _ = neutralize_charges(&mol);
            let _ = normalize_groups(&mol);
            let _ = reionize(&mol);
            let _ = uncharge(&mol);
            let _ = normalize_zwitterion(&mol);
            let _ = remove_isotopes(&mol);
            let _ = remove_stereo(&mol);
            let _ = largest_fragment(&mol);
            let _ = fragment_parent(&mol);
            let _ = charge_parent(&mol);
            let _ = isotope_parent(&mol);
            let _ = stereo_parent(&mol);
            let limits = TautomerLimits::default();
            let _ = super_parent(&mol, &limits);
            let _ = murcko_scaffold(&mol);
            let _ = generic_murcko_scaffold(&mol);
            let _ = scaffold_network(&mol);
            let _ = schuffenhauer_parents(&mol);
            let _ = brics_bonds(&mol);
            let _ = brics_fragments(&mol);
            let _ = recap_breakable_bond_count(&mol);
            let _ = recap_fragment(&mol);
            let _ = canonical_tautomer(&mol);
            let _ = enumerate_tautomers(&mol);
            let _ = tautomer_parent(&mol, &limits);
            let _ = pka_acid(&mol);
            let _ = pka_base(&mol);
            let _ = pka_both(&mol);
            let _ = predict_pka(&mol);
            let _ = admet_profile(&mol);
            let _ = boiled_egg(&mol);
        }
    }
    assert_eq!(exercised, 1_024);
}

#[test]
fn virtual_hydrogen_is_terminal_during_potential_stereocenter_ranking() {
    let phosphine = parse("[H]P([H])[H]").unwrap();
    assert!(potential_stereocenter_indices(&phosphine).is_empty());
}
