//! Panic-safety matrix for the public 3D APIs on unusual valid structures.

use chematic_3d::etkdg_knowledge::{
    TorsionKnowledgeConfig, TorsionOptimizationConfig, build_torsion_knowledge,
    evaluate_torsion_energy, macrocycle_14_bound_adjustments, optimize_torsions,
};
use chematic_3d::minimize::ForceFieldPolicy;
use chematic_3d::rdkit_embed::{RdkitEmbedOptions, rdkit_bounds_matrix, rdkit_embed_molecule};
use chematic_3d::{
    PipelineV2Config, align_coords, asphericity, assign_stereo_from_3d, autocorr_3d,
    build_constraints, calc_mol_sasa, eccentricity, generate_coords,
    generate_coords_connectivity_ordered, generate_coords_etkdg, getaway_descriptors, npr1, npr2,
    pharmacophore_fp_3d, plane_of_best_fit, pmi, pmi1, pmi2, pmi3, radius_of_gyration,
    rdf_descriptors, rmsd_no_align, sasa_descriptor, sasa_per_atom, sasa_per_element,
    satisfy_constraints, usr_descriptors, usr_similarity, whim_descriptors, whim_getaway_combined,
};
use chematic_smiles::parse;
use std::collections::HashSet;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn public_3d_surface_is_panic_free_on_diverse_stress_structures() {
    let mut exercised = 0usize;
    for line in STRESS.lines() {
        if exercised == 256 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 80 {
            continue;
        }
        exercised += 1;

        let coords = generate_coords(&mol);
        let ordered = generate_coords_connectivity_ordered(&mol);
        let etkdg = generate_coords_etkdg(&mol);
        assert_eq!(coords.atom_count(), mol.atom_count());
        assert_eq!(ordered.atom_count(), mol.atom_count());
        assert_eq!(etkdg.atom_count(), mol.atom_count());

        let constraints = build_constraints(&mol);
        let _ = satisfy_constraints(&coords, &mol, &constraints, 2);
        let raw: Vec<[f64; 3]> = coords.points.iter().map(|p| [p.x, p.y, p.z]).collect();
        let ordered_raw: Vec<[f64; 3]> = ordered.points.iter().map(|p| [p.x, p.y, p.z]).collect();
        let _ = align_coords(&raw, &ordered_raw);
        let _ = rmsd_no_align(&raw, &ordered_raw);
        let _ = assign_stereo_from_3d(&mol, &coords);
        let fp = pharmacophore_fp_3d(&mol, &coords);
        assert_eq!(chematic_3d::tanimoto_pharmacophore_3d(&fp, &fp), 1.0);

        let _ = pmi(&mol, &coords);
        let _ = pmi1(&mol, &coords);
        let _ = pmi2(&mol, &coords);
        let _ = pmi3(&mol, &coords);
        let _ = npr1(&mol, &coords);
        let _ = npr2(&mol, &coords);
        let _ = radius_of_gyration(&mol, &coords);
        let _ = asphericity(&mol, &coords);
        let _ = eccentricity(&mol, &coords);
        let _ = plane_of_best_fit(&mol, &coords);
        let _ = whim_descriptors(&mol, &coords);
        let _ = getaway_descriptors(&mol, &coords);
        let _ = whim_getaway_combined(&mol, &coords);
        let _ = rdf_descriptors(&mol, &coords);
        let _ = autocorr_3d(&mol, &coords);
        let _ = calc_mol_sasa(&mol, &coords);
        let _ = sasa_per_atom(&mol, &coords, 1.4, 32);
        let _ = sasa_descriptor(&mol, &coords);
        let _ = sasa_per_element(&mol, &coords);

        let usr = usr_descriptors(&raw);
        assert!((usr_similarity(&usr, &usr) - 1.0).abs() < 1e-12);

        if exercised <= 128 {
            let explicit_h = chematic_chem::add_hydrogens(&mol);
            for set_15 in [false, true] {
                for smooth in [false, true] {
                    let _ = rdkit_bounds_matrix(
                        &explicit_h,
                        set_15,
                        smooth,
                        exercised.is_multiple_of(2),
                    );
                }
            }
            if exercised <= 64 {
                for options in [
                    RdkitEmbedOptions::default(),
                    RdkitEmbedOptions {
                        use_small_ring_torsions: true,
                        use_macrocycle_torsions: true,
                        use_macrocycle14config: true,
                        ..RdkitEmbedOptions::default()
                    },
                ] {
                    let _ = rdkit_embed_molecule(&explicit_h, &options);
                }
            }
        }

        if exercised <= 12 {
            let policy = match exercised % 3 {
                0 => ForceFieldPolicy::UffOnly,
                1 => ForceFieldPolicy::Mmff94BondAngleStrict,
                _ => ForceFieldPolicy::Mmff94WithUffFallback,
            };
            let mut config = PipelineV2Config::minimal(policy);
            config.total_timeout_ms = Some(2_000);
            config.include_legacy_torsion_heuristic = exercised.is_multiple_of(2);
            config.gate_mmff94_torsion_oop = exercised % 3 == 1;
            config.gate_mmff94_stretch_bend = exercised % 4 == 1;
            let _ = chematic_3d::embed_pipeline_v2(&mol, &config);
        }
    }
    assert_eq!(exercised, 256);
}

#[test]
fn torsion_preferences_cover_every_path_in_the_special_chemistry_corpus() {
    let mut molecules = 0usize;
    let mut paths = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 160 {
            continue;
        }
        molecules += 1;

        let ring_set = chematic_perception::find_sssr(&mol);
        let ring_bonds: HashSet<(u32, u32)> = ring_set
            .rings()
            .iter()
            .flat_map(|ring| {
                (0..ring.len()).flat_map(move |index| {
                    let a = ring[index].0;
                    let b = ring[(index + 1) % ring.len()].0;
                    [(a, b), (b, a)]
                })
            })
            .collect();
        let _ = chematic_3d::etkdg_knowledge::build_smarts_torsion_map(&mol, &ring_bonds);

        for (b, _) in mol.atoms() {
            for (c, _) in mol.neighbors(b) {
                for (a, _) in mol.neighbors(b).filter(|(a, _)| *a != c) {
                    for (d, _) in mol.neighbors(c).filter(|(d, _)| *d != b && *d != a) {
                        paths += 1;
                        if let Some(preference) =
                            chematic_3d::etkdg_knowledge::get_torsion_preference(&mol, a, b, c, d)
                        {
                            assert!(preference.angle_deg.is_finite());
                            assert!(preference.penalty_per_degree.is_finite());
                            let score = chematic_3d::etkdg_knowledge::score_torsion(
                                preference.angle_deg + 360.0,
                                &preference,
                            );
                            assert!(score.abs() < 1e-12);
                        }
                    }
                }
            }
        }
    }
    assert!(molecules > 3_500);
    assert!(paths > 20_000);
}

#[test]
fn rdkit_bounds_cover_the_full_special_chemistry_corpus() {
    let mut molecules = 0usize;
    let mut successful = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 160 {
            continue;
        }
        molecules += 1;
        let explicit_h = chematic_chem::add_hydrogens(&mol);
        for set_15 in [false, true] {
            for scale_vdw in [false, true] {
                if let Ok(matrix) = rdkit_bounds_matrix(&explicit_h, set_15, true, scale_vdw) {
                    assert_eq!(matrix.len(), explicit_h.atom_count());
                    successful += 1;
                }
            }
        }
    }
    assert!(molecules > 3_500);
    assert!(successful > 10_000);
}

#[test]
fn v2_torsion_knowledge_covers_all_rule_tiers_on_special_chemistry() {
    let all_tiers = TorsionKnowledgeConfig {
        use_exp_torsions: true,
        use_small_ring_torsions: true,
        use_macrocycle_torsions: true,
        use_macrocycle_14_bounds: true,
        include_legacy_heuristic: true,
    };
    let mut molecules = 0usize;
    let mut potentials = 0usize;
    let mut evaluated = 0usize;
    let mut adjustments = 0usize;
    let mut optimized = 0usize;

    for line in STRESS.lines() {
        if molecules == 512 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() < 4 || mol.atom_count() > 100 {
            continue;
        }
        molecules += 1;

        // Default remains a true no-op; each enabled tier is then exercised
        // together so priority and same-tier conflict handling run on the
        // same unusual structures used by the cross-engine comparison.
        let disabled = build_torsion_knowledge(&mol, &TorsionKnowledgeConfig::default());
        assert!(disabled.potentials.is_empty());
        let report = build_torsion_knowledge(&mol, &all_tiers);
        potentials += report.potentials.len();
        adjustments += macrocycle_14_bound_adjustments(&mol, &all_tiers)
            .expect("valid stress molecule topology")
            .len();

        let coords = generate_coords_connectivity_ordered(&mol);
        let energy = evaluate_torsion_energy(&mol, &coords, &report.potentials)
            .expect("coordinate count follows molecule");
        assert!(energy.total_energy.is_finite());
        evaluated += report.potentials.len();

        if optimized < 12 && !report.potentials.is_empty() {
            let config = TorsionOptimizationConfig {
                max_iterations: 2,
                max_line_search_steps: 2,
                ..TorsionOptimizationConfig::default()
            };
            let (next, optimization) =
                optimize_torsions(&mol, &coords, &report.potentials, &config)
                    .expect("valid torsion optimization input");
            assert_eq!(next.atom_count(), mol.atom_count());
            assert!(optimization.energy_before.is_finite());
            assert!(optimization.energy_after.is_finite());
            optimized += 1;
        }
    }

    assert_eq!(molecules, 512);
    assert!(potentials > 100, "expected broad torsion-rule coverage");
    assert_eq!(evaluated, potentials);
    assert!(adjustments > 0, "expected at least one macrocycle 1-4 case");
    assert!(optimized > 0);
}
