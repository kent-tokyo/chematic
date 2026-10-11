//! Panic-safety matrix for the public 3D APIs on unusual valid structures.

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
