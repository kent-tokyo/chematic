//! Boundary tests for RDKit-compatible molecule embedding.

use super::*;

#[test]
fn embedding_rejects_empty_graphs_and_unsupported_modes_with_context() {
    let empty = chematic_core::MoleculeBuilder::new().build();
    let error = rdkit_embed_molecule(&empty, &RdkitEmbedOptions::default()).unwrap_err();
    assert_eq!(
        error,
        RdkitEmbedError::RdkitException("molecule has no atoms".into())
    );
    assert!(error.to_string().contains("molecule has no atoms"));
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let options = RdkitEmbedOptions {
        use_basic_knowledge: false,
        ..RdkitEmbedOptions::default()
    };
    let error = rdkit_embed_molecule(&molecule, &options).unwrap_err();
    assert!(matches!(error, RdkitEmbedError::Unsupported(_)));
    assert!(
        error
            .to_string()
            .contains("plain ETDG (no basic knowledge)")
    );
    assert_eq!(molecule.atom_count(), 9);
}

#[test]
fn embedding_error_variants_keep_actionable_context() {
    assert_eq!(
        RdkitEmbedError::Failed.to_string(),
        "RDKit embedding failed (EmbedMolecule returns -1)"
    );
    assert_eq!(
        RdkitEmbedError::from(BoundsInvariant("bad bounds")).to_string(),
        "RDKit embedding: RDKit raises: bad bounds"
    );
    assert!(
        RdkitEmbedError::from(EmbedError::BadDirection)
            .to_string()
            .contains("BadDirection")
    );
    assert!(
        RdkitEmbedError::from(InternalCoordinateConstraintError::TriangleSmoothingFailed)
            .to_string()
            .contains("invalid internal-coordinate constraint")
    );
    assert_eq!(
        RdkitEmbedMmffError::ForceField("missing parameters".into()).to_string(),
        "RDKit-compatible MMFF94: missing parameters"
    );
    assert_eq!(
        RdkitEmbedMmffError::from(RdkitEmbedError::Failed).to_string(),
        RdkitEmbedError::Failed.to_string()
    );
}

#[test]
fn ensemble_boundaries_do_not_construct_unrequested_conformers() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
    let options = RdkitEmbedOptions::default();
    assert!(
        rdkit_embed_multiple_confs(&molecule, 0, &options)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff(&molecule, 0, &options, 50)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_diagnostics(&molecule, 0, &options, 50)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_adaptive(&molecule, 0, &options, 5, 5, 1)
            .unwrap()
            .is_empty()
    );
    assert!(
        rdkit_embed_multiple_confs_mmff_rescue(&molecule, 0, &options, 5, 5, 1)
            .unwrap()
            .is_empty()
    );
    let error =
        rdkit_embed_multiple_confs_mmff_gradient_rescue(&molecule, 0, &options, 5, 5, 1, f64::NAN)
            .unwrap_err();
    assert!(error.to_string().contains("finite and non-negative"));
    let error =
        rdkit_embed_multiple_confs_mmff_gradient_rescue(&molecule, 0, &options, 5, 5, 1, -0.1);
    assert!(error.unwrap_err().to_string().contains("non-negative"));

    let negative_seed = RdkitEmbedOptions {
        random_seed: -1,
        ..options.clone()
    };
    assert!(
        rdkit_embed_molecule(&molecule, &negative_seed)
            .unwrap_err()
            .to_string()
            .contains("randomSeed")
    );
}

#[test]
fn plain_distance_geometry_and_stereo_discovery_are_explicit() {
    let molecule =
        chematic_chem::add_hydrogens(&chematic_smiles::parse("C[C@H](F)Cl.C=C=C").unwrap());
    let view = chematic_smiles::rdkit_mol_view(&molecule).unwrap();
    let (chiral, tetrahedral) = find_chiral_sets(&view);
    assert!(!chiral.is_empty());
    assert!(tetrahedral.is_empty());
    let (ends, stereo) = find_double_bonds(&view);
    assert!(!ends.is_empty());
    assert!(stereo.is_empty());

    let plain = RdkitEmbedOptions {
        use_exp_torsion_angle_prefs: false,
        use_basic_knowledge: false,
        ..RdkitEmbedOptions::default()
    };
    let methane = chematic_chem::add_hydrogens(&chematic_smiles::parse("C").unwrap());
    assert_eq!(rdkit_embed_molecule(&methane, &plain).unwrap().len(), 5);
}

#[test]
fn embedding_covers_stereo_rings_fragments_and_constraints() {
    let cases = [
        "N[C@@H](C)C(=O)O",
        "F/C=C(/Cl)Br",
        "C1CC2CCC1C2",
        "C1CCCC/C=C/CCC1",
        "CC(=O)NC",
        "CCO.[Na+]",
    ];
    for text in cases {
        let parsed = chematic_smiles::parse(text).unwrap_or_else(|error| panic!("{text}: {error}"));
        let molecule = chematic_chem::add_hydrogens(&parsed);
        let options = RdkitEmbedOptions {
            random_seed: 17,
            max_iterations: 2,
            ..RdkitEmbedOptions::default()
        };
        if let Ok(coords) = rdkit_embed_molecule(&molecule, &options) {
            assert_eq!(coords.len(), molecule.atom_count(), "{text}");
            assert!(coords.iter().flatten().all(|value| value.is_finite()));
        }
        let bounds = rdkit_bounds_matrix(&molecule, true, true, true)
            .unwrap_or_else(|error| panic!("bounds {text}: {error}"));
        assert_eq!(bounds.len(), molecule.atom_count());
        assert!(bounds.iter().all(|row| row.len() == molecule.atom_count()));
    }

    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCCO").unwrap());
    let constraints = InternalCoordinateConstraints {
        bond_lengths: vec![super::super::constraints::BondLengthConstraint {
            atom1: 0,
            atom2: 1,
            lower: 1.45,
            upper: 1.60,
        }],
        ..InternalCoordinateConstraints::default()
    };
    let matrix =
        rdkit_bounds_matrix_with_constraints(&molecule, &constraints, true, true, true).unwrap();
    assert!((1.45..=1.60).contains(&matrix[0][1]));
}

#[test]
fn quality_candidate_selection_is_finite_stable_and_bounded() {
    let candidate = |energy: f64, status: i32, marker: f64| RdkitMmffConformer {
        coords: vec![[marker, 0.0, 0.0]],
        status,
        energy,
    };
    assert!(select_quality_candidates(Vec::new(), Vec::new(), 0).is_empty());

    let non_finite = select_quality_candidates(
        vec![
            candidate(f64::NAN, 1, 1.0),
            candidate(f64::INFINITY, 0, 2.0),
        ],
        Vec::new(),
        1,
    );
    assert_eq!(non_finite.len(), 1);
    assert_eq!(non_finite[0].coords[0][0], 1.0);
    assert!(lowest_finite_energy(&non_finite).is_none());

    let selected = select_quality_candidates(
        vec![candidate(4.0, 1, 4.0), candidate(2.0, 1, 2.0)],
        vec![candidate(3.0, 0, 3.0), candidate(f64::NAN, 0, 9.0)],
        3,
    );
    assert_eq!(selected.len(), 3);
    assert_eq!(selected[0].energy, 2.0);
    assert_eq!(selected[1].energy, 3.0);
    assert_eq!(selected[2].energy, 4.0);
    assert_eq!(lowest_finite_energy(&selected), Some(2.0));
}
