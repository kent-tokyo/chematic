use chematic_3d::rdkit_embed::{
    BondLengthConstraint, InternalCoordinateConstraints, RdkitEmbedOptions,
    rdkit_embed_molecule_with_constraints,
};

#[test]
fn disconnected_fragments_keep_their_own_bond_length_constraints() {
    let molecule = chematic_chem::add_hydrogens(&chematic_smiles::parse("CC.CC").unwrap());
    let constraints = InternalCoordinateConstraints {
        bond_lengths: vec![
            BondLengthConstraint::around(0, 1, 1.5, 0.2),
            BondLengthConstraint::around(2, 3, 1.5, 0.2),
        ],
        ..Default::default()
    };
    let coords = rdkit_embed_molecule_with_constraints(
        &molecule,
        &RdkitEmbedOptions::default(),
        &constraints,
    )
    .unwrap();
    assert_eq!(coords.len(), molecule.atom_count());
    for (a, b) in [(0, 1), (2, 3)] {
        let distance = (0..3)
            .map(|i| (coords[a][i] - coords[b][i]).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!((1.3..=1.7).contains(&distance), "{a}-{b}: {distance}");
    }
}

#[test]
fn nonfinite_torsion_bounds_fail_before_embedding() {
    use chematic_3d::rdkit_embed::{TorsionAngleConstraint, rdkit_bounds_matrix_with_constraints};
    let molecule = chematic_smiles::parse("CCCC").unwrap();
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let constraints = InternalCoordinateConstraints {
            torsions: vec![TorsionAngleConstraint {
                atom1: 0,
                atom2: 1,
                atom3: 2,
                atom4: 3,
                lower: value,
                upper: 1.,
            }],
            ..Default::default()
        };
        let error =
            rdkit_bounds_matrix_with_constraints(&molecule, &constraints, true, true, false)
                .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("torsion constraint contains a non-finite value")
        );
    }
}
