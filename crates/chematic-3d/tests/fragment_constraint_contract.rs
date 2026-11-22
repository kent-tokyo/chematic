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
