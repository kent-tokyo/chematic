use chematic_mol::{MolMetadata, parse_mol_with_coords, write_sdf_record};

#[test]
fn sdf_with_supplied_coordinates_does_not_invent_unspecified_azo_stereo() {
    let molecule = chematic_smiles::parse("CN=NC").unwrap();
    let coords = [(-1.0, 1.0), (0.0, 0.0), (1.0, 0.0), (2.0, -1.0)];
    let text = write_sdf_record(
        &molecule,
        &MolMetadata::default(),
        &coords,
        &Default::default(),
    );
    assert!(text.contains("  2  3  2  3\n"), "{text}");
    let (_, _, observed_coords) = parse_mol_with_coords(&text).unwrap();
    assert_eq!(observed_coords, coords);
}

#[test]
fn sdf_keeps_declared_double_bond_stereo_and_supplied_coordinates() {
    let molecule = chematic_smiles::parse("C/C=C/C").unwrap();
    let coords = [(-1.0, 1.0), (0.0, 0.0), (1.0, 0.0), (2.0, -1.0)];
    let text = write_sdf_record(
        &molecule,
        &MolMetadata::default(),
        &coords,
        &Default::default(),
    );
    assert!(!text.contains("  2  3  2  3\n"), "{text}");
    let (parsed, _, observed_coords) = parse_mol_with_coords(&text).unwrap();
    assert_eq!(observed_coords, coords);
    assert_eq!(
        chematic_chem::assign_ez_bonds(&parsed),
        chematic_chem::assign_ez_bonds(&molecule)
    );
}
