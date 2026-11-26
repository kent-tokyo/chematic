use super::*;

#[test]
fn coincident_layout_retry_recovers_tetrahedral_stereo() {
    let mol = chematic_smiles::parse("F[C@](Cl)(Br)I").unwrap();
    let centres = tetrahedral_centres(&mol);
    let mut coords = layout_angstrom(&mol, false);
    coords[1] = coords[0];
    let (wedges, unexpressed_centres) = place_wedges(&mol, &coords, &centres, &[]);
    assert_eq!(unexpressed_centres, vec![AtomIdx(1)]);
    let original = coords.clone();
    let mut drawing = StereoDepiction {
        coords,
        wedges,
        unexpressed_centres,
        ..Default::default()
    };

    retry_coincident_centres(&mol, &mut drawing, &centres, &[], &HashMap::new());

    assert!(drawing.unexpressed_centres.is_empty());
    assert_eq!(drawing.wedges.len(), 1);
    assert_eq!(drawing.wedges.values().next().unwrap().start, AtomIdx(1));
    assert_ne!(drawing.coords, original);
}
