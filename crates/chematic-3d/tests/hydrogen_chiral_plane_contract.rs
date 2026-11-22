use chematic_3d::add_hydrogens_with_coords;
use chematic_core::{AtomIdx, Coords3D, Point3};
use chematic_smiles::parse;

#[test]
fn opposite_chiral_tags_place_hydrogen_on_opposite_sides_of_plane() {
    // Fixed-coordinate reference: RDKit 2025.09.3 Chem.AddHs(addCoords=True).
    let coords = Coords3D {
        points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(-0.5, 0.866, 0.0),
            Point3::new(-0.5, -0.866, 0.0),
        ],
    };
    for (smiles, expected_z) in [("N[C@H](F)Cl", 1.1), ("N[C@@H](F)Cl", -1.1)] {
        let mol = parse(smiles).unwrap();
        let (expanded, full) = add_hydrogens_with_coords(&mol, &coords).unwrap();
        let h = expanded
            .neighbors(AtomIdx(1))
            .map(|(idx, _)| idx)
            .find(|&idx| expanded.atom(idx).element.atomic_number() == 1)
            .expect("chiral center hydrogen");
        let point = full.points[h.0 as usize];
        assert!(point.x.abs() < 1e-9 && point.y.abs() < 1e-9);
        assert!((point.z - expected_z).abs() < 1e-9, "{smiles}: {point:?}");
        for i in 0..coords.points.len() {
            assert!((full.points[i].sub(&coords.points[i])).norm() < 1e-9);
        }
    }
}
