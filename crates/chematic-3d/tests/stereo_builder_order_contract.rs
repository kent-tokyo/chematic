use chematic_3d::{
    coords::{Coords3D, Point3},
    stereo_constraints::verify_stereo,
};
use chematic_core::{AtomIdx, MoleculeBuilder};

#[test]
fn builder_stereocenters_recover_the_same_neighbor_order_as_parser() {
    for input in ["[C@H](F)(Cl)Br", "F[C@H](Cl)Br"] {
        let parsed = chematic_smiles::parse(input).unwrap();
        let mut builder = MoleculeBuilder::new();
        for (_, atom) in parsed.atoms() {
            builder.add_atom(atom.clone());
        }
        for (_, bond) in parsed.bonds() {
            builder
                .add_bond(bond.atom1, bond.atom2, bond.order)
                .unwrap();
        }
        let rebuilt = builder.build();
        let center = rebuilt
            .atoms()
            .find(|(_, a)| a.chirality.is_tetrahedral())
            .unwrap()
            .0;
        assert!(rebuilt.stereo_neighbor_order(center).is_none());
        let mut coords = Coords3D::new_zeroed(parsed.atom_count());
        for (i, point) in [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
            .iter()
            .enumerate()
        {
            coords.set(AtomIdx(i as u32), Point3::new(point[0], point[1], point[2]));
        }
        let expected = verify_stereo(&parsed, &coords);
        let actual = verify_stereo(&rebuilt, &coords);
        assert_eq!(expected.tetrahedral.len(), 1, "{input}");
        assert!(
            expected.tetrahedral[0].status.is_satisfied()
                || expected.tetrahedral[0].status.is_violated(),
            "{input}"
        );
        assert_eq!(
            actual.tetrahedral[0].status, expected.tetrahedral[0].status,
            "{input}"
        );
    }
}
