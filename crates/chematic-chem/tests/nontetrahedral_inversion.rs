use chematic_chem::invert_stereocenter;
use chematic_core::{AtomIdx, Chirality, SquarePlanarPermutation};

#[test]
fn inversion_preserves_square_planar_and_mirrors_tb_and_octahedral_tags() {
    for (smiles, expected) in [
        (
            "[Pt@SP1](F)(Cl)(Br)I",
            Chirality::SquarePlanar(SquarePlanarPermutation::SP1),
        ),
        ("[P@TB1](F)(Cl)(Br)(I)N", Chirality::TrigonalBipyramidal(2)),
        ("[Co@OH1](F)(Cl)(Br)(I)(N)O", Chirality::Octahedral(2)),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let original = mol.atom(AtomIdx(0)).chirality;
        let inverse = invert_stereocenter(&mol, AtomIdx(0));
        assert_eq!(inverse.atom(AtomIdx(0)).chirality, expected, "{smiles}");
        assert_eq!(inverse.atom_count(), mol.atom_count());
        assert_eq!(inverse.bond_count(), mol.bond_count());
        let restored = invert_stereocenter(&inverse, AtomIdx(0));
        assert_eq!(restored.atom(AtomIdx(0)).chirality, original, "{smiles}");
        assert_eq!(
            restored.stereo_neighbor_order(AtomIdx(0)),
            mol.stereo_neighbor_order(AtomIdx(0))
        );
    }
}
