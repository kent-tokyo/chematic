use chematic_core::nontetrahedral::{NonTetrahedralClass as Class, across, permute, swap};

#[test]
fn invalid_permutations_and_ligand_positions_do_not_produce_stereo_assignments() {
    for class in [
        Class::SquarePlanar,
        Class::TrigonalBipyramidal,
        Class::Octahedral,
    ] {
        let outside = class.max_nbors();
        assert_eq!(swap(class, 1, 0, outside), 0);
        assert_eq!(swap(class, class.max_permutation() + 1, 0, 1), 0);
        assert_eq!(across(class, 0, 0), None);
        assert_eq!(across(class, class.max_permutation() + 1, 0), None);
        assert_eq!(across(class, 1, outside), None);
        assert_eq!(permute(class, 1, &[Some(1), Some(2)], &[Some(1)]), None);
    }
}
