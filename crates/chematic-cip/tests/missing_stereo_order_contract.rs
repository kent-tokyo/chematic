use chematic_cip::{CipBudget, SkipReason, assign_cip_accurate_experimental};
use chematic_core::{AtomIdx, MoleculeBuilder};
use chematic_smiles::parse;

#[test]
fn tagged_center_without_ligand_order_is_skipped_not_guessed() {
    let mol = parse("[C@](F)(Cl)(Br)I").unwrap();
    let ordinary = assign_cip_accurate_experimental(&mol, CipBudget::default_budget()).unwrap();
    assert_eq!(ordinary.assignments.len(), 1);
    let mut edited = MoleculeBuilder::from_molecule(&mol);
    edited.clear_stereo_neighbor_order(AtomIdx(0));
    let incomplete =
        assign_cip_accurate_experimental(&edited.build(), CipBudget::default_budget()).unwrap();
    assert!(incomplete.assignments.is_empty());
    assert_eq!(
        incomplete.skipped,
        [(AtomIdx(0), SkipReason::NotFourSubstituents)]
    );
}

#[test]
fn malformed_ligand_frame_is_skipped_instead_of_assigned() {
    let mol = parse("[C@](F)(Cl)(Br)I").unwrap();
    for order in [vec![1, 2, 3, 3], vec![0, 1, 2, 3]] {
        let mut edited = MoleculeBuilder::from_molecule(&mol);
        edited.set_stereo_neighbor_order(AtomIdx(0), order);
        let result =
            assign_cip_accurate_experimental(&edited.build(), CipBudget::default_budget()).unwrap();
        assert!(result.assignments.is_empty());
        assert_eq!(
            result.skipped,
            [(AtomIdx(0), SkipReason::NotFourSubstituents)]
        );
    }
}
