use chematic_core::{AtomIdx, Chirality};
use chematic_mol::read_mol_with_diagnostics;
use chematic_perception::StereoRejectionReason;

#[test]
fn sdf_frontier_row_745_retains_other_stereo_but_rejects_the_overlapping_center() {
    let report = read_mol_with_diagnostics(include_str!(
        "fixtures/sdf_overlapping_flat_neighbors_v2000.mol"
    ))
    .unwrap();
    assert_eq!(report.mol.atom(AtomIdx(7)).chirality, Chirality::None);
    assert!(report.mol.stereo_neighbor_order(AtomIdx(7)).is_none());
    assert_eq!(report.stereo_diagnostics.len(), 1);
    assert_eq!(report.stereo_diagnostics[0].atom, AtomIdx(7));
    assert_eq!(
        report.stereo_diagnostics[0].reason,
        StereoRejectionReason::DegenerateGeometry
    );
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&report.mol).unwrap(),
        "CCO[C@H]1O[C@@H]2OC3(C)CCC4[C@H](C)CCC([C@H]1C)[C@]42OO3"
    );
}
