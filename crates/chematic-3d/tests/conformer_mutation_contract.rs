use chematic_3d::{ConformerEnsemble, ConformerError, Coords3D, Point3};
use chematic_core::AtomIdx;
use chematic_smiles::parse;

#[test]
fn ensemble_validates_and_mutates_conformer_coordinates() {
    let molecule = parse("CC").unwrap();
    assert!(matches!(
        ConformerEnsemble::with_conformer(molecule.clone(), Coords3D::new_zeroed(1)),
        Err(ConformerError::AtomCountMismatch {
            expected: 2,
            got: 1
        })
    ));
    let mut ensemble =
        ConformerEnsemble::with_conformer(molecule, Coords3D::new_zeroed(2)).unwrap();
    assert_eq!(ensemble.add_conformer(Coords3D::new_zeroed(2)).unwrap(), 1);
    ensemble
        .get_conformer_mut(1)
        .unwrap()
        .set(AtomIdx(0), Point3::new(2.0, 0.0, 0.0));
    assert!((ensemble.conformer_rmsd_no_align(0, 1).unwrap() - 2.0_f64.sqrt()).abs() < 1e-12);
    assert!(ensemble.get_conformer_mut(2).is_none());
    assert!(ensemble.remove_conformer(0).is_some());
    assert_eq!(ensemble.get_conformer(0).unwrap().get(AtomIdx(0)).x, 2.0);
}
