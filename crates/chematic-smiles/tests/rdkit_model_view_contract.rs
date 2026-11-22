use chematic_core::{AtomIdx, BondIdx, BondOrder};
use chematic_smiles::{parse, rdkit_model_correct_view};

#[test]
fn correction_refuses_different_graph_shapes_and_preserves_inputs() {
    let original = parse("C=C").unwrap();
    assert!(rdkit_model_correct_view(&original, &parse("C").unwrap()).is_none());
    assert!(rdkit_model_correct_view(&original, &original).is_none());

    let mut changed_view = original.clone();
    changed_view.set_charge(AtomIdx(0), 1);
    changed_view.set_bond_order(BondIdx(0), BondOrder::Triple);
    let corrected = rdkit_model_correct_view(&original, &changed_view).unwrap();
    assert_eq!(corrected.atom(AtomIdx(0)).charge, 0);
    assert_eq!(corrected.bond(BondIdx(0)).order, BondOrder::Double);
    assert_eq!(changed_view.atom(AtomIdx(0)).charge, 1);
    assert_eq!(changed_view.bond(BondIdx(0)).order, BondOrder::Triple);
    assert_eq!(original.bond(BondIdx(0)).order, BondOrder::Double);

    for order in [BondOrder::Quadruple, BondOrder::QueryAny, BondOrder::Dative] {
        let mut view = original.clone();
        view.set_bond_order(BondIdx(0), order);
        let corrected = rdkit_model_correct_view(&original, &view).unwrap();
        assert_eq!(
            corrected.bond(BondIdx(0)).order,
            BondOrder::Double,
            "{order:?}"
        );
        assert_eq!(view.bond(BondIdx(0)).order, order);
    }
}
