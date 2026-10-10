use crate::*;

#[test]
fn every_bond_order_has_an_unambiguous_display_and_valence_contract() {
    for (order, token, value, integer, pi) in [
        (BondOrder::Single, "-", Some(1.0), 1, false),
        (BondOrder::Double, "=", Some(2.0), 2, true),
        (BondOrder::Triple, "#", Some(3.0), 3, true),
        (BondOrder::Quadruple, "$", Some(4.0), 4, true),
        (BondOrder::Aromatic, ":", None, 1, false),
        (BondOrder::Up, "/", Some(1.0), 1, false),
        (BondOrder::Down, "\\", Some(1.0), 1, false),
        (BondOrder::Zero, "~", Some(0.0), 0, false),
        (BondOrder::Dative, "->", Some(1.0), 1, false),
        (BondOrder::QueryAny, "~", None, 1, false),
        (BondOrder::QuerySingleOrDouble, "~", None, 1, false),
        (BondOrder::QuerySingleOrAromatic, "~", None, 1, false),
        (BondOrder::QueryDoubleOrAromatic, "~", None, 1, false),
    ] {
        assert_eq!(order.to_string(), token);
        assert_eq!(order.smiles_char(), token.chars().next().unwrap());
        assert_eq!(order.order_value(), value);
        assert_eq!(order.order_int(), integer);
        assert_eq!(crate::valence::is_pi_bond(order), pi);
    }
}

#[test]
fn coordination_chirality_tokens_and_mirror_permutations_are_reversible() {
    for class in [
        NonTetrahedralClass::SquarePlanar,
        NonTetrahedralClass::TrigonalBipyramidal,
        NonTetrahedralClass::Octahedral,
    ] {
        for permutation in 0..=class.max_permutation() {
            let tag = Chirality::from_nontetrahedral(class, permutation).unwrap();
            assert!(!tag.is_tetrahedral());
            assert_eq!(tag.nontetrahedral(), Some((class, permutation)));
            let token = tag.smiles_token();
            assert!(token.starts_with(&format!("@{}", class.token())));
            if permutation == 0 {
                assert_eq!(token, format!("@{}", class.token()));
            } else {
                assert_eq!(token, format!("@{}{permutation}", class.token()));
            }
            let mirror = crate::nontetrahedral::invert(class, permutation);
            assert!(mirror <= class.max_permutation());
            assert_eq!(crate::nontetrahedral::invert(class, mirror), permutation);
        }
        assert_eq!(
            Chirality::from_nontetrahedral(class, class.max_permutation() + 1),
            None
        );
    }
    assert_eq!(Chirality::None.smiles_token(), "");
    assert_eq!(Chirality::CounterClockwise.smiles_token(), "@");
    assert_eq!(Chirality::Clockwise.smiles_token(), "@@");
    let mut builder = MoleculeBuilder::new();
    let mut carbon = Atom::new(Element::C);
    carbon.hydrogen_count = Some(3);
    assert_eq!(carbon.explicit_hcount(), Some(3));
    let idx = builder.add_atom(carbon);
    let mol = builder.build();
    #[allow(deprecated)]
    let count = total_hcount(&mol, idx);
    assert_eq!(count, implicit_hcount(&mol, idx));
}
