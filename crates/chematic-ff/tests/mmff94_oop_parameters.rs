use chematic_ff::{OOP_SP2_TYPES, mmff94_energy::MMFF94_OOP, mmff94_oop};

#[test]
fn oop_equivalence_precedes_wildcard_for_every_outer_atom_order() {
    // RDKit 2026.03.6: benzamide-like carbonyl and conjugated iminium CH.
    for (center, outer, expected) in [(3, [7, 10, 37], 0.116), (3, [2, 5, 54], 0.081)] {
        let [a, b, c] = outer;
        for [i, k, l] in [
            [a, b, c],
            [a, c, b],
            [b, a, c],
            [b, c, a],
            [c, a, b],
            [c, b, a],
        ] {
            assert_eq!(mmff94_oop(center, i, k, l), Some(expected));
        }
    }
}

#[test]
fn oop_eligible_centers_include_all_parameterized_types() {
    let mut expected: Vec<_> = MMFF94_OOP.iter().map(|row| row.1).collect();
    expected.sort_unstable();
    expected.dedup();
    // Current main also retains eligible centers without a table row.
    // Coverage must include every parameterized center, without dropping those.
    for center in expected {
        assert!(
            OOP_SP2_TYPES.binary_search(&center).is_ok(),
            "missing center {center}"
        );
    }
    // Pyridine N-oxide and cationic conjugated nitrogen were omitted.
    assert!(OOP_SP2_TYPES.binary_search(&69).is_ok());
    assert!(OOP_SP2_TYPES.binary_search(&55).is_ok());
}

#[test]
fn unknown_outer_types_cannot_be_hidden_by_a_wildcard_hit() {
    assert_eq!(mmff94_oop(3, 250, 5, 7), None);
    // Literal wildcard queries remain available for table inspection.
    assert_eq!(mmff94_oop(3, 0, 0, 0), Some(0.13));
}
