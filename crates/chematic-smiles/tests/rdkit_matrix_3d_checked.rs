//! Numeric and allocation contracts extending the RDKit #9647 guard.
use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_smiles::{
    parse, rdkit_distance_matrix_3d, rdkit_distance_matrix_3d_checked as checked,
};

#[test]
fn ordinary_distances_and_atomic_weight_diagonals_are_exactly_compatible() {
    for smiles in ["CCO", "C.O.N", "c1ccccc1", "[H]CO"] {
        let mol = parse(smiles).unwrap();
        let n = if smiles == "[H]CO" {
            2
        } else {
            mol.atom_count()
        };
        let coords: Vec<_> = (0..n)
            .map(|i| [i as f64 * 1.25, (i * i) as f64 * 0.125, -0.25])
            .collect();
        for weights in [false, true] {
            assert_eq!(
                checked(&mol, &coords, weights).unwrap(),
                rdkit_distance_matrix_3d(&mol, &coords, weights).unwrap()
            );
        }
    }
}

#[test]
fn representable_large_and_tiny_distances_are_recovered() {
    let mol = parse("CC").unwrap();
    for length in [1e308, 1e-160, 1e-300, f64::from_bits(1)] {
        let coords = [[0., 0., 0.], [length, 0., 0.]];
        let matrix = checked(&mol, &coords, false).unwrap();
        assert_eq!(matrix[0][1], length);
        assert_eq!(matrix[1][0], length);
    }
    let diagonal = checked(&mol, &[[0., 0., 0.], [1e308, 1e308, 1e308]], false).unwrap()[0][1];
    assert!((diagonal / 1e308 - 3.0_f64.sqrt()).abs() < 1e-15);
}

#[test]
fn nonfinite_coordinates_refuse_even_for_single_atom_diagonals() {
    let mol = parse("[He]").unwrap();
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for axis in 0..3 {
            let mut xyz = [0.; 3];
            xyz[axis] = value;
            let err = checked(&mol, &[xyz], false).unwrap_err();
            assert!(err.to_string().contains("non-finite coordinate at atom 0"));
        }
    }
}

#[test]
fn truly_unrepresentable_distances_refuse_without_a_partial_result() {
    let mol = parse("CCC").unwrap();
    for coords in [
        [[0., 0., 0.], [1., 0., 0.], [f64::MAX, f64::MAX, 0.]],
        [[1e308, 0., 0.], [0., 0., 0.], [-1e308, 0., 0.]],
    ] {
        let err = checked(&mol, &coords, false).unwrap_err();
        assert!(
            err.to_string()
                .contains("non-finite 3D distance between atoms")
        );
    }
}

#[test]
fn weighted_dummy_diagonals_are_explicitly_undefined() {
    let mol = parse("*").unwrap();
    assert_eq!(checked(&mol, &[[0., 0., 0.]], false).unwrap(), [[0.]]);
    let err = checked(&mol, &[[0., 0., 0.]], true).unwrap_err();
    assert!(err.to_string().contains("atomic number zero"));
}

#[test]
fn coordinate_count_matches_prepared_atom_order_after_hydrogen_removal() {
    let mol = parse("[H]C([H])([H])[H]").unwrap();
    assert_eq!(checked(&mol, &[[3., 4., 5.]], false).unwrap(), [[0.]]);
    assert!(
        checked(&mol, &[[0.; 3]; 5], false)
            .unwrap_err()
            .to_string()
            .contains("5 coordinates for 1 atoms")
    );
    assert!(
        checked(&parse("CC").unwrap(), &[], false)
            .unwrap_err()
            .to_string()
            .contains("coordinates")
    );
    assert!(
        checked(&MoleculeBuilder::new().build(), &[], true)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn exact_atom_cap_is_accepted() {
    let mol = parse(&vec!["[He]"; 1024].join(".")).unwrap();
    let coords: Vec<_> = (0..1024).map(|i| [i as f64, 0., 0.]).collect();
    let matrix = checked(&mol, &coords, false).unwrap();
    assert_eq!(matrix.len(), 1024);
    assert_eq!(matrix[0][1023], 1023.);
}

#[test]
fn oversized_native_inputs_refuse_before_preparation_or_coordinate_inspection() {
    for n in [1025, 46341, 65535, 65536, 65537] {
        let mut builder = MoleculeBuilder::new();
        for _ in 0..n {
            builder.add_atom(Atom::new(Element::C));
        }
        let err = checked(&builder.build(), &[[f64::NAN; 3]], true).unwrap_err();
        assert!(err.to_string().contains("distance matrix resource limit"));
    }
}

#[test]
fn legacy_api_is_not_silently_changed() {
    let mol = parse("CC").unwrap();
    assert!(
        rdkit_distance_matrix_3d(&mol, &[[0.; 3], [1e308, 0., 0.]], false).unwrap()[0][1]
            .is_infinite()
    );
    assert_eq!(
        rdkit_distance_matrix_3d(&mol, &[[0.; 3], [1e-300, 0., 0.]], false).unwrap()[0][1],
        0.
    );
    assert!(
        rdkit_distance_matrix_3d(&mol, &[[0.; 3], [f64::NAN, 0., 0.]], false).unwrap()[0][1]
            .is_nan()
    );
}
