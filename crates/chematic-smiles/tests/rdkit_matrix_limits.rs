//! RDKit #9647: reject dense matrix dimensions before preparing large graphs.
use chematic_core::{Atom, Element, MoleculeBuilder};
use chematic_smiles::{RdkitSmilesError, parse, rdkit_distance_matrix, rdkit_distance_matrix_3d};

#[test]
fn excessive_native_graphs_fail_before_rdkit_preparation() {
    // No SMILES provenance: without the preflight this would fail in model
    // conversion instead. A matrix-size error proves the early guard ran.
    for n in [1025, 46_341, 65_535, 65_536, 65_537] {
        let mut builder = MoleculeBuilder::new();
        for _ in 0..n {
            builder.add_atom(Atom::new(Element::C));
        }
        let mol = builder.build();
        for (use_bo, use_atom_wts) in [(false, false), (true, false), (false, true), (true, true)] {
            let err = rdkit_distance_matrix(&mol, use_bo, use_atom_wts).unwrap_err();
            assert!(matches!(err, RdkitSmilesError::Unsupported(_)));
            assert!(
                err.to_string().contains("distance matrix resource limit"),
                "{err}"
            );
            assert!(err.to_string().contains(&format!("{n} atoms")), "{err}");
        }
        // The guard also precedes coordinate validation and allocation.
        let err = rdkit_distance_matrix_3d(&mol, &[], false).unwrap_err();
        assert!(
            err.to_string().contains("distance matrix resource limit"),
            "{err}"
        );
    }
}

#[test]
fn exact_input_cap_still_supports_cartesian_distances() {
    let input = vec!["[He]"; 1024].join(".");
    let mol = parse(&input).unwrap();
    let coords: Vec<_> = (0..1024).map(|i| [i as f64, 0., 0.]).collect();
    let matrix = rdkit_distance_matrix_3d(&mol, &coords, false).unwrap();
    assert_eq!(matrix.len(), 1024);
    assert!(matrix.iter().all(|row| row.len() == 1024));
    assert_eq!(matrix[0][1023], 1023.);
    assert_eq!(matrix[1023][0], 1023.);
    assert_eq!(matrix[512][512], 0.);
}

#[test]
fn input_cap_is_explicitly_checked_before_hydrogen_removal() {
    let mol = parse(&vec!["[H]C"; 513].join(".")).unwrap();
    assert_eq!(mol.atom_count(), 1026);
    let err = rdkit_distance_matrix(&mol, false, false).unwrap_err();
    assert!(err.to_string().contains("1026 atoms"));
}

#[test]
fn small_matrix_weights_fragment_sentinel_and_h_removal_are_unchanged() {
    let mol = parse("C=CO").unwrap();
    assert_eq!(
        rdkit_distance_matrix(&mol, false, false).unwrap(),
        [[0., 1., 2.], [1., 0., 1.], [2., 1., 0.]]
    );
    assert_eq!(
        rdkit_distance_matrix(&mol, true, false).unwrap(),
        [[0., 0.5, 1.5], [0.5, 0., 1.], [1.5, 1., 0.]]
    );
    assert_eq!(
        rdkit_distance_matrix(&mol, false, true).unwrap()[2][2],
        0.75
    );
    let coords = [[0., 0., 0.], [3., 0., 0.], [3., 4., 0.]];
    assert_eq!(
        rdkit_distance_matrix_3d(&mol, &coords, false).unwrap(),
        [[0., 3., 5.], [3., 0., 4.], [5., 4., 0.]]
    );
    assert!(
        rdkit_distance_matrix_3d(&mol, &coords[..2], false)
            .unwrap_err()
            .to_string()
            .contains("coordinates")
    );
    let fragments = parse("C.O").unwrap();
    assert_eq!(
        rdkit_distance_matrix(&fragments, false, false).unwrap()[0][1],
        1e8
    );
    let explicit = parse("[H]C([H])([H])[H]").unwrap();
    assert_eq!(
        rdkit_distance_matrix(&explicit, false, false).unwrap(),
        [[0.]]
    );
}
