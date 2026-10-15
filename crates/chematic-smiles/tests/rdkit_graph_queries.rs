//! Original MIT OR Apache-2.0 graph fixtures; pinned oracle/provenance/hashes
//! live in validation/rdkit_issues/graph-queries/fixtures.json.
use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
use chematic_smiles::{RdkitSmilesError, parse, rdkit_adjacency_matrix, rdkit_shortest_path};

#[test]
fn binary_and_bond_order_matrices_are_distinct() {
    let m = parse("C=CC#N").unwrap();
    assert_eq!(
        rdkit_adjacency_matrix(&m, false, 0).unwrap(),
        [
            [0., 1., 0., 0.],
            [1., 0., 1., 0.],
            [0., 1., 0., 1.],
            [0., 0., 1., 0.]
        ]
    );
    assert_eq!(
        rdkit_adjacency_matrix(&m, true, 0).unwrap(),
        [
            [0., 2., 0., 0.],
            [2., 0., 1., 0.],
            [0., 1., 0., 3.],
            [0., 0., 3., 0.]
        ]
    );
}

#[test]
fn dative_valence_is_asymmetric_while_paths_are_undirected() {
    let m = parse("N->[Cu]").unwrap();
    assert_eq!(
        rdkit_adjacency_matrix(&m, false, 0).unwrap(),
        [[0., 1.], [1., 0.]]
    );
    assert_eq!(
        rdkit_adjacency_matrix(&m, true, 0).unwrap(),
        [[0., 0.], [1., 0.]]
    );
    assert_eq!(rdkit_shortest_path(&m, 0, 1).unwrap(), [0, 1]);
    assert_eq!(rdkit_shortest_path(&m, 1, 0).unwrap(), [1, 0]);
}

#[test]
fn aromaticity_h_removal_and_empty_graphs_use_existing_model() {
    let m = parse("C1=CC=CC=C1").unwrap();
    let a = rdkit_adjacency_matrix(&m, true, 0).unwrap();
    assert_eq!(a[0][1], 1.5);
    assert_eq!(a[0][5], 1.5);
    let m = parse("[H]C([H])([H])[H]").unwrap();
    assert_eq!(rdkit_adjacency_matrix(&m, true, 0).unwrap(), [[0.]]);
    assert!(rdkit_shortest_path(&m, 0, 1).is_err());
    let m = MoleculeBuilder::new().build();
    assert!(rdkit_adjacency_matrix(&m, false, 0).unwrap().is_empty());
    assert!(rdkit_shortest_path(&m, 0, 1).is_err());
}

#[test]
fn tied_paths_follow_bond_insertion_order_in_both_directions() {
    let m = parse("C1CCC1").unwrap();
    assert_eq!(rdkit_shortest_path(&m, 0, 2).unwrap(), [0, 1, 2]);
    assert_eq!(rdkit_shortest_path(&m, 2, 0).unwrap(), [2, 1, 0]);
    assert_eq!(rdkit_shortest_path(&m, 1, 3).unwrap(), [1, 0, 3]);
    assert_eq!(rdkit_shortest_path(&m, 3, 1).unwrap(), [3, 2, 1]);
}

#[test]
fn disconnected_is_empty_invalid_or_equal_is_error() {
    let m = parse("CC.O").unwrap();
    assert!(rdkit_shortest_path(&m, 0, 2).unwrap().is_empty());
    for pair in [(0, 0), (3, 0), (0, usize::MAX)] {
        assert!(matches!(
            rdkit_shortest_path(&m, pair.0, pair.1),
            Err(RdkitSmilesError::Unsupported(_))
        ));
    }
    for value in [-1, 1, 127, i32::MAX] {
        assert!(rdkit_adjacency_matrix(&m, false, value).is_err());
    }
}

#[test]
fn input_limits_precede_preparation_and_hydrogen_removal() {
    let m = parse(&vec!["[H]C"; 513].join(".")).unwrap();
    for error in [
        rdkit_adjacency_matrix(&m, false, 0).unwrap_err(),
        rdkit_shortest_path(&m, 0, 1).unwrap_err(),
    ] {
        assert!(error.to_string().contains("graph query resource limit"));
        assert!(error.to_string().contains("1026 atoms"));
    }
    let mut builder = MoleculeBuilder::new();
    for _ in 0..65_536 {
        builder.add_atom(Atom::new(Element::C));
    }
    let m = builder.build();
    assert!(
        rdkit_adjacency_matrix(&m, false, 0)
            .unwrap_err()
            .to_string()
            .contains("resource limit")
    );
    // Dense invalid-valence graph: edge preflight must win over sanitization.
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<_> = (0..182)
        .map(|_| builder.add_atom(Atom::new(Element::C)))
        .collect();
    for i in 0..atoms.len() {
        for j in i + 1..atoms.len() {
            builder
                .add_bond(atoms[i], atoms[j], BondOrder::Single)
                .unwrap();
        }
    }
    let m = builder.build();
    let error = rdkit_shortest_path(&m, 0, 1).unwrap_err();
    assert!(error.to_string().contains("16471 bonds"));
    assert!(error.to_string().contains("graph query resource limit"));
}

#[test]
fn exact_cap_supports_dense_and_linear_queries() {
    let m = parse(&"C".repeat(1024)).unwrap();
    let matrix = rdkit_adjacency_matrix(&m, false, 0).unwrap();
    assert_eq!(matrix.len(), 1024);
    assert!(matrix.iter().all(|r| r.len() == 1024));
    assert_eq!(matrix[0][1023], 0.);
    assert_eq!(
        rdkit_shortest_path(&m, 0, 1023).unwrap(),
        (0..1024).collect::<Vec<_>>()
    );
}

#[test]
fn unsupported_bond_types_and_self_bonds_fail_explicitly() {
    assert!(rdkit_adjacency_matrix(&parse("C~N").unwrap(), true, 0).is_err());
    let mut b = MoleculeBuilder::new();
    let a = b.add_atom(Atom::new(Element::C));
    b.add_bond(a, a, BondOrder::Single).unwrap();
    assert!(
        rdkit_shortest_path(&b.build(), 0, 1)
            .unwrap_err()
            .to_string()
            .contains("self bonds")
    );
}
