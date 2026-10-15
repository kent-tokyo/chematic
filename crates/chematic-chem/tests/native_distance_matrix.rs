//! Allocation and native-graph contract for the #9647 sibling API.
use chematic_chem::{
    NativeDistanceMatrixError, NativeDistanceMatrixLimits, topological_distance_matrix,
    try_topological_distance_matrix, try_topological_distance_matrix_with_limits,
};
use chematic_core::{Atom, AtomIdx, BondOrder, Element, Molecule, MoleculeBuilder};

fn isolated_atoms(count: usize, element: Element) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..count {
        builder.add_atom(Atom::new(element));
    }
    builder.build()
}

#[test]
fn rejects_issue_boundaries_without_provenance_or_dense_allocation() {
    for n in [1025, 46_341, 65_535, 65_536, 65_537] {
        let mol = isolated_atoms(n, Element::C);
        assert_eq!(
            diagnostic(try_topological_distance_matrix(&mol)),
            Err(NativeDistanceMatrixError::ResourceLimit {
                resource: "atoms",
                requested: n,
                limit: 1024
            })
        );
    }
}

#[test]
fn exact_cap_accepts_connected_and_disconnected_graphs() {
    let mol = isolated_atoms(1024, Element::C);
    let matrix = diagnostic(try_topological_distance_matrix(&mol)).unwrap();
    assert_eq!(matrix.len(), 1024);
    for (i, row) in matrix.iter().enumerate() {
        assert_eq!(row.len(), 1024);
        for (j, &d) in row.iter().enumerate() {
            assert_eq!(d, if i == j { 0 } else { u32::MAX });
        }
    }
    let mut builder = MoleculeBuilder::from_molecule(&mol);
    for i in 1..1024 {
        builder
            .add_bond(AtomIdx(i - 1), AtomIdx(i), BondOrder::Single)
            .unwrap();
    }
    let matrix = diagnostic(try_topological_distance_matrix(&builder.build())).unwrap();
    for (i, row) in matrix.iter().enumerate() {
        for (j, &d) in row.iter().enumerate() {
            assert_eq!(d as usize, i.abs_diff(j));
        }
    }
}

#[test]
fn all_input_hydrogens_count_even_when_output_would_be_empty() {
    assert_eq!(
        diagnostic(try_topological_distance_matrix(&isolated_atoms(
            1024,
            Element::H
        )))
        .unwrap(),
        Vec::<Vec<u32>>::new()
    );
    assert!(matches!(
        diagnostic(try_topological_distance_matrix(&isolated_atoms(
            1025,
            Element::H
        ))),
        Err(NativeDistanceMatrixError::ResourceLimit {
            resource: "atoms",
            ..
        })
    ));
}

#[test]
fn hydrogens_are_excluded_as_endpoints_and_intermediates() {
    let mol = chematic_smiles::parse("C[H]O.*.[2H].[H][H]").unwrap();
    let expected = vec![
        vec![0, u32::MAX, u32::MAX],
        vec![u32::MAX, 0, u32::MAX],
        vec![u32::MAX, u32::MAX, 0],
    ];
    assert_eq!(
        diagnostic(try_topological_distance_matrix(&mol)).unwrap(),
        expected
    );
    assert_eq!(topological_distance_matrix(&mol), expected);
}

#[test]
fn tighter_limits_and_zero_budgets_are_exact() {
    let mol = chematic_smiles::parse("CCC").unwrap();
    for (limits, resource, requested, limit) in [
        (
            NativeDistanceMatrixLimits {
                max_atoms: 2,
                ..Default::default()
            },
            "atoms",
            3,
            2,
        ),
        (
            NativeDistanceMatrixLimits {
                max_bonds: 1,
                ..Default::default()
            },
            "bonds",
            2,
            1,
        ),
        (
            NativeDistanceMatrixLimits {
                max_entries: 8,
                ..Default::default()
            },
            "entries",
            9,
            8,
        ),
    ] {
        assert_eq!(
            diagnostic(try_topological_distance_matrix_with_limits(&mol, &limits)),
            Err(NativeDistanceMatrixError::ResourceLimit {
                resource,
                requested,
                limit
            })
        );
    }
    let exact = NativeDistanceMatrixLimits {
        max_atoms: 3,
        max_bonds: 2,
        max_entries: 9,
    };
    assert_eq!(
        diagnostic(try_topological_distance_matrix_with_limits(&mol, &exact)).unwrap(),
        topological_distance_matrix(&mol)
    );
    let zero = NativeDistanceMatrixLimits {
        max_atoms: 0,
        max_bonds: 0,
        max_entries: 0,
    };
    assert!(
        diagnostic(try_topological_distance_matrix_with_limits(
            &MoleculeBuilder::new().build(),
            &zero
        ))
        .unwrap()
        .is_empty()
    );
    assert!(
        diagnostic(try_topological_distance_matrix_with_limits(
            &isolated_atoms(1, Element::C),
            &zero
        ))
        .is_err()
    );
}

#[test]
fn callers_cannot_raise_hard_atom_or_entry_cap() {
    let limits = NativeDistanceMatrixLimits {
        max_atoms: usize::MAX,
        max_bonds: usize::MAX,
        max_entries: usize::MAX,
    };
    assert!(matches!(
        diagnostic(try_topological_distance_matrix_with_limits(
            &isolated_atoms(1025, Element::C),
            &limits
        )),
        Err(NativeDistanceMatrixError::ResourceLimit { limit: 1024, .. })
    ));
}

#[test]
fn dense_bond_cap_is_checked_including_hydrogens() {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..182 {
        builder.add_atom(Atom::new(Element::H));
    }
    for i in 0..182 {
        for j in i + 1..182 {
            builder
                .add_bond(AtomIdx(i), AtomIdx(j), BondOrder::Single)
                .unwrap();
        }
    }
    let mol = builder.build();
    assert_eq!(mol.bond_count(), 16_471);
    let limits = NativeDistanceMatrixLimits {
        max_bonds: usize::MAX,
        ..Default::default()
    };
    assert!(matches!(
        diagnostic(try_topological_distance_matrix_with_limits(&mol, &limits)),
        Err(NativeDistanceMatrixError::ResourceLimit {
            resource: "bonds",
            requested: 16_471,
            limit: 16_384
        })
    ));
}

#[test]
fn exact_bond_cap_is_accepted_and_next_bond_is_rejected() {
    let mut builder = MoleculeBuilder::new();
    for _ in 0..182 {
        builder.add_atom(Atom::new(Element::H));
    }
    let edges: Vec<_> = (0..182)
        .flat_map(|i| (i + 1..182).map(move |j| (i, j)))
        .collect();
    for &(i, j) in &edges[..16_384] {
        builder
            .add_bond(AtomIdx(i), AtomIdx(j), BondOrder::Single)
            .unwrap();
    }
    let mol = builder.build();
    assert!(
        diagnostic(try_topological_distance_matrix(&mol))
            .unwrap()
            .is_empty()
    );
    let mut builder = MoleculeBuilder::from_molecule(&mol);
    let (i, j) = edges[16_384];
    builder
        .add_bond(AtomIdx(i), AtomIdx(j), BondOrder::Single)
        .unwrap();
    assert!(matches!(
        diagnostic(try_topological_distance_matrix(&builder.build())),
        Err(NativeDistanceMatrixError::ResourceLimit {
            resource: "bonds",
            requested: 16_385,
            limit: 16_384
        })
    ));
}

#[test]
fn ordinary_native_graph_results_match_existing_api() {
    for smiles in [
        "CCC",
        "C=CO",
        "c1ccccc1",
        "C1CC2CCC1C2",
        "C.O",
        "[H]C([H])O",
        "*CC",
        "N->[Cu+2]",
        "[Na+].[Cl-]",
        "[2H]OC",
        "F/C=C\\Cl",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        assert_eq!(
            diagnostic(try_topological_distance_matrix(&mol)).unwrap(),
            topological_distance_matrix(&mol),
            "{smiles}"
        );
    }
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
