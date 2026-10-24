use chematic_smiles::parse;

#[test]
fn indexed_graph_relation_covers_structural_boundary_diagnostics() {
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    let methane = parse("C").unwrap();
    let ethane = parse("CC").unwrap();
    assert!(matches!(
        same_structure_under_identity(&methane, &ethane, mode),
        StructuralComparison::Distinct
    ));

    let complete = parse("F[C@](Cl)(Br)I").unwrap();
    let mut incomplete = complete.clone();
    incomplete.clear_stereo_neighbor_order(AtomIdx(1));
    assert!(matches!(
        same_structure_under_identity(&complete, &incomplete, mode),
        StructuralComparison::Inconclusive(
            IndexedGraphRelationDiagnostic::TetrahedralParityNotRepresented
        )
    ));

    incomplete.set_stereo_neighbor_order(AtomIdx(1), vec![0, 2, 3]);
    assert_eq!(canonical_tetrahedral_parity(&incomplete, AtomIdx(1)), None);
}

#[cfg(feature = "native-inchi")]
#[test]
fn native_inchi_corroboration_reports_invalid_input_without_guessing() {
    let empty = chematic_core::MoleculeBuilder::new().build();
    let carbon = parse("C").unwrap();
    let mut diagnostics = Vec::new();
    assert_eq!(
        shared_native_inchi_evidence(&empty, &carbon, &mut diagnostics),
        (None, None)
    );
    assert!(
        diagnostics.contains(&IndexedGraphRelationDiagnostic::InvalidMolecule),
        "{diagnostics:?}"
    );
}
