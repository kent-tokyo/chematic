use chematic_smiles::{
    RdkitSmilesParams, parse, rdkit_added_hs_molecule, rdkit_cx_smiles, rdkit_parsed_molecule,
};

#[test]
fn already_explicit_hydrogen_graph_is_not_projected_twice() {
    let molecule = parse("CO").unwrap();
    let expanded = rdkit_added_hs_molecule(&molecule).expect("ordinary SMILES can add H");
    assert!(expanded.atom_count() > molecule.atom_count());
    // RDKit 2025.09.3 CXSMILES on the in-memory AddHs graph.
    assert_eq!(
        rdkit_cx_smiles(&expanded, &RdkitSmilesParams::default()).unwrap(),
        "[H]OC([H])([H])[H]"
    );
    assert!(rdkit_added_hs_molecule(&expanded).is_none());
    assert!(rdkit_parsed_molecule(&expanded).is_none());
}
