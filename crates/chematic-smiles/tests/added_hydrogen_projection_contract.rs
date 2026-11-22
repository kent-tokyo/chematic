use chematic_smiles::{parse, rdkit_added_hs_molecule, rdkit_parsed_molecule};

#[test]
fn already_explicit_hydrogen_graph_is_not_projected_twice() {
    let molecule = parse("CO").unwrap();
    let expanded = rdkit_added_hs_molecule(&molecule).expect("ordinary SMILES can add H");
    assert!(expanded.atom_count() > molecule.atom_count());
    assert!(rdkit_added_hs_molecule(&expanded).is_none());
    assert!(rdkit_parsed_molecule(&expanded).is_none());
}
