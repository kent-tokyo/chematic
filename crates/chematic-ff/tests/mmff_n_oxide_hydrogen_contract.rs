use chematic_core::Element;
use chematic_ff::assign_mmff94_numeric_types;

#[test]
fn protonated_n_oxides_assign_rdkit_hydrogen_types() {
    // Independently checked with RDKit 2025.09.3 MMFFGetMoleculeProperties.
    for (smiles, hydrogen_count) in [
        ("C[NH+]([O-])C", 1),
        ("C[NH2+]([O-])", 2),
        ("[H][N+](C)(C)[O-]", 1),
    ] {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let (nitrogen, _) = mol.atoms().find(|(_, a)| a.element == Element::N).unwrap();
        assert_eq!(types[nitrogen.0 as usize], 68, "{smiles}");
        let hydrogens: Vec<_> = mol
            .neighbors(nitrogen)
            .filter_map(|(idx, _)| {
                (mol.atom(idx).element == Element::H).then_some(types[idx.0 as usize])
            })
            .collect();
        assert_eq!(hydrogens, vec![23; hydrogen_count], "{smiles}");
    }
}
