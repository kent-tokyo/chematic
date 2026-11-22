use chematic_perception::rdkit_organometallic_dative_bonds;
use chematic_smiles::parse;

#[test]
fn hypervalent_nitrogen_selects_one_metal_bond_for_dative_cleanup() {
    // Bond flags from RDKit 2025.09.3 MolFromSmiles cleanup.
    for (smiles, expected) in [
        ("C[N]([Cu])([Cu])[Cu]", [false, false, false, true]),
        ("C[N]([Cu])([Zn])[Zn]", [false, false, false, true]),
        ("C[N]([Cu])([Zn])C", [false, false, true, false]),
    ] {
        let molecule = parse(smiles).unwrap();
        assert_eq!(
            rdkit_organometallic_dative_bonds(&molecule),
            expected,
            "{smiles}"
        );
    }
}
