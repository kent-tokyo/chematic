use chematic_chem::{calc_mol_formula, cleanup};
use chematic_smiles::{parse, rdkit_cleanup_canonical_smiles};

#[test]
fn preserves_copper_hydrogen_without_changing_the_selected_organic_site() {
    let input = "CN1CCc2cc[c]3c4c2[C@H]1Cc1ccc[c](c1-4)[Cu@@H]3[CH3]";
    let molecule = parse(input).unwrap();
    let cleaned = cleanup(&molecule);
    let expected = "C.CN1CCc2cc[c]3c4c2[C@H]1[CH-]c1ccc[c](c1-4)[CuH+]3";
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
        Ok(expected)
    );
    assert_eq!(calc_mol_formula(&cleaned), calc_mol_formula(&molecule));
    assert_eq!(
        calc_mol_formula(&parse(expected).unwrap()),
        calc_mol_formula(&molecule)
    );
}
