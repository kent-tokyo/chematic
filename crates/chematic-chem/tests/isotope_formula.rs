use chematic_chem::{add_hydrogens, calc_mol_formula, calc_mol_formula_with_isotopes};

#[test]
fn authored_isotope_charge_and_hydrogen_formulas_match() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/openbabel-1464-2026-10-24/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let smiles = case["smiles"].as_str().unwrap();
        let mol = chematic_smiles::parse(smiles).unwrap();
        let legacy = calc_mol_formula(&mol);
        let got = calc_mol_formula_with_isotopes(&mol);
        assert_eq!(got, case["expected_formula"], "{smiles}");
        assert_eq!(
            calc_mol_formula(&mol),
            legacy,
            "legacy formula changed: {smiles}"
        );
        assert_eq!(calc_mol_formula_with_isotopes(&mol.clone()), got);
        // Adding ordinary H graph atoms must neither erase isotope labels nor
        // double count the explicit and bracket hydrogen representations.
        if !mol.atoms().any(|(_, a)| a.wildcard) {
            let expanded = add_hydrogens(&mol);
            assert_eq!(calc_mol_formula_with_isotopes(&expanded), got, "{smiles}");
        }
    }
}

#[test]
fn empty_formula_and_core_explicit_atom_helper_keep_their_contracts() {
    let empty = chematic_core::MoleculeBuilder::new().build();
    assert_eq!(calc_mol_formula_with_isotopes(&empty), "H0");
    let labelled = chematic_smiles::parse("[13CH4]").unwrap();
    assert_eq!(labelled.formula_with_isotopes(), "13C");
    assert_eq!(calc_mol_formula_with_isotopes(&labelled), "[13C]H4");
    assert_eq!(calc_mol_formula(&labelled), "CH4");
}
