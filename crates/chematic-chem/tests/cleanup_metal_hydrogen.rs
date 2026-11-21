use chematic_chem::{CleanupError, calc_mol_formula, cleanup, cleanup_checked};
use chematic_smiles::{parse, rdkit_cleanup_canonical_smiles};

#[test]
fn preserves_copper_hydrogen_without_changing_the_selected_organic_site() {
    let input = "CN1CCc2cc[c]3c4c2[C@H]1Cc1ccc[c](c1-4)[Cu@@H]3[CH3]";
    let molecule = parse(input).unwrap();
    let cleaned = cleanup(&molecule).unwrap();
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

#[test]
fn checked_cleanup_refuses_unresolved_hydrogen_and_keeps_the_input() {
    let input = "C/N=C(\\S)N/N=C/C(C)=[N]/[SnH]=[CH]/N";
    let molecule = parse(input).unwrap();
    let before = chematic_smiles::canonical_smiles(&molecule);
    assert!(matches!(
        cleanup(&molecule),
        Err(CleanupError::AmbiguousMetalHydrogen)
    ));
    assert!(matches!(
        cleanup_checked(&molecule),
        Err(CleanupError::AmbiguousMetalHydrogen)
    ));
    assert_eq!(chematic_smiles::canonical_smiles(&molecule), before);
}

#[test]
fn checked_cleanup_retains_ferrocene_and_resolved_metal_hydrogen() {
    let ferrocene = parse("C1=CC=CC1[Fe]C1C=CC=C1").unwrap();
    let cleaned = cleanup_checked(&ferrocene).unwrap();
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
        Ok("C1=C[CH-]C=C1.C1=C[CH-]C=C1.[Fe+2]")
    );
    let hydride = parse("C[PtH]C").unwrap();
    let cleaned = cleanup_checked(&hydride).unwrap();
    assert_eq!(calc_mol_formula(&cleaned), calc_mol_formula(&hydride));
}
