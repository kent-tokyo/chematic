use chematic_chem::{CleanupError, add_hydrogens, calc_mol_formula, cleanup, cleanup_checked};
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
    for input in [
        "C/N=C(\\S)N/N=C/C(C)=[N]/[SnH]=[CH]/N",
        "C/N=C(\\S)N/N=C/C(C)=[N]/[Sn]([H])=[CH]/N",
    ] {
        let molecule = parse(input).unwrap();
        let before = chematic_smiles::canonical_smiles(&molecule);
        assert!(
            matches!(
                cleanup(&molecule),
                Err(CleanupError::AmbiguousMetalHydrogen)
            ),
            "{input}"
        );
        assert!(
            matches!(
                cleanup_checked(&molecule),
                Err(CleanupError::AmbiguousMetalHydrogen)
            ),
            "{input}"
        );
        assert_eq!(chematic_smiles::canonical_smiles(&molecule), before);
    }
}

#[test]
fn checked_cleanup_retains_ferrocene_and_resolved_metal_hydrogen() {
    let ferrocene = parse("C1=CC=CC1[Fe]C1C=CC=C1").unwrap();
    let cleaned = cleanup_checked(&ferrocene).unwrap();
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&cleaned).as_deref(),
        Ok("C1=C[CH-]C=C1.C1=C[CH-]C=C1.[Fe+2]")
    );
    for input in ["C[PtH]C", "C[Pt]([H])C"] {
        let hydride = parse(input).unwrap();
        let cleaned = cleanup_checked(&hydride).unwrap();
        assert_eq!(
            calc_mol_formula(&cleaned),
            calc_mol_formula(&hydride),
            "{input}"
        );
    }
}

#[test]
fn explicit_ligand_hydrogens_do_not_gain_hydrogen_after_metal_cut() {
    for input in [
        "O=C1C[CH]([Pt]2=[CH]C=C(Cl)C=[CH]2)CN1",
        "O=C(O)C1=[C]2C=C(F)C=[CH][Fe]2=[N]C(C2=CC=C(C3=CC=CC=C3)C=C2)=C1",
        "C=[C](C1=CC=CC=C1OCCC1=CC=CC=C1)[Fe]1[CH]=CN=[CH]1",
    ] {
        let molecule = parse(input).unwrap();
        let expanded = add_hydrogens(&molecule);
        let compact_result = cleanup(&molecule).unwrap();
        let expanded_result = cleanup(&expanded).unwrap();
        assert_eq!(
            calc_mol_formula(&expanded),
            calc_mol_formula(&molecule),
            "{input}"
        );
        assert_eq!(
            calc_mol_formula(&expanded_result),
            calc_mol_formula(&molecule),
            "{input}"
        );
        assert_eq!(
            rdkit_cleanup_canonical_smiles(&expanded_result),
            rdkit_cleanup_canonical_smiles(&compact_result),
            "{input}"
        );
    }
}

#[test]
fn metal_disconnection_does_not_create_hydrogen() {
    for input in ["C[N]([SnH])C"] {
        let source = parse(input).unwrap();
        let expanded = add_hydrogens(&source);
        let compact = cleanup(&source).unwrap();
        let explicit = cleanup(&expanded).unwrap();
        assert_eq!(
            calc_mol_formula(&compact),
            calc_mol_formula(&source),
            "{input}"
        );
        assert_eq!(
            calc_mol_formula(&explicit),
            calc_mol_formula(&source),
            "{input}"
        );
        assert_eq!(
            rdkit_cleanup_canonical_smiles(&compact),
            rdkit_cleanup_canonical_smiles(&explicit),
            "{input}"
        );
    }
}

#[test]
fn refuses_detached_metals_that_would_change_on_repeat() {
    for input in ["O->[Mg](<-O)(<-O)<-O", "C(=[O][Li])O"] {
        let source = parse(input).unwrap();
        for variant in [&source, &add_hydrogens(&source)] {
            assert!(
                matches!(
                    cleanup(variant),
                    Err(CleanupError::UnstableMetalDisconnection)
                ),
                "{input}"
            );
        }
    }
}
