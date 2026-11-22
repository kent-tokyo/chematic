use chematic_chem::{
    CleanupError, add_hydrogens, calc_mol_formula, cleanup_checked, cleanup_checked as cleanup,
};
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
fn checked_cleanup_preserves_tin_hydrogen_and_keeps_the_input() {
    let mut previous = None;
    for input in [
        "C/N=C(\\S)N/N=C/C(C)=[N]/[SnH]=[CH]/N",
        "C/N=C(\\S)N/N=C/C(C)=[N]/[Sn]([H])=[CH]/N",
    ] {
        let molecule = parse(input).unwrap();
        let before = chematic_smiles::canonical_smiles(&molecule);
        let result = cleanup(&molecule).unwrap();
        let smiles = rdkit_cleanup_canonical_smiles(&result).unwrap();
        assert_eq!(calc_mol_formula(&result), calc_mol_formula(&molecule));
        assert_eq!(
            rdkit_cleanup_canonical_smiles(&cleanup_checked(&result).unwrap()).unwrap(),
            smiles
        );
        if let Some(first) = &previous {
            assert_eq!(first, &smiles);
        }
        previous = Some(smiles);
        assert_eq!(chematic_smiles::canonical_smiles(&molecule), before);
    }
}

#[test]
fn checked_cleanup_refuses_unresolved_metal_hydrogen() {
    let molecule = parse("[Na]N=[N+]=[N-]").unwrap();
    assert!(matches!(
        cleanup(&molecule),
        Err(CleanupError::AmbiguousMetalHydrogen)
    ));
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
    let input = "C[N]([SnH])C";
    let source = parse(input).unwrap();
    let expanded = add_hydrogens(&source);
    let compact = cleanup(&source).unwrap();
    let explicit = cleanup(&expanded).unwrap();
    assert_eq!(calc_mol_formula(&compact), calc_mol_formula(&source));
    assert_eq!(calc_mol_formula(&explicit), calc_mol_formula(&source));
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&compact),
        rdkit_cleanup_canonical_smiles(&explicit)
    );
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

#[test]
fn platinum_aromatic_cut_preserves_the_ligand_radical_and_formula() {
    let original = parse("Cc1c[c]([Pt])c(Oc2ccc(O)cc2CC(=O)O)s(Cl)c1").unwrap();
    let expanded = parse("C(C1=C([H])[C]([Pt])=C(Oc2c([H])c([H])c(O[H])c([H])c2C(C(=O)O[H])([H])[H])S(Cl)=C1[H])([H])([H])[H]").unwrap();
    let a = cleanup(&original).unwrap();
    let b = cleanup(&expanded).unwrap();
    assert_eq!(calc_mol_formula(&a), calc_mol_formula(&original));
    assert_eq!(calc_mol_formula(&b), calc_mol_formula(&original));
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&a),
        rdkit_cleanup_canonical_smiles(&b)
    );
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&a),
        rdkit_cleanup_canonical_smiles(&cleanup(&a).unwrap())
    );
}

#[test]
fn copper_aromatic_ligand_is_independent_of_kekule_and_explicit_hydrogens() {
    let original =
        parse("CN(C)C1=[C]([Cu][C@@H](CC2=CC=C(NC(=O)C3=C(Cl)C=NC=C3Cl)C=C2)C(=O)O)C(=O)C1=O")
            .unwrap();
    let explicit = parse("C(N(C([H])([H])[H])c1[c]([Cu][C@@](C(c2c([H])c([H])c(N(C(=O)c3c(Cl)c([H])nc([H])c3Cl)[H])c([H])c2[H])([H])[H])(C(=O)O[H])[H])c(=O)c1=O)([H])([H])[H]").unwrap();
    let a = cleanup(&original).unwrap();
    let b = cleanup(&explicit).unwrap();
    assert_eq!(calc_mol_formula(&a), calc_mol_formula(&original));
    assert_eq!(calc_mol_formula(&b), calc_mol_formula(&original));
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&a),
        rdkit_cleanup_canonical_smiles(&b)
    );
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&a),
        rdkit_cleanup_canonical_smiles(&cleanup(&a).unwrap())
    );
}

#[test]
fn metal_disconnection_does_not_invent_an_imine_hydrogen() {
    let source =
        parse("[Na]12=[C](C=C(F)C=[CH]1)N(C1=C[C+]=C(C=C1)Cl)C1C(=[N]2)C(=O)N(C)C(=O)N=1").unwrap();
    let expanded = add_hydrogens(&source);
    let compact = cleanup(&source).unwrap();
    let explicit = cleanup(&expanded).unwrap();
    assert_eq!(calc_mol_formula(&compact), calc_mol_formula(&source));
    assert_eq!(calc_mol_formula(&explicit), calc_mol_formula(&source));
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&compact),
        rdkit_cleanup_canonical_smiles(&explicit)
    );
    assert_eq!(
        rdkit_cleanup_canonical_smiles(&compact),
        rdkit_cleanup_canonical_smiles(&cleanup(&compact).unwrap())
    );
}

#[test]
fn refuses_unrepresentable_metal_tetrahedral_stereo_for_both_spellings() {
    for input in [
        "O=C(O)[C@@H]1([SnH])C[C@@H](/C=C/CP(=O)(O)O)CCN1",
        "O=C(O[H])[C@]1(->[Sn][H])([H])C([H])([H])[C@](/C(=C(/C(P(=O)(O[H])O[H])([H])[H])[H])[H])([H])C([H])([H])C([H])([H])N1[H]",
    ] {
        let source = parse(input).unwrap();
        assert!(
            matches!(cleanup(&source), Err(CleanupError::UnstableMetalStereo)),
            "{input}"
        );
    }
}
