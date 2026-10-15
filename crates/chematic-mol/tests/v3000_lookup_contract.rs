use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::read_mol_v3000_with_diagnostics;

fn block(ids: &[u32], bonds: &[(u32, u32, u32)]) -> String {
    let mut text = format!(
        "lookup\n  validation        2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n",
        ids.len(),
        bonds.len()
    );
    for (index, id) in ids.iter().enumerate() {
        text.push_str(&format!("M  V30 {id} C {index} 0 0 0\n"));
    }
    text.push_str("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for (id, a, b) in bonds {
        text.push_str(&format!("M  V30 {id} 1 {a} {b}\n"));
    }
    text.push_str("M  V30 END BOND\nM  V30 END CTAB\nM  END\n");
    text
}

#[test]
fn dense_sparse_permuted_and_zero_ids_resolve_to_source_order() {
    for ids in [[1, 2, 3], [10, 20, 30], [3, 1, 2], [0, 1, 2]] {
        let report = read_mol_v3000_with_diagnostics(&block(
            &ids,
            &[(1, ids[0], ids[2]), (2, ids[2], ids[1])],
        ))
        .unwrap();
        assert_eq!(report.mol.bond_count(), 2);
        assert!(report.mol.bond_between(AtomIdx(0), AtomIdx(2)).is_some());
        assert!(report.mol.bond_between(AtomIdx(2), AtomIdx(1)).is_some());
    }
}

#[test]
fn duplicate_atom_ids_keep_the_historical_first_match() {
    // This parser historically accepts duplicate atom ids. Its lookup must
    // continue to choose the first occurrence until that contract changes.
    let report = read_mol_v3000_with_diagnostics(&block(&[2, 2, 3], &[(1, 2, 3)])).unwrap();
    assert!(report.mol.bond_between(AtomIdx(0), AtomIdx(2)).is_some());
    assert!(report.mol.bond_between(AtomIdx(1), AtomIdx(2)).is_none());
}

#[test]
fn duplicate_bond_ids_and_missing_endpoints_still_fail() {
    for input in [
        block(&[1, 2, 3], &[(1, 1, 2), (1, 2, 3)]),
        block(&[1, 2, 3], &[(8, 1, 2), (8, 2, 3)]),
        block(&[1, 2, 3], &[(1, 0, 2)]),
        block(&[1, 2, 3], &[(1, 1, 99)]),
    ] {
        assert!(read_mol_v3000_with_diagnostics(&input).is_err());
    }
}

#[test]
fn continued_atom_properties_keep_typed_and_opaque_fields() {
    let input = block(&[1, 2], &[(1, 1, 2)]).replace(
        "M  V30 1 C 0 0 0 0",
        "M  V30 1 C 0 0 0 0 CHG=1 -\nM  V30 HCOUNT=2 FOO=(3 7 8 9)",
    );
    let report = read_mol_v3000_with_diagnostics(&input).unwrap();
    assert_eq!(report.mol.atom(AtomIdx(0)).charge, 1);
    assert_eq!(report.mol.atom(AtomIdx(0)).hydrogen_count, Some(2));
    assert_eq!(
        report.metadata.v3000_atom_properties,
        [(1, "FOO=(3 7 8 9)".into())]
    );
    assert_eq!(
        report.mol.bonds().next().unwrap().1.order,
        BondOrder::Single
    );
}

#[test]
fn total_valence_overrides_do_not_invent_implicit_hydrogens() {
    for (value, expected) in [
        (0, None),
        (-1, Some(0)),
        (15, Some(0)),
        (1, Some(0)),
        (4, Some(3)),
    ] {
        let input = block(&[1, 2], &[(1, 1, 2)]).replace(
            "M  V30 1 C 0 0 0 0",
            &format!("M  V30 1 C 0 0 0 0 VAL={value}"),
        );
        let report = read_mol_v3000_with_diagnostics(&input).unwrap();
        assert_eq!(
            report.mol.atom(AtomIdx(0)).hydrogen_count,
            expected,
            "VAL={value}"
        );
        assert!(
            report.metadata.v3000_atom_properties[0]
                .1
                .contains(&format!("VAL={value}"))
        );
    }
}

#[test]
fn sulfur_radical_total_valence_is_retained() {
    let report = read_mol_v3000_with_diagnostics(include_str!(
        "fixtures/v3000_total_valence_sulfur_radical.mol"
    ))
    .unwrap();
    assert_eq!(report.mol.atom(AtomIdx(11)).hydrogen_count, Some(0));
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&report.mol).unwrap(),
        "CCCCOc1ccc(NC[S](=O)=O)cn1"
    );
}

#[test]
fn iodine_radical_total_valence_is_retained() {
    let report = read_mol_v3000_with_diagnostics(include_str!(
        "fixtures/v3000_total_valence_iodine_radical.mol"
    ))
    .unwrap();
    assert_eq!(report.mol.atom(AtomIdx(0)).hydrogen_count, Some(0));
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&report.mol).unwrap(),
        "c1ccc2c(c1)[I]c1ccccc1-2"
    );
}

#[test]
fn unspecified_azo_stays_unspecified_after_v3000_rewriting() {
    let report =
        read_mol_v3000_with_diagnostics(include_str!("fixtures/v3000_unspecified_azo.mol"))
            .unwrap();
    let expected = chematic_smiles::rdkit_canonical_smiles(&report.mol).unwrap();
    let output = chematic_mol::write_mol_v3000(&report.mol, &report.metadata, &report.coords);
    assert!(output.contains("CFG=2"));
    let reread = read_mol_v3000_with_diagnostics(&output).unwrap();
    assert_eq!(
        chematic_smiles::rdkit_canonical_smiles(&reread.mol).unwrap(),
        expected
    );
}
