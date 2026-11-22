use chematic_mol::read_mol_v3000_with_diagnostics;

fn with_collection(line: &str) -> String {
    format!(
        "collection\n  validation        2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 END BOND\nM  V30 BEGIN COLLECTION\nM  V30 {line}\nM  V30 END COLLECTION\nM  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn stereo_collection_rejects_invalid_group_ids_and_members() {
    for (line, expected) in [
        (
            "MDLV30/STERELx ATOMS=(1 1)",
            "OR group id is not an integer",
        ),
        (
            "MDLV30/STERACx ATOMS=(1 1)",
            "AND group id is not an integer",
        ),
        (
            "MDLV30/STEBRELx BONDS=(1 1)",
            "OR group id is not an integer",
        ),
        (
            "MDLV30/STEBRACx BONDS=(1 1)",
            "AND group id is not an integer",
        ),
        ("MDLV30/STEABS ATOMS=(1 x)", "non-integer member id"),
        ("MDLV30/STEABS ATOMS=(1 9)", "unknown atom id 9"),
        ("MDLV30/STEBABS BONDS=(1 9)", "unknown bond id 9"),
        ("MDLV30/STEABS ATOMS=(0)", "no atom ids"),
        ("MDLV30/STEBABS BONDS=(0)", "no bond ids"),
        ("MDLV30/STEABS ATOMS=(2 1)", "declares 2 ids but contains 1"),
    ] {
        let error = read_mol_v3000_with_diagnostics(&with_collection(line))
            .err()
            .expect("invalid collection");
        assert!(error.to_string().contains(expected), "{line}: {error}");
    }
}

#[test]
fn repeated_absolute_collection_lines_merge_and_deduplicate_members() {
    let input = with_collection(
        "MDLV30/STEABS ATOMS=(2 2 1)\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 MDLV30/STEBABS BONDS=(1 1)\nM  V30 MDLV30/STEBABS BONDS=(1 1)",
    );
    let report = read_mol_v3000_with_diagnostics(&input).unwrap();
    assert_eq!(report.mol.stereo_groups().len(), 1);
    assert_eq!(
        report.mol.stereo_groups()[0].atom_indices,
        [chematic_core::AtomIdx(0), chematic_core::AtomIdx(1)]
    );
    assert_eq!(
        report.mol.stereo_group_bond_indices(0).unwrap(),
        [chematic_core::BondIdx(0)]
    );
}

#[test]
fn malformed_rgroup_lists_are_reported() {
    let input = with_collection("MDLV30/STEABS ATOMS=(1 1)");
    for (property, expected) in [("RGROUPS=(x 7)", "RGROUPS property lacks a valid count")] {
        let molecule = input.replace(
            "M  V30 1 C 0 0 0 0\n",
            &format!("M  V30 1 C 0 0 0 0 {property}\n"),
        );
        let error = read_mol_v3000_with_diagnostics(&molecule)
            .err()
            .expect("invalid RGROUPS");
        assert!(error.to_string().contains(expected), "{property}: {error}");
    }
}
