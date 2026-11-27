use chematic_mol::{MrvError, parse_mrv};

#[test]
fn mrv_radical_literals_are_validated() {
    for radical in [
        "monovalent",
        "divalent",
        "divalent1",
        "divalent3",
        "trivalent",
        "trivalent2",
        "trivalent4",
        "4",
    ] {
        let source = format!(
            r#"<cml><molecule><atomArray><atom id="a1" elementType="C" radical="{radical}"/></atomArray></molecule></cml>"#
        );
        let record = parse_mrv(&source).unwrap_or_else(|error| panic!("{radical}: {error}"));
        assert_eq!(record.mol.atom_count(), 1);
    }
    let source = r#"<cml><molecule><atomArray><atom id="a1" elementType="C" radical="pentavalent"/></atomArray></molecule></cml>"#;
    assert!(matches!(
        parse_mrv(source),
        Err(MrvError::InvalidRadical { .. })
    ));
}
