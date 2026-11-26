use chematic_mol::{MrvError, parse_mrv};

#[test]
fn repeated_bonds_are_rejected_in_either_direction() {
    for refs in ["a1 a2", "a2 a1"] {
        let text = format!(
            r#"<cml><molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="O"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="1"/><bond atomRefs2="{refs}" order="1"/></bondArray></molecule></cml>"#
        );
        assert!(parse_mrv(&text).is_err());
    }
}

#[test]
fn nested_unlabelled_molecule_is_not_silently_discarded() {
    let text = r#"<cml><molecule><atomArray><atom id="a1" elementType="C"/></atomArray><molecule/></molecule></cml>"#;
    assert!(matches!(
        parse_mrv(text),
        Err(MrvError::UnsupportedFeature { .. })
    ));
}

#[test]
fn xml_comments_do_not_change_the_parsed_molecule() {
    let base =
        r#"<cml><molecule><atomArray><atom id="a1" elementType="C"/></atomArray></molecule></cml>"#;
    let with_comment = base.replace("<atomArray>", "<atomArray><!-- carbon -->");
    let plain = parse_mrv(base).unwrap();
    let commented = parse_mrv(&with_comment).unwrap();
    assert_eq!(commented.mol.atom_count(), plain.mol.atom_count());
    assert_eq!(commented.mol.total_formula(), plain.mol.total_formula());
}
