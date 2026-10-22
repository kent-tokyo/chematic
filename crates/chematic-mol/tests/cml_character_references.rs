//! Original minimal XML lexical fixtures under the project license.
use chematic_core::{AtomIdx, BondOrder};
use chematic_mol::{parse_cml, parse_cml_strict, write_cml};

const SOURCE: &str = "<molecule>\n<atom id=\"a1\" elementType=\"N\" formalCharge=\"1\" isotopeNumber=\"15\" hydrogenCount=\"4\" x2=\"1.25\" y2=\"-2.5\"/>\n<atom id=\"h1\" elementType=\"H\" isotopeNumber=\"2\" hydrogenCount=\"0\" x2=\"2.25\" y2=\"-2.5\"/>\n<bond atomRefs2=\"a1 h1\" order=\"1\"/>\n</molecule>";

#[test]
fn literal_decimal_and_hex_spellings_have_the_same_chemistry() {
    let mut variants = vec![SOURCE.to_owned()];
    for (before, after) in [
        ("elementType=\"N\"", "elementType=\"&#78;\""),
        ("formalCharge=\"1\"", "formalCharge=\"&#x2b;&#49;\""),
        ("isotopeNumber=\"15\"", "isotopeNumber=\"&#49;&#x35;\""),
        ("hydrogenCount=\"4\"", "hydrogenCount=\"&#x34;\""),
        ("x2=\"1.25\"", "x2=\"&#49;.25\""),
        ("id=\"a1\"", "id=\"a&#x31;\""),
        ("atomRefs2=\"a1 h1\"", "atomRefs2=\"a1&#32;h1\""),
        ("order=\"1\"", "order=\"&#49;\""),
    ] {
        variants.push(SOURCE.replace(before, after));
    }
    for source in variants {
        for parser in [parse_cml, parse_cml_strict] {
            let (mol, coords) = parser(&source).unwrap();
            assert_eq!((mol.atom_count(), mol.bond_count()), (2, 1));
            assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
            assert_eq!(mol.atom(AtomIdx(0)).isotope, Some(15));
            assert_eq!(mol.atom(AtomIdx(0)).hydrogen_count, Some(3));
            assert_eq!(mol.atom(AtomIdx(1)).isotope, Some(2));
            assert_eq!(mol.bonds().next().unwrap().1.order, BondOrder::Single);
            assert_eq!(coords, [(1.25, -2.5), (2.25, -2.5)]);
            let output = write_cml(&mol, Some(&coords));
            let (after, after_coords) = parse_cml_strict(&output).unwrap();
            assert_eq!(after.atom(AtomIdx(0)).charge, 1);
            assert_eq!(after.atom(AtomIdx(0)).isotope, Some(15));
            assert_eq!(after.atom(AtomIdx(0)).hydrogen_count, Some(3));
            assert_eq!(after_coords, coords);
        }
    }
}

#[test]
fn malformed_or_out_of_range_chemistry_never_defaults_to_neutral_or_zero() {
    for (field, values) in [
        (
            "formalCharge",
            vec!["oops", "128", "-129", "", "&amp;#49;", "&#0;", "&#xD800;"],
        ),
        (
            "isotopeNumber",
            vec!["oops", "-1", "65536", "", "&#x110000;"],
        ),
        ("x2", vec!["oops", "", "NaN", "inf", "1e400"]),
        ("y2", vec!["oops", "", "NaN", "-inf"]),
    ] {
        for value in values {
            let original = match field {
                "formalCharge" => "1",
                "isotopeNumber" => "15",
                "x2" => "1.25",
                _ => "-2.5",
            };
            let source = SOURCE.replace(
                &format!("{field}=\"{original}\""),
                &format!("{field}=\"{value}\""),
            );
            assert!(parse_cml(&source).is_err(), "accepted {field}={value}");
            assert!(
                parse_cml_strict(&source).is_err(),
                "strict accepted {field}={value}"
            );
        }
    }
}

#[test]
fn duplicate_decoded_atom_ids_and_conflicting_isotope_aliases_refuse() {
    for id in ["a1", "a&#49;", "a&#x31;"] {
        let source = SOURCE.replace(
            "<atom id=\"h1\"",
            &format!("<atom id=\"{id}\" elementType=\"C\"/>\n<atom id=\"h1\""),
        );
        assert!(
            parse_cml(&source)
                .err()
                .unwrap()
                .to_string()
                .contains("duplicate CML atom id")
        );
    }
    let conflict = SOURCE.replace(
        "isotopeNumber=\"15\"",
        "isotope=\"14\" isotopeNumber=\"15\"",
    );
    assert!(parse_cml(&conflict).is_err());
    let same = SOURCE.replace(
        "isotopeNumber=\"15\"",
        "isotope=\"015\" isotopeNumber=\"&#49;5\"",
    );
    assert_eq!(
        parse_cml(&same).unwrap().0.atom(AtomIdx(0)).isotope,
        Some(15)
    );
}

#[test]
fn absent_fields_and_numeric_representation_bounds_remain_supported() {
    let (mol, coords) =
        parse_cml("<molecule>\n<atom id=\"a\" elementType=\"He\"/>\n</molecule>").unwrap();
    assert_eq!(mol.atom(AtomIdx(0)).charge, 0);
    assert_eq!(mol.atom(AtomIdx(0)).isotope, None);
    assert_eq!(coords, [(0.0, 0.0)]);
    // Representation bounds, not assertions of chemically valid charge states.
    for charge in i8::MIN..=i8::MAX {
        let source = SOURCE.replace("formalCharge=\"1\"", &format!("formalCharge=\"{charge}\""));
        assert_eq!(
            parse_cml(&source).unwrap().0.atom(AtomIdx(0)).charge,
            charge
        );
    }
    for isotope in [0, 1, 65535] {
        let source = SOURCE.replace(
            "isotopeNumber=\"15\"",
            &format!("isotopeNumber=\"{isotope}\""),
        );
        assert_eq!(
            parse_cml(&source).unwrap().0.atom(AtomIdx(0)).isotope,
            (isotope != 0).then_some(isotope)
        );
    }
}

#[test]
fn shared_cdxml_attribute_reader_does_not_decode_literal_entities_twice() {
    let source = "<CDXML Name=\"α &amp;lt; β &amp;gt; γ &#x1f642;\"><page id=\"1\"/></CDXML>";
    let doc = chematic_mol::CdxmlDocument::parse(source).unwrap();
    assert_eq!(
        doc.document_attributes["Name"],
        serde_json::json!("α &lt; β &gt; γ 🙂")
    );
    assert_eq!(doc.write(), source);
}
