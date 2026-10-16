//! Open Babel #2333: uppercase aliases belong only to XYZ species columns.
use chematic_core::{AtomIdx, Element};
use chematic_mol::{
    ExtxyzFileReader, ExtxyzReader, XyzError, XyzFileReader, XyzReader, XyzValue, parse_extxyz,
    parse_xyz, parse_xyz_all, write_extxyz, write_xyz,
};
use std::io::Cursor;

fn periodic_table(uppercase: bool, extended: bool) -> String {
    let comment = if extended {
        "Properties=species:S:1:pos:R:3:name:S:1 label=CL"
    } else {
        "CL BR CA: opaque comment"
    };
    let mut text = format!("118\n{comment}\n");
    for n in 1..=118 {
        let e = Element::from_atomic_number(n).unwrap();
        let symbol = if uppercase {
            e.symbol().to_ascii_uppercase()
        } else {
            e.symbol().to_owned()
        };
        text.push_str(&format!(
            "{symbol} {n} -2.5 1.25{}\n",
            if extended { " CA" } else { "" }
        ));
    }
    text
}

#[test]
fn all_118_elements_have_the_same_identity_in_both_cases_and_native_readers() {
    for uppercase in [false, true] {
        let text = periodic_table(uppercase, false);
        let plain = parse_xyz(&text).unwrap();
        let (mol, coords) = chematic_3d::parse_xyz(&text).unwrap();
        let extended = parse_extxyz(&periodic_table(uppercase, true)).unwrap();
        assert_eq!((mol.atom_count(), mol.bond_count()), (118, 0));
        for n in 1..=118u8 {
            let i = usize::from(n - 1);
            assert_eq!(plain.atoms[i].element.atomic_number(), n);
            assert_eq!(extended.atoms[i].element.atomic_number(), n);
            assert_eq!(mol.atom(AtomIdx(i as u32)).element.atomic_number(), n);
            assert_eq!(plain.coords()[i], (f64::from(n), -2.5, 1.25));
            let p = coords.get(AtomIdx(i as u32));
            assert_eq!((p.x, p.y, p.z), plain.coords()[i]);
        }
    }
}

#[test]
fn trajectory_and_buffered_readers_share_the_species_rules() {
    let text = "1\nfirst\nCL 1 2 3\n1\nsecond\nBR -1 -2 -3\n";
    let expected = parse_xyz_all(text).unwrap();
    assert_eq!(
        XyzReader::new(text).collect::<Result<Vec<_>, _>>().unwrap(),
        expected
    );
    assert_eq!(
        XyzFileReader::new(Cursor::new(text))
            .collect::<Result<Vec<_>, _>>()
            .unwrap(),
        expected
    );
    assert_eq!(
        ExtxyzReader::new(text)
            .collect::<Result<Vec<_>, _>>()
            .unwrap(),
        expected
    );
    assert_eq!(
        ExtxyzFileReader::new(Cursor::new(text))
            .collect::<Result<Vec<_>, _>>()
            .unwrap(),
        expected
    );
    assert_eq!(
        expected
            .iter()
            .map(|f| f.atoms[0].element.symbol())
            .collect::<Vec<_>>(),
        ["Cl", "Br"]
    );
}

#[test]
fn only_species_is_normalized_not_comments_or_string_properties() {
    let text = periodic_table(true, true);
    let frame = parse_extxyz(&text).unwrap();
    assert!(frame.comment.contains("label=CL"));
    assert_eq!(frame.info, vec![("label".into(), "CL".into())]);
    assert!(
        frame.properties[0]
            .values
            .iter()
            .all(|row| row == &vec![XyzValue::Str("CA".into())])
    );
    let out = write_extxyz(&frame).unwrap();
    let again = parse_extxyz(&out).unwrap();
    assert_eq!(again.atoms, frame.atoms);
    assert_eq!(again.properties, frame.properties);
    assert_eq!(again.info, frame.info);
}

#[test]
fn writers_emit_canonical_symbols_and_preserve_coordinate_order() {
    let text = periodic_table(true, false);
    let frame = parse_xyz(&text).unwrap();
    let out = write_xyz(&frame);
    assert_eq!(parse_xyz(&out).unwrap(), frame);
    for (n, row) in (1..=118).zip(out.lines().skip(2)) {
        assert_eq!(
            row.split_whitespace().next().unwrap(),
            Element::from_atomic_number(n).unwrap().symbol()
        );
    }
    let (mol, coords) = chematic_3d::parse_xyz(&text).unwrap();
    let out = chematic_3d::write_xyz(&mol, &coords, "canonical");
    assert_eq!(parse_xyz(&out).unwrap().atoms, frame.atoms);
}

#[test]
fn invalid_tokens_are_not_prefix_guessed_or_changed_to_dummy_atoms() {
    for token in [
        "cl", "cL", "c", "ZZ", "CL1", "CLO", "17", "D", "T", "*", "ＣＬ", "Å", "СL",
    ] {
        let text = format!("1\ninvalid\n{token} 0 0 0\n");
        for error in [
            parse_xyz(&text).unwrap_err(),
            parse_extxyz(&text).unwrap_err(),
        ] {
            assert!(
                matches!(error,XyzError::UnknownElement {symbol,line:3} if symbol==token),
                "{token}"
            );
        }
        assert!(
            matches!(chematic_3d::parse_xyz(&text),Err(chematic_3d::XyzError::UnknownElement(symbol)) if symbol==token)
        );
    }
}

#[test]
fn core_symbols_and_smiles_are_still_case_sensitive() {
    for token in ["CL", "BR", "CA", "FE", "CO", "cl"] {
        assert_eq!(Element::from_symbol(token), None);
    }
    assert_eq!(Element::from_symbol("Ca"), Some(Element::CA));
    assert!(chematic_smiles::parse("CL").is_err());
    assert!(chematic_smiles::parse("CCl").is_ok());
    assert!(
        chematic_smiles::parse("c1ccccc1")
            .unwrap()
            .atoms()
            .all(|(_, a)| a.aromatic)
    );
}

#[test]
fn explicit_hydrogens_and_missing_connectivity_keep_existing_semantics() {
    let text = include_str!(
        "../../../validation/openbabel_issues/xyz-elements/chloromethane-uppercase.xyz"
    );
    let frame = parse_xyz(text).unwrap();
    let (mol, coords) = chematic_3d::parse_xyz(text).unwrap();
    assert_eq!(
        frame
            .atoms
            .iter()
            .map(|a| a.element.atomic_number())
            .collect::<Vec<_>>(),
        [6, 17, 1, 1, 1]
    );
    assert_eq!(
        (mol.atom_count(), mol.bond_count(), coords.points.len()),
        (5, 0, 5)
    );
    assert_eq!(frame.to_molecule().bond_count(), 0);
}
