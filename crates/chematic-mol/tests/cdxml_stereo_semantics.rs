//! Independent CDXML screen-coordinate and reference-end expectations.
use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, CipCode};
use chematic_mol::{CdxmlParseOptions, parse_cdxml, parse_cdxml_with_options};

const DRAWING: &str =
    include_str!("../../../validation/rdkit_issues/cdxml/nondirectional-stereo.cdxml");

#[test]
fn screen_y_and_solid_hashed_sense_match_directional_reference_cip() {
    let options = CdxmlParseOptions {
        infer_nondirectional_stereo: true,
    };
    for (display, expected) in [
        ("Bold", CipCode::S),
        ("Hash", CipCode::R),
        ("Dash", CipCode::R),
        ("WedgeBegin", CipCode::S),
        ("WedgedHashBegin", CipCode::R),
    ] {
        let (mol, coords) =
            parse_cdxml_with_options(&DRAWING.replace("Bold", display), &options).unwrap();
        assert_eq!(
            chematic_chem::assign_cip(&mol).get(AtomIdx(0)),
            Some(expected),
            "{display}"
        );
        assert_eq!(
            coords[1],
            (-1., -1.),
            "source Y-down coordinates must be retained"
        );
    }
}

#[test]
fn directional_end_glyphs_belong_to_the_narrow_endpoint() {
    for (display, order, expected) in [
        ("WedgeEnd", BondOrder::Up, CipCode::S),
        ("WedgedHashEnd", BondOrder::Down, CipCode::R),
    ] {
        let text = DRAWING.replace("Bold", display);
        let (mol, _) = parse_cdxml(&text).unwrap();
        assert_eq!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
        assert_eq!(mol.bond(BondIdx(0)).atom1, AtomIdx(1));
        assert_eq!(mol.bond(BondIdx(0)).order, order);
        let reversed = text.replace("B=\"1\" E=\"2\" Display=", "B=\"2\" E=\"1\" Display=");
        let (mol, _) = parse_cdxml(&reversed).unwrap();
        assert_eq!(
            chematic_chem::assign_cip(&mol).get(AtomIdx(0)),
            Some(expected)
        );
    }
}

#[test]
fn directional_bond_does_not_assign_its_wide_endpoint() {
    let text =
        include_str!("../../../validation/rdkit_issues/cdxml/ambiguous-nondirectional.cdxml")
            .replace("Bold", "WedgeBegin");
    let (mol, _) = parse_cdxml(&text).unwrap();
    assert_ne!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
    assert_eq!(mol.atom(AtomIdx(1)).chirality, Chirality::None);
}

#[test]
fn bold_on_a_multiple_bond_does_not_become_a_stereo_marker() {
    let text = DRAWING.replace(
        "B=\"1\" E=\"2\" Display=",
        "B=\"1\" E=\"2\" Order=\"2\" Display=",
    );
    for infer in [false, true] {
        let (mol, _) = parse_cdxml_with_options(
            &text,
            &CdxmlParseOptions {
                infer_nondirectional_stereo: infer,
            },
        )
        .unwrap();
        assert_eq!(mol.bond(BondIdx(0)).order, BondOrder::Double);
        assert_eq!(mol.atom(AtomIdx(0)).chirality, Chirality::None);
    }
}
