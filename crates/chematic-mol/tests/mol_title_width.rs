//! Open Babel #2112 regression: header names must fit one 80-byte line.
use chematic_core::{Coords3D, Point3};
use chematic_mol::mol2000::*;
use chematic_mol::mol3000::*;
use chematic_smiles::parse;
use std::collections::HashMap;

fn outputs(name: &str) -> Vec<String> {
    let mol = parse("CO").unwrap();
    let meta = MolMetadata::default().with_name(name);
    let before = meta.clone();
    let coords = [(0.0, 0.0), (1.4, 0.2)];
    let conf = Coords3D {
        points: vec![Point3::new(0.0, 0.0, 0.5), Point3::new(1.4, 0.2, 0.8)],
    };
    let props = HashMap::from([("ORIGINAL_TITLE".into(), "retained SD field".into())]);
    let mut reused = "old buffer".to_owned();
    write_mol_with_coords_into(&mut reused, &mol, &meta, &coords);
    let mut out = vec![
        write_mol(&mol, &meta),
        write_mol_with_coords(&mol, &meta, &coords),
        reused,
        write_mol_with_stereo_report(&mol, &meta, &coords).0,
        write_laid_out_mol(&mol, &meta).0,
        write_mol_with_conformer(&mol, &meta, &conf),
        write_mol_with_conformer_checked(&mol, &meta, &conf).unwrap(),
        write_mol_v3000(&mol, &meta, &coords),
        write_mol_v3000_with_stereo_report(&mol, &meta, &coords).0,
        write_mol_v3000_with_conformer(&mol, &meta, &conf),
        write_mol_v3000_with_conformer_checked(&mol, &meta, &conf).unwrap(),
        write_sdf_record(&mol, &meta, &coords, &props),
        write_sdf_record_v3000(&mol, &meta, &coords, &props),
        write_sdf_record_with_conformer(&mol, &meta, &conf, &props),
        write_sdf_record_with_conformer_checked(&mol, &meta, &conf, &props).unwrap(),
        write_sdf(&[(&mol, &meta, &coords)]),
        write_sdf_with_charges(&[(&mol, &meta, &coords, &[0.1, -0.1])]),
    ];
    let mut sdf_buffer = String::from("stale");
    write_sdf_record_into(&mut sdf_buffer, &mol, &meta, &coords, &props);
    out.push(sdf_buffer);
    assert_eq!(
        meta, before,
        "serialization must not change source metadata"
    );
    out
}

fn check(name: &str, expected: &str) {
    for block in outputs(name) {
        assert!(
            block.lines().next().unwrap().len() <= 80,
            "title exceeds 80 bytes"
        );
        assert_eq!(block.lines().next().unwrap(), expected);
        let mol_block = block.split_once("M  END\n").unwrap().0.to_owned() + "M  END\n";
        let (mol, meta) = if block.lines().nth(3).unwrap().contains("V3000") {
            parse_mol_v3000(&mol_block).unwrap()
        } else {
            parse_mol(&mol_block).unwrap()
        };
        assert_eq!(meta.name, expected.trim());
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }
}

#[test]
fn empty_title_keeps_blank_line() {
    check("", "");
}
#[test]
fn ordinary_title_unchanged() {
    check("methanol", "methanol");
}
#[test]
fn eighty_ascii_bytes_unchanged() {
    check(&"X".repeat(80), &"X".repeat(80));
}
#[test]
fn long_ascii_title_is_bounded() {
    check(&"0123456789".repeat(10), &"0123456789".repeat(8));
}
#[test]
fn utf8_is_not_split() {
    check(&"炭".repeat(30), &"炭".repeat(26));
}
#[test]
fn utf8_boundary_keeps_complete_prefix() {
    check(&("X".repeat(79) + "炭Z"), &"X".repeat(79));
}
#[test]
fn controls_cannot_inject_header_lines() {
    check("a\r\nb\tc\0d\u{7f}e", "a  b c d e");
}
#[test]
fn unicode_line_separators_are_spaces() {
    check("a\u{85}b\u{2028}c\u{2029}d", "a b c d");
}
#[test]
fn very_long_input_is_bounded() {
    check(&"Z".repeat(1_000_000), &"Z".repeat(80));
}
#[test]
fn reused_buffer_has_no_previous_title_tail() {
    let mol = parse("C").unwrap();
    let mut out = String::new();
    write_mol_with_coords_into(
        &mut out,
        &mol,
        &MolMetadata::default().with_name(&"X".repeat(100)),
        &[],
    );
    write_mol_with_coords_into(
        &mut out,
        &mol,
        &MolMetadata::default().with_name("new"),
        &[],
    );
    assert_eq!(out.lines().next().unwrap(), "new");
}

#[test]
fn emit_independent_reader_cases() {
    let cases = [
        ("plain".to_owned(), "plain".to_owned()),
        ("".to_owned(), "".to_owned()),
        ("0".repeat(80), "0".repeat(80)),
        ("0123456789".repeat(10), "0123456789".repeat(8)),
        ("炭".repeat(30), "炭".repeat(26)),
        ("X".repeat(79) + "炭Z", "X".repeat(79)),
        ("a\r\nb\tc\0d".to_owned(), "a  b c d".to_owned()),
    ];
    for (source, expected) in cases {
        for (writer, block) in outputs(&source).into_iter().enumerate() {
            println!(
                "TITLE_ORACLE {}",
                serde_json::json!({"writer": writer, "source_title": source, "expected_title": expected, "block": block})
            );
        }
    }
}
