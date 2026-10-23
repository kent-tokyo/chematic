use chematic_mol::{
    MolParseError, parse_mol_with_coords, parse_mol_with_coords_3d, read_mol_with_diagnostics,
};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/v2000_coordinate_projection/cases.json"
    ))
    .unwrap()
}

#[test]
fn every_source_triplet_survives_independently_of_perception() {
    for case in packet()["cases"].as_array().unwrap() {
        let input = case["input"].as_str().unwrap();
        for text in [
            input.to_owned(),
            input.replace('\n', "\r\n"),
            input.trim_end_matches('\n').to_owned(),
        ] {
            let (mol, meta, coords) = parse_mol_with_coords_3d(&text).unwrap();
            assert_eq!(meta.name, case["id"]);
            let expected: Vec<Vec<f64>> =
                serde_json::from_value(case["expected"]["coordinates"].clone()).unwrap();
            assert_eq!(coords.len(), mol.atom_count());
            for ((x, y, z), row) in coords.iter().zip(&expected) {
                assert_eq!(
                    [x.to_bits(), y.to_bits(), z.to_bits()],
                    [row[0].to_bits(), row[1].to_bits(), row[2].to_bits()],
                    "{}",
                    meta.name
                );
            }
            let (legacy, old_meta, xy) = parse_mol_with_coords(&text).unwrap();
            assert_eq!(old_meta, meta);
            assert_eq!(
                chematic_smiles::canonical_smiles(&legacy),
                chematic_smiles::canonical_smiles(&mol)
            );
            for ((x, y), (a, b, _)) in xy.iter().zip(&coords) {
                assert_eq!([x.to_bits(), y.to_bits()], [a.to_bits(), b.to_bits()]);
            }
        }
    }
    let tiny = packet()["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"] == "tiny-z-declared-2D")
        .unwrap()["input"]
        .as_str()
        .unwrap()
        .to_owned();
    assert!(
        read_mol_with_diagnostics(&tiny)
            .unwrap()
            .conformer
            .is_none()
    );
    assert_ne!(parse_mol_with_coords_3d(&tiny).unwrap().2[0].2, 0.);
}

#[test]
fn malformed_or_nonfinite_fields_are_not_replaced_by_zero() {
    let source = packet()["cases"][0]["input"].as_str().unwrap().to_owned();
    for axis in 0..3 {
        for raw in ["", "NaN", "inf", "-inf", "1e999", "not-a-num"] {
            let mut lines: Vec<_> = source.lines().map(str::to_owned).collect();
            lines[4].replace_range(axis * 10..axis * 10 + 10, &format!("{raw:>10}"));
            let error = match parse_mol_with_coords_3d(&lines.join("\n")) {
                Err(error) => error,
                Ok(_) => panic!("malformed coordinate was accepted"),
            };
            assert!(matches!(
                error,
                MolParseError::InvalidAtomLine { line: 5, .. }
            ));
        }
    }
}

#[test]
fn empty_truncated_wrong_version_and_oversized_inputs_remain_bounded() {
    for input in [
        "",
        "a\nb\nc\n",
        "a\nb\nc\n  2  0  0  0  0  0            999 V2000\n",
    ] {
        assert!(parse_mol_with_coords_3d(input).is_err());
    }
    let source = packet()["cases"][0]["input"].as_str().unwrap().to_owned();
    assert!(parse_mol_with_coords_3d(&source.replace("V2000", "V3000")).is_err());
    assert!(matches!(
        parse_mol_with_coords_3d(&" ".repeat(16 * 1024 * 1024 + 1)),
        Err(MolParseError::ResourceLimit { .. })
    ));
    let (_, _, coords) =
        parse_mol_with_coords_3d("empty\n\n\n  0  0  0  0  0  0            999 V2000\nM  END\n")
            .unwrap();
    assert!(coords.is_empty());
}
