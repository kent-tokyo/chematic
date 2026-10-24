use chematic_core::{Coords3D, Point3};
use chematic_mol::{
    parse_mol_v3000_with_coords_3d, write_mol_v3000_with_conformer,
    write_mol_v3000_with_conformer_checked,
};
use serde_json::Value;

#[test]
fn finite_coordinates_survive_both_native_3d_writers_and_three_regenerations() {
    let packet: Value = serde_json::from_str(include_str!(
        "../../../validation/v3000_3d_numeric/cases.json"
    ))
    .unwrap();
    let mut continuation_seen = false;
    for c in packet["cases"].as_array().unwrap() {
        let expected: Vec<[u64; 3]> = c["expected_coordinates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| {
                std::array::from_fn(|i| v[i].as_str().unwrap().parse::<f64>().unwrap().to_bits())
            })
            .collect();
        let mut source = c["input"].as_str().unwrap().to_string();
        for _ in 0..3 {
            let (mol, meta, coords) = parse_mol_v3000_with_coords_3d(&source).unwrap();
            let bits: Vec<_> = coords
                .iter()
                .map(|&(x, y, z)| [x.to_bits(), y.to_bits(), z.to_bits()])
                .collect();
            assert_eq!(bits, expected, "{}", c["id"]);
            let conformer = Coords3D {
                points: coords
                    .into_iter()
                    .map(|(x, y, z)| Point3::new(x, y, z))
                    .collect(),
            };
            source = write_mol_v3000_with_conformer(&mol, &meta, &conformer);
            assert_eq!(
                source,
                write_mol_v3000_with_conformer_checked(&mol, &meta, &conformer).unwrap()
            );
            assert!(
                source
                    .lines()
                    .filter(|l| l.starts_with("M  V30 "))
                    .all(|l| l.len() <= 80),
                "{}",
                c["id"]
            );
            continuation_seen |= source.contains(" -\nM  V30 ");
            let (_, _, reparsed) = parse_mol_v3000_with_coords_3d(&source).unwrap();
            let actual: Vec<_> = reparsed
                .iter()
                .map(|&(x, y, z)| [x.to_bits(), y.to_bits(), z.to_bits()])
                .collect();
            assert_eq!(actual, expected);
        }
    }
    assert!(continuation_seen);
}

#[test]
fn exact_short_four_decimal_coordinates_keep_the_old_atom_row() {
    let mol = chematic_smiles::parse("C").unwrap();
    let conformer = Coords3D {
        points: vec![Point3::new(0.125, -2.5, 3.0001)],
    };
    let output = write_mol_v3000_with_conformer(&mol, &Default::default(), &conformer);
    assert!(output.contains("M  V30 1 C 0.1250 -2.5000 3.0001 0\n"));
}
