use chematic_core::AtomIdx;
use chematic_mol::{
    MolParseError, parse_cjson, parse_mol_v3000_with_coords, parse_mol_v3000_with_coords_3d,
    read_mol_v3000_with_diagnostics, write_cjson,
};
use serde_json::Value;

const PACKET: &str = include_str!("../../../validation/v3000_coordinate_projection/cases.json");
fn cases() -> Vec<Value> {
    serde_json::from_str::<Value>(PACKET).unwrap()["cases"]
        .as_array()
        .unwrap()
        .clone()
}

#[test]
fn source_xyz_and_indexed_graph_survive_cjson_projection() {
    for case in cases() {
        let source = case["input"].as_str().unwrap();
        let expected = &case["expected"];
        for text in [
            source.to_owned(),
            source.replace('\n', "\r\n"),
            source.trim_end_matches('\n').to_owned(),
        ] {
            let (mol, metadata, coords) = parse_mol_v3000_with_coords_3d(&text).unwrap();
            assert_eq!(metadata.name, case["id"].as_str().unwrap());
            assert_eq!(coords.len(), mol.atom_count());
            for (i, &(x, y, z)) in coords.iter().enumerate() {
                for (axis, v) in [x, y, z].into_iter().enumerate() {
                    let wanted: f64 = expected["coordinates"][i][axis]
                        .as_str()
                        .unwrap()
                        .parse()
                        .unwrap();
                    assert_eq!(
                        v.to_bits(),
                        wanted.to_bits(),
                        "{} row {i} axis {axis}",
                        case["id"]
                    );
                }
                let atom = mol.atom(AtomIdx(i as u32));
                assert_eq!(
                    u64::from(atom.element.atomic_number()),
                    expected["numbers"][i].as_u64().unwrap()
                );
                assert_eq!(
                    i64::from(atom.charge),
                    expected["formal_charges"][i].as_i64().unwrap()
                );
                assert_eq!(
                    u64::from(atom.isotope.unwrap_or(0)),
                    expected["isotopes"][i].as_u64().unwrap()
                );
            }
            let output = write_cjson(&mol, &coords);
            let value: Value = serde_json::from_str(&output).unwrap();
            assert_eq!(value["atoms"]["elements"]["number"], expected["numbers"]);
            let mut bonds = Vec::new();
            for (i, pair) in value["bonds"]["connections"]["index"]
                .as_array()
                .unwrap()
                .as_chunks::<2>()
                .0
                .iter()
                .enumerate()
            {
                bonds.push(serde_json::json!([
                    pair[0],
                    pair[1],
                    value["bonds"]["order"][i]
                ]));
            }
            assert_eq!(Value::Array(bonds), expected["bonds"]);
            let (_, reread) = parse_cjson(&output).unwrap();
            for (a, b) in coords.iter().zip(reread) {
                assert_eq!(
                    [a.0.to_bits(), a.1.to_bits(), a.2.to_bits()],
                    [b.0.to_bits(), b.1.to_bits(), b.2.to_bits()]
                );
            }
            let (_, _, xy) = parse_mol_v3000_with_coords(&text).unwrap();
            for ((x, y), xyz) in xy.iter().zip(&coords) {
                assert_eq!(
                    [x.to_bits(), y.to_bits()],
                    [xyz.0.to_bits(), xyz.1.to_bits()]
                );
            }
        }
    }
}

#[test]
fn source_z_is_independent_of_perception_conformer() {
    let packet = cases();
    let input = packet[0]["input"].as_str().unwrap();
    let report = read_mol_v3000_with_diagnostics(input).unwrap();
    assert!(report.conformer.is_none());
    let (_, _, xyz) = parse_mol_v3000_with_coords_3d(input).unwrap();
    assert_ne!(xyz[0].2, 0.0);
}

#[test]
fn malformed_coordinates_and_input_overruns_fail() {
    let packet = cases();
    let input = packet[0]["input"].as_str().unwrap();
    for bad in [
        input.replace("M  END", ""),
        input.replace("0.123456789012345", "NaN"),
        input.replace("0.123456789012345", "inf"),
        input.replace("31 1 1 2", "31 1 1 99"),
    ] {
        assert!(parse_mol_v3000_with_coords_3d(&bad).is_err());
    }
    assert!(matches!(
        parse_mol_v3000_with_coords_3d(&" ".repeat(16 * 1024 * 1024 + 1)),
        Err(MolParseError::ResourceLimit { .. })
    ));
}
