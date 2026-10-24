use chematic_mol::*;
use serde_json::Value;
use std::io::Cursor;

fn packet() -> Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/xyz_numeric_read/cases.json"
    )))
    .unwrap()
}

fn check(result: Result<XyzFrame, XyzError>, case: &Value, field: usize) {
    if let Some(bits) = case["bits"].as_str() {
        let frame = result.unwrap();
        let value = match field {
            0..=2 => [frame.atoms[0].x, frame.atoms[0].y, frame.atoms[0].z][field],
            3 => frame.lattice.unwrap()[4],
            _ => match frame.properties[0].values[0][1] {
                XyzValue::Real(v) => v,
                _ => panic!("real column"),
            },
        };
        assert_eq!(value.to_bits(), u64::from_str_radix(bits, 16).unwrap());
    } else {
        assert!(result.is_err(), "{} field {field}", case["id"]);
    }
}

fn source(token: &str, field: usize) -> String {
    let mut coords = ["0", "0", "0"];
    if field < 3 {
        coords[field] = token;
    }
    let header = match field {
        3 => format!("Lattice=\"1 0 0 0 {token} 0 0 0 1\""),
        4 => "Properties=species:S:1:pos:R:3:force:R:3".into(),
        _ => "number".into(),
    };
    let extra = if field == 4 {
        format!(" 0 {token} 0")
    } else {
        String::new()
    };
    format!("1\n{header}\nC {}{extra}\n", coords.join(" "))
}

#[test]
fn xyz_memory_file_and_batch_paths_agree() {
    let packet = packet();
    for group in ["valid", "invalid"] {
        for case in packet[group].as_array().unwrap() {
            for axis in 0..3 {
                let text = source(case["token"].as_str().unwrap(), axis);
                check(parse_xyz(&text), case, axis);
                check(parse_xyz_all(&text).map(|mut f| f.remove(0)), case, axis);
                check(XyzReader::new(&text).next().unwrap(), case, axis);
                check(
                    XyzFileReader::new(Cursor::new(&text)).next().unwrap(),
                    case,
                    axis,
                );
                let batch = XyzBatchReader::new(Cursor::new(&text), 2).next().unwrap();
                check(batch.frames.into_iter().next().unwrap(), case, axis);
            }
        }
    }
}

#[test]
fn extxyz_coordinates_lattice_and_real_properties_agree_across_readers() {
    let packet = packet();
    for group in ["valid", "invalid"] {
        for case in packet[group].as_array().unwrap() {
            for field in 0..5 {
                let text = source(case["token"].as_str().unwrap(), field);
                check(parse_extxyz(&text), case, field);
                check(
                    parse_extxyz_all(&text).map(|mut f| f.remove(0)),
                    case,
                    field,
                );
                check(ExtxyzReader::new(&text).next().unwrap(), case, field);
                check(
                    ExtxyzFileReader::new(Cursor::new(&text)).next().unwrap(),
                    case,
                    field,
                );
                let batch = ExtxyzBatchReader::new(Cursor::new(&text), 2)
                    .next()
                    .unwrap();
                check(batch.frames.into_iter().next().unwrap(), case, field);
            }
        }
    }
}

#[test]
fn errors_retain_field_and_source_context() {
    assert!(matches!(
        parse_xyz(&source("NaN", 1)),
        Err(XyzError::NonFiniteCoordinate {
            line: 3,
            axis: "y",
            ..
        })
    ));
    assert!(matches!(parse_xyz(&source("1e-9999", 1)),
        Err(XyzError::InvalidAtomLine { line: 3, detail }) if detail.contains("underflows") && detail.contains("y")));
    assert!(matches!(parse_extxyz(&source("1e-9999", 3)),
        Err(XyzError::InvalidLattice { line: 2, detail }) if detail.contains("underflows")));
    assert!(matches!(parse_extxyz(&source("1e-9999", 4)),
        Err(XyzError::InvalidPropertyValue { line: 3, property, raw }) if property == "force" && raw == "1e-9999"));
}

#[test]
fn native_trajectory_readers_stop_on_the_bad_frame() {
    let input = source("1", 0) + &source("1e-9999", 0) + &source("2", 0);
    let mut xyz = XyzReader::new(&input);
    assert!(xyz.next().unwrap().is_ok());
    assert!(xyz.next().unwrap().is_err());
    assert!(xyz.next().is_none());
    let mut ext = ExtxyzFileReader::new(Cursor::new(input));
    assert!(ext.next().unwrap().is_ok());
    assert!(ext.next().unwrap().is_err());
    assert!(ext.next().is_none());
}

#[test]
fn strings_and_raw_info_are_not_reinterpreted_as_reals() {
    let text = "1\nProperties=species:S:1:pos:R:3:label:S:1 energy=1e-9999\nC 0 0 0 1e-9999\n";
    let frame = parse_extxyz(text).unwrap();
    assert_eq!(
        frame.properties[0].values[0][0],
        XyzValue::Str("1e-9999".into())
    );
    assert_eq!(frame.info, [("energy".into(), "1e-9999".into())]);
    let output = write_extxyz(&frame).unwrap();
    assert_eq!(parse_extxyz(&output).unwrap().properties, frame.properties);
}
