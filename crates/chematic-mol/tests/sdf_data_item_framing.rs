//! Original unique-field packet; see validation/sdf_data_item_framing/cases.json.
use chematic_mol::{
    SdfFileReader, SdfRecord, SdfRecordReader, write_sdf_record, write_sdf_record_v3000,
};
use std::collections::HashMap;
use std::io::{BufReader, Cursor};

fn check(record: &SdfRecord, expected: &HashMap<String, String>) {
    assert_eq!((record.mol.atom_count(), record.mol.bond_count()), (2, 1));
    assert_eq!(&record.properties, expected);
}

#[test]
fn string_stream_fast_and_rewrite_paths_preserve_authored_field_values() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/sdf_data_item_framing/cases.json"
    ))
    .unwrap();
    for (_, ctab) in packet["ctabs"].as_object().unwrap() {
        for case in packet["cases"].as_array().unwrap() {
            let source = format!(
                "{}{}$$$$\n",
                ctab.as_str().unwrap(),
                case["data"].as_str().unwrap()
            );
            let expected: HashMap<String, String> =
                serde_json::from_value(case["expected_properties"].clone()).unwrap();
            for source in [source.clone(), source.replace('\n', "\r\n")] {
                let mut strings = SdfRecordReader::new(&source);
                let record = strings.next().unwrap().unwrap();
                check(&record, &expected);
                assert!(strings.next().is_none());
                for capacity in [1, 17, 8192] {
                    let reader =
                        || BufReader::with_capacity(capacity, Cursor::new(source.as_bytes()));
                    for mut stream in [SdfFileReader::new(reader()), SdfFileReader::fast(reader())]
                    {
                        check(&stream.next().unwrap().unwrap(), &expected);
                        assert!(stream.next().is_none());
                    }
                }
                for output in [
                    write_sdf_record(
                        &record.mol,
                        &record.meta,
                        &record.coords,
                        &record.properties,
                    ),
                    write_sdf_record_v3000(
                        &record.mol,
                        &record.meta,
                        &record.coords,
                        &record.properties,
                    ),
                ] {
                    check(
                        &SdfRecordReader::new(&output).next().unwrap().unwrap(),
                        &expected,
                    );
                }
            }
        }
    }
}

#[test]
fn empty_fields_and_genuine_repeated_keys_keep_the_legacy_map_projection() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/sdf_data_item_framing/cases.json"
    ))
    .unwrap();
    let source = format!(
        "{}> <EMPTY>\n\n> <REPEAT>\nfirst\n\n> <REPEAT>\nlast\n\n$$$$\n",
        packet["ctabs"]["v2000"].as_str().unwrap()
    );
    let record = SdfRecordReader::new(&source).next().unwrap().unwrap();
    assert_eq!(
        record.properties,
        HashMap::from([
            ("EMPTY".into(), "".into()),
            ("REPEAT".into(), "last".into())
        ])
    );
}
