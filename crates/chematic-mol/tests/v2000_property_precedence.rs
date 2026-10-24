use chematic_core::Molecule;
use chematic_mol::{
    SdfFileReader, SdfReader, parse_mol, parse_mol_with_coords, read_mol_with_diagnostics,
};
use serde_json::{Value, json};
use std::io::{BufReader, Cursor};

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/v2000_property_precedence/cases.json"
    ))
    .unwrap()
}

fn fields(mol: &Molecule) -> Value {
    let mut bonds: Vec<_> = mol
        .bonds()
        .map(|(_, b)| {
            [
                b.atom1.0.min(b.atom2.0),
                b.atom1.0.max(b.atom2.0),
                b.order.order_value().unwrap() as u32,
            ]
        })
        .collect();
    bonds.sort();
    json!({
        "numbers": mol.atoms().map(|(_, a)| a.element.atomic_number()).collect::<Vec<_>>(),
        "charges": mol.atoms().map(|(_, a)| a.charge).collect::<Vec<_>>(),
        "isotopes": mol.atoms().map(|(_, a)| a.isotope.unwrap_or(0)).collect::<Vec<_>>(),
        "atom_maps": mol.atoms().map(|(_, a)| a.atom_map.unwrap_or(0)).collect::<Vec<_>>(),
        "bonds": bonds,
    })
}

#[test]
fn property_channels_replace_atom_block_values_in_all_native_reader_modes() {
    for case in packet()["cases"].as_array().unwrap() {
        let source = case["input"].as_str().unwrap();
        for input in [source.to_string(), source.replace('\n', "\r\n")] {
            let sdf = format!("{input}$$$$\n");
            let ordinary = SdfFileReader::new(BufReader::new(Cursor::new(sdf.as_bytes())))
                .next()
                .unwrap()
                .unwrap();
            let fast = SdfFileReader::fast(BufReader::new(Cursor::new(sdf.as_bytes())))
                .next()
                .unwrap()
                .unwrap();
            let string_record = SdfReader::new(&sdf).next().unwrap().unwrap().0;
            for mol in [
                parse_mol(&input).unwrap().0,
                parse_mol_with_coords(&input).unwrap().0,
                read_mol_with_diagnostics(&input).unwrap().mol,
                ordinary.mol,
                fast.mol,
                string_record,
            ] {
                assert_eq!(fields(&mol), case["expected"], "{}", case["id"]);
            }
        }
    }
}

#[test]
fn property_presence_never_leaks_between_sdf_records_or_sd_data() {
    let packet = packet();
    let cases = packet["cases"].as_array().unwrap();
    let sdf: String = cases
        .iter()
        .chain(cases.iter().rev())
        .map(|c| {
            // Apparent property records after M END are opaque SD field data.
            format!(
                "{}> <note>\nM  CHG  1   1   3\nM  ISO  1   1  99\n\n$$$$\n",
                c["input"].as_str().unwrap()
            )
        })
        .collect();
    let expected: Vec<_> = cases
        .iter()
        .chain(cases.iter().rev())
        .map(|c| c["expected"].clone())
        .collect();
    let normal: Vec<_> = SdfFileReader::new(BufReader::new(Cursor::new(sdf.as_bytes())))
        .map(|r| fields(&r.unwrap().mol))
        .collect();
    let fast: Vec<_> = SdfFileReader::fast(BufReader::new(Cursor::new(sdf.as_bytes())))
        .map(|r| fields(&r.unwrap().mol))
        .collect();
    assert_eq!(normal, expected);
    assert_eq!(fast, expected);
}

#[test]
fn accepted_legacy_empty_property_records_still_have_channel_presence() {
    let source = packet()["cases"][1]["input"].as_str().unwrap().to_string();
    for (record, charges, isotopes) in [
        ("M  CHG  0", vec![0, 0], vec![15, 18]),
        ("M  ISO  0", vec![1, -1], vec![0, 0]),
    ] {
        let input = source.replace("M  END", &format!("{record}\nM  END"));
        let actual = fields(&parse_mol(&input).unwrap().0);
        assert_eq!(actual["charges"], json!(charges));
        assert_eq!(actual["isotopes"], json!(isotopes));
    }
}
