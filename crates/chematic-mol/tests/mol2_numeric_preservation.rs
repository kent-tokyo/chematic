use chematic_mol::mol2_source::write_mol2_source;
use chematic_mol::{
    Mol2SectionOrder, parse_mol2_record, write_mol2, write_mol2_record,
    write_mol2_record_with_section_order,
};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/mol2_numeric_preservation/cases.json"
    ))
    .unwrap()
}

fn coordinate_bits(values: &[(f64, f64, f64)]) -> Vec<[u64; 3]> {
    values
        .iter()
        .map(|&(x, y, z)| [x.to_bits(), y.to_bits(), z.to_bits()])
        .collect()
}

#[test]
fn every_native_mol2_writer_retains_authored_numeric_fields() {
    for case in packet()["cases"].as_array().unwrap() {
        let record = parse_mol2_record(case["input"].as_str().unwrap()).unwrap();
        let expected: Vec<[u64; 3]> = case["expected_coordinates"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| {
                std::array::from_fn(|i| row[i].as_str().unwrap().parse::<f64>().unwrap().to_bits())
            })
            .collect();
        let charges: Vec<_> = case["expected_partial_charges"]
            .as_array()
            .unwrap()
            .iter()
            .map(|value| value.as_str().unwrap().parse::<f64>().unwrap().to_bits())
            .collect();
        assert_eq!(coordinate_bits(&record.coords), expected);
        let outputs = [
            write_mol2_record(&record),
            write_mol2_record_with_section_order(&record, Mol2SectionOrder::AtomBondAdjacent),
            write_mol2(&record.molecule, &record.coords),
            write_mol2_source(&record.molecule, &record.coords),
        ];
        for (i, output) in outputs.iter().enumerate() {
            let parsed = parse_mol2_record(output).unwrap();
            assert_eq!(
                coordinate_bits(&parsed.coords),
                expected,
                "{}: {i}",
                case["id"]
            );
            if i < 2 {
                let actual: Vec<_> = parsed
                    .atoms
                    .iter()
                    .map(|a| a.partial_charge.unwrap().to_bits())
                    .collect();
                assert_eq!(actual, charges, "{}: {i}", case["id"]);
            }
        }
    }
}

#[test]
fn edited_values_and_stratified_finite_bits_survive_repeated_regeneration() {
    let fixture = packet();
    let mut record = parse_mol2_record(fixture["cases"][0]["input"].as_str().unwrap()).unwrap();
    // Cover both signs and many exponents/mantissas without another dependency.
    let mut state = 0x517c_c1b7_2722_0a95_u64;
    for _ in 0..512 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let value = f64::from_bits(state);
        if !value.is_finite() {
            continue;
        }
        record.coords[0] = (value, -value, -0.0);
        record.atoms[0].partial_charge = Some(value);
        let expected = coordinate_bits(&record.coords);
        for _ in 0..3 {
            record = parse_mol2_record(&write_mol2_record(&record)).unwrap();
            assert_eq!(coordinate_bits(&record.coords), expected);
            assert_eq!(
                record.atoms[0].partial_charge.unwrap().to_bits(),
                value.to_bits()
            );
        }
    }
}

#[test]
fn exactly_represented_four_decimal_fields_keep_legacy_spelling() {
    let fixture = packet();
    let record = parse_mol2_record(fixture["cases"][0]["input"].as_str().unwrap()).unwrap();
    let output = write_mol2_record(&record);
    let fields: Vec<_> = output
        .lines()
        .find(|line| line.contains("C.3"))
        .unwrap()
        .split_whitespace()
        .collect();
    assert_eq!(&fields[2..5], &["0.1250", "-2.5000", "3.0001"]);
    assert_eq!(fields[8], "0.2500");
}
