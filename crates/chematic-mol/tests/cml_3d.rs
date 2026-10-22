use chematic_core::{AtomIdx, Element};
use chematic_mol::{CmlParseLimits, parse_cml_with_3d, parse_cml_with_3d_with_limits};
use serde_json::Value;
fn packet() -> Value {
    serde_json::from_str(include_str!("../../../validation/cml_3d/cases.json")).unwrap()
}
fn layouts(input: &str) -> [String; 3] {
    [
        input.to_owned(),
        input.replace('\n', "\r\n"),
        input.lines().map(str::trim).collect(),
    ]
}
#[test]
fn source_order_and_cartesian_bits_survive_all_layouts() {
    for case in packet()["valid"].as_array().unwrap() {
        for input in layouts(case["input"].as_str().unwrap()) {
            let (mol, coords) = parse_cml_with_3d(&input).unwrap();
            let coords = coords.unwrap();
            for (actual, expected) in coords.iter().flatten().zip(
                case["coords"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .flat_map(|r| r.as_array().unwrap()),
            ) {
                assert_eq!(
                    actual.to_bits(),
                    expected.as_f64().unwrap().to_bits(),
                    "{}",
                    case["id"]
                );
            }
            assert_eq!(coords.len(), case["coords"].as_array().unwrap().len());
            let numbers: Vec<_> = mol
                .atoms()
                .map(|(_, a)| a.element.atomic_number())
                .collect();
            assert_eq!(serde_json::to_value(numbers).unwrap(), case["numbers"]);
            let indices: Vec<_> = mol
                .bonds()
                .flat_map(|(_, b)| [b.atom1.0, b.atom2.0])
                .collect();
            assert_eq!(serde_json::to_value(indices).unwrap(), case["bond_indices"]);
            let output = chematic_mol::write_cjson(
                &mol,
                &coords
                    .iter()
                    .map(|p| (p[0], p[1], p[2]))
                    .collect::<Vec<_>>(),
            );
            let (_, back) = chematic_mol::parse_cjson(&output).unwrap();
            for (actual, expected) in back
                .iter()
                .flat_map(|&(x, y, z)| [x, y, z])
                .zip(coords.iter().flatten())
            {
                assert_eq!(actual.to_bits(), expected.to_bits());
            }
        }
    }
}
#[test]
fn absent_geometry_stays_absent_and_legacy_2d_stays_compatible() {
    for case in packet()["absent"].as_array().unwrap() {
        for input in layouts(case["input"].as_str().unwrap()) {
            assert!(parse_cml_with_3d(&input).unwrap().1.is_none());
        }
    }
    let input = packet()["absent"][1]["input"].as_str().unwrap().to_owned();
    assert_eq!(chematic_mol::parse_cml(&input).unwrap().1, vec![(1.5, 2.5)]);
}
#[test]
fn invalid_geometry_and_ambiguous_correspondence_refuse() {
    for case in packet()["invalid"].as_array().unwrap() {
        for input in layouts(case["input"].as_str().unwrap()) {
            assert!(parse_cml_with_3d(&input).is_err(), "{}", case["id"]);
        }
    }
}
#[test]
fn named_writer_three_dimensional_output_can_be_read_back() {
    let mol = chematic_smiles::parse("[13CH3][NH3+]").unwrap();
    let coords = [[-0.0, 1.4175225513503917, 2.125], [1.25, 2.5, -3.75]];
    let text = chematic_mol::write_cml_named_checked(
        &mol,
        &["carbon".into(), "nitrogen".into()],
        Some(&coords),
    )
    .unwrap();
    let (back, positions) = parse_cml_with_3d(&text).unwrap();
    for (a, b) in coords
        .iter()
        .flatten()
        .zip(positions.unwrap().iter().flatten())
    {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(back.atom(AtomIdx(0)).isotope, Some(13));
    assert_eq!(back.atom(AtomIdx(1)).element, Element::N);
    assert_eq!(back.atom(AtomIdx(1)).charge, 1);
}
#[test]
fn bounds_are_inclusive_on_physical_input_and_normalized_projection() {
    let input = packet()["valid"][0]["input"].as_str().unwrap().to_owned();
    let exact = CmlParseLimits {
        max_input_bytes: input.len(),
        max_line_bytes: input.lines().map(str::len).max().unwrap(),
        max_atoms: 3,
        max_bonds: 2,
        ..Default::default()
    };
    assert!(parse_cml_with_3d_with_limits(&input, &exact).is_ok());
    let mut bad = Vec::new();
    let mut limit = exact;
    limit.max_input_bytes -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.max_line_bytes -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.max_atoms -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.max_bonds -= 1;
    bad.push(limit);
    let mut limit = exact;
    limit.max_lines = 1;
    bad.push(limit);
    let mut limit = exact;
    limit.max_xml_elements = 1;
    bad.push(limit);
    for limit in bad {
        assert!(matches!(
            parse_cml_with_3d_with_limits(&input, &limit),
            Err(chematic_mol::CmlError::ResourceLimit { .. })
        ));
    }
    let compact = input.lines().map(str::trim).collect::<String>();
    for limit in [
        CmlParseLimits {
            max_lines: 1,
            ..Default::default()
        },
        CmlParseLimits {
            max_xml_elements: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            parse_cml_with_3d_with_limits(&compact, &limit),
            Err(chematic_mol::CmlError::ResourceLimit { .. })
        ));
    }
}
