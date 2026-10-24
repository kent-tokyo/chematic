use chematic_core::{Coords3D, Molecule, Point3};
use chematic_mol::*;
use serde_json::{Value, json};
use std::collections::HashMap;

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/v2000_isotope_output/cases.json"
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
fn absolute_labels_survive_all_native_v2000_and_sdf_writers() {
    for c in packet()["cases"].as_array().unwrap() {
        let mol = chematic_smiles::parse(c["smiles"].as_str().unwrap()).unwrap();
        assert_eq!(fields(&mol), c["expected"], "input {}", c["id"]);
        let meta = MolMetadata::default();
        let coords: Vec<_> = (0..mol.atom_count()).map(|i| (i as f64, 0.0)).collect();
        let conformer = Coords3D {
            points: coords
                .iter()
                .map(|&(x, y)| Point3::new(x, y, 0.5))
                .collect(),
        };
        let mut buffer = "old contents".to_string();
        chematic_mol::mol2000::write_mol_with_coords_into(&mut buffer, &mol, &meta, &coords);
        for output in [
            write_mol(&mol, &meta),
            write_mol_with_coords(&mol, &meta, &coords),
            buffer,
            write_mol_with_stereo_report(&mol, &meta, &coords).0,
            write_mol_with_conformer(&mol, &meta, &conformer),
            write_mol_with_conformer_checked(&mol, &meta, &conformer).unwrap(),
        ] {
            let (after, _) = parse_mol(&output).unwrap();
            assert_eq!(fields(&after), c["expected"], "{}\n{output}", c["id"]);
            for line in output.lines().skip(4).take(mol.atom_count()) {
                let difference: i16 = line[34..36].trim().parse().unwrap();
                assert!((-3..=4).contains(&difference));
            }
        }
        let props = HashMap::from([("marker".to_string(), "kept".to_string())]);
        for sdf in [
            write_sdf_record(&mol, &meta, &coords, &props),
            write_sdf_record_with_conformer(&mol, &meta, &conformer, &props),
            write_sdf_record_with_conformer_checked(&mol, &meta, &conformer, &props).unwrap(),
            write_sdf(&[(&mol, &meta, &coords)]),
        ] {
            let records = parse_sdf_with_coords(&sdf).unwrap();
            assert_eq!(records.len(), 1);
            assert_eq!(fields(&records[0].0), c["expected"], "{}", c["id"]);
        }
    }
}

#[test]
fn iso_records_include_all_labels_and_split_after_eight_pairs() {
    for c in packet()["cases"].as_array().unwrap() {
        let mol = chematic_smiles::parse(c["smiles"].as_str().unwrap()).unwrap();
        let output = write_mol(&mol, &MolMetadata::default());
        let actual: Vec<[u16; 2]> = output
            .lines()
            .filter(|s| s.starts_with("M  ISO"))
            .flat_map(|line| {
                let count: usize = line[6..9].trim().parse().unwrap();
                assert!((1..=8).contains(&count));
                assert_eq!(line.len(), 9 + count * 8);
                (0..count)
                    .map(|i| {
                        let start = 9 + 8 * i;
                        [
                            line[start..start + 4].trim().parse().unwrap(),
                            line[start + 4..start + 8].trim().parse().unwrap(),
                        ]
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        let expected: Vec<_> = mol
            .atoms()
            .filter_map(|(idx, a)| a.isotope.map(|mass| [idx.0 as u16 + 1, mass]))
            .collect();
        assert_eq!(actual, expected, "{}", c["id"]);
    }
}
