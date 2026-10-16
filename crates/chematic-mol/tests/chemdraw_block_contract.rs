use chematic_core::BondOrder;
use chematic_mol::{
    ChemDrawError, ChemDrawFormat as Format, ChemDrawLimits, ChemDrawOptions,
    parse_chemdraw_block as read, parse_chemdraw_block_with_limits,
};

fn bytes(hex: &str) -> Vec<u8> {
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}
fn cases() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../validation/rdkit_issues/chemdraw-block/cases.json"
    ))
    .unwrap()
}
fn prop(tag: u16, data: &[u8]) -> Vec<u8> {
    [
        tag.to_le_bytes().as_slice(),
        (data.len() as u16).to_le_bytes().as_slice(),
        data,
    ]
    .concat()
}
fn obj(tag: u16, id: u32, data: &[u8]) -> Vec<u8> {
    [
        tag.to_le_bytes().as_slice(),
        id.to_le_bytes().as_slice(),
        data,
        &[0, 0],
    ]
    .concat()
}
fn simple(node: &[u8], bond: &[u8]) -> Vec<u8> {
    let pos = prop(0x0200, &[0; 8]);
    let a = obj(0x8004, 10, &[pos.as_slice(), node].concat());
    let b = obj(0x8004, 20, &pos);
    let e = obj(
        0x8005,
        30,
        &[
            prop(0x0604, &10u32.to_le_bytes()).as_slice(),
            prop(0x0605, &20u32.to_le_bytes()).as_slice(),
            bond,
        ]
        .concat(),
    );
    let f = obj(0x8003, 3, &[a, b, e].concat());
    let p = obj(0x8001, 2, &f);
    let d = obj(0x8000, 1, &p);
    [
        b"VjCD0100\x04\x03\x02\x01".as_slice(),
        &[0; 16],
        &d,
        &[0, 0],
    ]
    .concat()
}
#[test]
fn all_authored_formats_preserve_indexed_graphs_and_coordinates() {
    let data = cases();
    assert_eq!(data["cases"].as_array().unwrap().len(), 56);
    for row in data["cases"].as_array().unwrap() {
        let b = bytes(row["cdx_hex"].as_str().unwrap());
        let xml = row["cdxml"].as_str().unwrap();
        let options = ChemDrawOptions {
            allow_metadata_loss: row["allow_metadata_loss"].as_bool().unwrap(),
        };
        let x = read(xml.as_bytes(), Format::Cdxml, &options).unwrap_or_else(|e| panic!("{e}"));
        for format in [Format::Auto, Format::Cdx] {
            let result =
                read(&b, format, &options).unwrap_or_else(|e| panic!("{}: {e}", row["id"]));
            assert_eq!(result.format, Format::Cdx);
            assert_eq!(result.source, b);
            assert_eq!(result.molecules.len(), x.molecules.len());
            assert_eq!(
                result.omissions.len(),
                usize::from(options.allow_metadata_loss)
            );
            for ((a, ac), (e, ec)) in result.molecules.iter().zip(&x.molecules) {
                assert_eq!(ac, ec);
                assert_eq!(
                    a.atoms()
                        .map(|(_, v)| (v.element, v.charge, v.isotope, v.hydrogen_count))
                        .collect::<Vec<_>>(),
                    e.atoms()
                        .map(|(_, v)| (v.element, v.charge, v.isotope, v.hydrogen_count))
                        .collect::<Vec<_>>()
                );
                assert_eq!(
                    a.bonds()
                        .map(|(_, v)| (v.atom1, v.atom2, v.order))
                        .collect::<Vec<_>>(),
                    e.bonds()
                        .map(|(_, v)| (v.atom1, v.atom2, v.order))
                        .collect::<Vec<_>>()
                );
                for (_, edge) in a.bonds().filter(|(_, v)| v.order == BondOrder::Aromatic) {
                    assert!(a.atom(edge.atom1).aromatic && a.atom(edge.atom2).aromatic);
                }
            }
        }
        let result = read(xml.as_bytes(), Format::Auto, &options).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(result.format, Format::Cdxml);
        assert_eq!(result.projected_cdxml, xml);
        assert_eq!(result.source, xml.as_bytes());
    }
}
#[test]
fn metadata_is_never_silently_omitted() {
    let data = cases();
    let row = data["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|r| r["allow_metadata_loss"] == true)
        .unwrap();
    let b = bytes(row["cdx_hex"].as_str().unwrap());
    assert!(matches!(
        read(&b, Format::Auto, &ChemDrawOptions::default()),
        Err(ChemDrawError::MetadataLoss { tag: 0x0805, .. })
    ));
    let result = read(
        &b,
        Format::Auto,
        &ChemDrawOptions {
            allow_metadata_loss: true,
        },
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result.omissions[0].tag, 0x0805);
    assert_eq!(result.source, b);
}
#[test]
fn unsupported_chemistry_still_refuses_with_metadata_consent() {
    for (node, bond) in [
        (prop(0x0400, &2u16.to_le_bytes()), vec![]),
        (prop(0x0422, &[2]), vec![]),
        (prop(0x0446, &[1]), vec![]),
        (prop(0x4000, &[1, 2]), vec![]),
        (vec![], prop(0x0600, &3u16.to_le_bytes())),
        (vec![], prop(0x0601, &6u16.to_le_bytes())),
        (vec![], prop(0x0602, &8u16.to_le_bytes())),
    ] {
        let input = simple(&node, &bond);
        assert!(matches!(
            read(
                &input,
                Format::Cdx,
                &ChemDrawOptions {
                    allow_metadata_loss: true
                }
            ),
            Err(ChemDrawError::Unsupported { .. })
        ));
    }
}
#[test]
fn malformed_records_and_every_truncation_are_errors() {
    let input = simple(&[], &[]);
    for n in 0..input.len() {
        assert!(
            read(&input[..n], Format::Cdx, &ChemDrawOptions::default()).is_err(),
            "prefix {n}"
        );
    }
    for (n, b) in [(0, b'X'), (8, 0), (12, 1), (input.len() - 1, 1)] {
        let mut bad = input.clone();
        bad[n] = b;
        assert!(read(&bad, Format::Cdx, &ChemDrawOptions::default()).is_err());
    }
    let mut trailing = input.clone();
    trailing.push(0);
    assert!(read(&trailing, Format::Cdx, &ChemDrawOptions::default()).is_err());
    for extra in [
        prop(0x0200, &[0; 8]),
        prop(0x0421, &[1, 2]),
        prop(0x0420, &u16::MAX.to_le_bytes()),
        prop(0x0402, &119u16.to_le_bytes()),
        prop(0x042b, &256u16.to_le_bytes()),
    ] {
        assert!(
            read(
                &simple(&extra, &[]),
                Format::Cdx,
                &ChemDrawOptions::default()
            )
            .is_err()
        );
    }
    assert!(read(b"<CDXML/>", Format::Cdx, &ChemDrawOptions::default()).is_err());
    assert!(read(&input, Format::Cdxml, &ChemDrawOptions::default()).is_err());
}
#[test]
fn limits_cover_binary_and_delegated_xml_resources() {
    let data = cases();
    let input = bytes(data["cases"][0]["cdx_hex"].as_str().unwrap());
    let defaults = ChemDrawLimits::default();
    let mut limits = Vec::new();
    let mut l = defaults;
    l.max_input_bytes = input.len() - 1;
    limits.push(l);
    let mut l = defaults;
    l.max_items = 1;
    limits.push(l);
    let mut l = defaults;
    l.max_depth = 2;
    limits.push(l);
    let mut l = defaults;
    l.max_property_bytes = 7;
    limits.push(l);
    let mut l = defaults;
    l.max_output_bytes = 1;
    limits.push(l);
    let mut l = defaults;
    l.cdxml.max_atoms = 2;
    limits.push(l);
    let mut l = defaults;
    l.cdxml.max_bonds = 1;
    limits.push(l);
    let mut l = defaults;
    l.cdxml.max_fragments = 0;
    limits.push(l);
    for l in limits {
        assert!(
            parse_chemdraw_block_with_limits(&input, Format::Auto, &ChemDrawOptions::default(), &l)
                .is_err()
        );
    }
}
#[test]
fn absent_element_defaults_to_carbon_and_format_names_are_explicit() {
    let result = read(&simple(&[], &[]), Format::Auto, &ChemDrawOptions::default())
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result.molecules[0].0.atom_count(), 2);
    assert!(
        result.molecules[0]
            .0
            .atoms()
            .all(|(_, a)| a.element.atomic_number() == 6)
    );
    assert_eq!(Format::from_name("auto"), Some(Format::Auto));
    assert_eq!(Format::from_name("cdx"), Some(Format::Cdx));
    assert_eq!(Format::from_name("mol"), None);
}

#[test]
fn extended_lengths_references_and_xml_structure_are_checked() {
    let input = simple(&[], &[]);
    let needle = prop(0x0200, &[0; 8]);
    let index = input
        .windows(needle.len())
        .position(|w| w == needle)
        .unwrap();
    let mut extended = input.clone();
    extended.splice(index + 2..index + 4, [255, 255, 8, 0, 0, 0]);
    assert!(read(&extended, Format::Cdx, &ChemDrawOptions::default()).is_ok());
    for value in [10u32, 999] {
        let mut bad = input.clone();
        let p = prop(0x0605, &20u32.to_le_bytes());
        let at = bad.windows(p.len()).position(|w| w == p).unwrap();
        bad[at + 4..at + 8].copy_from_slice(&value.to_le_bytes());
        assert!(read(&bad, Format::Cdx, &ChemDrawOptions::default()).is_err());
    }
    let mut duplicate = input.clone();
    let id = [0x04, 0x80, 20, 0, 0, 0];
    let at = duplicate.windows(6).position(|w| w == id).unwrap();
    duplicate[at + 2] = 10;
    assert!(read(&duplicate, Format::Cdx, &ChemDrawOptions::default()).is_err());
    for xml in [
        "not xml",
        "<CDXML><page></CDXML>",
        "<CDXML/><CDXML/>",
        "<other/>",
        "<CDXML>",
    ] {
        assert!(read(xml.as_bytes(), Format::Auto, &ChemDrawOptions::default()).is_err());
    }
    for l in [
        ChemDrawLimits {
            max_depth: 1,
            ..Default::default()
        },
        ChemDrawLimits {
            max_items: 1,
            ..Default::default()
        },
    ] {
        assert!(
            parse_chemdraw_block_with_limits(
                b"<CDXML><page/></CDXML>",
                Format::Auto,
                &ChemDrawOptions::default(),
                &l
            )
            .is_err()
        );
    }
}

#[test]
fn new_projection_retains_alkene_direction_metadata() {
    let data = cases();
    let row = data["cases"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["id"].as_str().unwrap().starts_with("alkene"))
        .unwrap();
    let expected = chematic_smiles::canonical_smiles(&chematic_smiles::parse("F/C=C/F").unwrap());
    for (input, format) in [
        (bytes(row["cdx_hex"].as_str().unwrap()), Format::Cdx),
        (
            row["cdxml"].as_str().unwrap().as_bytes().to_vec(),
            Format::Cdxml,
        ),
    ] {
        let block =
            read(&input, format, &ChemDrawOptions::default()).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(
            chematic_smiles::canonical_smiles(&block.molecules[0].0),
            expected
        );
    }
}

#[test]
fn groups_empty_documents_required_fields_and_text_omissions() {
    let wrap = |body: Vec<u8>| {
        [
            b"VjCD0100\x04\x03\x02\x01".as_slice(),
            &[0; 16],
            &obj(0x8000, 1, &body),
            &[0, 0],
        ]
        .concat()
    };
    let empty = read(&wrap(vec![]), Format::Cdx, &ChemDrawOptions::default())
        .unwrap_or_else(|e| panic!("{e}"));
    assert!(empty.molecules.is_empty());
    for body in [
        obj(0x8004, 10, &[]),
        obj(0x8005, 10, &[]),
        obj(0x8010, 10, &[]),
    ] {
        assert!(
            read(
                &wrap(obj(0x8001, 2, &obj(0x8003, 3, &body))),
                Format::Cdx,
                &ChemDrawOptions {
                    allow_metadata_loss: true
                }
            )
            .is_err()
        );
    }
    let fragment = obj(0x8003, 3, &obj(0x8004, 10, &prop(0x0200, &[0; 8])));
    let text = obj(0x8006, 4, &prop(0x0700, &[0, 0, b'C']));
    let grouped = wrap(obj(0x8001, 2, &obj(0x8002, 5, &[fragment, text].concat())));
    assert!(matches!(
        read(&grouped, Format::Cdx, &ChemDrawOptions::default()),
        Err(ChemDrawError::MetadataLoss { .. })
    ));
    let result = read(
        &grouped,
        Format::Cdx,
        &ChemDrawOptions {
            allow_metadata_loss: true,
        },
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(result.molecules.len(), 1);
    assert_eq!(
        result.omissions.iter().map(|o| o.tag).collect::<Vec<_>>(),
        [0x8006, 0x0700]
    );
    let hidden = wrap(obj(
        0x8001,
        2,
        &obj(0x8006, 4, &obj(0x8004, 10, &prop(0x0200, &[0; 8]))),
    ));
    assert!(
        read(
            &hidden,
            Format::Cdx,
            &ChemDrawOptions {
                allow_metadata_loss: true
            }
        )
        .is_err()
    );
}
